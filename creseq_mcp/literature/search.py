"""
creseq_mcp/literature/library.py
===========================

Stats and interpretation tools for CRE-seq / MPRA activity data.

This module operates AFTER library QC. It takes barcode-level or element-level
DNA/RNA count tables, computes log2 RNA/DNA activity, calls active elements,
ranks CRE candidates, and prepares top hits/motifs for downstream RAG-style
literature interpretation.

All public functions return:
    tuple[pd.DataFrame, dict]

This matches the existing QC tool convention.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import re
import time
import requests
import xml.etree.ElementTree as ET


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def _read_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    if not path.exists():
        raise ValueError(f"Path does not exist: {path}")

    name = path.name.lower()

    if name.endswith(".tsv") or name.endswith(".txt"):
        return pd.read_csv(path, sep="\t")
    if name.endswith(".tsv.gz") or name.endswith(".txt.gz"):
        return pd.read_csv(path, sep="\t", compression="gzip")
    if name.endswith(".csv"):
        return pd.read_csv(path)
    if name.endswith(".csv.gz"):
        return pd.read_csv(path, compression="gzip")

    raise ValueError(f"Unsupported file type: {path.suffix}. Use CSV/TSV/TXT, optionally gzipped.")

def _check_cols(df: pd.DataFrame, required: set[str], source: str) -> None:
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"{source}: missing required columns {missing}. "
            f"Found columns: {list(df.columns)}"
        )


_TF_SYNONYMS: dict[str, list[str]] = {
    "AP1": ["AP1", "AP-1", "FOS", "JUN", "FOSL1", "FOSL2", "JUNB", "JUND"],
    "NFE2L2": ["NFE2L2", "NRF2", "NF-E2-related factor 2"],
    "NRF2": ["NRF2", "NFE2L2", "NF-E2-related factor 2"],
    "HNF4A": ["HNF4A", "HNF4 alpha", "hepatocyte nuclear factor 4 alpha"],
    "HNF1A": ["HNF1A", "HNF1 alpha", "hepatocyte nuclear factor 1 alpha"],
    "CEBPA": ["CEBPA", "C/EBP alpha", "CCAAT enhancer binding protein alpha"],
    "CEBPB": ["CEBPB", "C/EBP beta", "CCAAT enhancer binding protein beta"],
    "FOXA1": ["FOXA1", "forkhead box A1"],
    "FOXA2": ["FOXA2", "forkhead box A2", "HNF3 beta"],
    "GATA1": ["GATA1", "GATA binding protein 1"],
    "GATA4": ["GATA4", "GATA binding protein 4"],
    "SP1": ["SP1", "specificity protein 1"],
    "CTCF": ["CTCF", "CCCTC-binding factor"],
    "YY1": ["YY1", "Yin Yang 1"],
    "KLF4": ["KLF4", "Kruppel-like factor 4"],
    "TEAD": ["TEAD", "TEAD1", "TEAD2", "TEAD3", "TEAD4"],
    "ETS": ["ETS", "ELK1", "ELK4", "ERG", "ETV1", "ETV4"],
    "RUNX": ["RUNX", "RUNX1", "RUNX2", "RUNX3"],
    "STAT": ["STAT", "STAT1", "STAT3", "STAT5A", "STAT5B"],
    "E2F": ["E2F", "E2F1", "E2F3", "E2F4"],
    "NFY": ["NFY", "NFYA", "NFYB", "NFYC"],
    "IRF": ["IRF", "IRF1", "IRF3", "IRF4", "IRF8"],
    "MAF": ["MAF", "MAFF", "MAFG", "MAFK"],
    "BACH": ["BACH", "BACH1", "BACH2"],
}

_REGULATORY_TERMS = [
    "enhancer",
    "promoter",
    "cis-regulatory element",
    "regulatory element",
    "regulatory sequence",
    "transcriptional regulation",
]

_ASSAY_TERMS = [
    "MPRA",
    "massively parallel reporter assay",
    "CRE-seq",
    "lentiMPRA",
    "STARR-seq",
    "reporter assay",
]

_CAUSAL_REGULATORY_TERMS = [
    "activates",
    "activation",
    "represses",
    "repression",
    "required",
    "drives",
    "necessary",
    "controls",
    "regulates",
]

_CELL_TYPE_SYNONYMS: dict[str, list[str]] = {
    "HEPG2": ["HepG2", "hepatocyte", "hepatocytes", "liver cell", "hepatic cell", "hepatocellular carcinoma"],
    "K562": ["K562", "erythroid", "erythroleukemia", "leukemia cell"],
    "WTC11": ["WTC11", "iPSC", "induced pluripotent stem cell"],
    "HEK293": ["HEK293", "293T", "human embryonic kidney"],
    "HCT116": ["HCT116", "colon cancer", "colorectal carcinoma"],
}

_EVIDENCE_STRENGTH_MAP = {
    "MPRA": 1.0,
    "Reporter_assay": 0.9,
    "Perturbation": 0.8,
    "ChIP_binding": 0.6,
    "Motif_analysis": 0.4,
    "Review": 0.2,
    "Database": 0.3,
    "Unknown": 0.1,
}

_MAX_SCORE = 2.0 + 1.5 + 2.5 + 1.2 + 1.0 + 0.5


def _unique_terms(terms: list[str]) -> list[str]:
    seen = set()
    out = []
    for term in terms:
        key = term.strip().lower()
        if key and key not in seen:
            seen.add(key)
            out.append(term.strip())
    return out


def expand_tf_terms(tf_name: str) -> list[str]:
    """Return query synonyms for a TF/motif name, preserving the input term first."""
    tf = str(tf_name).strip()
    synonyms = _TF_SYNONYMS.get(tf.upper(), [])
    return _unique_terms([tf, *synonyms])


def expand_tf_terms_detailed(tf_name: str) -> dict[str, Any]:
    """
    Return weighted TF synonym groups.

    Exact symbols get the highest weight, aliases lower weight, and family-like
    terms the lowest.  The simple ``expand_tf_terms`` list API is preserved for
    existing callers.
    """
    tf = str(tf_name).strip()
    aliases = [term for term in _TF_SYNONYMS.get(tf.upper(), []) if term.strip().lower() != tf.lower()]
    family_terms = [
        term for term in aliases
        if tf.upper() in {"AP1", "TEAD", "ETS", "RUNX", "STAT", "E2F", "NFY", "IRF", "MAF", "BACH"}
        and _term_key(term) == _term_key(tf)
    ]
    alias_terms = [term for term in aliases if term not in family_terms]
    all_terms = _unique_terms([tf, *alias_terms, *family_terms])
    weights = {tf: 1.0}
    weights.update({term: 0.8 for term in alias_terms})
    weights.update({term: 0.5 for term in family_terms})
    return {
        "symbols": [tf] if tf else [],
        "aliases": alias_terms,
        "family_terms": family_terms,
        "all_terms": all_terms,
        "weights": weights,
    }


def expand_cell_type_terms(cell_type: str | None) -> dict[str, Any]:
    """Return lightweight cell-type synonyms and weights."""
    if not cell_type:
        return {"symbols": [], "aliases": [], "all_terms": [], "weights": {}}
    cell = str(cell_type).strip()
    aliases = [term for term in _CELL_TYPE_SYNONYMS.get(cell.upper(), []) if term.strip().lower() != cell.lower()]
    all_terms = _unique_terms([cell, *aliases])
    weights = {cell: 1.0}
    weights.update({term: 0.7 for term in aliases})
    return {
        "symbols": [cell],
        "aliases": aliases,
        "all_terms": all_terms,
        "weights": weights,
    }


def _pubmed_field_group(terms: list[str], field: str = "Title/Abstract") -> str:
    escaped = []
    for term in _unique_terms(terms):
        if " " in term or "-" in term or "/" in term:
            escaped.append(f'"{term}"[{field}]')
        else:
            escaped.append(f"{term}[{field}]")
    return "(" + " OR ".join(escaped) + ")"


def build_queries(
    tf: str,
    cell_type: str | None = None,
    species: str = "human",
) -> dict[str, str]:
    """
    Build multi-intent PubMed queries for TF regulatory evidence.

    The returned queries use OR within concept groups and AND across concepts.
    They intentionally keep TF terms in Title/Abstract for precision.
    """
    tf_group = _pubmed_field_group(expand_tf_terms_detailed(tf)["all_terms"])
    regulatory_group = _pubmed_field_group(["enhancer", "cis-regulatory", "promoter"])
    species_terms = [species] if species else ["human"]
    if str(species or "").lower() in {"human", "homo sapiens", "9606"}:
        species_terms = ["human", "Homo sapiens"]
    species_group = _pubmed_field_group(species_terms)
    exclusion = "NOT " + _pubmed_field_group(["yeast", "drosophila"])

    concept_groups = [tf_group]
    cell_terms = expand_cell_type_terms(cell_type)["all_terms"]
    if cell_terms:
        concept_groups.append(_pubmed_field_group(cell_terms))
    concept_groups.extend([regulatory_group, species_group])

    intent_filters = {
        "mpra": _pubmed_field_group(["MPRA", "massively parallel reporter assay", "luciferase"]),
        "binding": _pubmed_field_group(["ChIP-seq", "CUT&RUN", "CUT&Tag"]),
        "motif": _pubmed_field_group(["motif", "binding motif", "transcription factor binding site"]),
        "perturbation": _pubmed_field_group(["knockdown", "knockout", "CRISPR", "overexpression"]),
    }

    base = " AND ".join(concept_groups)
    return {
        intent: f"{base} AND {intent_filter} {exclusion}"
        for intent, intent_filter in intent_filters.items()
    }


def build_pubmed_query(
    tf_name: str,
    target_cell_type: str | None = None,
    off_target_cell_type: str | None = None,
    species: str = "human",
) -> str:
    """
    Build a synonym-aware PubMed query for TF/motif evidence.

    The query keeps the TF requirement strict but broadens the biology terms so
    relevant enhancer/reporter papers are not missed just because they omit MPRA
    or a specific cell type from the title/abstract.
    """

    query = build_queries(tf_name, target_cell_type, species=species)["mpra"]

    if off_target_cell_type:
        query = f"{query} NOT {_pubmed_field_group([off_target_cell_type])}"

    return query


def _text_has_any(text: str, terms: list[str]) -> bool:
    haystack = str(text or "").lower()
    return any(term.lower() in haystack for term in terms if term)


def _term_key(term: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(term or "").lower())


def _is_tf_name_match(candidate: Any, tf_name: str) -> bool:
    candidate_key = _term_key(candidate)
    if not candidate_key:
        return False
    for term in expand_tf_terms(tf_name):
        term_key = _term_key(term)
        if term_key and (candidate_key == term_key or candidate_key.startswith(term_key)):
            return True
    return False


def _extract_pub_year(pubdate: Any) -> float | None:
    match = re.search(r"\b(19|20)\d{2}\b", str(pubdate or ""))
    return float(match.group(0)) if match else None


def _record_text(record: pd.Series | dict[str, Any]) -> str:
    get = record.get
    return ". ".join(
        _clean_space(get(field))
        for field in ("title", "abstract")
        if _clean_space(get(field))
    )


def classify_evidence(record: pd.Series | dict[str, Any]) -> str:
    """
    Classify the biological evidence type from title/abstract keywords.

    The order is intentional: MPRA and reporter assays are most direct for
    CRE-seq interpretation, followed by perturbation and binding evidence.
    """
    text = _record_text(record).lower()
    source = str(record.get("source") or "")
    if "review" in text:
        return "Review"
    if _text_has_any(text, ["mpra", "massively parallel reporter assay", "mpractive", "lentimpra"]):
        return "MPRA"
    if _text_has_any(text, ["luciferase", "reporter assay", "dual reporter", "enhancer assay"]):
        return "Reporter_assay"
    if _text_has_any(text, ["knockdown", "knockout", "crispr", "overexpression", "sirna", "shrna"]):
        return "Perturbation"
    if _text_has_any(text, ["chip-seq", "chip seq", "cut&run", "cut&tag", "cut and run", "cut and tag"]):
        return "ChIP_binding"
    if _text_has_any(text, ["motif", "binding site", "transcription factor binding site", "position weight matrix"]):
        return "Motif_analysis"
    if source in {"JASPAR", "ENCODE"}:
        return "Database"
    return "Unknown"


def _confidence_from_score(normalized_score: float) -> str:
    if normalized_score > 0.85:
        return "high"
    if normalized_score >= 0.6:
        return "medium"
    return "low"


def compute_encode_support(
    tf: str,
    region: str | None = None,
    cell_type: str | None = None,
    record: pd.Series | dict[str, Any] | None = None,
    motif_present: bool = False,
) -> dict[str, Any]:
    """
    Compute lightweight ENCODE support from available metadata.

    This is conservative because this project does not currently pass genomic
    regions or peak files into the literature layer.  If future rows include
    ``peak_overlap`` or ``signal_strength``, those fields are used.
    """
    if record is None:
        record = {}
    target = record.get("target_label") or record.get("tf_name") or ""
    biosample = record.get("biosample_term_name") or record.get("cell_type_query") or ""
    direct_binding = _is_tf_name_match(target, tf)

    cell_terms = expand_cell_type_terms(cell_type)["all_terms"] if cell_type else []
    if cell_type and str(biosample).lower() == str(cell_type).lower():
        cell_match = "exact"
    elif cell_terms and _text_has_any(str(biosample), cell_terms):
        cell_match = "partial"
    else:
        cell_match = "none"

    peak_overlap = bool(record.get("peak_overlap", False))
    if region and str(record.get("region", "")) == str(region):
        peak_overlap = True

    signal_strength = pd.to_numeric(record.get("signal_strength", record.get("signal", 0.0)), errors="coerce")
    signal_strength = float(0.0 if pd.isna(signal_strength) else signal_strength)

    return {
        "direct_binding": bool(direct_binding),
        "peak_overlap": bool(peak_overlap),
        "signal_strength": signal_strength,
        "cell_type_match": cell_match,
        "motif_peak_support": bool(motif_present and peak_overlap),
    }


def score_evidence_records(
    evidence_df: pd.DataFrame,
    target_cell_type: str | None = None,
) -> pd.DataFrame:
    """
    Add simple evidence-strength flags and an evidence_score column.

    Scores are intentionally transparent rather than statistical: they help the
    agent sort retrieved records by relevance without hiding why a record ranked
    highly.
    """

    if evidence_df.empty:
        out = evidence_df.copy()
        for col in (
            "tf_match",
            "cell_type_match",
            "regulatory_keyword_match",
            "assay_keyword_match",
            "recency_bonus",
            "evidence_score",
        ):
            out[col] = pd.Series(dtype="float64" if col.endswith(("bonus", "score")) else "bool")
        return out

    out = evidence_df.copy()
    text_cols = [
        c for c in (
            "title",
            "abstract",
            "journal",
            "authors",
            "name",
            "target_label",
            "biosample_term_name",
            "assay_title",
        )
        if c in out.columns
    ]
    combined_text = (
        out[text_cols].fillna("").astype(str).agg(" ".join, axis=1)
        if text_cols
        else pd.Series("", index=out.index)
    )

    tf_matches = []
    for _, row in out.iterrows():
        motif = str(row.get("motif") or row.get("tf_name") or "")
        row_text = " ".join(str(row.get(c, "")) for c in text_cols)
        tf_matches.append(_text_has_any(row_text, expand_tf_terms(motif)))

    out["tf_match"] = tf_matches
    out["cell_type_match"] = (
        combined_text.map(lambda text: _text_has_any(text, expand_cell_type_terms(target_cell_type)["all_terms"]))
        if target_cell_type
        else False
    )
    out["regulatory_keyword_match"] = combined_text.map(lambda text: _text_has_any(text, _REGULATORY_TERMS))
    out["assay_keyword_match"] = combined_text.map(lambda text: _text_has_any(text, _ASSAY_TERMS))
    out["causal_regulatory_keyword_match"] = combined_text.map(
        lambda text: _text_has_any(text, _CAUSAL_REGULATORY_TERMS)
    )
    out["evidence_class"] = out.apply(classify_evidence, axis=1)
    out["evidence_strength"] = out["evidence_class"].map(_EVIDENCE_STRENGTH_MAP).fillna(0.1)

    years = out.get("pubdate", pd.Series("", index=out.index)).map(_extract_pub_year)
    out["recency_bonus"] = years.map(lambda y: 0.5 if y and y >= 2020 else 0.0)
    out["evidence_score"] = (
        out["tf_match"].astype(float) * 2.0
        + out["cell_type_match"].astype(float) * 1.5
        + out["evidence_strength"].astype(float) * 2.5
        + out["causal_regulatory_keyword_match"].astype(float) * 1.2
        + out["assay_keyword_match"].astype(float) * 1.0
        + out["recency_bonus"]
    )
    if "source" in out.columns:
        is_encode = out["source"].astype(str).eq("ENCODE")
        if is_encode.any():
            supports = [
                compute_encode_support(
                    str(row.get("motif") or row.get("tf_name") or ""),
                    cell_type=target_cell_type,
                    record=row,
                    motif_present=bool(row.get("motif")),
                )
                if str(row.get("source")) == "ENCODE"
                else None
                for _, row in out.iterrows()
            ]
            out["encode_support"] = supports
            out.loc[
                is_encode & out["encode_support"].map(lambda x: bool(x and x.get("motif_peak_support"))),
                "evidence_score",
            ] += 1.0

    out["evidence_score_normalized"] = (out["evidence_score"] / _MAX_SCORE).clip(0.0, 1.0)
    out["confidence"] = out["evidence_score_normalized"].map(_confidence_from_score)

    sort_cols = ["evidence_score", "motif"] if "motif" in out.columns else ["evidence_score"]
    ascending = [False, True] if "motif" in out.columns else [False]
    return out.sort_values(sort_cols, ascending=ascending).reset_index(drop=True)


def _clean_space(text: Any) -> str:
    if text is None or pd.isna(text):
        return ""
    return re.sub(r"\s+", " ", str(text or "")).strip()


def _clean_identifier(value: Any) -> str:
    text = _clean_space(value)
    if re.fullmatch(r"\d+\.0", text):
        return text[:-2]
    return text


def _source_id(row: pd.Series) -> str:
    source = str(row.get("source") or "")
    if source == "PubMed" and pd.notna(row.get("pmid")):
        return f"PMID:{_clean_identifier(row.get('pmid'))}"
    if source == "JASPAR" and pd.notna(row.get("matrix_id")):
        return f"JASPAR:{_clean_identifier(row.get('matrix_id'))}"
    if source == "ENCODE" and pd.notna(row.get("accession")):
        return f"ENCODE:{_clean_identifier(row.get('accession'))}"
    return f"{source or 'Evidence'}:{int(row.name) + 1}"


def _citation(row: pd.Series) -> str:
    source = str(row.get("source") or "")
    if source == "PubMed" and pd.notna(row.get("pmid")):
        return f"PMID {_clean_identifier(row.get('pmid'))}"
    if source == "JASPAR" and pd.notna(row.get("matrix_id")):
        return f"JASPAR {_clean_identifier(row.get('matrix_id'))}"
    if source == "ENCODE" and pd.notna(row.get("accession")):
        return f"ENCODE {_clean_identifier(row.get('accession'))}"
    return source or "Evidence record"


def _snippet_from_text(text: str, terms: list[str], max_chars: int) -> str:
    text = _clean_space(text)
    if len(text) <= max_chars:
        return text

    lowered = text.lower()
    match_positions = [
        lowered.find(term.lower())
        for term in terms
        if term and lowered.find(term.lower()) >= 0
    ]
    center = min(match_positions) if match_positions else 0
    start = max(0, center - max_chars // 3)
    end = min(len(text), start + max_chars)
    start = max(0, end - max_chars)
    snippet = text[start:end].strip()
    if start > 0:
        snippet = "... " + snippet
    if end < len(text):
        snippet = snippet + " ..."
    return snippet


def _rag_context_text(row: pd.Series, max_chars: int) -> str:
    source = str(row.get("source") or "")
    title = _clean_space(row.get("title"))
    abstract = _clean_space(row.get("abstract"))
    motif = str(row.get("motif") or row.get("tf_name") or "")

    if source == "PubMed":
        base = abstract or title
        terms = expand_tf_terms(motif) + _REGULATORY_TERMS + _ASSAY_TERMS
        return _snippet_from_text(base, terms, max_chars)

    if source == "JASPAR":
        parts = [
            f"JASPAR motif profile {row.get('matrix_id')}",
            f"name {row.get('name')}" if pd.notna(row.get("name")) else "",
            f"collection {row.get('collection')}" if pd.notna(row.get("collection")) else "",
            f"family {row.get('family')}" if pd.notna(row.get("family")) else "",
        ]
        return _snippet_from_text(". ".join(p for p in parts if p), expand_tf_terms(motif), max_chars)

    if source == "ENCODE":
        parts = [
            f"ENCODE experiment {row.get('accession')}",
            f"target {row.get('target_label')}" if pd.notna(row.get("target_label")) else "",
            f"biosample {row.get('biosample_term_name')}" if pd.notna(row.get("biosample_term_name")) else "",
            f"assay {row.get('assay_title')}" if pd.notna(row.get("assay_title")) else "",
            f"status {row.get('status')}" if pd.notna(row.get("status")) else "",
        ]
        return _snippet_from_text(". ".join(p for p in parts if p), expand_tf_terms(motif), max_chars)

    return _snippet_from_text(" ".join(_clean_space(v) for v in row.dropna().tolist()), [], max_chars)


def extract_claim(record: pd.Series | dict[str, Any]) -> dict[str, Any]:
    """Extract a compact rule-based biological claim from a literature record."""
    text = _record_text(record)
    tf = str(record.get("motif") or record.get("tf_name") or "")
    terms = expand_tf_terms(tf)
    sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]
    selected = ""
    for sentence in sentences:
        if _text_has_any(sentence, terms) and _text_has_any(sentence, _CAUSAL_REGULATORY_TERMS + _REGULATORY_TERMS):
            selected = sentence
            break
    if not selected and sentences:
        selected = sentences[0]
    if not selected:
        source = str(record.get("source") or "")
        if source == "ENCODE":
            biosample = record.get("biosample_term_name")
            selected = (
                f"ENCODE reports {record.get('assay_title') or 'functional genomics'} "
                f"for {record.get('target_label') or tf}"
                f"{' in ' + str(biosample) if pd.notna(biosample) else ''}."
            )
        elif source == "JASPAR":
            selected = (
                f"JASPAR provides motif profile {record.get('matrix_id') or ''} "
                f"for {record.get('name') or tf}."
            )
        else:
            selected = f"{tf} has retrieved regulatory evidence."

    lower = selected.lower()
    if _text_has_any(lower, ["represses", "repression", "suppresses", "inhibits", "decreases", "silences"]):
        direction = "repression"
    elif _text_has_any(lower, ["activates", "activation", "drives", "enhances", "increases", "upregulates"]):
        direction = "activation"
    else:
        direction = "unknown"

    evidence_class = str(record.get("evidence_class") or classify_evidence(record))
    species = "human" if _text_has_any(text, ["human", "homo sapiens", "hepg2", "k562"]) else "unknown"
    return {
        "claim": selected or f"{tf} has retrieved regulatory evidence.",
        "direction": direction,
        "assay": evidence_class,
        "context": selected,
        "species": species,
    }


def _why_relevant(row: pd.Series) -> list[str]:
    reasons = []
    if bool(row.get("tf_match", False)):
        reasons.append("matches TF/synonym")
    if bool(row.get("cell_type_match", False)):
        reasons.append("matches target cell type")
    if bool(row.get("regulatory_keyword_match", False)):
        reasons.append("mentions regulatory context")
    if bool(row.get("assay_keyword_match", False)):
        reasons.append("mentions reporter/MPRA-style assay context")
    if pd.notna(row.get("evidence_class")) and row.get("evidence_class") not in {None, "Unknown"}:
        reasons.append(f"classified as {row.get('evidence_class')}")
    if pd.notna(row.get("query_scope")) and row.get("query_scope") == "tf_regulatory_fallback":
        reasons.append("broader TF/regulatory fallback evidence")
    return reasons if reasons else ["database record for this motif/TF"]


def _claim_key(claim: Any) -> str:
    return re.sub(r"[^a-z0-9]", "", str(claim or "").lower())[:160]


def _balance_rag_records(df: pd.DataFrame, max_records: int, max_per_tf: int = 5) -> pd.DataFrame:
    """Select diverse RAG rows across evidence classes while capping per TF."""
    if df.empty:
        return df

    df = df.sort_values("evidence_score", ascending=False).copy()
    selected_indices: list[int] = []
    selected_keys: set[str] = set()
    per_tf: dict[str, int] = {}

    def add_row(idx: int) -> None:
        if len(selected_indices) >= max_records or idx in selected_indices:
            return
        row = df.loc[idx]
        tf = str(row.get("motif") or row.get("tf_name") or "unknown")
        if per_tf.get(tf, 0) >= max_per_tf:
            return
        key = _claim_key(row.get("claim"))
        if key and key in selected_keys:
            return
        selected_indices.append(idx)
        selected_keys.add(key)
        per_tf[tf] = per_tf.get(tf, 0) + 1

    priority_groups = [
        {"MPRA", "Reporter_assay"},
        {"ChIP_binding"},
        {"Motif_analysis", "Perturbation"},
    ]
    for group in priority_groups:
        candidates = df[df["evidence_class"].isin(group)] if "evidence_class" in df.columns else pd.DataFrame()
        if len(candidates):
            add_row(int(candidates.index[0]))

    for idx in df.index:
        add_row(int(idx))
        if len(selected_indices) >= max_records:
            break

    return df.loc[selected_indices].reset_index(drop=True)



# ---------------------------------------------------------------------------
# rank_cre_candidates
# ---------------------------------------------------------------------------

def rank_cre_candidates(
    activity_table_path: str | Path,
    top_n: int = 20,
    activity_col: str = "log2_ratio",
    q_col: str = "fdr",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Rank CREs by activity strength and statistical confidence.

    Score logic:
        high activity is good
        low q-value is good
        low DNA coverage is penalized if present
    """

    df = _read_table(activity_table_path)
    _check_cols(df, {activity_col}, str(activity_table_path))

    result = df.copy()

    if q_col not in result.columns:
        result[q_col] = 1.0

    result[activity_col] = pd.to_numeric(result[activity_col], errors="coerce").fillna(0)
    result[q_col] = pd.to_numeric(result[q_col], errors="coerce").fillna(1.0).clip(1e-12, 1)

    result["confidence_score"] = -np.log10(result[q_col])
    result["rank_score"] = result[activity_col] + 0.25 * result["confidence_score"]

    if "low_dna_coverage" in result.columns:
        result.loc[result["low_dna_coverage"] == True, "rank_score"] -= 1.0

    result = result.sort_values("rank_score", ascending=False).reset_index(drop=True)
    result["rank"] = np.arange(1, len(result) + 1)

    top = result.head(top_n).copy()

    summary = {
        "n_ranked": int(len(result)),
        "top_n": int(top_n),
        "top_element": str(top.iloc[0].get("element_id", top.iloc[0].get("oligo_id", "unknown")))
        if len(top)
        else None,
        "median_top_activity": float(top[activity_col].median()) if len(top) else None,
        "warnings": [],
        "pass": True,
    }

    return top, summary


# ---------------------------------------------------------------------------
# Tool 4: motif_enrichment_summary
# ---------------------------------------------------------------------------

def motif_enrichment_summary(
    activity_table_path: str | Path,
    motif_col: str = "top_motif",
    active_col: str = "active",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Lightweight motif enrichment summary.

    This does not run MEME/HOMER/FIMO. It summarizes motifs already annotated
    in the table and computes active-vs-inactive enrichment ratios.
    """

    df = _read_table(activity_table_path)
    _check_cols(df, {motif_col, active_col}, str(activity_table_path))

    rows = []
    motifs = sorted(m for m in df[motif_col].dropna().unique() if str(m).lower() != "none")

    active_df = df[df[active_col] == True]
    inactive_df = df[df[active_col] == False]

    for motif in motifs:
        active_rate = float((active_df[motif_col] == motif).mean()) if len(active_df) else 0.0
        inactive_rate = float((inactive_df[motif_col] == motif).mean()) if len(inactive_df) else 0.0

        enrichment = (active_rate + 1e-6) / (inactive_rate + 1e-6)

        rows.append(
            {
                "motif": motif,
                "active_fraction": active_rate,
                "inactive_fraction": inactive_rate,
                "enrichment_ratio": enrichment,
                "n_active_with_motif": int((active_df[motif_col] == motif).sum()),
                "n_inactive_with_motif": int((inactive_df[motif_col] == motif).sum()),
            }
        )

    result = pd.DataFrame(rows).sort_values("enrichment_ratio", ascending=False)

    summary = {
        "n_motifs_tested": int(len(result)),
        "top_enriched_motif": str(result.iloc[0]["motif"]) if len(result) else None,
        "warnings": [],
        "pass": True,
    }

    return result, summary


# ---------------------------------------------------------------------------
# Tool 5: prepare_rag_context
# ---------------------------------------------------------------------------

def prepare_rag_context(
    ranked_table_path: str | Path,
    top_n: int = 10,
    motif_col: str = "top_motif",
    target_cell_type: str | None = None,
    off_target_cell_type: str | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Prepare a compact table that an LLM/RAG layer can use for literature search.

    This does not call external APIs directly. It creates suggested search queries
    for PubMed/JASPAR/ENCODE-style lookups.
    """

    df = _read_table(ranked_table_path)
    top = df.head(top_n).copy()

    queries = []

    if motif_col in top.columns:
        motifs = [m for m in top[motif_col].dropna().unique() if str(m).lower() != "none"]
        for motif in motifs:
            if target_cell_type and off_target_cell_type:
                queries.append(
                    f"{motif} transcription factor {target_cell_type} enhancer activity "
                    f"off target {off_target_cell_type}"
                )
            elif target_cell_type:
                queries.append(f"{motif} transcription factor {target_cell_type} enhancer activity")
            else:
                queries.append(f"{motif} transcription factor enhancer MPRA CRE-seq")

    top["rag_search_terms"] = "; ".join(queries[:5]) if queries else ""

    summary = {
        "n_top_elements": int(len(top)),
        "suggested_queries": queries[:10],
        "target_cell_type": target_cell_type,
        "off_target_cell_type": off_target_cell_type,
        "warnings": [],
        "pass": True,
    }

    return top, summary



# ---------------------------------------------------------------------------
# Literature/API retrieval helpers
# ---------------------------------------------------------------------------

def _safe_get_json(
    url: str,
    params: dict[str, Any] | None = None,
    timeout: int = 15,
    headers: dict[str, str] | None = None,
) -> dict[str, Any]:
    """
    Small helper for GET-based JSON APIs.

    Returns a dictionary with either:
        {"ok": True, "data": ...}
    or:
        {"ok": False, "error": "...", "url": "..."}
    """

    try:
        response = requests.get(
            url,
            params=params,
            timeout=timeout,
            headers=headers or {"Accept": "application/json"},
        )
        response.raise_for_status()
        return {"ok": True, "data": response.json()}
    except Exception as exc:
        return {
            "ok": False,
            "error": str(exc),
            "url": url,
            "params": params or {},
        }


def _safe_get_text(
    url: str,
    params: dict[str, Any] | None = None,
    timeout: int = 15,
    headers: dict[str, str] | None = None,
) -> dict[str, Any]:
    try:
        response = requests.get(
            url,
            params=params,
            timeout=timeout,
            headers=headers or {"Accept": "application/xml"},
        )
        response.raise_for_status()
        return {"ok": True, "data": response.text}
    except Exception as exc:
        return {
            "ok": False,
            "error": str(exc),
            "url": url,
            "params": params or {},
        }


def _parse_pubmed_abstracts(xml_text: str) -> dict[str, str]:
    abstracts: dict[str, str] = {}
    root = ET.fromstring(xml_text)
    for article in root.findall(".//PubmedArticle"):
        pmid_el = article.find(".//PMID")
        if pmid_el is None or not pmid_el.text:
            continue
        parts = []
        for abstract_el in article.findall(".//AbstractText"):
            text = "".join(abstract_el.itertext()).strip()
            label = abstract_el.attrib.get("Label")
            if text and label:
                parts.append(f"{label}: {text}")
            elif text:
                parts.append(text)
        abstracts[pmid_el.text.strip()] = " ".join(parts)
    return abstracts


def _fetch_pubmed_abstracts(
    ids: list[str],
    email: str | None = None,
    api_key: str | None = None,
) -> tuple[dict[str, str], list[str]]:
    if not ids:
        return {}, []

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    params: dict[str, Any] = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "xml",
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key

    response = _safe_get_text(f"{base}/efetch.fcgi", params=params)
    if not response["ok"]:
        return {}, [f"PubMed EFetch failed: {response['error']}"]

    try:
        return _parse_pubmed_abstracts(response["data"]), []
    except Exception as exc:
        return {}, [f"PubMed abstract parsing failed: {exc}"]


def search_pubmed(
    query: str,
    max_results: int = 5,
    email: str | None = None,
    api_key: str | None = None,
    include_abstracts: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Search PubMed using NCBI E-utilities.

    API flow:
        1. ESearch finds PubMed IDs for the query.
        2. ESummary retrieves title/journal/date metadata for those IDs.

    Notes:
        - email is recommended by NCBI for responsible API use.
        - api_key is optional but helps with rate limits.
    """

    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

    esearch_params: dict[str, Any] = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": max_results,
        "sort": "relevance",
    }

    if email:
        esearch_params["email"] = email
    if api_key:
        esearch_params["api_key"] = api_key

    esearch = _safe_get_json(f"{base}/esearch.fcgi", esearch_params)

    if not esearch["ok"]:
        return pd.DataFrame(), {
            "query": query,
            "n_results": 0,
            "warnings": [f"PubMed ESearch failed: {esearch['error']}"],
            "pass": False,
        }

    ids = esearch["data"].get("esearchresult", {}).get("idlist", [])

    if not ids:
        return pd.DataFrame(), {
            "query": query,
            "n_results": 0,
            "warnings": ["No PubMed results found."],
            "pass": True,
        }

    # Be polite between NCBI requests.
    time.sleep(0.34)

    esummary_params: dict[str, Any] = {
        "db": "pubmed",
        "id": ",".join(ids),
        "retmode": "json",
    }

    if email:
        esummary_params["email"] = email
    if api_key:
        esummary_params["api_key"] = api_key

    esummary = _safe_get_json(f"{base}/esummary.fcgi", esummary_params)

    if not esummary["ok"]:
        return pd.DataFrame({"pmid": ids}), {
            "query": query,
            "n_results": len(ids),
            "warnings": [f"PubMed ESummary failed: {esummary['error']}"],
            "pass": False,
        }

    result_obj = esummary["data"].get("result", {})
    abstract_map: dict[str, str] = {}
    warnings_list = []
    if include_abstracts:
        time.sleep(0.34)
        abstract_map, warnings_list = _fetch_pubmed_abstracts(ids, email=email, api_key=api_key)

    rows = []

    for pmid in ids:
        item = result_obj.get(pmid, {})
        rows.append(
            {
                "source": "PubMed",
                "query": query,
                "pmid": pmid,
                "title": item.get("title"),
                "journal": item.get("fulljournalname"),
                "pubdate": item.get("pubdate"),
                "abstract": abstract_map.get(pmid),
                "authors": ", ".join(
                    author.get("name", "")
                    for author in item.get("authors", [])[:5]
                    if isinstance(author, dict)
                ),
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            }
        )

    df = pd.DataFrame(rows)

    summary = {
        "query": query,
        "n_results": int(len(df)),
        "warnings": warnings_list,
        "pass": True,
    }

    return df, summary


def search_jaspar_motif(
    tf_name: str,
    species: int = 9606,
    collection: str = "CORE",
    max_results: int = 5,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Search JASPAR for transcription factor motif matrix profiles.

    species=9606 means human.
    """

    url = "https://jaspar.elixir.no/api/v1/matrix/"

    params = {
        "search": tf_name,
        "species": species,
        "collection": collection,
        "format": "json",
        "page_size": max_results,
    }

    response = _safe_get_json(url, params=params)

    if not response["ok"]:
        return pd.DataFrame(), {
            "tf_name": tf_name,
            "n_results": 0,
            "warnings": [f"JASPAR search failed: {response['error']}"],
            "pass": False,
        }

    data = response["data"]
    results = data.get("results", []) if isinstance(data, dict) else []
    matched_results = [
        item for item in results
        if _is_tf_name_match(item.get("name"), tf_name)
    ]

    rows = []
    for item in matched_results[:max_results]:
        matrix_id = item.get("matrix_id") or item.get("matrix_id_base")
        rows.append(
            {
                "source": "JASPAR",
                "tf_name": tf_name,
                "matrix_id": matrix_id,
                "name": item.get("name"),
                "collection": item.get("collection"),
                "tax_group": item.get("tax_group"),
                "species": item.get("species"),
                "class": item.get("class"),
                "family": item.get("family"),
                "url": f"https://jaspar.elixir.no/matrix/{matrix_id}/" if matrix_id else None,
            }
        )

    df = pd.DataFrame(rows)

    summary = {
        "tf_name": tf_name,
        "n_results": int(len(df)),
        "warnings": [],
        "pass": True,
    }

    if len(df) == 0:
        summary["warnings"].append("No JASPAR motif records found for exact TF/synonym match.")

    return df, summary


def search_encode_tf(
    tf_name: str,
    cell_type: str | None = None,
    max_results: int = 5,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Search ENCODE for TF-related experiment records.

    This is intentionally broad because ENCODE metadata terms vary.
    It works best for TF names like GATA1, TAL1, CTCF, SPI1, HNF4A.
    """

    url = "https://www.encodeproject.org/search/"

    params: dict[str, Any] = {
        "type": "Experiment",
        "searchTerm": tf_name,
        "format": "json",
        "limit": max_results,
    }

    if cell_type:
        params["searchTerm"] = f"{tf_name} {cell_type}"

    response = _safe_get_json(url, params=params)

    if not response["ok"]:
        if "404" in response["error"]:
            return pd.DataFrame(), {
                "tf_name": tf_name,
                "cell_type": cell_type,
                "n_results": 0,
                "warnings": ["No ENCODE experiment records found."],
                "pass": True,
            }
        return pd.DataFrame(), {
            "tf_name": tf_name,
            "cell_type": cell_type,
            "n_results": 0,
            "warnings": [f"ENCODE search failed: {response['error']}"],
            "pass": False,
        }

    graph = response["data"].get("@graph", [])

    rows = []
    for item in graph[:max_results]:
        accession = item.get("accession")
        biosample = item.get("biosample_ontology", {}) or {}
        target = item.get("target", {}) or {}
        target_label = target.get("label") if isinstance(target, dict) else None
        if not _is_tf_name_match(target_label, tf_name):
            continue

        rows.append(
            {
                "source": "ENCODE",
                "tf_name": tf_name,
                "cell_type_query": cell_type,
                "accession": accession,
                "assay_title": item.get("assay_title"),
                "target_label": target_label,
                "biosample_term_name": biosample.get("term_name") if isinstance(biosample, dict) else None,
                "status": item.get("status"),
                "url": f"https://www.encodeproject.org/experiments/{accession}/"
                if accession
                else None,
            }
        )

    df = pd.DataFrame(rows)

    summary = {
        "tf_name": tf_name,
        "cell_type": cell_type,
        "n_results": int(len(df)),
        "warnings": [],
        "pass": True,
    }

    if len(df) == 0:
        summary["warnings"].append("No ENCODE experiment records found for exact TF/synonym match.")

    return df, summary


def literature_search_for_motifs(
    motif_table_path: str | Path,
    motif_col: str = "motif",
    target_cell_type: str | None = None,
    off_target_cell_type: str | None = None,
    species: str = "human",
    top_n_motifs: int = 5,
    max_pubmed_results_per_motif: int = 3,
    max_database_results_per_motif: int = 3,
    email: str | None = None,
    ncbi_api_key: str | None = None,
    output_path: str | Path | None = None,
    multi_intent_queries: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Run actual API-backed literature/database retrieval for enriched TF motifs.

    Inputs:
        motif_table_path:
            Output from motif_enrichment_summary(), or any table with a motif column.

    For each top motif, this queries:
        - PubMed through NCBI E-utilities
        - JASPAR REST API
        - ENCODE REST API

    Output:
        One combined evidence table that can be passed to an LLM for interpretation.
    """

    motif_df = _read_table(motif_table_path)
    if motif_col not in motif_df.columns and motif_col == "motif" and "tf_name" in motif_df.columns:
        motif_col = "tf_name"
    _check_cols(motif_df, {motif_col}, str(motif_table_path))

    if "enrichment_ratio" in motif_df.columns:
        motif_df = motif_df.sort_values("enrichment_ratio", ascending=False)
    elif "odds_ratio" in motif_df.columns:
        motif_df = motif_df.sort_values("odds_ratio", ascending=False)
    elif "fdr" in motif_df.columns:
        motif_df = motif_df.sort_values("fdr", ascending=True)

    motifs = [
        str(m)
        for m in motif_df[motif_col].dropna().head(top_n_motifs).tolist()
        if str(m).strip() and str(m).lower() != "none"
    ]

    all_rows = []
    warnings_list = []

    for motif in motifs:
        query_map = (
            build_queries(motif, target_cell_type, species=species)
            if multi_intent_queries
            else {"mpra": build_pubmed_query(
                motif,
                target_cell_type=target_cell_type,
                off_target_cell_type=off_target_cell_type,
                species=species,
            )}
        )
        if off_target_cell_type and multi_intent_queries:
            query_map = {
                intent: f"{query} NOT {_pubmed_field_group([off_target_cell_type])}"
                for intent, query in query_map.items()
            }

        pubmed_frames = []
        pubmed_any = False
        for intent, pubmed_query in query_map.items():
            pubmed_df, pubmed_summary = search_pubmed(
                query=pubmed_query,
                max_results=max_pubmed_results_per_motif,
                email=email,
                api_key=ncbi_api_key,
            )
            if len(pubmed_df):
                pubmed_any = True
                pubmed_df["motif"] = motif
                pubmed_df["tf_synonyms"] = "; ".join(expand_tf_terms(motif))
                pubmed_df["query_scope"] = "target_cell" if target_cell_type else "tf_regulatory"
                pubmed_df["query_intent"] = intent
                pubmed_df["evidence_type"] = "literature"
                pubmed_frames.append(pubmed_df)
            warnings_list.extend(pubmed_summary.get("warnings", []))

        if not pubmed_any and target_cell_type:
            fallback_query_map = (
                build_queries(motif, None, species=species)
                if multi_intent_queries
                else {"mpra": build_pubmed_query(
                    motif,
                    target_cell_type=None,
                    off_target_cell_type=off_target_cell_type,
                    species=species,
                )}
            )
            if off_target_cell_type and multi_intent_queries:
                fallback_query_map = {
                    intent: f"{query} NOT {_pubmed_field_group([off_target_cell_type])}"
                    for intent, query in fallback_query_map.items()
                }
            for intent, fallback_query in fallback_query_map.items():
                pubmed_df, fallback_summary = search_pubmed(
                    query=fallback_query,
                    max_results=max_pubmed_results_per_motif,
                    email=email,
                    api_key=ncbi_api_key,
                )
                if len(pubmed_df):
                    pubmed_df["motif"] = motif
                    pubmed_df["tf_synonyms"] = "; ".join(expand_tf_terms(motif))
                    pubmed_df["query_scope"] = "tf_regulatory_fallback"
                    pubmed_df["query_intent"] = intent
                    pubmed_df["evidence_type"] = "literature"
                    pubmed_frames.append(pubmed_df)
                warnings_list.extend(fallback_summary.get("warnings", []))
            if pubmed_frames:
                warnings_list.append(
                    f"No target-cell PubMed records for {motif}; used broader TF/regulatory query."
                )

        if pubmed_frames:
            all_rows.append(pd.concat(pubmed_frames, ignore_index=True, sort=False))

        jaspar_df, jaspar_summary = search_jaspar_motif(
            tf_name=motif,
            max_results=max_database_results_per_motif,
        )
        if len(jaspar_df):
            jaspar_df["motif"] = motif
            jaspar_df["tf_synonyms"] = "; ".join(expand_tf_terms(motif))
            jaspar_df["evidence_type"] = "motif_database"
            all_rows.append(jaspar_df)
        warnings_list.extend(jaspar_summary.get("warnings", []))

        encode_df, encode_summary = search_encode_tf(
            tf_name=motif,
            cell_type=target_cell_type,
            max_results=max_database_results_per_motif,
        )
        if len(encode_df):
            encode_df["motif"] = motif
            encode_df["tf_synonyms"] = "; ".join(expand_tf_terms(motif))
            encode_df["evidence_type"] = "functional_genomics"
            all_rows.append(encode_df)
        warnings_list.extend(encode_summary.get("warnings", []))

        # Avoid hammering APIs.
        time.sleep(0.25)

    if all_rows:
        combined = pd.concat(all_rows, ignore_index=True, sort=False)
        if "pmid" in combined.columns:
            dedup_subset = [c for c in ("source", "motif", "pmid") if c in combined.columns]
            pubmed_part = combined[combined["source"].astype(str).eq("PubMed")].drop_duplicates(
                subset=dedup_subset,
                keep="first",
            )
            non_pubmed = combined[~combined["source"].astype(str).eq("PubMed")]
            combined = pd.concat([pubmed_part, non_pubmed], ignore_index=True, sort=False)
    else:
        combined = pd.DataFrame(
            columns=[
                "source",
                "motif",
                "tf_synonyms",
                "evidence_type",
                "evidence_class",
                "query_intent",
                "title",
                "journal",
                "pubdate",
                "matrix_id",
                "accession",
                "url",
            ]
        )

    combined = score_evidence_records(combined, target_cell_type=target_cell_type)

    output_path_str = None
    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        combined.to_csv(output, sep="\t", index=False)
        output_path_str = str(output)

    summary = {
        "motifs_searched": motifs,
        "n_motifs": int(len(motifs)),
        "n_evidence_records": int(len(combined)),
        "target_cell_type": target_cell_type,
        "off_target_cell_type": off_target_cell_type,
        "species": species,
        "output_path": output_path_str,
        "warnings": warnings_list,
        "pass": True,
    }

    return combined, summary


def interpret_literature_evidence(
    evidence_table_path: str | Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Convert retrieved API records into a simple interpretation-ready summary.

    This does not use an LLM. It creates a structured summary that the MCP agent
    or frontend can display.
    """

    df = _read_table(evidence_table_path)
    if "evidence_score" not in df.columns:
        df = score_evidence_records(df)

    if df.empty:
        return df, {
            "n_sources": 0,
            "interpretation": "No literature or database evidence was retrieved.",
            "warnings": ["Evidence table is empty."],
            "pass": False,
        }

    source_counts = df["source"].value_counts().to_dict() if "source" in df.columns else {}

    motif_counts = df["motif"].value_counts().to_dict() if "motif" in df.columns else {}
    motif_scores = (
        df.groupby("motif")["evidence_score"].max().sort_values(ascending=False).to_dict()
        if "motif" in df.columns and "evidence_score" in df.columns
        else {}
    )

    interpretation_parts = []

    if motif_scores:
        top_motif = max(motif_scores, key=motif_scores.get)
        interpretation_parts.append(
            f"The strongest retrieved evidence is for {top_motif} "
            f"(top evidence score {motif_scores[top_motif]:.1f}, "
            f"{motif_counts.get(top_motif, 0)} supporting records)."
        )
    elif motif_counts:
        top_motif = max(motif_counts, key=motif_counts.get)
        interpretation_parts.append(
            f"The retrieved evidence has the most records for {top_motif} "
            f"({motif_counts[top_motif]} supporting records)."
        )

    if "PubMed" in source_counts:
        interpretation_parts.append(
            f"PubMed returned {source_counts['PubMed']} literature records."
        )

    if "JASPAR" in source_counts:
        interpretation_parts.append(
            f"JASPAR returned {source_counts['JASPAR']} motif profile records."
        )

    if "ENCODE" in source_counts:
        interpretation_parts.append(
            f"ENCODE returned {source_counts['ENCODE']} functional genomics records."
        )

    result = pd.DataFrame(
        {
            "metric": ["source_counts", "motif_counts", "motif_top_scores"],
            "value": [str(source_counts), str(motif_counts), str(motif_scores)],
        }
    )

    summary = {
        "n_sources": int(len(source_counts)),
        "source_counts": source_counts,
        "motif_counts": motif_counts,
        "motif_top_scores": motif_scores,
        "interpretation": " ".join(interpretation_parts),
        "warnings": [],
        "pass": True,
    }

    return result, summary


def prepare_literature_rag_context(
    evidence_table_path: str | Path,
    max_records: int = 8,
    min_score: float = 4.0,
    max_context_chars: int = 700,
    output_path: str | Path | None = None,
    max_per_tf: int = 5,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """
    Convert scored evidence records into citation-ready RAG context chunks.

    This function is deterministic: it does not ask an LLM to decide relevance.
    It filters by evidence_score, keeps the best records, and emits compact
    context strings plus source IDs that an LLM can cite.
    """

    df = _read_table(evidence_table_path)
    if "evidence_score" not in df.columns:
        df = score_evidence_records(df)
    elif "evidence_class" not in df.columns:
        df["evidence_class"] = df.apply(classify_evidence, axis=1)
    if "evidence_score_normalized" not in df.columns and "evidence_score" in df.columns:
        df["evidence_score_normalized"] = (
            pd.to_numeric(df["evidence_score"], errors="coerce").fillna(0.0) / _MAX_SCORE
        ).clip(0.0, 1.0)
    if "confidence" not in df.columns and "evidence_score_normalized" in df.columns:
        df["confidence"] = df["evidence_score_normalized"].map(_confidence_from_score)

    if "evidence_score" in df.columns:
        df["evidence_score"] = pd.to_numeric(df["evidence_score"], errors="coerce").fillna(0.0)
        filtered = df[df["evidence_score"] >= float(min_score)].copy()
    else:
        filtered = df.copy()

    if filtered.empty:
        empty = pd.DataFrame(
            columns=[
                "source_id",
                "motif",
                "tf",
                "motif_id",
                "source",
                "evidence_type",
                "source_evidence_type",
                "title",
                "url",
                "citation",
                "claim",
                "direction",
                "cell_type",
                "assay",
                "encode_support",
                "evidence_score",
                "confidence",
                "context",
                "why_relevant",
                "contradicts",
            ]
        )
        output_path_str = None
        if output_path is not None:
            output = Path(output_path)
            output.parent.mkdir(parents=True, exist_ok=True)
            empty.to_csv(output, sep="\t", index=False)
            output_path_str = str(output)

        summary = {
            "n_context_records": 0,
            "n_input_records": int(len(df)),
            "min_score": float(min_score),
            "max_records": int(max_records),
            "max_per_tf": int(max_per_tf),
            "output_path": output_path_str,
            "warnings": ["No evidence records met the RAG context score threshold."],
            "pass": False,
        }
        return empty, summary

    filtered = filtered.sort_values("evidence_score", ascending=False).copy()

    rows = []
    for _, row in filtered.iterrows():
        context = _rag_context_text(row, max_context_chars)
        claim = extract_claim(row)
        tf = row.get("motif") or row.get("tf_name")
        encode_support = row.get("encode_support")
        if not isinstance(encode_support, dict) and str(row.get("source") or "") == "ENCODE":
            encode_support = compute_encode_support(str(tf or ""), cell_type=row.get("cell_type_query"), record=row)
        elif not isinstance(encode_support, dict):
            encode_support = {}
        rows.append(
            {
                "source_id": _source_id(row),
                "motif": row.get("motif"),
                "tf": tf,
                "motif_id": row.get("matrix_id"),
                "source": row.get("source"),
                "evidence_type": row.get("evidence_class") or row.get("evidence_type"),
                "source_evidence_type": row.get("evidence_type"),
                "title": row.get("title") if pd.notna(row.get("title")) else row.get("name"),
                "url": row.get("url"),
                "citation": _citation(row),
                "claim": claim["claim"],
                "direction": claim["direction"],
                "cell_type": row.get("biosample_term_name") or row.get("cell_type_query"),
                "assay": claim["assay"],
                "encode_support": encode_support,
                "evidence_score": float(row.get("evidence_score", 0.0)),
                "evidence_score_normalized": float(row.get("evidence_score_normalized", 0.0)),
                "confidence": row.get("confidence") or _confidence_from_score(float(row.get("evidence_score_normalized", 0.0))),
                "context": context,
                "why_relevant": _why_relevant(row),
                "contradicts": False,
            }
        )

    rag_df = _balance_rag_records(pd.DataFrame(rows), max_records=max_records, max_per_tf=max_per_tf)

    output_path_str = None
    if output_path is not None:
        output = Path(output_path)
        output.parent.mkdir(parents=True, exist_ok=True)
        rag_df.to_csv(output, sep="\t", index=False)
        output_path_str = str(output)

    summary = {
        "n_context_records": int(len(rag_df)),
        "n_input_records": int(len(df)),
        "min_score": float(min_score),
        "max_records": int(max_records),
        "max_per_tf": int(max_per_tf),
        "motifs": sorted(str(m) for m in rag_df["motif"].dropna().unique()),
        "sources": rag_df["source"].value_counts().to_dict() if "source" in rag_df.columns else {},
        "evidence_types": rag_df["evidence_type"].value_counts().to_dict() if "evidence_type" in rag_df.columns else {},
        "output_path": output_path_str,
        "warnings": [],
        "pass": True,
    }

    return rag_df, summary
