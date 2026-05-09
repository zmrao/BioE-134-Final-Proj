# Literature Search + RAG Contributions

This document summarizes my contributions to the CRE-seq MCP final project. My section focused on the literature-search and retrieval-augmented generation (RAG) layer that helps interpret enriched transcription factor motifs from CRE-seq / MPRA results using external biological evidence.

## Project Area

My work is centered in:

- `creseq_mcp/literature/search.py`
- `creseq_mcp/server.py`
- `frontend/agent.py`
- `tests/stats/test_library.py`
- `pyproject.toml`

The literature/RAG pipeline sits after activity calling and motif enrichment. It takes enriched TF motifs, searches biological evidence sources, scores and classifies the evidence, and produces citation-ready context rows for an LLM to use when explaining possible CRE regulators.

## MCP Tools Added Or Extended

I contributed the following MCP-facing tools:

- `tool_rank_cre_candidates`
- `tool_motif_enrichment_summary`
- `tool_prepare_rag_context`
- `tool_search_pubmed`
- `tool_search_jaspar_motif`
- `tool_search_encode_tf`
- `tool_literature_search_for_motifs`
- `tool_interpret_literature_evidence`
- `tool_prepare_literature_rag_context`

These tools expose the literature and interpretation layer through the MCP server so the frontend agent can call them automatically.

## Pipeline Workflow

The literature/RAG workflow is:

1. Rank active CRE candidates from activity results.
2. Summarize motifs enriched among active CREs.
3. Search for supporting evidence for top TF motifs.
4. Retrieve evidence from PubMed, JASPAR, and ENCODE.
5. Classify the evidence type.
6. Score each evidence record for biological relevance.
7. Convert the best records into citation-ready RAG context.
8. Return structured rows that an LLM can cite without inventing unsupported claims.

The main output files are:

- `literature_evidence.tsv`
- `literature_rag_context.tsv`

## Query Generation

I implemented synonym-aware and intent-aware query generation for PubMed.

The query builder expands TF names and aliases, for example:

- `NRF2`
- `NFE2L2`
- `NF-E2-related factor 2`

It also supports lightweight cell-type expansion, for example:

- `HepG2`
- `hepatocyte`
- `liver cell`
- `hepatic cell`
- `hepatocellular carcinoma`

The upgraded query system builds multiple search intents:

- `mpra`: MPRA, massively parallel reporter assay, luciferase evidence
- `binding`: ChIP-seq, CUT&RUN, CUT&Tag evidence
- `motif`: motif and TF binding-site evidence
- `perturbation`: knockdown, knockout, CRISPR, overexpression evidence

This improves retrieval quality because a single PubMed query often misses useful papers. The multi-intent approach separates direct reporter-assay evidence from binding and perturbation evidence.

## External Evidence Sources

I added or extended API-backed searches for:

- PubMed through NCBI E-utilities
- JASPAR motif matrix profiles
- ENCODE functional genomics records

For PubMed, the tool retrieves:

- PMID
- title
- journal
- publication date
- authors
- abstract
- URL

For JASPAR, the tool retrieves:

- matrix ID
- motif/TF name
- collection
- taxonomic group
- family/class metadata
- motif URL

For ENCODE, the tool retrieves:

- accession
- assay title
- target TF label
- biosample/cell type
- release status
- experiment URL

## Evidence Classification

I added rule-based evidence classification so the RAG layer can distinguish stronger biological evidence from weaker database-only evidence.

Evidence records are classified as:

- `MPRA`
- `Reporter_assay`
- `ChIP_binding`
- `Motif_analysis`
- `Perturbation`
- `Review`
- `Database`
- `Unknown`

This matters because direct MPRA or reporter-assay evidence should count more strongly for CRE-seq interpretation than a generic motif database hit.

## Evidence Scoring

I rewrote the evidence scoring logic to prioritize biological relevance.

The score uses:

- TF/synonym match
- target cell-type match
- evidence type strength
- regulatory keywords
- assay keywords
- recency

Evidence type strength is weighted so direct CRE/MPRA evidence ranks highest:

- `MPRA`: strongest
- `Reporter_assay`: strong
- `Perturbation`: useful functional evidence
- `ChIP_binding`: useful binding evidence
- `Motif_analysis`: mechanistic but less direct
- `Review`: background support
- `Database`: useful but lower confidence alone

The scoring output includes:

- `evidence_score`
- `evidence_score_normalized`
- `confidence`

Confidence is mapped to:

- `high`
- `medium`
- `low`

## RAG Context Formatting

I implemented citation-ready RAG rows through `prepare_literature_rag_context`.

Each RAG row includes:

- `source_id`
- `tf`
- `motif_id`
- `source`
- `evidence_type`
- `claim`
- `direction`
- `cell_type`
- `assay`
- `encode_support`
- `evidence_score`
- `confidence`
- `citation`
- `url`
- `context`
- `why_relevant`
- `contradicts`

This gives the LLM enough structured context to explain findings while citing specific evidence records.

## Claim Extraction

I added simple rule-based claim extraction.

The claim extractor looks for sentences containing:

- the TF or a synonym
- regulatory language such as activates, represses, required, drives, necessary, controls, or regulates

It also estimates direction:

- `activation`
- `repression`
- `unknown`

For database records without paper abstracts, it creates clean fallback claims, for example:

- `ENCODE reports TF ChIP-seq for NFE2L2 in HepG2.`
- `JASPAR provides motif profile MA0150.2 for NRF2.`

This prevents the LLM from receiving empty or unhelpful context.

## RAG Balancing

I added balancing logic so the final RAG context is not dominated by one evidence type or one TF.

The formatter tries to include diverse evidence when available:

- at least one MPRA/reporter record
- at least one ChIP/binding record
- at least one motif/mechanistic record

It also:

- caps records per TF
- deduplicates highly similar claims
- sorts records by evidence quality

This makes the final context more useful for biological interpretation.

## Frontend Agent Integration

I updated the frontend agent prompt so the chat assistant knows the correct pipeline order.

The prompt instructs the agent to:

- run literature search after motif enrichment
- prepare RAG context from `literature_evidence.tsv`
- ground answers in RAG context rows
- cite `source_id`, `citation`, or `url`
- avoid unsupported literature claims

This helps make the LLM behavior safer and more scientifically grounded.

## Testing

I added and extended tests in `tests/stats/test_library.py`.

The test coverage includes:

- CRE candidate ranking
- motif enrichment summary
- synonym-aware PubMed query construction
- multi-intent query generation
- evidence classification
- evidence scoring and sorting
- PubMed API parsing with mocked responses
- JASPAR API parsing with mocked responses
- ENCODE target filtering with mocked responses
- combined literature search across sources
- writing `literature_evidence.tsv`
- RAG context filtering
- citation/source ID generation
- claim extraction
- database fallback context
- RAG output schema

The focused literature/RAG test suite passes:

```bash
python3 -m pytest tests/stats/test_library.py
```

Expected result:

```text
23 passed
```

## Dependencies

I added `requests` as a project dependency because the literature module calls external APIs:

- PubMed / NCBI E-utilities
- JASPAR REST API
- ENCODE REST API

## Why This Matters

The goal of this section is to turn motif enrichment results into interpretable biological evidence.

Without this layer, the pipeline can say that a motif is enriched, but it cannot explain whether that TF has known support in:

- MPRA or reporter assays
- TF binding experiments
- motif databases
- cell-type-specific ENCODE datasets
- perturbation studies

My contribution adds that interpretive bridge. It lets the MCP answer questions like:

- Which TFs might explain the active CREs?
- Is there literature support for this TF in this cell type?
- Is the support direct MPRA evidence or weaker motif/database evidence?
- Which citations support the interpretation?
- What evidence should the LLM use when writing a biological explanation?

## Example Output

A final RAG row can look conceptually like:

```text
source_id: PMID:111
tf: NRF2
evidence_type: MPRA
claim: NFE2L2 NRF2 drives enhancer activity in HepG2 by MPRA.
direction: activation
assay: MPRA
confidence: high
citation: PMID 111
why_relevant:
  - matches TF/synonym
  - matches target cell type
  - mentions regulatory context
  - mentions reporter/MPRA-style assay context
```

This is the kind of structured context the LLM can safely use to write a grounded interpretation.
