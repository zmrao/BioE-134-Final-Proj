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

- `tool_rank_cre_candidates`: ranks CREs by activity and statistical confidence.
- `tool_motif_enrichment_summary`: summarizes TF motifs enriched among active CREs.
- `tool_prepare_rag_context`: prepares top CREs and motif search terms.
- `tool_search_pubmed`: searches PubMed with NCBI E-utilities.
- `tool_search_jaspar_motif`: searches JASPAR motif matrix profiles.
- `tool_search_encode_tf`: searches ENCODE TF/cell-type functional genomics records.
- `tool_literature_search_for_motifs`: runs PubMed, JASPAR, and ENCODE searches for top motifs and writes `literature_evidence.tsv`.
- `tool_interpret_literature_evidence`: summarizes retrieved evidence by motif and source.
- `tool_prepare_literature_rag_context`: converts scored evidence into citation-ready RAG rows and writes `literature_rag_context.tsv`.

These tools expose the literature and interpretation layer through the MCP server so the frontend agent can call them automatically.

Each MCP wrapper includes typed inputs, default parameters, a docstring description, and JSON-serializable output. The main tool inputs include paths to motif/evidence tables, target/off-target cell type, species, maximum result counts, output paths, and RAG filtering parameters such as `min_score`, `max_records`, and `max_per_tf`.

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

## Test Prompts

These are realistic prompts that trigger my literature/RAG functionality through the MCP-backed chat agent:

```text
Run literature search for the top enriched motifs in my motif_enrichment.tsv for HepG2 and prepare citation-ready RAG context.
```

Expected behavior:

- Calls `tool_literature_search_for_motifs`
- Searches PubMed, JASPAR, and ENCODE
- Writes `literature_evidence.tsv`
- Calls `tool_prepare_literature_rag_context`
- Returns citation-ready rows with source IDs, claims, evidence types, and confidence

```text
For the enriched NRF2 and HNF4A motifs, find PubMed and ENCODE evidence in HepG2 and summarize which TF has stronger support.
```

Expected behavior:

- Uses synonym expansion, such as `NRF2` and `NFE2L2`
- Searches for MPRA/reporter, binding, motif, and perturbation evidence
- Scores records by TF match, cell-type match, assay relevance, regulatory keywords, and recency
- Summarizes strongest evidence by motif

```text
Use the literature evidence table to answer: are the top motifs supported by direct MPRA evidence or only by motif/database evidence?
```

Expected behavior:

- Calls `tool_prepare_literature_rag_context` or `tool_interpret_literature_evidence`
- Distinguishes `MPRA`, `Reporter_assay`, `ChIP_binding`, `Motif_analysis`, and `Database` evidence
- Grounds the answer in `source_id`, `citation`, or `url`

```text
Prepare RAG context from literature_evidence.tsv with at most 5 records per TF and only include records with evidence score above 4.
```

Expected behavior:

- Calls `tool_prepare_literature_rag_context`
- Applies score filtering
- Balances evidence types
- Caps records per TF
- Writes `literature_rag_context.tsv`

## Rubric Coverage

### Function Code

The code solves the defined problem: converting enriched CRE-seq motifs into biologically meaningful, citation-backed evidence for interpretation. It uses structured query generation, external API retrieval, evidence classification, scoring, claim extraction, and balanced RAG formatting. The implementation stays within the existing MCP architecture and does not add unnecessary dependencies.

### Documentation

The functions include docstrings explaining purpose, inputs, outputs, and behavior. This README explains the scope of my section and my individual contribution.

### MCP Wrappers

The MCP wrappers are registered in `creseq_mcp/server.py` with typed inputs, defaults, docstring descriptions, and JSON-safe outputs. The literature tools expose search, scoring, interpretation, and RAG formatting through the same structure as the rest of the MCP.

### Test Prompts

The prompt examples above are aligned with the literature/RAG scope and trigger realistic tool usage through the frontend chat agent.

### Pytests

The tests in `tests/stats/test_library.py` cover typical and edge cases with mocked APIs. They verify query generation, evidence classification, scoring, PubMed/JASPAR/ENCODE parsing, combined search behavior, file output, empty-result handling, RAG context formatting, citation IDs, and output schema.

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
