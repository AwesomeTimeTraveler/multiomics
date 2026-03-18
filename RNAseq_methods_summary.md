
# RNA-seq aging pipeline: code-focused methods summary

This document summarizes the current analysis logic in a manuscript-style methods format, with emphasis on how objects are filtered, merged, classified, and reused across downstream analyses. It is intentionally written to be easy to expand later with citations and statistical justification.

## Overview

The pipeline integrates three DESeq2-based analysis spaces:

1. **Acomys-only analysis** using a species-specific count matrix and metadata.
2. **Mus-only analysis** using a species-specific count matrix and metadata.
3. **Dual-species ortholog analysis** using an ortholog-collapsed matrix aligned to a mouse-centric ortholog key.

Across these spaces, the workflow produces a shared annotation scaffold, variance-stabilized matrices for QC and downstream systems analyses, differential expression statistics, a unified gene-program classification table, pathway enrichment, TF activity inference, WGCNA modules, and publication-ready plots.

## Input objects and assumptions

The pipeline assumes the following R objects are already available:

- `dds_acomys.rds`
- `dds_mus.rds`
- `dds_dual_ortholog.rds`
- sample metadata embedded in each DESeq2 object
- ortholog mappings linking Acomys and Mus genes
- optional annotation linking ortholog IDs to gene symbols and Entrez IDs
- optional curated gene-set overlays such as secreted genes

The ortholog-centric analyses are anchored to the **Mus Ensembl gene identifier** wherever possible, so all joins and downstream annotations can be performed consistently with mouse annotation resources.

## Ortholog preprocessing and annotation

Ortholog relationships are loaded from the OrthoFinder-derived object and reformatted into explicit Acomys and Mus gene ID columns. Transcript/version suffixes are removed so that identifiers can be matched consistently across DESeq2 results, VST matrices, and annotation resources.

Mouse annotation is then added using `org.Mm.eg.db`, typically including:

- `gene_id`
- `Symbol`
- `description`
- `entrez`

This produces a one-row-per-ortholog scaffold that serves as the master annotation table for downstream merging.

### Processing details

- Acomys and Mus IDs are first parsed from raw ortholog strings.
- Version suffixes are stripped using regex.
- Duplicate ortholog rows are collapsed conservatively to a single representative row.
- Mouse-centric symbol and Entrez annotation are added to the Mus ortholog key.
- The annotation scaffold is saved and reused rather than rebuilt later.

## QC and variance-stabilized expression matrices

Each DESeq2 object is normalized for metadata consistency by explicitly coercing age and species variables into a fixed biological order. Variance-stabilizing transformation is then applied using `vst(..., blind = FALSE)`.

The transformed matrices are used for:

- PCA
- sample-to-sample distance heatmaps
- group mean summaries
- divergent trajectory heatmaps
- stable-gene analyses
- TF activity inference
- WGCNA

### Processing details

- `ageGroup` is coerced to the ordered factor: `embryonic -> young -> old`.
- `species` is coerced to a fixed factor order.
- VST objects are written to disk and reused across all later modules.
- PCA is generated from transformed data rather than raw counts.
- sample distance matrices are computed from Euclidean distances on the transformed expression matrix.

## Differential expression modeling

### Single-species analyses

For each species separately, the pipeline fits:

- a **likelihood ratio test (LRT)** model to detect genes with any age-dependent expression change across the full series
- a **Wald model** for pairwise age contrasts

The primary pairwise contrasts currently emphasized in the classification layer are:

- `young vs embryonic`
- `old vs embryonic`
- `old vs young`

Where relevant, `lfcShrink(..., type = "apeglm")` is used to obtain shrunken log2 fold changes for downstream ranking and interpretation.

### Dual-species analysis

The ortholog matrix is modeled with a species-by-age interaction design:

`~ species_model + ageGroup + species_model:ageGroup`

This is used to generate:

- an **interaction LRT** testing whether trajectories differ between species
- Wald coefficients or linear combinations for species effects at specific ages
- Wald contrasts representing age effects within each species

### Processing details

- The interaction LRT is treated as the primary test of divergence between species.
- Single-species shrunken `old vs embryonic` results are used for the master classification table.
- Results are exported as plain tables and merged back onto the ortholog scaffold by parsed gene ID.

## Construction of the master classification table

A central `master` table is created by joining the ortholog scaffold to:

- Acomys shrunken age-contrast results
- Mus shrunken age-contrast results
- dual-species interaction results
- VST-based age-group mean summaries for each species

This table is the central organizing object for nearly all downstream analyses.

### VST-derived summaries used in classification

For each species separately, the transformed expression matrix is collapsed to age-group means:

- embryonic mean
- young mean
- old mean

The variance across these three means is then computed per gene and used as a stability metric.

### Pattern calling

Each gene is assigned a coarse temporal pattern independently within each species using the three age-group means:

- `stable`
- `monotonic_up`
- `monotonic_down`
- `peak`
- `trough`
- `mixed`

### Significance and stability flags

Per-gene logical flags are then defined, including:

- significant age responsiveness in Acomys
- significant age responsiveness in Mus
- significant species-by-age interaction
- stable expression in Acomys
- stable expression in Mus

Stability is defined by combining:

1. absence of significant age-response under the chosen adjusted P-value cutoff, and
2. low variance across age-group means relative to a lower-tail variance quantile threshold.

### Main program classes

The pipeline then assigns each ortholog to a major class such as:

- `divergent_trajectory`
- `shared_age_responsive`
- `acomys_specific_age_responsive`
- `mus_specific_age_responsive`
- `shared_stable`
- `acomys_uniquely_stable`
- `mus_uniquely_stable_reference`
- `other`

A secondary subclass records matched or mismatched temporal patterns between species and also highlights special subclasses such as embryonic-like retention in Acomys.

## Divergent trajectory analysis

Genes classified as `divergent_trajectory` are further summarized using dual-species group means across six ordered sample groups:

- Acomys_embryonic
- Acomys_young
- Acomys_old
- Mus_embryonic
- Mus_young
- Mus_old

These values are row-scaled and hierarchically clustered for heatmap visualization. Cluster assignments can then be summarized by average trajectory shape across species and age.

### Processing details

- Only genes already classified as divergent are retained for this module.
- Group means are computed directly from the dual-species VST matrix.
- Genes are joined to the master table by the Mus ortholog key.
- Heatmaps are generated from scaled means rather than raw sample-level values.

## Stable and embryonic-retention analysis

Genes classified as stable in Acomys, or specifically as embryonic-like retained programs, are re-examined using transformed expression profiles.

This module evaluates whether adult Acomys expression remains closer to embryonic Acomys than the corresponding Mus profile does.

### Processing details

- The embryonic-retention score is derived from VST group means.
- Stable-gene heatmaps are generated from age-group means.
- Optional sample-level correlation to an embryonic centroid can be used as a higher-level retention score.

## Functional enrichment

The enrichment layer uses the master classification table joined to Entrez annotation and performs:

- GO Biological Process over-representation analysis
- Reactome pathway enrichment
- Hallmark `fgsea` on ranked interaction statistics

### Processing details

- ORA is applied to discrete classes or subclasses.
- `fgsea` uses the interaction statistic as a continuous ranking variable.
- Hallmark gene sets are retrieved with `msigdbr` for mouse.
- Results are saved both as tables and dot plots/bar plots.

## TF analysis with decoupleR

The current decoupleR workflow infers **TF regulatory activity from target-gene expression**, not from TF gene expression alone and not from direct binding occupancy.

The dual-species VST matrix is converted to a Mus-keyed expression matrix, duplicated Mus IDs are collapsed by mean, and a TF–target prior network such as CollecTRI is used with `run_ulm()`.

### What is currently being inferred

This analysis estimates whether the expression of known target genes is consistent with activation or repression of a transcription factor. It is therefore a **regulon activity** analysis.

### What it is not doing

It is not, by itself, testing:

- whether a TF transcript is differentially expressed
- whether a TF physically binds DNA in the analyzed condition

### Recommended separation of TF analyses

The codebase should treat the following as distinct layers:

1. **TF transcript abundance**  
   Differential expression of TF genes themselves.

2. **TF regulatory activity**  
   decoupleR/ULM or similar inference from target-gene behavior.

3. **TF binding potential or occupancy**  
   Separate analysis based on motif databases, promoter/enhancer target predictions, ATAC-seq motif deviation, ChIP-seq, CUT&RUN, or other binding-derived priors.

### Processing details

- Rows are converted to Mus gene IDs for compatibility with mouse regulatory priors.
- Duplicate Mus IDs are collapsed by mean.
- `run_ulm()` is run with a signed TF–target network.
- Downstream plotting should use the inferred numeric activity column, not the method label column.

## WGCNA

WGCNA is used as a systems-level complementary analysis on the dual-species transformed expression matrix.

### Processing details

- Expression is converted to a sample-by-gene matrix.
- Genes can be variance-filtered before network construction.
- A soft threshold is chosen using scale-free topology diagnostics.
- Signed adjacency and topological overlap are computed.
- Dynamic tree cutting identifies modules.
- Closely related modules are merged.
- Module eigengenes are correlated with encoded age, species, and interaction traits.
- Module membership, eigengene values, and hub genes are exported.

## Secretome overlay

A secreted-gene layer should be added as a **parallel annotation axis** on top of the master classification table rather than replacing the existing classes.

Recommended columns:

- `is_secreted`
- `secretome_source`
- `class_secretome`

This allows each gene to retain its original program class while also being labeled as part of the secretome.

### Processing details

- A curated secreted-gene reference is loaded and normalized to the same Mus ortholog namespace as `master$gene_id`.
- A boolean `is_secreted` flag is added by join.
- A secondary `class_secretome` label combines secretome status with `class_main`.
- Secreted-only subsets can then be used for heatmaps, enrichment, candidate ranking, and figure panels.

## Redundancies that should be removed in code

The current QMD contains multiple repeated blocks that should be consolidated:

- repeated helper definitions for age and species normalization
- repeated ID parsing helpers
- repeated enrichment sections
- repeated decoupleR sections
- inconsistent use of `paths$annot` versus `paths$annotation`

These should be centralized into a single config/helper section or sourced helper script.

## Immediate code priorities

1. Centralize all helper functions and thresholds in one config block.
2. Standardize path names and annotation loading.
3. Fix Entrez mapping so `ENTREZID` is actually retrieved.
4. Fix decoupleR downstream code to use the numeric activity score.
5. Add secreted-gene overlay columns directly into the master table.
6. Filter and summarize downstream analyses by `is_secreted` where relevant.
7. Keep one canonical version of each module rather than maintaining duplicated V2/V3 blocks.

## Recommended outputs to maintain

The cleaned pipeline should continue to write:

- annotation scaffold
- VST objects
- DESeq2 result tables
- master classification table
- divergent and stable gene plots
- enrichment tables and figures
- TF activity tables and figures
- WGCNA module tables and figures
- secretome overlay tables and figures

## Suggested manuscript wording for the TF section

A concise statement that matches the code would be:

> Transcription factor activity was inferred from the expression behavior of known TF target genes using a signed TF–target regulatory prior network and univariate linear modeling. This analysis estimates regulator activity from regulon expression patterns and is distinct from both transcription factor mRNA abundance and direct measurements of DNA binding occupancy.

## Suggested manuscript wording for the secretome layer

> To prioritize extracellular signaling programs, we overlaid a curated secreted-gene reference onto the ortholog-level master classification table. Secretome membership was treated as an independent annotation axis, allowing genes to be stratified simultaneously by temporal/species expression program and predicted secreted status.

