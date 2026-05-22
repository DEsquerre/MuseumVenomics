# Data and code for: *Museum venomics reveals an untapped biochemical archive in natural history collections*

This repository accompanies the manuscript:

> Esquerré, D., Keogh, J. S., Dashevsky, D., Boileau, J., Carroll, A.,
> Dunstan, N., & Mikheyev, A. S. *Museum venomics reveals an untapped
> biochemical archive in natural history collections.*

Corresponding author: Damien Esquerré (desquerre@uow.edu.au), Environmental
Futures Research Centre, School of Science, University of Wollongong, NSW
2522, Australia.

---

## Overview

We applied quantitative proteomic mass spectrometry (LC–MS/MS) to preserved
venom glands from 64 specimens (37 venomous snake species, 32 elapids and 5
viperids) housed in Australian natural history collections, alongside fresh
milked venom from 25 of the same species. Specimens spanned 0–57 years of
formalin fixation followed by ethanol storage. We compared preserved-gland
profiles to fresh venom profiles, tested whether specimen age degrades data
quality, and integrated the results with a curated database of published
snake venom proteomes to place museum venomics in the broader context of
elapid and viperid venom diversity.

This repository provides the processed data needed to reproduce all
non-instrumental analyses and figures in the paper, plus the R script that
generates them.

---

## Files

| File | Description |
| --- | --- |
| `MassSpecData.csv` | Per-protein iBAQ values from MaxQuant for every sample (raw output, joined with sample metadata and a manual toxin-family classification). |
| `SummaryResults_VenomMassSpecExperiment.csv` | Per-sample percentage of each major toxin family. Derived from `MassSpecData.csv` by summing iBAQ values per family and re-expressing as a proportion of the total venom proteome of that sample. This is the main input to the analyses in the paper. |
| `PublishedProteomes.csv` | Curated database of published venom proteomes from the literature, harmonised into the same toxin-family categorisation as the experimental data. Used for the cross-study PCA and comparison analyses. |
| `AnalyseCleanData_publication_tidy.R` | R script that reproduces every figure and statistical analysis in the paper from the three CSVs above. |
| `squamates_Title_Science2024_ultrametric_constrained.tre` | Pruned squamate phylogeny (Newick, ultrametric) used for the phylogeny-aligned dot plot. Sourced from Title et al. (2024) *Science*. **Required by the script** but not described as part of the proteomic dataset. |

---

## File details

### `MassSpecData.csv`

iBAQ values for every protein identified by MaxQuant in any sample. The file
is wide and stores both protein-level annotation and per-sample abundance:

* **Columns:** 1 row-label column followed by 7,648 protein columns (one per
  protein group), each headed by the FASTA description string used by
  MaxQuant.
* **Rows 1–7 (annotation rows, identified by the value in column 1):**
  `Protein`, `Family`, `Group`, `Subgroup`, `Function`, `Notes`,
  `UniProt_ID`. These rows describe each protein column: the toxin family
  it was assigned to (used in the summary table), its sub-classification
  where applicable, putative function, and primary UniProt accession.
* **Rows 8 onwards:** one row per sample. Cell values are MaxQuant iBAQ
  intensities for that protein in that sample. Empty cells indicate the
  protein was not detected. Row labels follow the pattern
  `iBAQ <DEVC_id>_<initials>_<form>` where `form` is one of `Venom`
  (milked venom), `Gland` (preserved gland), or other variants for the
  same-individual *Pseudechis australis* validation pair.

This file is provided so others can re-derive family-level summaries with
alternative groupings or run their own protein-level analyses; the figures
in the paper are produced from the summary tables below.

### `SummaryResults_VenomMassSpecExperiment.csv`

One row per sample (n = 93 after sample-quality filtering). The first 11
columns are metadata; columns 12 onwards are the percentage of each major
toxin family in that sample's venom proteome.

Metadata columns:

| Column | Description |
| --- | --- |
| `Sample` | Sample label as it appears in the MaxQuant output (`iBAQ ...`). |
| `D.E.Venom.Collection.Number` | Internal collection identifier (`DEVC_NNN`); cross-references Table S1 of the manuscript. |
| `Family`, `Genus`, `Species` | Taxonomy. |
| `Species2` | `Genus species` concatenation, used as the species key throughout the analysis. |
| `Form_short` | Sample type: `MilkedVenom`, `PreservedFixedGland` (formalin-fixed, ethanol-stored), or `FrozenEthanolGland` (frozen and/or ethanol-only, no formalin). |
| `Sex` | `M`, `F`, or empty when unknown. |
| `Date.collected` | Year of specimen collection. Empty for milked venom from Venom Supplies stock and for captive-bred animals without a recorded date. |
| `SVL`, `Tail.length` | Snout–vent length and tail length when recorded; empty otherwise. |

Toxin-family columns (values are percentages, summing to ~100 per row):
`3FTx`, `PLA2`, `SVMP`, `SVSP`, `CTL`, `CRiSP`, `Kun`, `LAAO`, `5N`,
`AChE`, `AmPep`, `Cys`, `Dis`, `Hyal`, `NGF`, `NP`, `Oha-Vesp`, `PDE`,
`PLB`, `PLC`, `VEGF`, `Venom factor`, `Venom peroxiredoxin`, `Waprin`,
`Other`. See Table S2 of the manuscript and `MassSpecData.csv` row 2
(`Family`) for full descriptions and protein-level membership of each
family.

`Other` aggregates rare toxins not represented elsewhere. The
analysis script re-normalises the percentages to 100 after dropping
`Other` and the published-proteome `Unidentified` column, so all
comparisons in the paper are over the 24 shared toxin families.

Missing values: empty cells indicate either "not measured" (metadata) or
"not detected" (toxin family proportions; treated as zero for the
analyses).

### `PublishedProteomes.csv`

Curated database of published venom proteomes used for the cross-study
PCA (Fig. 3) and the within/within-genus/outside-genus correlation
controls. One row per published proteome.

The first 11 columns mirror the metadata columns of
`SummaryResults_VenomMassSpecExperiment.csv` (with an additional
`Locality` column). The remaining columns are toxin-family percentages,
including the same 24 families used in the experimental data plus a few
additional categories that appear in some published studies (`DEF`,
`Mpi`, `Waglerin`, `Endopeptidase`, `Trypsinogen`, `Protease inhibitor`,
`Peptides`, `Unidentified`). The `TOTAL` column gives the sum the entry
was originally reported with; `Reference` cites the source publication.

The starting points for the database were Tasoulis & Isbister (2017)
*Toxins* 9:290 and Oliveira et al. (2022) *Frontiers in Ecology and
Evolution* 10:1066144, supplemented with subsequent literature.

### `AnalyseCleanData_publication_tidy.R`

R script that reproduces every figure and table-style result reported in
the paper. The analyses are organised into 21 numbered sections:

1. Package check and load.
2. Path resolution (works under `Rscript`, RStudio, and interactive
   sessions; outputs go to `./outputs/`).
3. Load and prepare data; drop `Other`/`Unidentified` and re-normalise to
   100% across the 24 shared toxin families.
4. Build a combined experimental + published data frame and per-species
   averages; produce the toxin matrix used for PCA.
5. Define the toxin-family colour palette used throughout.
6. Per-sample pie charts, full and filtered to families ≥1%.
7. Per-species-average pie charts.
8. Multi-pie PDF pages grouped by family.
9. Toxin-legend grid (visual key).
10. Per-species gland-vs-venom regression plots.
11. Within-species pairwise correlation table
    (`CorrelationResultsVenomVsGlands.csv`) using linear, Spearman, and
    beta regression. This table underlies Fig. 1A.
12. Within-species pairwise correlations (alternative formatting used
    for the comparison-category histograms).
13. Within-genus, between-species comparisons.
14. Outside-genus comparisons.
15. Histogram and density plots of correlation metrics by category
    (Fig. 1A) plus a CSV of category means.
16. Effect of specimen age on gland–venom correlation (Fig. 2A); fits
    LMs and a species-random-intercept mixed model.
17. Phylogeny import and imputation of missing tips (used in the
    phylogeny-aligned dot plot).
18. Per-species LM R² dot plot aligned with phylogeny (Fig. 2B).
19. Same-individual *Pseudechis australis* validation (Fig. 1B).
20. PCA of mean toxin profiles by species (Fig. 3) plus an interactive
    Plotly version saved as HTML.
21. `sessionInfo.txt` written to `outputs/`.

#### Software environment

The script was developed under R ≥ 4.4 and uses the following CRAN
packages. They are checked at startup; if any are missing the script
prints the full `install.packages()` line:

`dplyr`, `tidyr`, `reshape2`, `ggplot2`, `ggrepel`, `ggfortify`,
`patchwork`, `cowplot`, `plotly`, `RColorBrewer`, `grDevices`, `pals`,
`colorspace`, `scico`, `grid`, `betareg`, `ape`, `phytools`, `lme4`,
`MuMIn`. `htmlwidgets` and `rmarkdown` are optional (used only to save
the interactive PCA as a self-contained HTML).

`sessionInfo.txt` is written to the outputs folder on every run as a
provenance record.

#### How to run

Place all four input files (the three CSVs and the `.tre` file) and the
script in the same directory, then either:

* run `Rscript AnalyseCleanData_publication_tidy.R` from that directory, or
* open the script in RStudio and source it.

The script auto-detects its own directory (via `here::here()` →
`rstudioapi` → `Rscript --file=` → `getwd()`); no `setwd()` or path
editing is needed. All outputs are written to `./outputs/`, including:

* `outputs/PieCharts/` and `PieCharts_filteredOver1percent/` — per-sample
  pies (full and ≥1% filtered).
* `outputs/PieCharts_AveragePerSpecies/` (+ filtered) — per-species
  averages.
* `outputs/PieCharts_4x6_byFamily/` — multi-pie PDF pages by family.
* `outputs/RegressionPlots/`, `RegressionPlots_filtered/`,
  `RegressionPlots_combined/` — per-species gland-vs-venom regressions.
* `outputs/Tables/` — `CorrelationResultsVenomVsGlands.csv`,
  `results_within_genus_comparisons.csv`,
  `results_outside_genus_comparisons.csv`,
  `ComparisonCategory_means.csv`, `GlandAge_models.txt`.
* Top-level PDFs: `ToxinLegendGrid.pdf`,
  `WithinSpecies_metric_distributions.pdf`,
  `ComparisonCategories_histograms.pdf`,
  `ComparisonCategories_densities.pdf`, `GlandAge_vs_metrics.pdf`,
  `PhylogeneticTree_imputed.pdf`, `DotPlot_LMR2_by_species.pdf`,
  `PCA_meanToxinProfiles.pdf`, `PCA_interactive.html`.
* `sessionInfo.txt`.

The script catches numerical edge cases internally (sparse or perfectly
collinear pairs that destabilise beta regression are reported via
`message()` and skipped); these are warnings, not errors.

---

## Methods (brief)

Sample preparation followed a modified filter-aided sample preparation
(FASP) protocol for formalin-fixed tissues. Tryptic peptides were
analysed by data-dependent LC–MS/MS on a Thermo Orbitrap Fusion ETD
coupled to a Dionex UltiMate 3000 RSLCnano. Raw files were searched in
MaxQuant v2.0.3.0 against UniProt restricted to taxon = Serpentes, with
a 1% FDR at peptide-spectrum, protein, and site levels. Quantification
used intensity-based absolute quantification (iBAQ). Each protein was
manually assigned to one of 24 toxin families (Table S2 of the paper);
relative abundance of each toxin family in each sample was computed as
the sum of iBAQ values for that family divided by the sum of iBAQ values
across all proteins identified in the sample. Nine gland samples with
<15% venom proteins were removed prior to downstream analyses; in the
retained samples all non-venom proteins were excluded before computing
proportions. Full methodological detail is in the manuscript and its
Supplementary Materials.

Raw mass spectrometry data are deposited in the PRIDE proteomics
repository (accession listed in the published manuscript).

---

## License

Data and code are released under CC0 1.0 (public domain dedication) for
maximum re-use. If you use these data or code in your work, please cite
the manuscript above.

---

## Contact

For questions about the data, methods, or analysis pipeline please
contact Damien Esquerré (desquerre@uow.edu.au).
