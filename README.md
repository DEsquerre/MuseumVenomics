# Data and code for: *Unlocking the venom vault: Museum venomics reveals an untapped biochemical archive in natural history collections*

This repository accompanies the manuscript:

> Esquerré, D., Keogh, J. S., Dashevsky, D., Boileau, J., Carroll, A.,
> Dunstan, N., Balchan, N. R., & Mikheyev, A. S. *Unlocking the venom vault:
> Museum venomics reveals an untapped biochemical archive in natural history
> collections.*

Corresponding author: Damien Esquerré (desquerre@uow.edu.au), Environmental
Futures Research Centre, School of Science, University of Wollongong, NSW
2522, Australia.

Archived version: https://doi.org/10.5281/zenodo.XXXXXXX

---

## Overview

We used quantitative proteomic mass spectrometry (LC–MS/MS) on venom glands
from 61 preserved snakes in Australian natural history collections, stored for
up to 57 years, and compared them with lyophilised milked venom from 25 of the
same species. We also compared a formalin-fixed gland with milked venom from
the same individual for five species (*Pseudechis australis*, *Agkistrodon
contortrix*, *A. piscivorus*, *Crotalus atrox* and *Sistrurus miliarius*), and
placed the results in the context of a curated database of published snake
venom proteomes.

This repository holds the processed data and the R script that reproduces the
statistical analyses and figures in the paper. Raw mass spectrometry data are
deposited in PRIDE/ProteomeXchange under accession **PXDXXXXXX**.

---

## Files

| File | Description |
| --- | --- |
| `MuseumVenomics_analysis.R` | R script that reproduces the analyses and figures (see below). |
| `SummaryResults_ThisStudy.csv` | Toxin-family composition of every sample from this study used in the paper (92 samples: museum glands, milked venoms and the five same-individual gland–venom pairs). |
| `PublishedProteomes.csv` | The 413 published snake venom proteomes used in the paper, harmonised to the same toxin families. |
| `squamates_Title_Science2024_ultrametric_constrained.tre` | Ultrametric squamate phylogeny (Newick) from Title et al. (2024) *Science*, used to order species in Fig. 3B. |
| `MassSpecData.csv` | Per-protein MaxQuant iBAQ values with the manual toxin-family classification. Not read by the script; provided for protein-level re-analysis. |
| `LICENSE` | CC0 1.0. |

---

## Data files

### `SummaryResults_ThisStudy.csv`

One row per sample: 11 metadata columns followed by the percentage of each
toxin family in the sample's venom proteome.

| Column | Description |
| --- | --- |
| `Sample` | Sample label as in the MaxQuant output (`iBAQ DEVC_NNN...`). |
| `D.E.Venom.Collection.Number` | Collection identifier (`DEVC NNN`); cross-references the supplementary sample tables. |
| `Family`, `Genus`, `Species` | Taxonomy (`Species` is the epithet). |
| `Species2` | `Genus species`, used as the species key. |
| `Form_short` | `MilkedVenom`, `PreservedFixedGland` (formalin-fixed, ethanol-stored) or `FrozenEthanolGland` (frozen and/or ethanol only). |
| `Sex` | Where recorded. |
| `Date.collected` | Year the specimen was collected (`NA` for most milked venoms). |
| `SVL`, `Tail.length` | Snout–vent and tail length (mm) where recorded. |

Toxin-family columns: `3FTx`, `PLA2`, `SVMP`, `SVSP`, `CTL`, `CRiSP`, `Kun`,
`LAAO`, `5N`, `AChE`, `AmPep`, `Cys`, `Dis`, `Hyal`, `NGF`, `NP`, `Oha-Vesp`,
`PDE`, `PLB`, `PLC`, `VEGF`, `Venom factor`, `Venom peroxiredoxin`, `Waprin`,
`Other`. Each is the summed iBAQ of the family's proteins divided by the summed
iBAQ of all venom proteins in the sample (×100). A value of 0 means the family
was not detected.

Nine samples that failed quality control are not included: DEVC_015
(labelled *Crotalus unicolor* but with a cobra-like profile, a labelling
error), DEVC_016 (*Agkistrodon bilineatus*), DEVC_048 and DEVC_054
(*Pseudonaja*), DEVC_087 (*Denisonia devisi*), DEVC_095 (*Drysdalia
coronoides*) and DEVC_096 (*Echiopsis curta*), all with low venom-protein
content and/or highly divergent profiles; and DEVC_020_g and DEVC_021_g, gland
tissue from specimens for which venom extracted from the same gland was
analysed instead (DEVC_020_v, DEVC_021_v).

### `PublishedProteomes.csv`

One row per published proteome. Metadata columns as above, plus `Locality`,
`Technique`, `Reference` and `DOI`. Toxin families are reported as percentages
of the total venom proteome; the 24 families used in the analyses already
include any subclasses that some studies report separately (e.g. `L-3FTx`,
`SV-GII-PLA2`), which are retained for reference along with extra categories
(`Other`, `Unidentified`, `Peptides`, etc.). `TOTAL` is the sum as originally
reported. Empty cells mean the family was not reported. The database was built
from Tasoulis & Isbister (2017) *Toxins* 9:290 and Oliveira et al. (2022)
*Nature Reviews Chemistry*, plus later literature. One proteome (Vanuopadath
et al. 2022) is used in the pie charts but excluded from the PCA in Fig. 4.

### Processing applied by the script

`Other` and `Unidentified` are dropped and the 24 shared toxin families are
re-normalised to 100% in every sample, so all comparisons use the same
families.

---

## Running the script

Requires R (developed under R 4.4.2) and the CRAN packages `dplyr`, `tidyr`,
`ggplot2`, `ggrepel`, `patchwork`, `betareg`, `ape`, `phytools`, `scatterpie`
and `ggforce`. The script checks for them at start-up and prints the
`install.packages()` command for any that are missing.

Put the script and the three input files (two CSVs and the `.tre` file) in
the same folder, then either run

```
Rscript MuseumVenomics_analysis.R
```

from that folder, or open the script in RStudio and source it. No paths need
editing. All outputs are written to `./outputs/`, and `outputs/sessionInfo.txt`
records the package versions used.

Beta regressions that fail to converge for an individual comparison (sparse
profiles) are reported with a message and recorded as `NA`. These are not
errors.

---

## Outputs and where they appear in the paper

| Paper | Output (in `outputs/`) |
| --- | --- |
| Fig. 1 | `Fig1_ComparisonCategories.pdf` |
| Results: mean metrics per comparison category | `Tables/ComparisonCategory_means.csv` (all pairs: `Tables/Comparisons_*.csv`) |
| Fig. 2 | `Fig2_SameSpecimen_pies.pdf`, `Fig2_SameSpecimen_regressions.pdf` |
| Results: same-specimen statistics | `Tables/SameSpecimen_correlations.csv` |
| Fig. 3A and age models | `Fig3A_SpecimenAge.pdf`, `Tables/GlandAge_models.txt` |
| Fig. 3B | `Fig3B_Species.pdf`, `Fig3B_Phylogeny.pdf` |
| Fig. 4 | `Fig4_PCA_log.pdf`; loadings and variance in `Tables/Fig4_PCA_log_*.csv` (raw-proportion PCA: `Fig4_PCA_raw.pdf`) |
| Fig. S2 | `FigS2_ToxinFamilyRichness.pdf`, `Tables/ToxinFamilyRichness.csv` |
| Fig. S3 | `FigS3_IntraspecificPCA.pdf` |
| Appendix 1 | `Appendix1_PieCharts/`, `Appendix1_Regressions/` |

Fig. S1 (venom-protein content of each sample) was produced from the
protein-level data and is not generated by this script.

---

## Methods (brief)

Gland tissue was processed with a filter-aided sample preparation (FASP)
protocol for formalin-fixed tissue, digested with Trypsin/LysC and analysed by
data-dependent LC–MS/MS on a Thermo Orbitrap Fusion ETD coupled to a Dionex
UltiMate 3000 RSLCnano. Raw files were searched in MaxQuant v2.0.3.0 against
UniProt (taxon = Serpentes) at 1% FDR, and proteins were quantified by iBAQ.
Each venom protein was assigned manually to one of 24 toxin families. Gland and
venom profiles were compared with linear regression, Spearman rank correlation
and beta regression. Full details are in the manuscript and its supplementary
material.

---

## License

Data and code are released under CC0 1.0 (public domain dedication). If you
use them, please cite the paper above.

## Contact

Damien Esquerré (desquerre@uow.edu.au)
