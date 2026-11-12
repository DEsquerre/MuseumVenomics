# =============================================================================
# Unlocking the venom vault: Museum venomics reveals an untapped biochemical
# archive in natural history collections
#
# Esquerré D, Keogh JS, Dashevsky D, Boileau J, Carroll A, Dunstan N,
# Balchan NR, Mikheyev AS
#
# Reproduces the statistical analyses and figures in the paper (except
# Fig. S1) from the processed toxin-family tables.
#
# Inputs (same folder as this script):
#   SummaryResults_ThisStudy.csv  - toxin-family % for every sample analysed
#   PublishedProteomes.csv        - published venom proteomes used in the paper
#   squamates_Title_Science2024_ultrametric_constrained.tre - Title et al. 2024
#
# Outputs are written to ./outputs/. See README.md for which output file
# corresponds to which figure or result.
#
# Developed under R 4.4.2.
# =============================================================================


# ---- 0. Settings ------------------------------------------------------------

# Folder holding the inputs: detected automatically when run with Rscript or
# sourced in RStudio; otherwise the current working directory.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", args[grep("^--file=", args)])
  if (length(f) == 1) return(dirname(normalizePath(f)))
  p <- tryCatch(
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable())
      rstudioapi::getSourceEditorContext()$path else "",
    error = function(e) "")
  if (nzchar(p)) return(dirname(p))
  getwd()
}
base_dir    <- get_script_dir()
outputs_dir <- file.path(base_dir, "outputs")

# Gland and milked venom from the same individual (Fig. 2)
same_specimen_pairs <- data.frame(
  Species = c("Pseudechis australis", "Agkistrodon contortrix",
              "Agkistrodon piscivorus", "Crotalus atrox", "Sistrurus miliarius"),
  Venom   = c("DEVC_145", "DEVC_331", "DEVC_319", "DEVC_328", "DEVC_310"),
  Gland   = c("DEVC_146", "DEVC_333", "DEVC_321", "DEVC_330", "DEVC_312")
)

# Specimen age = year the museum glands were dissected minus collection year
dissection_year <- 2022

# The four viper same-specimen glands were fixed in 2024, so their collection
# year does not give a storage age comparable with the museum specimens; they
# are left out of the age analysis (Fig. 3).
age_excluded_glands <- c("DEVC_333", "DEVC_321", "DEVC_330", "DEVC_312")

# Published proteome excluded from the PCA (Fig. 4)
excluded_published_pca <- c("Vanuopadath et al. 2022")

# Beta regression requires values strictly inside (0, 1): proportions of 0 or
# 1 are replaced by eps or 1 - eps. The comparison-category analysis (Fig. 1)
# used eps = 1e-11; all other analyses used eps = 1e-9.
eps_categories <- 1e-11
eps_default    <- 1e-9

# Toxin-family richness threshold, % abundance (Fig. S2)
richness_threshold <- 0.1

# Minimum number of proteomes per species for the within-species PCAs (Fig. S3)
min_samples_intraspecific <- 3


# ---- 1. Packages ------------------------------------------------------------

required_pkgs <- c("dplyr", "tidyr", "ggplot2", "ggrepel", "patchwork",
                   "betareg", "ape", "phytools", "scatterpie", "ggforce")
missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nInstall with: install.packages(c(\"",
       paste(missing_pkgs, collapse = "\", \""), "\"))")
}
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(ggrepel)
  library(patchwork); library(betareg); library(ape); library(phytools)
  library(scatterpie); library(ggforce)
})

out <- function(...) file.path(outputs_dir, ...)
for (d in c("Tables", "Appendix1_PieCharts", "Appendix1_Regressions")) {
  dir.create(out(d), showWarnings = FALSE, recursive = TRUE)
}


# ---- 2. Data ----------------------------------------------------------------

input_files <- c("SummaryResults_ThisStudy.csv", "PublishedProteomes.csv",
                 "squamates_Title_Science2024_ultrametric_constrained.tre")
missing_inputs <- input_files[!file.exists(file.path(base_dir, input_files))]
if (length(missing_inputs) > 0) {
  stop("Input file(s) not found in ", base_dir, ": ", paste(missing_inputs, collapse = ", "))
}

# The 24 toxin families compared throughout ("Other" and "Unidentified" are
# excluded). Names follow R's check.names conversion (3FTx -> X3FTx, etc.).
toxins <- c("X3FTx", "PLA2", "SVMP", "SVSP", "CTL", "CRiSP", "Kun", "LAAO",
            "X5N", "AChE", "AmPep", "Cys", "Dis", "Hyal", "NGF", "NP",
            "Oha.Vesp", "PDE", "PLB", "PLC", "VEGF", "Venom.factor",
            "Venom.peroxiredoxin", "Waprin")
toxin_labels <- setNames(sub("^X", "", gsub("\\.", " ", toxins)), toxins)
toxin_labels["Oha.Vesp"] <- "Oha-Vesp"

# Empty cells (not detected) become 0, and the 24 families are re-normalised
# to sum to 100% in each sample.
normalise_toxins <- function(df) {
  df[toxins] <- lapply(df[toxins], function(x) { x <- as.numeric(x); x[is.na(x)] <- 0; x })
  df[toxins] <- df[toxins] / rowSums(df[toxins]) * 100
  df
}

experimental <- normalise_toxins(read.csv(file.path(base_dir, "SummaryResults_ThisStudy.csv")))
published    <- normalise_toxins(read.csv(file.path(base_dir, "PublishedProteomes.csv")))
message("Samples from this study: ", nrow(experimental),
        " | published proteomes: ", nrow(published))

color_mapping <- c(
  "X3FTx" = "#F47B5B", "PLA2" = "#B5EFB5", "SVMP" = "#FBE426", "SVSP" = "#2ED9FF",
  "CTL" = "#AA0DFE", "CRiSP" = "#AAF400", "Kun" = "#3283FE", "LAAO" = "#F7E1A0",
  "X5N" = "#F8A19F", "AChE" = "#F6222E", "AmPep" = "#DD6F91", "Cys" = "#F45366",
  "Dis" = "#1C7F93", "Hyal" = "#1C8356", "NGF" = "#B10DA1", "NP" = "#1CBE4F",
  "Oha.Vesp" = "#7ED7D1", "PDE" = "#C075A6", "PLB" = "#FC1CBF", "PLC" = "#FA0087",
  "VEGF" = "#BDCDFF", "Venom.factor" = "#822E1C", "Venom.peroxiredoxin" = "#C4451C",
  "Waprin" = "#782AB6"
)

safe_name <- function(x) gsub("[/\\\\:*?\"<>|]", "_", x)

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[[which(grepl("guide-box", g$layout$name) &
                   !sapply(g$grobs, inherits, "zeroGrob"))[1]]]
}


# ---- 3. Correlation metrics -------------------------------------------------

to_beta <- function(x, eps) pmin(pmax(x / 100, eps), 1 - eps)

# Compare a milked venom profile with another profile using LM adjusted R2,
# Spearman correlation and beta-regression pseudo-R2. Families absent from
# both samples are dropped.
compare_profiles <- function(venom, other, eps) {
  keep <- !(venom == 0 & other == 0)
  d <- data.frame(MilkedVenom = venom[keep], Gland = other[keep])
  lm_s  <- summary(lm(MilkedVenom ~ Gland, data = d))
  spear <- suppressWarnings(cor.test(d$Gland, d$MilkedVenom, method = "spearman"))
  br <- tryCatch({
    s <- suppressWarnings(summary(betareg(MilkedVenom ~ Gland,
           data = data.frame(MilkedVenom = to_beta(d$MilkedVenom, eps),
                             Gland = to_beta(d$Gland, eps)))))
    c(s$pseudo.r.squared, s$coefficients$mean["Gland", "Pr(>|z|)"])
  }, error = function(e) c(NA, NA))
  data.frame(n_families    = nrow(d),
             LM_R_squared  = round(lm_s$adj.r.squared, 3),
             LM_p          = round(lm_s$coefficients["Gland", "Pr(>|t|)"], 3),
             Spearman_R    = round(unname(spear$estimate), 3),
             Spearman_p    = round(spear$p.value, 3),
             GLM_R_squared = round(br[1], 3),
             GLM_p         = round(br[2], 3))
}

# Compare every milked venom with each sample returned by `candidates(venom)`
run_comparisons <- function(candidates, eps) {
  venoms <- experimental[experimental$Form_short == "MilkedVenom", ]
  rows <- list()
  for (i in seq_len(nrow(venoms))) {
    v <- venoms[i, ]
    pool <- candidates(v)
    for (j in seq_len(nrow(pool))) {
      o <- pool[j, ]
      rows[[length(rows) + 1]] <- cbind(
        data.frame(Species1 = v$Species2, Sample1 = v$Sample,
                   Species2 = o$Species2, Sample2 = o$Sample,
                   Form_Sample2 = o$Form_short,
                   Date_Collected_Sample2 = o$Date.collected),
        compare_profiles(unlist(v[toxins]), unlist(o[toxins]), eps))
    }
  }
  do.call(rbind, rows)
}

same_species  <- function(v) experimental[experimental$Species2 == v$Species2 &
                                             experimental$Sample != v$Sample, ]
same_genus    <- function(v) experimental[experimental$Genus == v$Genus &
                                             experimental$Species2 != v$Species2, ]
other_genera  <- function(v) experimental[experimental$Genus != v$Genus, ]


# ---- 4. Fig. 1: within species vs within genus vs outside genus ----------

within_species <- run_comparisons(same_species, eps_categories)
within_genus   <- run_comparisons(same_genus,   eps_categories)
outside_genus  <- run_comparisons(other_genera, eps_categories)

write.csv(within_species, out("Tables", "Comparisons_within_species.csv"), row.names = FALSE)
write.csv(within_genus,   out("Tables", "Comparisons_within_genus.csv"),   row.names = FALSE)
write.csv(outside_genus,  out("Tables", "Comparisons_outside_genus.csv"),  row.names = FALSE)

categories <- c("Within species", "Within genus", "Outside genus")
all_comp <- bind_rows(mutate(within_species, Category = categories[1]),
                      mutate(within_genus,   Category = categories[2]),
                      mutate(outside_genus,  Category = categories[3])) %>%
  mutate(Category = factor(Category, levels = categories))

# Means per category (reported in Results)
category_means <- all_comp %>%
  group_by(Category) %>%
  summarise(n_comparisons = n(),
            across(c(LM_R_squared, LM_p, Spearman_R, Spearman_p, GLM_R_squared, GLM_p),
                   ~ mean(.x, na.rm = TRUE), .names = "{.col}_mean"),
            .groups = "drop")
write.csv(category_means, out("Tables", "ComparisonCategory_means.csv"), row.names = FALSE)
print(as.data.frame(category_means))

metric_names <- c(LM_R_squared = "LM R²", Spearman_R = "Spearman R",
                  GLM_R_squared = "Beta regression pseudo-R²")
category_colors <- c("Within species" = "#E69F00", "Within genus" = "#009E73",
                     "Outside genus" = "#0072B2")

p_fig1 <- all_comp %>%
  select(Category, all_of(names(metric_names))) %>%
  pivot_longer(-Category, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(metric_names[Metric], levels = metric_names)) %>%
  ggplot(aes(Value, fill = Category, colour = Category)) +
  geom_density(alpha = 0.5, adjust = 1.2, na.rm = TRUE) +
  facet_grid(Metric ~ ., scales = "free_y") +
  scale_fill_manual(values = category_colors) +
  scale_colour_manual(values = category_colors) +
  labs(x = "Metric value", y = "Density", fill = "Comparison", colour = "Comparison") +
  theme_minimal(base_size = 12)
ggsave(out("Fig1_ComparisonCategories.pdf"), p_fig1, width = 8, height = 10)


# ---- 5. Fig. 2: gland vs venom from the same individual -------------------

get_profile <- function(id) {
  r <- experimental[grepl(paste0("^iBAQ ", id, "(_|$)"), experimental$Sample), ]
  if (nrow(r) != 1) stop("Expected one sample for ", id, ", found ", nrow(r))
  unlist(r[toxins])
}

same_specimen_results <- list()
same_specimen_plots   <- list()
for (k in seq_len(nrow(same_specimen_pairs))) {
  sp <- same_specimen_pairs$Species[k]
  d  <- data.frame(Toxin = toxins,
                   Venom = get_profile(same_specimen_pairs$Venom[k]),
                   Gland = get_profile(same_specimen_pairs$Gland[k]))
  d  <- d[!(d$Venom == 0 & d$Gland == 0), ]

  lm_s  <- summary(lm(Venom ~ Gland, data = d))
  spear <- cor.test(d$Venom, d$Gland, method = "spearman", exact = FALSE)
  br    <- summary(betareg(Venom ~ Gland,
                           data = data.frame(Venom = to_beta(d$Venom, eps_default),
                                             Gland = to_beta(d$Gland, eps_default))))

  same_specimen_results[[k]] <- data.frame(
    Species = sp, Venom_sample = same_specimen_pairs$Venom[k],
    Gland_sample = same_specimen_pairs$Gland[k], n_toxin_families = nrow(d),
    LM_R_squared = round(lm_s$adj.r.squared, 3),
    LM_p = signif(lm_s$coefficients["Gland", "Pr(>|t|)"], 3),
    Spearman_R = round(unname(spear$estimate), 3), Spearman_p = signif(spear$p.value, 3),
    GLM_R_squared = round(br$pseudo.r.squared, 3),
    GLM_p = signif(br$coefficients$mean["Gland", "Pr(>|z|)"], 3))

  same_specimen_plots[[sp]] <- ggplot(d, aes(Gland, Venom, colour = Toxin)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = toxin_labels[Toxin]), box.padding = 0.6,
                    max.overlaps = Inf, size = 3) +
    geom_smooth(inherit.aes = FALSE, aes(Gland, Venom), method = "lm", formula = y ~ x) +
    scale_x_sqrt() + scale_y_sqrt() +
    scale_colour_manual(values = color_mapping) +
    annotate("text", x = 0, y = Inf, hjust = 0, vjust = 1.3, size = 3.5, fontface = "italic",
             label = sprintf("LM R² = %.2f\nSpearman R = %.2f\nBeta pseudo-R² = %.2f",
                             lm_s$adj.r.squared, spear$estimate, br$pseudo.r.squared)) +
    labs(title = sp, x = "Gland (%)", y = "Venom (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = "none")
}
same_specimen_results <- do.call(rbind, same_specimen_results)
write.csv(same_specimen_results, out("Tables", "SameSpecimen_correlations.csv"), row.names = FALSE)
print(same_specimen_results)

ids <- unlist(same_specimen_pairs[, c("Venom", "Gland")])
toxins_present <- toxins[rowSums(sapply(ids, get_profile)) > 0]

ss_legend <- extract_legend(
  ggplot(data.frame(Toxin = factor(toxins_present, levels = toxins_present), x = 1, y = 1),
         aes(x, y, colour = Toxin)) +
    geom_point(size = 4) +
    scale_colour_manual(values = color_mapping, labels = toxin_labels, name = NULL) +
    guides(colour = guide_legend(nrow = 9, byrow = TRUE)) +
    theme_void() +
    theme(legend.position = "bottom", legend.text = element_text(size = 10),
          legend.key.size = unit(0.45, "cm")))
ggsave(out("Fig2_SameSpecimen_regressions.pdf"),
       wrap_plots(same_specimen_plots, ncol = 3) / wrap_elements(full = ss_legend) +
         plot_layout(heights = c(10, 2.5)),
       width = 15, height = 11.5)

pie_df <- do.call(rbind, lapply(seq_len(nrow(same_specimen_pairs)), function(k) {
  rbind(data.frame(Species = same_specimen_pairs$Species[k], Sample = "Milked venom",
                   Toxin = toxins_present,
                   Pct = get_profile(same_specimen_pairs$Venom[k])[toxins_present]),
        data.frame(Species = same_specimen_pairs$Species[k], Sample = "Preserved gland",
                   Toxin = toxins_present,
                   Pct = get_profile(same_specimen_pairs$Gland[k])[toxins_present]))
}))
pie_df <- pie_df[pie_df$Pct > 0, ] %>%
  mutate(Toxin   = factor(Toxin, levels = toxins_present),
         Species = factor(Species, levels = same_specimen_pairs$Species),
         Sample  = factor(Sample, levels = c("Milked venom", "Preserved gland")))

p_ss_pies <- ggplot(pie_df, aes(x = 1, y = Pct, fill = Toxin)) +
  geom_col(width = 1, colour = NA) +
  coord_polar(theta = "y") +
  facet_grid(Sample ~ Species, switch = "y") +
  scale_fill_manual(values = color_mapping, labels = toxin_labels, name = NULL, drop = FALSE) +
  guides(fill = guide_legend(nrow = 3, byrow = TRUE)) +
  theme_void() +
  theme(strip.text.x = element_text(face = "italic", size = 11, margin = margin(b = 4)),
        strip.text.y.left = element_text(size = 11, angle = 90, margin = margin(r = 4)),
        legend.position = "bottom", legend.text = element_text(size = 9),
        legend.key.size = unit(0.4, "cm"))
ggsave(out("Fig2_SameSpecimen_pies.pdf"), p_ss_pies, width = 12, height = 6)


# ---- 6. Fig. 3A: specimen age vs gland-venom correlation ----------------

age_excluded <- grepl(paste0("^iBAQ (", paste(age_excluded_glands, collapse = "|"), ")(_|$)"),
                      experimental$Sample)
within_species_age <- run_comparisons(same_species, eps_default)
results_dated <- within_species_age %>%
  filter(!is.na(Date_Collected_Sample2),
         !Sample2 %in% experimental$Sample[age_excluded]) %>%
  mutate(age = dissection_year - Date_Collected_Sample2)
message("Comparisons in the age analysis: ", nrow(results_dated))

sink(out("Tables", "GlandAge_models.txt"))
for (m in names(metric_names)) {
  cat("\n=====", metric_names[[m]], "~ specimen age =====\n")
  print(summary(lm(reformulate("age", m), data = results_dated)))
}
sink()

metric_colors <- setNames(c("#F43151", "#8362DD", "#BAEA75"), metric_names)
p_fig3a <- results_dated %>%
  select(age, all_of(names(metric_names))) %>%
  pivot_longer(-age, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(metric_names[Metric], levels = metric_names)) %>%
  ggplot(aes(age, Value, colour = Metric)) +
  geom_point(size = 3, alpha = 0.7, na.rm = TRUE) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linetype = "dashed",
              linewidth = 0.8, na.rm = TRUE) +
  scale_colour_manual(values = metric_colors) +
  labs(x = "Specimen age (years)", y = "Metric value", colour = NULL) +
  theme_minimal(base_size = 12)
ggsave(out("Fig3A_SpecimenAge.pdf"), p_fig3a, width = 8, height = 6)


# ---- 7. Fig. 3B: per-species results ordered by phylogeny ---------------

tree <- read.tree(file.path(base_dir, "squamates_Title_Science2024_ultrametric_constrained.tre"))
tree$tip.label <- gsub("_", " ", tree$tip.label)
tree$edge.length[is.na(tree$edge.length)] <- 0.01

# Species missing from the tree are grafted onto the MRCA of their congeners,
# or as sister to a single congener.
species_in_data <- unique(within_species_age$Species1)
absent_species  <- setdiff(species_in_data, tree$tip.label)
median_bl <- median(tree$edge.length)
for (sp in absent_species) {
  congeners <- grep(paste0("^", strsplit(sp, " ")[[1]][1], " "), tree$tip.label, value = TRUE)
  if (length(congeners) > 1) {
    mrca <- getMRCA(tree, congeners)
    bl   <- max(max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)[mrca], 1e-5)
    tree <- bind.tip(tree, tip.label = sp, where = mrca, edge.length = bl)
  } else if (length(congeners) == 1) {
    sister <- which(tree$tip.label == congeners)
    pos    <- min(median_bl / 2, 0.5 * tree$edge.length[tree$edge[, 2] == sister])
    tree   <- bind.tip(tree, tip.label = sp, where = sister, position = pos, edge.length = pos)
  } else {
    stop("No congeners of ", sp, " in the tree.")
  }
}
if (length(absent_species) > 0)
  message("Species grafted onto the tree: ", paste(absent_species, collapse = ", "))
trimmed_tree <- drop.tip(tree, setdiff(tree$tip.label, species_in_data))

pdf(out("Fig3B_Phylogeny.pdf"), width = 8, height = 12)
plot(trimmed_tree)
invisible(dev.off())

p_fig3b <- within_species_age %>%
  left_join(select(results_dated, Sample1, Sample2, age), by = c("Sample1", "Sample2")) %>%
  mutate(Species1 = factor(Species1, levels = trimmed_tree$tip.label)) %>%
  ggplot(aes(LM_R_squared, Species1, colour = age)) +
  geom_point(size = 3) +
  scale_colour_viridis_c(option = "cividis", na.value = "grey70") +
  labs(x = "LM R² (gland vs venom)", y = NULL, colour = "Specimen\nage (years)") +
  theme_minimal() +
  theme(axis.text.y = element_text(face = "italic"))
ggsave(out("Fig3B_Species.pdf"), p_fig3b, width = 8, height = 12)


# ---- 8. Fig. 4: PCA of this study + published proteomes ------------------

study_genera <- unique(experimental$Genus)
pca_data <- bind_rows(
  mutate(experimental, Source = "This study"),
  mutate(published, Source = ifelse(Genus %in% study_genera, "PublishedShared", "PublishedOther"))
) %>% filter(!Sample %in% excluded_published_pca)

clade_genera <- list(
  "Vipers" = c("Bitis", "Crotalus", "Agkistrodon", "Atropoides", "Bothriechis",
               "Bothriopsis", "Bothrocophias", "Bothrops", "Calloselasma", "Cerastes",
               "Cerrophidion", "Craspedocephalus", "Daboia", "Deinagkistrodon", "Echis",
               "Gloydius", "Hypnale", "Lachesis", "Macrovipera", "Metlapilcoatlus",
               "Mixcoatlus", "Montivipera", "Ophryacus", "Ovophis", "Photobothrops",
               "Porthidium", "Protobothrops", "Sistrurus", "Trimeresurus",
               "Tropidolaemus", "Vipera"),
  "Coral snakes" = c("Calliophis", "Sinomicrurus", "Micrurus"),
  "Afro-Asian elapids" = c("Aspidelaps", "Dendroaspis", "Ophiophagus", "Walterinnesia",
                           "Hemachatus", "Naja", "Bungarus"),
  "Terrestrial hydrophiines" = c("Acanthophis", "Austrelaps", "Denisonia", "Hoplocephalus",
                                 "Drysdalia", "Echiopsis", "Micropechis", "Notechis",
                                 "Oxyuranus", "Pseudechis", "Toxicocalamus", "Tropidechis",
                                 "Vermicella", "Pseudonaja", "Suta"),
  "Sea snakes" = c("Hydrophis", "Aipysurus"),
  "Sea kraits" = c("Laticauda")
)
group_colors <- c("Vipers" = "#FFD700", "Coral snakes" = "#984EA3",
                  "Afro-Asian elapids" = "#FF8C00", "Terrestrial hydrophiines" = "#4DAF4A",
                  "Sea snakes" = "#377EB8", "Sea kraits" = "#5BC5F4", "Other" = "grey70")
genus_to_group <- setNames(rep(names(clade_genera), lengths(clade_genera)), unlist(clade_genera))
pca_data$Group <- ifelse(pca_data$Genus %in% names(genus_to_group),
                         genus_to_group[pca_data$Genus], "Other")
magenta <- "#D81B60"

# PCA (centred, unscaled) and plot. Shared-genus published proteomes are drawn
# in full colour, other genera faded; samples from this study have magenta
# outlines (circle = milked venom, diamond = gland). Sea snakes and sea kraits
# are drawn last so they are not hidden.
run_pca_plot <- function(X, file_stub, label) {
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  ve  <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  s   <- cbind(pca_data[, c("Source", "Form_short", "Group")], pca$x[, 1:2])
  s$Color <- group_colors[s$Group]
  s$Shape <- ifelse(s$Source == "This study" & s$Form_short != "MilkedVenom", 23, 21)
  sea <- s$Group %in% c("Sea snakes", "Sea kraits")
  layer <- function(i, border, cex, lwd)
    points(s$PC1[i], s$PC2[i], pch = s$Shape[i], col = border, bg = s$Color[i], cex = cex, lwd = lwd)

  pdf(out(paste0(file_stub, ".pdf")), width = 9, height = 8)
  plot(s$PC1, s$PC2, type = "n", main = paste("PCA of toxin-family composition,", label),
       xlab = sprintf("PC1 (%.1f%%)", ve[1]), ylab = sprintf("PC2 (%.1f%%)", ve[2]))
  o <- s$Source == "PublishedOther"
  points(s$PC1[o], s$PC2[o], pch = s$Shape[o], col = NA,
         bg = adjustcolor(s$Color[o], alpha.f = 0.25), cex = 1.5)
  layer(s$Source == "PublishedShared" & !sea, "black", 1.6, 1.1)
  layer(s$Source == "This study"      & !sea, magenta, 1.8, 1.4)
  layer(s$Source == "PublishedShared" &  sea, "black", 1.6, 1.1)
  layer(s$Source == "This study"      &  sea, magenta, 1.8, 1.4)
  legend("topleft", title = "Sample type",
         legend = c("Published", "This study: milked venom", "This study: preserved gland"),
         pt.bg = "white", col = c("black", magenta, magenta), pch = c(21, 21, 23),
         pt.cex = 1.1, pt.lwd = 1.2, cex = 0.85)
  legend("bottomleft", title = "Clade", legend = names(group_colors),
         fill = group_colors, border = "black", cex = 0.85, ncol = 2)
  invisible(dev.off())

  write.csv(data.frame(PC = paste0("PC", seq_along(ve)), Variance_explained_pct = round(ve, 2)),
            out("Tables", paste0(file_stub, "_variance.csv")), row.names = FALSE)
  write.csv(round(pca$rotation, 3), out("Tables", paste0(file_stub, "_loadings.csv")))
}

# log(x+1) is reported in the paper; raw proportions gave consistent results.
run_pca_plot(log1p(as.matrix(pca_data[toxins])), "Fig4_PCA_log", "log(x+1)")
run_pca_plot(as.matrix(pca_data[toxins]), "Fig4_PCA_raw", "raw proportions")


# ---- 9. Fig. S3: within-species PCAs with pie charts --------------------

type_cols <- c("Preserved gland (this study)" = "black",
               "Milked venom (this study)"    = "#D55E00",
               "Published (milked venom)"     = "grey60")

dat_s3 <- bind_rows(mutate(experimental, Source = "This study"),
                    mutate(published, Source = "Published")) %>%
  mutate(Taxon = paste(Genus, trimws(Species)),
         Type  = case_when(Source == "Published"       ~ "Published (milked venom)",
                           Form_short == "MilkedVenom" ~ "Milked venom (this study)",
                           TRUE                        ~ "Preserved gland (this study)"))

keep_taxa <- dat_s3 %>%
  group_by(Taxon) %>%
  summarise(n = n(), own = any(Source == "This study"), .groups = "drop") %>%
  filter(n >= min_samples_intraspecific, own) %>%
  pull(Taxon)
message("Species in the within-species PCAs (Fig. S3): ", length(keep_taxa))

plot_species_pca <- function(taxon, pie_scale = 0.07) {
  d <- as.data.frame(dat_s3[dat_s3$Taxon == taxon, ])
  X <- log1p(as.matrix(d[toxins]))
  X <- X[, apply(X, 2, var) > 0, drop = FALSE]
  if (ncol(X) < 2) return(NULL)
  pca <- prcomp(X, center = TRUE, scale. = FALSE)
  ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  d$PC1 <- pca$x[, 1]
  d$PC2 <- if (ncol(pca$x) > 1) pca$x[, 2] else 0
  rng   <- max(diff(range(d$PC1)), diff(range(d$PC2)), 1e-6)
  d$r   <- pie_scale * rng
  half  <- rng / 2 + 1.5 * d$r[1]
  mx <- mean(range(d$PC1)); my <- mean(range(d$PC2))

  ggplot() +
    geom_scatterpie(data = d, aes(x = PC1, y = PC2, r = r),
                    cols = toxins[colSums(d[toxins]) > 0], color = NA) +
    geom_circle(data = d, aes(x0 = PC1, y0 = PC2, r = r, colour = Type),
                linewidth = 0.8, inherit.aes = FALSE) +
    scale_fill_manual(values = color_mapping, limits = toxins, drop = FALSE) +
    scale_colour_manual(values = type_cols, limits = names(type_cols), drop = FALSE) +
    coord_equal(xlim = mx + c(-half, half), ylim = my + c(-half, half)) +
    labs(title = taxon, x = paste0("PC1 (", ve[1], "%)"),
         y = paste0("PC2 (", ifelse(length(ve) > 1, ve[2], 0), "%)")) +
    theme_bw(base_size = 10) +
    theme(plot.title = element_text(face = "italic", size = 14),
          axis.title = element_text(size = 13), axis.text = element_text(size = 11),
          panel.grid = element_blank(), legend.position = "none")
}
s3_plots <- Filter(Negate(is.null), setNames(lapply(keep_taxa, plot_species_pca), keep_taxa))

s3_legend <- extract_legend(
  ggplot(data.frame(x = 1, y = 1, Toxin = factor(toxins, levels = toxins),
                    Type = factor(rep(names(type_cols), length.out = length(toxins)),
                                  levels = names(type_cols)))) +
    geom_point(aes(x, y, fill = Toxin), shape = 22, size = 7, colour = NA) +
    geom_point(aes(x, y, colour = Type), shape = 1, size = 7, stroke = 1.4) +
    scale_fill_manual(values = color_mapping, labels = toxin_labels, name = "Toxin family") +
    scale_colour_manual(values = type_cols, name = "Sample") +
    guides(colour = guide_legend(order = 1), fill = guide_legend(order = 2, ncol = 1)) +
    theme_void(base_size = 12) +
    theme(legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14), legend.key.size = unit(0.8, "cm")))

ggsave(out("FigS3_IntraspecificPCA.pdf"),
       (wrap_plots(s3_plots, ncol = 4) | wrap_elements(full = s3_legend)) +
         plot_layout(widths = c(4, 1)),
       width = 18, height = 4 * ceiling(length(s3_plots) / 4), limitsize = FALSE)


# ---- 10. Fig. S2: toxin-family richness, gland vs milked venom ----------

fam <- dat_s3 %>%
  filter(Source == "This study") %>%
  mutate(Sample_type = ifelse(Form_short == "MilkedVenom", "Venom", "Gland"),
         Families    = rowSums(across(all_of(toxins)) >= richness_threshold),
         ID          = sub("^iBAQ (DEVC_\\d+).*", "\\1", Sample))

by_species <- fam %>%
  group_by(Taxon, Sample_type) %>%
  summarise(Families = mean(Families), .groups = "drop") %>%
  pivot_wider(names_from = Sample_type, values_from = Families) %>%
  drop_na() %>%
  rename(Unit = Taxon)

by_specimen <- fam %>%
  inner_join(data.frame(Unit = rep(same_specimen_pairs$Species, 2),
                        ID   = c(same_specimen_pairs$Venom, same_specimen_pairs$Gland)),
             by = "ID") %>%
  select(Unit, Sample_type, Families) %>%
  pivot_wider(names_from = Sample_type, values_from = Families)

run_wilcox <- function(x, label) {
  wt <- wilcox.test(x$Gland, x$Venom, paired = TRUE, exact = FALSE)
  cat(sprintf("\n%s: n = %d | gland < venom: %d | gland > venom: %d | equal: %d\n",
              label, nrow(x), sum(x$Gland < x$Venom), sum(x$Gland > x$Venom),
              sum(x$Gland == x$Venom)))
  cat(sprintf("Median families: gland %.1f | venom %.1f | paired Wilcoxon V = %.0f, p = %.3f\n",
              median(x$Gland), median(x$Venom), wt$statistic, wt$p.value))
  mutate(x, Comparison = sprintf("%s (n = %d)\nV = %.0f, p = %.3f",
                                 label, nrow(x), wt$statistic, wt$p.value))
}
richness <- bind_rows(run_wilcox(by_species, "All species"),
                      run_wilcox(by_specimen, "Same specimen"))
write.csv(richness, out("Tables", "ToxinFamilyRichness.csv"), row.names = FALSE)

p_s2 <- richness %>%
  pivot_longer(c(Gland, Venom), names_to = "Sample", values_to = "Families") %>%
  mutate(Sample = factor(ifelse(Sample == "Gland", "Preserved gland", "Milked venom"),
                         levels = c("Preserved gland", "Milked venom")),
         Comparison = factor(Comparison, levels = unique(richness$Comparison))) %>%
  ggplot(aes(Sample, Families, group = Unit)) +
  geom_line(colour = "grey70") +
  geom_point(aes(colour = Sample), size = 2.5) +
  facet_wrap(~ Comparison) +
  scale_colour_manual(values = c("Preserved gland" = "black", "Milked venom" = "#D55E00"),
                      guide = "none") +
  labs(x = NULL, y = bquote("Toxin families (" >= .(richness_threshold) * "%)")) +
  theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        strip.text = element_text(size = 11))
ggsave(out("FigS2_ToxinFamilyRichness.pdf"), p_s2, width = 8, height = 4.5)


# ---- 11. Appendix 1: pie charts and species-level regressions ------------

# Pie charts of all samples and published proteomes, 24 per page, by snake
# family. Red frame = this study, blue frame = published. Families <1% omitted.
pies <- bind_rows(mutate(experimental, Source = "This study"),
                  mutate(published, Source = "Published")) %>%
  arrange(Family, Species2)
for (fam_name in sort(unique(pies$Family))) {
  fd <- pies[pies$Family == fam_name, ]
  for (page in seq_len(ceiling(nrow(fd) / 24))) {
    pdf(out("Appendix1_PieCharts", paste0(fam_name, "_page", page, ".pdf")),
        width = 11, height = 8.5)
    par(mfrow = c(4, 6), mar = c(1, 1, 2, 1))
    for (i in ((page - 1) * 24 + 1):min(page * 24, nrow(fd))) {
      r <- fd[i, ]
      v <- unlist(r[toxins]); v <- v[v >= 1]
      if (length(v) == 0) next
      pie(v, labels = NA, col = color_mapping[names(v)], cex.main = 0.65,
          main = paste(paste(r$Genus, r$Species), r$Form_short,
                       if (r$Source == "This study")
                         paste(r$D.E.Venom.Collection.Number, r$Date.collected) else r$Sample,
                       sep = "\n"))
      box(col = ifelse(r$Source == "This study", "red", "blue"), lwd = 2)
    }
    invisible(dev.off())
  }
}

# Colour key for the pie charts
p_key <- expand.grid(row = 1:4, col = 1:6) %>%
  arrange(col, row) %>%
  mutate(Toxin = toxins) %>%
  ggplot(aes(x = col, y = 5 - row)) +
  geom_point(aes(colour = Toxin), size = 10) +
  geom_text(aes(label = toxin_labels[Toxin]), hjust = 0, nudge_x = 0.35, size = 3.3) +
  scale_colour_manual(values = color_mapping, guide = "none") +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0.5, 4.6), expand = FALSE) +
  annotate("text", x = 1.5, y = 4.55, label = "Major toxin families", fontface = "bold", size = 4.2) +
  annotate("text", x = 4.5, y = 4.55, label = "Minor toxin families", fontface = "bold", size = 4.2) +
  theme_void()
ggsave(out("Appendix1_PieCharts", "ColourKey.pdf"), p_key, width = 9, height = 4)

# For each species, mean of all glands vs mean of all milked venoms
species_form_means <- experimental %>%
  mutate(Form = ifelse(Form_short == "MilkedVenom", "MilkedVenom", "Gland")) %>%
  group_by(Species2, Form) %>%
  summarise(across(all_of(toxins), mean), .groups = "drop")

for (sp in unique(species_form_means$Species2)) {
  m <- species_form_means[species_form_means$Species2 == sp, ]
  if (!all(c("Gland", "MilkedVenom") %in% m$Form)) next
  d <- data.frame(Toxin = toxins,
                  Gland = unlist(m[m$Form == "Gland", toxins]),
                  MilkedVenom = unlist(m[m$Form == "MilkedVenom", toxins]))
  d <- d[!(d$Gland == 0 & d$MilkedVenom == 0), ]
  res <- compare_profiles(d$MilkedVenom, d$Gland, eps_default)
  p <- ggplot(d, aes(Gland, MilkedVenom, colour = Toxin)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = toxin_labels[Toxin]), box.padding = 1, max.overlaps = Inf) +
    geom_smooth(inherit.aes = FALSE, aes(Gland, MilkedVenom), method = "lm", formula = y ~ x) +
    scale_x_sqrt() + scale_y_sqrt() +
    scale_colour_manual(values = color_mapping) +
    annotate("text", x = Inf, y = Inf, hjust = 1, vjust = 1.3, size = 5, fontface = "italic",
             label = sprintf("LM R² = %.2f\nSpearman R = %.2f\nBeta pseudo-R² = %.2f",
                             res$LM_R_squared, res$Spearman_R, res$GLM_R_squared)) +
    labs(title = sp, x = "Gland (%)", y = "Venom (%)") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"), legend.position = "none")
  ggsave(out("Appendix1_Regressions", paste0(safe_name(sp), ".pdf")), p, width = 8, height = 6)
}


# ---- 12. Session info -------------------------------------------------------

writeLines(capture.output(sessionInfo()), out("sessionInfo.txt"))
message("Done. Outputs written to ", outputs_dir)
