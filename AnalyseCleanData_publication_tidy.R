# =============================================================================
# 
# Companion script for Museum Venomics paper by Esquerré et al.
#
# Inputs:
#   - SummaryResults_VenomMassSpecExperiment.csv  : experimental venom comp.
#   - PublishedProteomes.csv                      : published venom proteomes
#   - squamates_Title_Science2024_ultrametric_constrained.tre : phylogeny from Title et al. 2024
#
# Outputs are written to <base_dir>/outputs/ (created automatically).
# =============================================================================


# ---- USER CONFIGURATION -----------------------------------------------------
# Set this to the folder containing your three input files.
# Use a full path, e.g.:
#   base_dir <- "/Users/yourname/Documents/venomics"   # Mac / Linux
#   base_dir <- "C:/Users/yourname/Documents/venomics" # Windows

base_dir <- "."   # <- CHANGE THIS

# Output folder (defaults to a subfolder inside base_dir; change if needed).
outputs_dir <- file.path(base_dir, "outputs")
# -----------------------------------------------------------------------------


# ---- 1. Packages ------------------------------------------------------------

required_pkgs <- c(
  "dplyr", "tidyr", "reshape2",
  "ggplot2", "ggrepel", "ggfortify", "patchwork", "cowplot", "plotly",
  "RColorBrewer", "grDevices", "pals", "colorspace", "scico", "grid",
  "betareg", "ape", "phytools",
  "lme4", "MuMIn"
)

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing_pkgs) > 0) {
  stop(
    "Missing packages: ", paste(missing_pkgs, collapse = ", "),
    "\nInstall with: install.packages(c(\"",
    paste(missing_pkgs, collapse = "\", \""), "\"))"
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(ggplot2)
  library(ggrepel)
  library(ggfortify)
  library(patchwork)
  library(cowplot)
  library(plotly)
  library(RColorBrewer)
  library(grDevices)
  library(pals)
  library(colorspace)
  library(scico)
  library(grid)
  library(betareg)
  library(ape)
  library(phytools)
  library(lme4)
  library(MuMIn)
})


# ---- 2. Output directories --------------------------------------------------

dir_pies_indiv         <- file.path(outputs_dir, "PieCharts")
dir_pies_indiv_filt    <- file.path(outputs_dir, "PieCharts_filteredOver1percent")
dir_pies_spavg         <- file.path(outputs_dir, "PieCharts_AveragePerSpecies")
dir_pies_spavg_filt    <- file.path(outputs_dir, "PieCharts_AveragePerSpecies_FilteredOnePercent")
dir_pies_byfamily      <- file.path(outputs_dir, "PieCharts_4x6_byFamily")
dir_regression         <- file.path(outputs_dir, "RegressionPlots")
dir_regression_filt    <- file.path(outputs_dir, "RegressionPlots_filtered")
dir_regression_combo   <- file.path(outputs_dir, "RegressionPlots_combined")
dir_tables             <- file.path(outputs_dir, "Tables")

for (d in c(outputs_dir, dir_pies_indiv, dir_pies_indiv_filt,
            dir_pies_spavg, dir_pies_spavg_filt, dir_pies_byfamily,
            dir_regression, dir_regression_filt, dir_regression_combo,
            dir_tables)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}


# ---- 3. Load and prepare data ----------------------------------------------

# Strip a leading row-number column if present (e.g. "X" added by write.csv).
read_proteome_csv <- function(path) {
  x <- read.csv(path, check.names = TRUE)
  drop <- names(x)[1] %in% c("X", "X.", "X...", "")
  if (drop) x <- x[, -1, drop = FALSE]
  x
}

# Experimental venom composition per sample
summary_group  <- read_proteome_csv(file.path(base_dir, "SummaryResults_VenomMassSpecExperiment.csv"))
summary_group3 <- summary_group

# Published proteomes (NA -> 0 in toxin columns starting at col 12)
published <- read_proteome_csv(file.path(base_dir, "PublishedProteomes.csv"))
published[, 12:ncol(published)] <- lapply(
  published[, 12:ncol(published)],
  function(x) replace(x, is.na(x), 0)
)

# Align columns. Add empty Unidentified column to the experimental table to
# match the published one. If the published table carries a
# `Keep_for_coanalysis` curation flag, filter on it; otherwise treat the file
# as already curated.
summary_group3$Unidentified <- 0
if ("Keep_for_coanalysis" %in% names(published)) {
  published_trimmed <- published[published$Keep_for_coanalysis == "Yes", ]
} else {
  published_trimmed <- published
}
published_trimmed <- published_trimmed[, names(published_trimmed) %in% names(summary_group3)]

# Versions with Other & Unidentified (last 2 cols) removed and percentages
# re-normalised so each sample sums to 100. Rows whose retained toxin columns
# sum to zero (no abundance in any of the 24 shared families) are dropped, as
# they carry no comparable composition information.
drop_other_unidentified <- function(df, label = "table") {
  toxin_start <- 12
  to_drop <- which(names(df) %in% c("Other", "Unidentified"))
  out <- if (length(to_drop) > 0) df[, -to_drop] else df
  toxin_end <- ncol(out)

  rs <- rowSums(out[, toxin_start:toxin_end])
  empty <- rs == 0 | is.na(rs)
  if (any(empty)) {
    message("Dropping ", sum(empty), " rows from ", label,
            " with zero abundance in shared toxin families.")
    out <- out[!empty, , drop = FALSE]
    rs  <- rs[!empty]
  }
  out[, toxin_start:toxin_end] <- out[, toxin_start:toxin_end] / rs * 100
  out
}

summary_group3_no_others <- drop_other_unidentified(summary_group3,  "summary_group3")
published_no_others      <- drop_other_unidentified(published_trimmed, "published_trimmed")

# Sanity checks: every retained row sums to 100 (within float tolerance).
stopifnot(
  all(abs(rowSums(summary_group3_no_others[, 12:ncol(summary_group3_no_others)]) - 100) < 1e-6),
  all(abs(rowSums(published_no_others[, 12:ncol(published_no_others)])           - 100) < 1e-6)
)


# ---- 4. Combined dataset and species averages ------------------------------

combined_df <- bind_rows(published_no_others, summary_group3_no_others)

meta_cols  <- c("Sample", "D.E.Venom.Collection.Number", "Family", "Genus",
                "Species", "Species2", "Form_short", "Sex", "Date.collected",
                "SVL", "Tail.length")
toxin_cols <- setdiff(names(combined_df), meta_cols)

combined_df[toxin_cols] <- lapply(combined_df[toxin_cols],
                                  function(x) as.numeric(as.character(x)))

species_avg <- combined_df %>%
  group_by(Species, Genus, Family) %>%
  summarise(across(all_of(toxin_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

# PCA-ready toxin matrix (drop zero-variance columns)
toxin_matrix <- species_avg %>%
  select(all_of(toxin_cols)) %>%
  select(where(~ var(.x, na.rm = TRUE) > 0))

# Quick listing of species being added vs the published set
new_species <- setdiff(summary_group3$Species2, published_no_others$Species2)
message("Species added (not in published set): ",
        if (length(new_species) == 0) "(none)" else paste(new_species, collapse = ", "))


# ---- 5. Colour palette ------------------------------------------------------

color_mapping <- c(
  "X3FTx" = "#F47B5B", "PLA2"  = "#B5EFB5", "SVMP"  = "#FBE426",
  "SVSP"  = "#2ED9FF", "CTL"   = "#AA0DFE", "CRiSP" = "#AAF400",
  "Kun"   = "#3283FE", "LAAO"  = "#F7E1A0", "X5N"   = "#F8A19F",
  "AChE"  = "#F6222E", "AmPep" = "#DD6F91", "Cys"   = "#F45366",
  "Dis"   = "#1C7F93", "Hyal"  = "#1C8356", "NGF"   = "#B10DA1",
  "NP"    = "#1CBE4F", "Oha.Vesp" = "#7ED7D1", "PDE" = "#C075A6",
  "PLB"   = "#FC1CBF", "PLC"   = "#FA0087", "VEGF"  = "#BDCDFF",
  "Venom.factor" = "#822E1C", "Venom.peroxiredoxin" = "#C4451C",
  "Waprin" = "#782AB6", "Other" = "#A2A2A5", "Unidentified" = "#D6D6D6"
)


# ---- 6. Per-sample pie charts (one PDF each) -------------------------------

# `df`         : data frame with metadata in cols 1:11, toxins in cols 12+
# `dir_full`   : directory for un-filtered charts
# `dir_filt`   : directory for charts with toxins >=1%
# `name_fn`    : function(row) returning the sample identifier used for the file
plot_pie_per_row <- function(df, dir_full, dir_filt,
                             name_fn = function(r) paste(r[4], r[5], r[7], r[2], r[9])) {
  toxin_idx <- 12:ncol(df)
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    sample_name <- name_fn(row)

    values <- as.numeric(row[, toxin_idx])
    names(values) <- names(df)[toxin_idx]
    used_colors <- color_mapping[names(values)]

    # Full
    pdf(file = file.path(dir_full, paste0(sample_name, "_pie_chart.pdf")))
    pie(values,
        labels = sprintf("%s: %1.2f", names(values), values),
        col = used_colors, main = sample_name)
    dev.off()

    # Filtered (>=1%)
    keep <- values >= 1
    if (any(keep)) {
      fv <- values[keep]
      pdf(file = file.path(dir_filt, paste0(sample_name, "_filtered_pie_chart.pdf")))
      pie(fv,
          labels = sprintf("%s: %1.2f", names(fv), fv),
          col = color_mapping[names(fv)],
          main = paste(sample_name, "(>=1%)"))
      dev.off()
    }
  }
}

plot_pie_per_row(summary_group3_no_others, dir_pies_indiv, dir_pies_indiv_filt)
plot_pie_per_row(published_no_others,      dir_pies_indiv, dir_pies_indiv_filt)


# ---- 7. Species-average pie charts -----------------------------------------

# species_avg has metadata in cols 1:3 (Species, Genus, Family) and toxins after.
plot_pie_species_avg <- function(df, dir_full, dir_filt) {
  toxin_idx <- 4:ncol(df)
  for (i in seq_len(nrow(df))) {
    row <- df[i, ]
    sample_name <- paste(row[[2]], row[[1]])  # Genus + Species

    values <- as.numeric(row[, toxin_idx])
    names(values) <- names(df)[toxin_idx]
    used_colors <- color_mapping[names(values)]

    pdf(file = file.path(dir_full, paste0(sample_name, "_pie_chart.pdf")))
    pie(values,
        labels = sprintf("%s: %1.2f", names(values), values),
        col = used_colors, main = sample_name)
    dev.off()

    keep <- values >= 1
    if (any(keep)) {
      fv <- values[keep]
      pdf(file = file.path(dir_filt, paste0(sample_name, "_filtered_pie_chart.pdf")))
      pie(fv,
          labels = sprintf("%s: %1.2f", names(fv), fv),
          col = color_mapping[names(fv)],
          main = paste(sample_name, "(>=1%)"))
      dev.off()
    }
  }
}

plot_pie_species_avg(species_avg, dir_pies_spavg, dir_pies_spavg_filt)


# ---- 8. Per-family multi-pie PDFs ------------------------------------------

summary_group3_no_others$Source <- "Experimental"
published_no_others$Source      <- "Published"

combined_no_others <- bind_rows(summary_group3_no_others, published_no_others) %>%
  arrange(Family, Species2)

plot_pies_by_family <- function(df, family_name, output_dir) {
  family_df  <- df[df$Family == family_name, , drop = FALSE]
  n_samples  <- nrow(family_df)
  n_per_page <- 24
  n_pages    <- ceiling(n_samples / n_per_page)

  toxin_names <- setdiff(names(df), c(meta_cols, "Source"))

  for (page in seq_len(n_pages)) {
    file_name <- file.path(output_dir, paste0(family_name, "_page", page, ".pdf"))
    pdf(file = file_name, width = 11, height = 8.5)
    par(mfrow = c(4, 6), mar = c(1, 1, 2, 1))

    start_idx <- (page - 1) * n_per_page + 1
    end_idx   <- min(start_idx + n_per_page - 1, n_samples)

    for (i in start_idx:end_idx) {
      row <- family_df[i, ]
      sample_title <- paste(
        paste(row$Genus, row$Species),
        row$Form_short,
        paste(row$D.E.Venom.Collection.Number, row$Date.collected),
        sep = "\n"
      )

      values <- as.numeric(row[, toxin_names])
      names(values) <- toxin_names
      values <- values[!is.na(values) & values >= 1]
      if (length(values) == 0) next

      pie(values, labels = NA, col = color_mapping[names(values)],
          main = sample_title, cex.main = 0.65)
      box(col = ifelse(row$Source == "Experimental", "red", "blue"), lwd = 2)
    }
    dev.off()
  }
}

for (fam in sort(unique(combined_no_others$Family))) {
  plot_pies_by_family(combined_no_others, fam, dir_pies_byfamily)
}


# ---- 9. Toxin legend grid (visual key) -------------------------------------

order_vec <- c(
  "X3FTx", "PLA2", "SVMP", "SVSP", "CTL", "CRiSP", "Kun", "LAAO",  # Major (8)
  "X5N", "AChE", "AmPep", "Cys", "Dis", "Hyal", "NGF", "NP",       # Minor (16)
  "Oha.Vesp", "PDE", "PLB", "PLC", "VEGF", "Venom.factor",
  "Venom.peroxiredoxin", "Waprin"
)
major <- order_vec[1:8]
minor <- order_vec[9:24]

legend_grid <- expand.grid(row = 1:4, col = 1:6) %>%
  arrange(col, row) %>%
  mutate(
    Family = c(major, minor),
    Group  = ifelse(col <= 2, "Major toxin families", "Minor toxin families")
  )

legend_plot_grid <- ggplot(legend_grid, aes(x = col, y = 5 - row)) +
  geom_tile(aes(width = 0.98, height = 0.98), fill = NA,
            color = "grey90", linewidth = 0.2) +
  geom_point(aes(color = Family), size = 10) +
  geom_text(aes(label = Family), hjust = 0, nudge_x = 0.35, size = 3.3) +
  scale_color_manual(values = color_mapping, guide = "none") +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0.5, 4.6), expand = FALSE) +
  theme_void() +
  annotate("text", x = 1.5, y = 4.55,
           label = "Major toxin families", fontface = "bold", size = 4.2) +
  annotate("text", x = 4.5, y = 4.55,
           label = "Minor toxin families", fontface = "bold", size = 4.2)

ggsave(file.path(outputs_dir, "ToxinLegendGrid.pdf"),
       legend_plot_grid, width = 9, height = 4)


# ---- 10. Regression: gland vs venom (per species) --------------------------

species_sample_m <- summary_group3_no_others %>%
  group_by(Species2, Form_short) %>%
  summarise(across(all_of(toxin_cols), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")

species_sample_m$Form_short[
  species_sample_m$Form_short %in% c("PreservedFixedGland", "FrozenEthanolGland")
] <- "Gland"

valid_species <- species_sample_m %>%
  group_by(Species2) %>%
  summarise(
    has_gland       = any(Form_short == "Gland"),
    has_milkedvenom = any(Form_short == "MilkedVenom"),
    .groups = "drop"
  ) %>%
  filter(has_gland & has_milkedvenom) %>%
  pull(Species2)

# Beta-regression-safe rescaling: clamp into (0, 1).
clamp_for_beta <- function(x) {
  x <- x / 100
  x[x <= 0] <- 1e-9
  x[x >= 1] <- 1 - 1e-9
  x
}

for (species in valid_species) {
  rows <- species_sample_m[species_sample_m$Species2 == species, ]
  specx <- as.data.frame(t(rows))
  speciesname <- specx[1, 1]
  colnames(specx) <- specx[2, ]
  specx <- as.data.frame(specx[-(1:2), ])

  specx$Gland       <- as.numeric(specx$Gland)
  specx$MilkedVenom <- as.numeric(specx$MilkedVenom)
  specx <- specx[!(specx$MilkedVenom == 0 & specx$Gland == 0), ]
  if (nrow(specx) == 0) next

  lm_model     <- lm(MilkedVenom ~ Gland, data = specx)
  r_squared    <- summary(lm_model)$adj.r.squared
  spearman_cor <- cor.test(specx$Gland, specx$MilkedVenom, method = "spearman")

  glm_r_squared <- NA_real_
  tryCatch({
    specx_beta <- specx
    specx_beta$Gland       <- clamp_for_beta(specx_beta$Gland)
    specx_beta$MilkedVenom <- clamp_for_beta(specx_beta$MilkedVenom)
    glm_model     <- betareg(MilkedVenom ~ Gland, data = specx_beta)
    glm_r_squared <- summary(glm_model)$pseudo.r.squared
  }, error = function(e) message("GLM failed for species ", species, ": ", e$message))

  specx$Toxin <- rownames(specx)
  base_plot <- ggplot(specx, aes(x = Gland, y = MilkedVenom, color = Toxin)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = Toxin), box.padding = 1, max.overlaps = Inf) +
    scale_x_sqrt() + scale_y_sqrt() +
    geom_smooth(inherit.aes = FALSE, method = "lm",
                aes(x = Gland, y = MilkedVenom)) +
    labs(title = species, x = "Gland", y = "Venom", color = "Toxin Group") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          legend.position = "none") +
    scale_color_manual(values = color_mapping) +
    annotate("text", x = Inf, y = Inf,
             label = paste("LM R^2 = ", round(r_squared, 2)),
             hjust = 1, vjust = 2.5, size = 5, fontface = "italic") +
    annotate("text", x = Inf, y = Inf,
             label = paste("Spearman R = ", round(spearman_cor$estimate, 2)),
             hjust = 1, vjust = 4, size = 5, fontface = "italic") +
    annotate("text", x = Inf, y = Inf,
             label = paste("GLM R^2 = ", round(glm_r_squared, 2)),
             hjust = 1, vjust = 5.5, size = 5, fontface = "italic")

  ggsave(file.path(dir_regression, paste0(species, "_plot_with_metrics.pdf")),
         base_plot, width = 8, height = 6)

  # Filtered (>=1% in either)
  specx_filtered <- specx %>% filter(!(Gland < 1 & MilkedVenom < 1))
  if (nrow(specx_filtered) == 0) next

  lm_model_f     <- lm(MilkedVenom ~ Gland, data = specx_filtered)
  r_squared_f    <- summary(lm_model_f)$adj.r.squared
  spearman_cor_f <- cor.test(specx_filtered$Gland, specx_filtered$MilkedVenom,
                             method = "spearman")

  glm_r_squared_f <- NA_real_
  tryCatch({
    specx_beta <- specx_filtered
    specx_beta$Gland       <- clamp_for_beta(specx_beta$Gland)
    specx_beta$MilkedVenom <- clamp_for_beta(specx_beta$MilkedVenom)
    glm_model_f     <- betareg(MilkedVenom ~ Gland, data = specx_beta)
    glm_r_squared_f <- summary(glm_model_f)$pseudo.r.squared
  }, error = function(e) message("Filtered GLM failed for species ", species, ": ", e$message))

  filt_plot <- ggplot(specx_filtered, aes(x = Gland, y = MilkedVenom, color = Toxin)) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = Toxin), box.padding = 1, max.overlaps = Inf) +
    scale_x_sqrt() + scale_y_sqrt() +
    geom_smooth(inherit.aes = FALSE, method = "lm",
                aes(x = Gland, y = MilkedVenom)) +
    labs(title = paste(species, "(filtered >=1%)"),
         x = "Gland", y = "Venom", color = "Toxin Group") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          legend.position = "none") +
    scale_color_manual(values = color_mapping) +
    annotate("text", x = Inf, y = Inf,
             label = paste("LM R^2 = ", round(r_squared_f, 2)),
             hjust = 1, vjust = 2.5, size = 5, fontface = "italic") +
    annotate("text", x = Inf, y = Inf,
             label = paste("Spearman R = ", round(spearman_cor_f$estimate, 2)),
             hjust = 1, vjust = 4, size = 5, fontface = "italic") +
    annotate("text", x = Inf, y = Inf,
             label = paste("GLM R^2 = ", round(glm_r_squared_f, 2)),
             hjust = 1, vjust = 5.5, size = 5, fontface = "italic")

  ggsave(file.path(dir_regression_filt, paste0(species, "_filtered_plot_with_metrics.pdf")),
         filt_plot, width = 8, height = 6)

  combined_plot <- base_plot + filt_plot + patchwork::plot_layout(ncol = 2)
  ggsave(file.path(dir_regression_combo, paste0(species, "_combined_plot.pdf")),
         combined_plot, width = 16, height = 6)
}


# ---- 11. Within-species correlation table ----------------------------------

# All MilkedVenom samples are paired against every other sample of the same
# species (typically gland samples) and a battery of correlation metrics is
# computed. Pairs where both values are zero are excluded (they destabilise
# beta regression).

results_table <- data.frame(
  Genus = character(), Species = character(),
  Sample1 = character(), Sample2 = character(),
  Date_Collected_Sample2 = character(),
  LM_R_squared = numeric(), LM_p_value = numeric(),
  Pearson_R = numeric(), Pearson_p = numeric(),
  Spearman_R = numeric(), Spearman_p = numeric(),
  GLM_R_squared = numeric(),
  GLM_r_squared_byhand = numeric(),
  GLM_r_squared_byhand_linear = numeric(),
  GLM_AIC = numeric(), GLM_p = numeric(),
  stringsAsFactors = FALSE
)

# Toxin column indices in summary_group3_no_others (start at 12).
sg_toxin_idx <- which(names(summary_group3_no_others) %in% toxin_cols)
sg_toxin_idx <- sg_toxin_idx[sg_toxin_idx >= 12]  # drop any meta hits

unique_species <- unique(summary_group3_no_others$Species2)

for (species in unique_species) {
  subset_species       <- subset(summary_group3_no_others, Species2 == species)
  milked_venom_samples <- subset(subset_species, Form_short == "MilkedVenom")
  if (nrow(milked_venom_samples) == 0) next

  for (i in seq_len(nrow(milked_venom_samples))) {
    sample1      <- milked_venom_samples[i, ]
    sample1_data <- as.numeric(sample1[, sg_toxin_idx])
    sample1_name <- sample1$Sample

    for (j in seq_len(nrow(subset_species))) {
      sample2 <- subset_species[j, ]
      if (sample1$Sample == sample2$Sample) next

      sample2_data           <- as.numeric(sample2[, sg_toxin_idx])
      sample2_name           <- sample2$Sample
      date_collected_sample2 <- sample2$Date.collected

      valid <- !is.na(sample1_data) & !is.na(sample2_data)
      a <- sample1_data[valid]; b <- sample2_data[valid]
      keep <- !(a == 0 & b == 0)
      a <- a[keep]; b <- b[keep]
      if (length(a) < 2) next

      specx <- data.frame(Gland = b, MilkedVenom = a)

      lm_model     <- lm(MilkedVenom ~ Gland, data = specx)
      r_squared    <- summary(lm_model)$adj.r.squared
      lm_p_value   <- summary(lm_model)$coefficients["Gland", "Pr(>|t|)"]
      pearson_cor  <- cor.test(specx$Gland, specx$MilkedVenom, method = "pearson")
      spearman_cor <- cor.test(specx$Gland, specx$MilkedVenom, method = "spearman")

      glm_r_squared               <- NA_real_
      glm_r_squared_byhand        <- NA_real_
      glm_r_squared_byhand_linear <- NA_real_
      glm_aic                     <- NA_real_
      glm_p                       <- NA_real_

      tryCatch({
        specx_beta <- specx
        specx_beta$Gland       <- clamp_for_beta(specx_beta$Gland)
        specx_beta$MilkedVenom <- clamp_for_beta(specx_beta$MilkedVenom)

        glm_model            <- betareg(MilkedVenom ~ Gland, data = specx_beta)
        glm_r_squared        <- summary(glm_model)$pseudo.r.squared
        glm_r_squared_byhand <- cor(
          log(specx_beta$MilkedVenom / (1 - specx_beta$MilkedVenom)),
          predict(glm_model, type = "link")
        ) ^ 2
        glm_r_squared_byhand_linear <- cor(
          specx_beta$MilkedVenom,
          predict(glm_model, type = "response")
        ) ^ 2
        glm_aic <- AIC(glm_model)
        glm_p   <- summary(glm_model)$coefficients$mean["Gland", "Pr(>|z|)"]
      }, error = function(e) message("GLM failed for species ", species, ": ", e$message))

      results_table <- rbind(results_table, data.frame(
        Genus = sample1$Genus, Species = species,
        Sample1 = sample1_name, Sample2 = sample2_name,
        Date_Collected_Sample2 = date_collected_sample2,
        LM_R_squared = round(r_squared, 3), LM_p_value = round(lm_p_value, 3),
        Pearson_R = round(pearson_cor$estimate, 3),
        Pearson_p = round(pearson_cor$p.value, 3),
        Spearman_R = round(spearman_cor$estimate, 3),
        Spearman_p = round(spearman_cor$p.value, 3),
        GLM_R_squared = round(glm_r_squared, 3),
        GLM_r_squared_byhand = round(glm_r_squared_byhand, 3),
        GLM_r_squared_byhand_linear = round(glm_r_squared_byhand_linear, 3),
        GLM_AIC = round(glm_aic, 3), GLM_p = round(glm_p, 3)
      ))
    }
  }
}

write.csv(results_table,
          file.path(dir_tables, "CorrelationResultsVenomVsGlands.csv"),
          row.names = FALSE)

# Quick distribution plot of the three R-like metrics
plot_data <- dplyr::select(results_table, LM_R_squared, Spearman_R, GLM_R_squared) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value")

dist_plot <- ggplot(plot_data, aes(x = Value, fill = Metric)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30,
                 alpha = 0.5, position = "identity") +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("LM_R_squared"  = "#F43151",
                               "Spearman_R"    = "#8362DD",
                               "GLM_R_squared" = "#BAEA75")) +
  labs(title = "Within-species metric distributions",
       x = "Metric value", y = "Density", fill = "Metric") +
  theme_minimal() +
  theme(plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text   = element_text(size = 12),
        axis.title  = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10))

ggsave(file.path(outputs_dir, "WithinSpecies_metric_distributions.pdf"),
       dist_plot, width = 8, height = 6)


# ---- 12. Within-species pairwise correlations (alternative table) ----------

results_within_species <- data.frame(
  Genus = character(), Species = character(),
  Sample1 = character(), Sample2 = character(),
  Pearson_R = numeric(), Pearson_p = numeric(),
  Spearman_R = numeric(), Spearman_p = numeric(),
  LM_R_squared = numeric(), LM_p = numeric(),
  GLM_R_squared = numeric(),
  GLM_r_squared_byhand = numeric(),
  GLM_r_squared_byhand_linear = numeric(),
  GLM_AIC = numeric(), GLM_p = numeric(),
  stringsAsFactors = FALSE
)

for (species in unique_species) {
  subset_species       <- subset(summary_group3_no_others, Species2 == species)
  milked_venom_samples <- subset(subset_species, Form_short == "MilkedVenom")
  if (nrow(milked_venom_samples) == 0) next

  for (i in seq_len(nrow(milked_venom_samples))) {
    sample1      <- milked_venom_samples[i, ]
    sample1_data <- as.numeric(sample1[, sg_toxin_idx])
    sample1_name <- sample1$Sample

    for (j in seq_len(nrow(subset_species))) {
      sample2 <- subset_species[j, ]
      if (sample1$Sample == sample2$Sample) next

      sample2_data <- as.numeric(sample2[, sg_toxin_idx])
      sample2_name <- sample2$Sample

      valid <- !is.na(sample1_data) & !is.na(sample2_data)
      a <- sample1_data[valid]; b <- sample2_data[valid]
      keep <- !(a == 0 & b == 0)
      a <- a[keep]; b <- b[keep]
      if (length(a) < 2) next

      specx <- data.frame(Gland = b, MilkedVenom = a)

      lm_model     <- lm(MilkedVenom ~ Gland, data = specx)
      lm_r_squared <- summary(lm_model)$adj.r.squared
      lm_p         <- summary(lm_model)$coefficients["Gland", "Pr(>|t|)"]
      pearson_cor  <- cor.test(specx$Gland, specx$MilkedVenom, method = "pearson")
      spearman_cor <- cor.test(specx$Gland, specx$MilkedVenom, method = "spearman")

      specx_beta <- specx
      specx_beta$Gland       <- clamp_for_beta(specx_beta$Gland)
      specx_beta$MilkedVenom <- clamp_for_beta(specx_beta$MilkedVenom)

      glm_r_squared               <- NA_real_
      GLM_r_squared_byhand        <- NA_real_
      GLM_r_squared_byhand_linear <- NA_real_
      glm_aic                     <- NA_real_
      glm_p                       <- NA_real_

      if (all(specx_beta$MilkedVenom > 0 & specx_beta$MilkedVenom < 1)) {
        tryCatch({
          glm_model            <- betareg(MilkedVenom ~ Gland, data = specx_beta)
          glm_r_squared        <- summary(glm_model)$pseudo.r.squared
          GLM_r_squared_byhand <- cor(
            log(specx_beta$MilkedVenom / (1 - specx_beta$MilkedVenom)),
            predict(glm_model, type = "link")
          ) ^ 2
          GLM_r_squared_byhand_linear <- cor(
            specx_beta$MilkedVenom,
            predict(glm_model, type = "response")
          ) ^ 2
          glm_aic <- AIC(glm_model)
          glm_p   <- summary(glm_model)$coefficients$mean["Gland", "Pr(>|z|)"]
        }, error = function(e) message("GLM failed for one pair: ", e$message))
      }

      results_within_species <- rbind(results_within_species, data.frame(
        Genus = sample1$Genus, Species = species,
        Sample1 = sample1_name, Sample2 = sample2_name,
        Pearson_R = round(pearson_cor$estimate, 3),
        Pearson_p = round(pearson_cor$p.value, 3),
        Spearman_R = round(spearman_cor$estimate, 3),
        Spearman_p = round(spearman_cor$p.value, 3),
        LM_R_squared = round(lm_r_squared, 3), LM_p = round(lm_p, 3),
        GLM_R_squared = round(glm_r_squared, 3),
        GLM_r_squared_byhand = round(GLM_r_squared_byhand, 3),
        GLM_r_squared_byhand_linear = round(GLM_r_squared_byhand_linear, 3),
        GLM_AIC = round(glm_aic, 3), GLM_p = round(glm_p, 3)
      ))
    }
  }
}


# ---- 13. Within-genus, between-species comparisons -------------------------

results_within_genus <- data.frame(
  Genus1 = character(), Species1 = character(), Sample1 = character(),
  Genus2 = character(), Species2 = character(), Sample2 = character(),
  Pearson_R = numeric(), Pearson_p = numeric(),
  Spearman_R = numeric(), Spearman_p = numeric(),
  LM_R_squared = numeric(), LM_p = numeric(),
  GLM_R_squared = numeric(),
  GLM_r_squared_byhand = numeric(),
  GLM_r_squared_byhand_linear = numeric(),
  GLM_AIC = numeric(), GLM_p = numeric(),
  stringsAsFactors = FALSE
)

unique_genera <- unique(summary_group3_no_others$Genus)

for (genus in unique_genera) {
  subset_genus         <- subset(summary_group3_no_others, Genus == genus)
  milked_venom_samples <- subset(subset_genus, Form_short == "MilkedVenom")
  if (nrow(milked_venom_samples) == 0) next

  for (k in seq_len(nrow(milked_venom_samples))) {
    milked_sample     <- milked_venom_samples[k, ]
    milked_venom_data <- as.numeric(milked_sample[, sg_toxin_idx])
    species1          <- milked_sample$Species2
    sample1           <- milked_sample$Sample

    for (i in seq_len(nrow(subset_genus))) {
      current_sample <- subset_genus[i, ]
      if (current_sample$Sample == sample1 || current_sample$Species2 == species1) next

      current_sample_data <- as.numeric(current_sample[, sg_toxin_idx])
      sample2  <- current_sample$Sample
      species2 <- current_sample$Species2

      valid <- !is.na(milked_venom_data) & !is.na(current_sample_data)
      a <- milked_venom_data[valid]; b <- current_sample_data[valid]
      keep <- !(a == 0 & b == 0)
      a <- a[keep]; b <- b[keep]
      if (length(a) < 2) next

      specx <- data.frame(Gland = b, MilkedVenom = a)

      lm_model     <- lm(MilkedVenom ~ Gland, data = specx)
      lm_r_squared <- summary(lm_model)$adj.r.squared
      lm_p         <- summary(lm_model)$coefficients["Gland", "Pr(>|t|)"]
      pearson_cor  <- cor.test(specx$Gland, specx$MilkedVenom, method = "pearson")
      spearman_cor <- cor.test(specx$Gland, specx$MilkedVenom, method = "spearman")

      specx_beta <- specx
      specx_beta$Gland       <- clamp_for_beta(specx_beta$Gland)
      specx_beta$MilkedVenom <- clamp_for_beta(specx_beta$MilkedVenom)

      glm_r_squared               <- NA_real_
      glm_r_squared_byhand        <- NA_real_
      glm_r_squared_byhand_linear <- NA_real_
      glm_aic                     <- NA_real_
      glm_p                       <- NA_real_

      if (all(specx_beta$MilkedVenom > 0 & specx_beta$MilkedVenom < 1)) {
        tryCatch({
          glm_model            <- betareg(MilkedVenom ~ Gland, data = specx_beta)
          glm_r_squared        <- summary(glm_model)$pseudo.r.squared
          glm_r_squared_byhand <- cor(
            log(specx_beta$MilkedVenom / (1 - specx_beta$MilkedVenom)),
            predict(glm_model, type = "link")
          ) ^ 2
          glm_r_squared_byhand_linear <- cor(
            specx_beta$MilkedVenom,
            predict(glm_model, type = "response")
          ) ^ 2
          glm_aic <- AIC(glm_model)
          glm_p   <- summary(glm_model)$coefficients$mean["Gland", "Pr(>|z|)"]
        }, error = function(e) message("GLM failed for one pair: ", e$message))
      }

      results_within_genus <- rbind(results_within_genus, data.frame(
        Genus1 = genus, Species1 = species1, Sample1 = sample1,
        Genus2 = genus, Species2 = species2, Sample2 = sample2,
        Pearson_R = round(pearson_cor$estimate, 3),
        Pearson_p = round(pearson_cor$p.value, 3),
        Spearman_R = round(spearman_cor$estimate, 3),
        Spearman_p = round(spearman_cor$p.value, 3),
        LM_R_squared = round(lm_r_squared, 3), LM_p = round(lm_p, 3),
        GLM_R_squared = round(glm_r_squared, 3),
        GLM_r_squared_byhand = round(glm_r_squared_byhand, 3),
        GLM_r_squared_byhand_linear = round(glm_r_squared_byhand_linear, 3),
        GLM_AIC = round(glm_aic, 3), GLM_p = round(glm_p, 3)
      ))
    }
  }
}

write.csv(results_within_genus,
          file.path(dir_tables, "results_within_genus_comparisons.csv"),
          row.names = FALSE)


# ---- 14. Outside-genus comparisons -----------------------------------------

results_outside_genus <- data.frame(
  Genus1 = character(), Species1 = character(),
  Genus2 = character(), Species2 = character(),
  Pearson_R = numeric(), Pearson_p = numeric(),
  Spearman_R = numeric(), Spearman_p = numeric(),
  LM_R_squared = numeric(), LM_p = numeric(),
  GLM_R_squared = numeric(),
  GLM_r_squared_byhand = numeric(),
  GLM_r_squared_byhand_linear = numeric(),
  GLM_AIC = numeric(), GLM_p = numeric(),
  stringsAsFactors = FALSE
)

for (genus1 in unique_genera) {
  subset_genus1        <- subset(summary_group3_no_others, Genus == genus1)
  milked_venom_samples <- subset(subset_genus1, Form_short == "MilkedVenom")
  if (nrow(milked_venom_samples) == 0) next

  for (i in seq_len(nrow(milked_venom_samples))) {
    sample1  <- milked_venom_samples[i, ]
    species1 <- sample1$Species2

    for (genus2 in setdiff(unique_genera, genus1)) {
      subset_genus2 <- subset(summary_group3_no_others, Genus == genus2)

      for (j in seq_len(nrow(subset_genus2))) {
        sample2  <- subset_genus2[j, ]
        species2 <- sample2$Species2

        a <- as.numeric(sample1[, sg_toxin_idx])
        b <- as.numeric(sample2[, sg_toxin_idx])
        valid <- !is.na(a) & !is.na(b)
        a <- a[valid]; b <- b[valid]
        keep <- !(a == 0 & b == 0)
        a <- a[keep]; b <- b[keep]
        if (length(a) < 2) next

        specx <- data.frame(Gland = b, MilkedVenom = a)

        lm_model     <- lm(MilkedVenom ~ Gland, data = specx)
        lm_r_squared <- summary(lm_model)$adj.r.squared
        lm_p         <- summary(lm_model)$coefficients["Gland", "Pr(>|t|)"]
        pearson_cor  <- cor.test(specx$Gland, specx$MilkedVenom, method = "pearson")
        spearman_cor <- cor.test(specx$Gland, specx$MilkedVenom, method = "spearman")

        specx_beta <- specx
        specx_beta$Gland       <- clamp_for_beta(specx_beta$Gland)
        specx_beta$MilkedVenom <- clamp_for_beta(specx_beta$MilkedVenom)

        glm_r_squared               <- NA_real_
        glm_r_squared_byhand        <- NA_real_
        glm_r_squared_byhand_linear <- NA_real_
        glm_aic                     <- NA_real_
        glm_p                       <- NA_real_

        if (all(specx_beta$MilkedVenom > 0 & specx_beta$MilkedVenom < 1)) {
          tryCatch({
            glm_model            <- betareg(MilkedVenom ~ Gland, data = specx_beta)
            glm_r_squared        <- summary(glm_model)$pseudo.r.squared
            glm_r_squared_byhand <- cor(
              log(specx_beta$MilkedVenom / (1 - specx_beta$MilkedVenom)),
              predict(glm_model, type = "link")
            ) ^ 2
            glm_r_squared_byhand_linear <- cor(
              specx_beta$MilkedVenom,
              predict(glm_model, type = "response")
            ) ^ 2
            glm_aic <- AIC(glm_model)
            glm_p   <- summary(glm_model)$coefficients$mean["Gland", "Pr(>|z|)"]
          }, error = function(e) message("GLM failed for one pair: ", e$message))
        }

        results_outside_genus <- rbind(results_outside_genus, data.frame(
          Genus1 = genus1, Species1 = species1,
          Genus2 = genus2, Species2 = species2,
          Pearson_R = round(pearson_cor$estimate, 3),
          Pearson_p = round(pearson_cor$p.value, 3),
          Spearman_R = round(spearman_cor$estimate, 3),
          Spearman_p = round(spearman_cor$p.value, 3),
          LM_R_squared = round(lm_r_squared, 3), LM_p = round(lm_p, 3),
          GLM_R_squared = round(glm_r_squared, 3),
          GLM_r_squared_byhand = round(glm_r_squared_byhand, 3),
          GLM_r_squared_byhand_linear = round(glm_r_squared_byhand_linear, 3),
          GLM_AIC = round(glm_aic, 3), GLM_p = round(glm_p, 3)
        ))
      }
    }
  }
}

write.csv(results_outside_genus,
          file.path(dir_tables, "results_outside_genus_comparisons.csv"),
          row.names = FALSE)


# ---- 15. Histograms / densities by comparison category ---------------------

results_within_species$Category <- "Within Species"
results_within_genus$Category   <- "Within Genus"
results_outside_genus$Category  <- "Outside Genus"

shared_cols <- c("Category", "LM_R_squared", "Spearman_R",
                 "GLM_R_squared", "GLM_r_squared_byhand_linear")

all_results <- rbind(
  results_within_species[, shared_cols],
  results_within_genus[, shared_cols],
  results_outside_genus[, shared_cols]
)
colnames(all_results) <- c("Category", "LM R^2", "Spearman", "GLM R^2", "Linear GLM R^2")

all_results_melted <- reshape2::melt(
  all_results, id.vars = "Category",
  variable.name = "Test", value.name = "Value"
)
all_results_melted$Category <- factor(
  all_results_melted$Category,
  levels = c("Within Species", "Within Genus", "Outside Genus")
)

custom_colors <- c("Within Species" = "#E69F00",
                   "Within Genus"   = "#009E73",
                   "Outside Genus"  = "#0072B2")

hist_plot <- ggplot(all_results_melted, aes(x = Value, fill = Category)) +
  geom_histogram(aes(y = after_stat(density)),
                 alpha = 0.6, position = "identity",
                 bins = 30, color = "black") +
  facet_grid(Test ~ ., scales = "free_y") +
  scale_fill_manual(values = custom_colors) +
  labs(title = "Relative frequency of metrics by test type",
       x = "Metric value", y = "Density (relative frequency)",
       fill = "Comparison type") +
  theme_minimal() +
  theme(plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10),
        strip.text   = element_text(size = 12))

ggsave(file.path(outputs_dir, "ComparisonCategories_histograms.pdf"),
       hist_plot, width = 8, height = 10)

density_plot <- ggplot(all_results_melted,
                       aes(x = Value, fill = Category, color = Category)) +
  geom_density(alpha = 0.5, adjust = 1.2) +
  facet_grid(Test ~ ., scales = "free_y") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  labs(title = "Kernel density of metrics by test type",
       x = "Metric value", y = "Density",
       fill = "Comparison type", color = "Comparison type") +
  theme_minimal() +
  theme(plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10),
        strip.text   = element_text(size = 12))

ggsave(file.path(outputs_dir, "ComparisonCategories_densities.pdf"),
       density_plot, width = 8, height = 10)

# Summary means written to a single CSV for the manuscript
summary_means <- data.frame(
  Category = c("Within species", "Within genus", "Outside genus"),
  LM_R_squared_mean = c(
    mean(results_within_species$LM_R_squared,  na.rm = TRUE),
    mean(results_within_genus$LM_R_squared,    na.rm = TRUE),
    mean(results_outside_genus$LM_R_squared,   na.rm = TRUE)
  ),
  LM_p_mean = c(
    mean(results_within_species$LM_p,  na.rm = TRUE),
    mean(results_within_genus$LM_p,    na.rm = TRUE),
    mean(results_outside_genus$LM_p,   na.rm = TRUE)
  ),
  Spearman_R_mean = c(
    mean(results_within_species$Spearman_R, na.rm = TRUE),
    mean(results_within_genus$Spearman_R,   na.rm = TRUE),
    mean(results_outside_genus$Spearman_R,  na.rm = TRUE)
  ),
  Spearman_p_mean = c(
    mean(results_within_species$Spearman_p, na.rm = TRUE),
    mean(results_within_genus$Spearman_p,   na.rm = TRUE),
    mean(results_outside_genus$Spearman_p,  na.rm = TRUE)
  ),
  GLM_R_squared_mean = c(
    mean(results_within_species$GLM_R_squared, na.rm = TRUE),
    mean(results_within_genus$GLM_R_squared,   na.rm = TRUE),
    mean(results_outside_genus$GLM_R_squared,  na.rm = TRUE)
  ),
  GLM_p_mean = c(
    mean(results_within_species$GLM_p, na.rm = TRUE),
    mean(results_within_genus$GLM_p,   na.rm = TRUE),
    mean(results_outside_genus$GLM_p,  na.rm = TRUE)
  )
)
write.csv(summary_means,
          file.path(dir_tables, "ComparisonCategory_means.csv"),
          row.names = FALSE)


# ---- 16. Effect of gland age on correlation --------------------------------

# Glands were dissected in 2022; age = 2022 - collection year of the gland.
results_dated     <- results_table[!is.na(results_table$Date_Collected_Sample2), ]
results_dated$age <- 2022 - as.numeric(results_dated$Date_Collected_Sample2)

results_dated_melted <- reshape2::melt(
  results_dated,
  id.vars = c("Genus", "Species", "Date_Collected_Sample2", "age"),
  measure.vars  = c("Spearman_R", "LM_R_squared", "GLM_R_squared"),
  variable.name = "Test",
  value.name    = "Metric"
)

metric_colors <- c("Spearman_R"    = "#8362DD",
                   "GLM_R_squared" = "#BAEA75",
                   "LM_R_squared"  = "#F43151")

age_plot <- ggplot(results_dated_melted, aes(x = age, y = Metric, color = Test)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = metric_colors) +
  labs(title = "Gland age vs metrics",
       x = "Gland age (years)", y = "Metric value", color = "Test") +
  theme_minimal() +
  theme(plot.title  = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text   = element_text(size = 12),
        axis.title  = element_text(size = 14),
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10))

ggsave(file.path(outputs_dir, "GlandAge_vs_metrics.pdf"),
       age_plot, width = 8, height = 6)

# Simple LM tests (gland age does not predict any of the correlation metrics)
agevsR2       <- lm(LM_R_squared  ~ age, data = results_dated)
agevsspearman <- lm(Spearman_R    ~ age, data = results_dated)
agevsglm      <- lm(GLM_R_squared ~ age, data = results_dated)

# Mixed model with random intercept for species
agevsR2_mixed <- lmer(LM_R_squared ~ age + (1 | Species), data = results_dated)

sink(file.path(dir_tables, "GlandAge_models.txt"))
cat("=== LM: LM_R_squared ~ age ===\n")
print(summary(agevsR2))
cat("\n=== LM: Spearman_R ~ age ===\n")
print(summary(agevsspearman))
cat("\n=== LM: GLM_R_squared ~ age ===\n")
print(summary(agevsglm))
cat("\n=== Mixed: LM_R_squared ~ age + (1 | Species) ===\n")
print(summary(agevsR2_mixed))
sink()


# ---- 17. Phylogeny: import, impute missing tips, plot ---------------------

results_table2 <- results_table %>%
  left_join(
    results_dated %>% dplyr::select(Sample1, Sample2, age),
    by = c("Sample1", "Sample2")
  )

tree <- read.tree(file.path(base_dir,
                            "squamates_Title_Science2024_ultrametric_constrained.tre"))
tree$tip.label <- gsub("_", " ", tree$tip.label)

if (any(is.na(tree$edge.length))) {
  message("NA branch lengths detected. Assigning default value of 0.01.")
  tree$edge.length[is.na(tree$edge.length)] <- 0.01
}

species_in_data <- unique(results_table$Species)
absent_species  <- setdiff(species_in_data, tree$tip.label)

if (length(absent_species) > 0) {
  median_branch_length <- median(tree$edge.length, na.rm = TRUE)
  if (is.na(median_branch_length) || median_branch_length == 0) {
    median_branch_length <- 0.01
  }

  for (species in absent_species) {
    genus     <- strsplit(species, " ")[[1]][1]
    congeners <- grep(paste0("^", genus, " "), tree$tip.label, value = TRUE)

    if (length(congeners) > 1) {
      mrca_node     <- getMRCA(tree, congeners)
      mrca_depth    <- node.depth.edgelength(tree)[mrca_node]
      desired_depth <- max(node.depth.edgelength(tree))
      new_branch    <- max(desired_depth - mrca_depth, 1e-5)
      tree <- bind.tip(tree, tip.label = species,
                       where = mrca_node, edge.length = new_branch)

    } else if (length(congeners) == 1) {
      sister_species <- congeners[1]
      new_branch     <- median_branch_length / 2
      sister_node    <- which(tree$tip.label == sister_species)
      tree <- bind.tip(tree, tip.label = species, where = sister_node,
                       position = new_branch, edge.length = new_branch)
      parent_edge <- which(tree$edge[, 2] == sister_node)
      tree$edge.length[parent_edge] <- new_branch

    } else {
      tree <- bind.tip(tree, tip.label = species, where = "root",
                       edge.length = 0.01)
      message("No congeners found for ", species, " - added to root.")
    }
  }
}

common_species <- intersect(tree$tip.label, species_in_data)
trimmed_tree   <- drop.tip(tree, setdiff(tree$tip.label, common_species))

message("Trimmed tree is ultrametric: ", is.ultrametric(trimmed_tree))

pdf(file.path(outputs_dir, "PhylogeneticTree_imputed.pdf"),
    width = 8, height = 12)
plot(trimmed_tree, main = "Phylogenetic tree with imputed species")
dev.off()


# ---- 18. Per-species LM R^2 dot plot aligned with phylogeny ----------------

results_table2$Species <- factor(results_table2$Species,
                                 levels = trimmed_tree$tip.label)

dot_plot <- ggplot(results_table2,
                   aes(x = LM_R_squared, y = Species, color = age)) +
  geom_point(size = 3) +
  labs(title = "LM R^2 by species (phylogeny-aligned)",
       x = "Linear regression R^2",
       y = "Species (aligned with phylogeny)",
       color = "Age") +
  theme_minimal() +
  theme(plot.title  = element_text(hjust = 0.5),
        axis.text.y = element_text(face = "italic")) +
  scale_color_viridis_c(option = "cividis")

ggsave(file.path(outputs_dir, "DotPlot_LMR2_by_species.pdf"),
       dot_plot, width = 8, height = 12)


# ---- 19. Pseudechis australis paired-specimen experiment -------------------

pseudechis_samples <- c("iBAQ DEVC_145_Y", "iBAQ DEVC_146_X")
pseudechis_rows    <- summary_group3_no_others$Sample %in% pseudechis_samples

if (sum(pseudechis_rows) == 2) {
  toxin_idx_sg <- which(names(summary_group3_no_others) %in% toxin_cols)
  pseudechisexp <- as.data.frame(t(summary_group3_no_others[pseudechis_rows, toxin_idx_sg]))
  colnames(pseudechisexp) <- c("Venom", "Gland")
  pseudechisexp$Venom <- as.numeric(pseudechisexp$Venom)
  pseudechisexp$Gland <- as.numeric(pseudechisexp$Gland)
  pseudechisexp <- pseudechisexp[!(pseudechisexp$Venom == 0 & pseudechisexp$Gland == 0), ]

  lm_model     <- lm(Venom ~ Gland, data = pseudechisexp)
  r_squared    <- summary(lm_model)$adj.r.squared
  spearman_cor <- cor.test(pseudechisexp$Venom, pseudechisexp$Gland,
                           method = "spearman")

  glm_data <- pseudechisexp
  glm_data$Venom <- clamp_for_beta(glm_data$Venom)
  glm_data$Gland <- clamp_for_beta(glm_data$Gland)

  glm_r_squared <- NA_real_
  tryCatch({
    glm_model     <- betareg(Venom ~ Gland, data = glm_data)
    glm_r_squared <- summary(glm_model)$pseudo.r.squared
  }, error = function(e) message("Pseudechis GLM failed: ", e$message))

  pseudechis_plot <- ggplot(pseudechisexp,
                            aes(x = Gland, y = Venom, color = rownames(pseudechisexp))) +
    geom_point(size = 5) +
    geom_text_repel(aes(label = rownames(pseudechisexp)),
                    box.padding = 1, max.overlaps = Inf) +
    scale_x_sqrt() + scale_y_sqrt() +
    geom_smooth(inherit.aes = FALSE, method = "lm",
                aes(x = Gland, y = Venom)) +
    labs(title = "Pseudechis australis experiment (data from same specimen)",
         x = "Gland", y = "Venom") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5, face = "italic"),
          legend.position = "none") +
    scale_color_manual(values = color_mapping) +
    annotate("text",
             x = min(pseudechisexp$Gland) * 1.4,
             y = max(pseudechisexp$Venom) * 0.9,
             label = paste("LM R^2 = ", round(r_squared, 2)),
             hjust = 0, vjust = 1, size = 5, fontface = "italic") +
    annotate("text",
             x = min(pseudechisexp$Gland) * 1.4,
             y = max(pseudechisexp$Venom) * 0.75,
             label = paste("Spearman R = ", round(spearman_cor$estimate, 2)),
             hjust = 0, vjust = 1, size = 5, fontface = "italic") +
    annotate("text",
             x = min(pseudechisexp$Gland) * 1.4,
             y = max(pseudechisexp$Venom) * 0.6,
             label = paste("GLM R^2 = ", round(glm_r_squared, 2)),
             hjust = 0, vjust = 1, size = 5, fontface = "italic")

  ggsave(file.path(dir_regression, "Pseudechis_experiment_plot.pdf"),
         pseudechis_plot, width = 8, height = 6)
} else {
  message("Pseudechis paired samples not found in data; skipping section 19.")
}


# ---- 20. PCA of all species' mean toxin profiles ---------------------------

# Run un-scaled PCA (the analysis used in the paper). Re-run with scale.=TRUE
# if you want to weight all toxin axes equally.
pca_result <- prcomp(toxin_matrix, center = TRUE, scale. = FALSE)

pca_scores <- as.data.frame(pca_result$x[, 1:2])
pca_scores$Genus   <- species_avg$Genus
pca_scores$Family  <- species_avg$Family
pca_scores$Species <- species_avg$Species

genus_levels <- sort(unique(pca_scores$Genus))
genus_colors <- scico(n = length(genus_levels), palette = "batlow")
names(genus_colors) <- genus_levels

pca_plot <- ggplot(pca_scores, aes(x = PC1, y = PC2,
                                   colour = Family, shape = Family)) +
  geom_point(size = 4, alpha = 0.9) +
  theme_minimal(base_size = 14) +
  ggtitle("PCA of mean toxin profiles by species") +
  xlab(paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1),
              "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1),
              "% variance)")) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

ggsave(file.path(outputs_dir, "PCA_meanToxinProfiles.pdf"),
       pca_plot, width = 8, height = 6)

# Interactive PCA (HTML widget)
interactive_pca <- plot_ly(
  data = pca_scores,
  x = ~PC1, y = ~PC2,
  type = "scatter", mode = "markers",
  color = ~Genus, colors = genus_colors,
  symbol = ~Family,
  marker = list(size = 10, opacity = 0.8),
  text = ~paste("<b>Species:</b>", Species,
                "<br><b>Genus:</b>",  Genus,
                "<br><b>Family:</b>", Family),
  hoverinfo = "text"
) %>%
  layout(
    title = "Interactive PCA of mean toxin profiles by species",
    xaxis = list(title = paste0("PC1 (",
                  round(summary(pca_result)$importance[2, 1] * 100, 1),
                  "% variance)")),
    yaxis = list(title = paste0("PC2 (",
                  round(summary(pca_result)$importance[2, 2] * 100, 1),
                  "% variance)")),
    legend = list(orientation = "v")
  )

if (requireNamespace("htmlwidgets", quietly = TRUE)) {
  # selfcontained=TRUE requires pandoc; fall back to a folder-based widget
  # if pandoc is unavailable.
  has_pandoc <- requireNamespace("rmarkdown", quietly = TRUE) &&
                rmarkdown::pandoc_available()
  tryCatch(
    htmlwidgets::saveWidget(
      interactive_pca,
      file.path(outputs_dir, "PCA_interactive.html"),
      selfcontained = isTRUE(has_pandoc)
    ),
    error = function(e) message("Could not save interactive PCA: ", e$message)
  )
} else {
  message("Install 'htmlwidgets' to save the interactive PCA as HTML.")
}


