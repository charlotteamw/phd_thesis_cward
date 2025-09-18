# Load libraries
library(devtools)
library(MixSIAR)
library(rjags)
library(R2jags)
library(tidyverse)

graphics.off()

# File paths
mix_file <- "//shiner_mix_liver.csv"
source_file <- "//source_data.csv"
discr_file  <- "//TEF_2source.csv"


# ------------------------------------------------------------------------------
# MODEL 1: NULL
# ------------------------------------------------------------------------------
mix1 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = NULL,
                      fac_random = NULL,
                      fac_nested = NULL,
                      cont_effects = NULL)

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix1)
discr <- load_discr_data(filename = discr_file, mix1)

plot_data("isospace_model1", FALSE, FALSE, mix1, source, discr)
plot_prior(alpha.prior = c(1), source)

write_JAGS_model("MixSIAR_model_1.txt", TRUE, TRUE, mix1, source)
jags.1 <- run_model(run = "test", mix1, source, discr, "MixSIAR_model_1.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options <- list(summary_save = TRUE, summary_name = "summary_model1", 
                       sup_post = TRUE, plot_post_save_pdf = FALSE, plot_post_name = "posterior_model1",
                       sup_pairs = TRUE, plot_pairs_save_pdf = TRUE, plot_pairs_name = "pairs_model1",
                       sup_xy = TRUE, plot_xy_save_pdf = TRUE, plot_xy_name = "xy_model1",
                       gelman = TRUE, heidel = FALSE, geweke = TRUE, diag_save = TRUE, diag_name = "diagnostics_model1",
                       indiv_effect = FALSE, return_obj = TRUE)

output_posteriors(jags.1, mix1, source, output_options)
output_diagnostics(jags.1, mix1, source, output_options)
output_stats(jags.1, mix1, source, output_options)

# ------------------------------------------------------------------------------
# MODEL 2: month + location
# ------------------------------------------------------------------------------
mix2 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("month", "location"),
                      fac_random = c(FALSE, FALSE),
                      fac_nested = c(FALSE, FALSE),
                      cont_effects = NULL)

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix2)
discr <- load_discr_data(filename = discr_file, mix2)


write_JAGS_model("MixSIAR_model_2.txt", TRUE, TRUE, mix2, source)
jags.2 <- run_model(run = "test", mix2, source, discr, "MixSIAR_model_2.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$plot_post_name <- "posterior_model2"
output_options$plot_pairs_name <- "pairs_model2"
output_options$plot_xy_name <- "xy_model2"
output_options$diag_name <- "diagnostics_model2"

output_posteriors(jags.2, mix2, source, output_options)
output_diagnostics(jags.2, mix2, source, output_options)
output_stats(jags.2, mix2, source, output_options)

graphics.off()

# ------------------------------------------------------------------------------
# MODEL 3: month + location + tl
# ------------------------------------------------------------------------------
mix3 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("month", "location"),
                      fac_random = c(FALSE, FALSE),
                      fac_nested = c(FALSE, FALSE),
                      cont_effects = "tl")

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix3)
discr <- load_discr_data(filename = discr_file, mix3)


plot_data("isospace_model3", FALSE, FALSE, mix3, source, discr)

write_JAGS_model("MixSIAR_model_3.txt", TRUE, TRUE, mix3, source)
jags.3 <- run_model(run = "test", mix3, source, discr, "MixSIAR_model_3.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model3"
output_options$plot_post_name <- "posterior_model3"
output_options$plot_pairs_name <- "pairs_model3"
output_options$plot_xy_name <- "xy_model3"
output_options$diag_name <- "diagnostics_model3"

output_posteriors(jags.3, mix3, source, output_options)
output_diagnostics(jags.3, mix3, source, output_options)
output_stats(jags.3, mix3, source, output_options)

# ------------------------------------------------------------------------------
# MODEL 4: location + tl
# ------------------------------------------------------------------------------
mix4 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("location"),
                      fac_random = c(FALSE),
                      fac_nested = c(FALSE),
                      cont_effects = "tl")

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix4)
discr <- load_discr_data(filename = discr_file, mix4)


write_JAGS_model("MixSIAR_model_4.txt", TRUE, TRUE, mix4, source)
jags.4 <- run_model(run = "test", mix4, source, discr, "MixSIAR_model_4.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model4"
output_options$plot_post_name <- "posterior_model4"
output_options$plot_pairs_name <- "pairs_model4"
output_options$plot_xy_name <- "xy_model4"
output_options$diag_name <- "diagnostics_model4"

output_posteriors(jags.4, mix4, source, output_options)
output_diagnostics(jags.4, mix4, source, output_options)
output_stats(jags.4, mix4, source, output_options)

# ------------------------------------------------------------------------------
# MODEL 5: month + tl
# ------------------------------------------------------------------------------
mix5 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("month"),
                      fac_random = c(FALSE),
                      fac_nested = c(FALSE),
                      cont_effects = "tl")

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix5)
discr <- load_discr_data(filename = discr_file, mix5)



write_JAGS_model("MixSIAR_model_5.txt", TRUE, TRUE, mix5, source)
jags.5 <- run_model(run = "test", mix5, source, discr, "MixSIAR_model_5.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

noutput_options$summary_name <- "summary_model5"
output_options$plot_post_name <- "posterior_model5"
output_options$plot_pairs_name <- "pairs_model5"
output_options$plot_xy_name <- "xy_model5"
output_options$diag_name <- "diagnostics_model5"

output_posteriors(jags.5, mix5, source, output_options)
output_diagnostics(jags.5, mix5, source, output_options)
output_stats(jags.5, mix5, source, output_options)
graphics.off()

# ------------------------------------------------------------------------------
# MODEL 6: tl
# ------------------------------------------------------------------------------
mix6 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = NULL,
                      fac_random = NULL,
                      fac_nested = NULL,
                      cont_effects = "tl")

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix6)
discr <- load_discr_data(filename = discr_file, mix6)


write_JAGS_model("MixSIAR_model_6.txt", TRUE, TRUE, mix6, source)
jags.6 <- run_model(run = "test", mix6, source, discr, "MixSIAR_model_6.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model6"
output_options$plot_post_name <- "posterior_model6"
output_options$plot_pairs_name <- "pairs_model6"
output_options$plot_xy_name <- "xy_model6"
output_options$diag_name <- "diagnostics_model6"

output_posteriors(jags.6, mix6, source, output_options)
output_diagnostics(jags.6, mix6, source, output_options)
output_stats(jags.6, mix6, source, output_options)

# ------------------------------------------------------------------------------
# MODEL 7: year + location + tl
# ------------------------------------------------------------------------------
mix7 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("year", "location"),
                      fac_random = c(TRUE, FALSE),
                      fac_nested = c(FALSE, FALSE),
                      cont_effects = "tl")

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix7)
discr <- load_discr_data(filename = discr_file, mix7)


write_JAGS_model("MixSIAR_model_7.txt", TRUE, TRUE, mix7, source)
jags.7 <- run_model(run = "test", mix7, source, discr, "MixSIAR_model_7.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model7"
output_options$plot_post_name <- "posterior_model7"
output_options$plot_pairs_name <- "pairs_model7"
output_options$plot_xy_name <- "xy_model7"
output_options$diag_name <- "diagnostics_model7"

output_posteriors(jags.7, mix7, source, output_options)
output_diagnostics(jags.7, mix7, source, output_options)
output_stats(jags.7, mix7, source, output_options)


# ------------------------------------------------------------------------------
# MODEL 8: year + tl
# ------------------------------------------------------------------------------
mix8 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("year"),
                      fac_random = c(TRUE),
                      fac_nested = c(FALSE),
                      cont_effects = "tl")


source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix8)
discr <- load_discr_data(filename = discr_file, mix8)

write_JAGS_model("MixSIAR_model_8.txt", TRUE, TRUE, mix8, source)
jags.8 <- run_model(run = "test", mix8, source, discr, "MixSIAR_model_8.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model8"
output_options$plot_post_name <- "posterior_model8"
output_options$plot_pairs_name <- "pairs_model8"
output_options$plot_xy_name <- "xy_model8"
output_options$diag_name <- "diagnostics_model8"

output_posteriors(jags.8, mix8, source, output_options)
output_diagnostics(jags.8, mix8, source, output_options)
output_stats(jags.8, mix8, source, output_options)

# ------------------------------------------------------------------------------
# MODEL 9: year + month + tl
# ------------------------------------------------------------------------------
mix9 <- load_mix_data(filename = mix_file,
                      iso_names = c("d13C","d15N"),
                      factors = c("year", "month"),
                      fac_random = c(TRUE, FALSE),
                      fac_nested = c(FALSE, TRUE),
                      cont_effects = "tl")

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix9)
discr <- load_discr_data(filename = discr_file, mix9)

write_JAGS_model("MixSIAR_model_9.txt", TRUE, TRUE, mix9, source)
jags.9 <- run_model(run = "test", mix9, source, discr, "MixSIAR_model_9.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model9"
output_options$plot_post_name <- "posterior_model9"
output_options$plot_pairs_name <- "pairs_model9"
output_options$plot_xy_name <- "xy_model9"
output_options$diag_name <- "diagnostics_model9"

output_posteriors(jags.9, mix9, source, output_options)
output_diagnostics(jags.9, mix9, source, output_options)
output_stats(jags.9, mix9, source, output_options)


# ------------------------------------------------------------------------------
# MODEL 10: month
# ------------------------------------------------------------------------------
mix10 <- load_mix_data(filename = mix_file,
                       iso_names = c("d13C","d15N"),
                       factors = c("month"),
                       fac_random = c(FALSE),
                       fac_nested = c(FALSE),
                       cont_effects = NULL)

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix10)
discr <- load_discr_data(filename = discr_file, mix10)



write_JAGS_model("MixSIAR_model_10.txt", TRUE, TRUE, mix10, source)
jags.10 <- run_model(run = "test", mix10, source, discr, "MixSIAR_model_10.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model10"
output_options$plot_post_name <- "posterior_model10"
output_options$plot_pairs_name <- "pairs_model10"
output_options$plot_xy_name <- "xy_model10"
output_options$diag_name <- "diagnostics_model10"

output_posteriors(jags.10, mix10, source, output_options)
output_diagnostics(jags.10, mix10, source, output_options)
output_stats(jags.10, mix10, source, output_options)


# ------------------------------------------------------------------------------
# MODEL 11: location
# ------------------------------------------------------------------------------
mix11 <- load_mix_data(filename = mix_file,
                       iso_names = c("d13C","d15N"),
                       factors = c("location"),
                       fac_random = c(FALSE),
                       fac_nested = c(FALSE),
                       cont_effects = NULL)

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix11)
discr <- load_discr_data(filename = discr_file, mix11)

plot_data("isospace_model11", FALSE, FALSE, mix11, source, discr)
plot_prior(alpha.prior = c(1), source)

write_JAGS_model("MixSIAR_model_11.txt", TRUE, TRUE, mix11, source)
jags.11 <- run_model(run = "test", mix11, source, discr, "MixSIAR_model_11.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model11"
output_options$plot_post_name <- "posterior_model11"
output_options$plot_pairs_name <- "pairs_model11"
output_options$plot_xy_name <- "xy_model11"
output_options$diag_name <- "diagnostics_model11"

output_posteriors(jags.11, mix11, source, output_options)
output_diagnostics(jags.11, mix11, source, output_options)
output_stats(jags.11, mix11, source, output_options)


# ------------------------------------------------------------------------------
# MODEL 12: year
# ------------------------------------------------------------------------------
mix12 <- load_mix_data(filename = mix_file,
                       iso_names = c("d13C","d15N"),
                       factors = c("year"),
                       fac_random = c(TRUE),
                       fac_nested = c(FALSE),
                       cont_effects = NULL)

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix12)
discr <- load_discr_data(filename = discr_file, mix12)

plot_data("isospace_model12", FALSE, FALSE, mix12, source, discr)
plot_prior(alpha.prior = c(1), source)

write_JAGS_model("MixSIAR_model_12.txt", TRUE, TRUE, mix12, source)
jags.12 <- run_model(run = "test", mix12, source, discr, "MixSIAR_model_12.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model12"
output_options$plot_post_name <- "posterior_model12"
output_options$plot_pairs_name <- "pairs_model12"
output_options$plot_xy_name <- "xy_model12"
output_options$diag_name <- "diagnostics_model12"

output_posteriors(jags.12, mix12, source, output_options)
output_diagnostics(jags.12, mix12, source, output_options)
output_stats(jags.12, mix12, source, output_options)


# ------------------------------------------------------------------------------
# MODEL 13: year + month nested
# ------------------------------------------------------------------------------
mix13 <- load_mix_data(filename = mix_file,
                       iso_names = c("d13C","d15N"),
                       factors = c("year", "month"),
                       fac_random = c(TRUE, FALSE),
                       fac_nested = c(FALSE, TRUE),
                       cont_effects = NULL)

source <- load_source_data(filename = source_file, source_factors = NULL, conc_dep = FALSE, data_type = "raw", mix13)
discr <- load_discr_data(filename = discr_file, mix13)

plot_data("isospace_model13", FALSE, FALSE, mix13, source, discr)
plot_prior(alpha.prior = c(1), source)

write_JAGS_model("MixSIAR_model_13.txt", TRUE, TRUE, mix13, source)
jags.13 <- run_model(run = "test", mix13, source, discr, "MixSIAR_model_13.txt", alpha.prior = 1, resid_err = TRUE, process_err = TRUE)

output_options$summary_name <- "summary_model13"
output_options$plot_post_name <- "posterior_model13"
output_options$plot_pairs_name <- "pairs_model13"
output_options$plot_xy_name <- "xy_model13"
output_options$diag_name <- "diagnostics_model13"

output_posteriors(jags.13, mix13, source, output_options)
output_diagnostics(jags.13, mix13, source, output_options)
output_stats(jags.13, mix13, source, output_options)

# ------------------------------------------------------------------------------
# Compare Models 1-13
# ------------------------------------------------------------------------------
comparison.table <- compare_models(list(
  model1 = jags.1,
  model2 = jags.2,
  model3 = jags.3,
  model4 = jags.4,
  model5 = jags.5,
  model6 = jags.6,
  model7 = jags.7,
  model8 = jags.8,
  model9 = jags.9,
  model10 = jags.10,
  model11 = jags.11, 
  model12 = jags.12,
  model13 = jags.13
))

xi.C <- sapply(list(jags.1, jags.2, jags.3, jags.4, jags.5, jags.6, jags.7, jags.8, jags.9, jags.10, jags.11, jags.12, jags.13), function(m) round(median(m$BUGSoutput$sims.list$resid.prop[, 1]), 2))
xi.N <- sapply(list(jags.1, jags.2, jags.3, jags.4, jags.5, jags.6, jags.7, jags.8, jags.9, jags.10, jags.11, jags.12, jags.13), function(m) round(median(m$BUGSoutput$sims.list$resid.prop[, 2]), 2))


min_looic <- min(comparison.table$LOOic)
se_min <- comparison.table$se_LOOic[which.min(comparison.table$LOOic)]

comparison.table <- comparison.table %>%
  mutate(`ΔLOOic` = LOOic - min_looic,
         `SE (ΔLOOic)` = sqrt((LOOic - min_looic)^2 + se_min^2),
         `ξC` = xi.C,
         `ξN` = xi.N) %>%
  select(Model, LOOic, se_LOOic, `ΔLOOic`, `SE (ΔLOOic)`, weight, `ξC`, `ξN`) %>%
  rename(`SE (LOOic)` = se_LOOic, wi = weight)

print(comparison.table)


write_csv(comparison.table, "/Users/charlotteward/Documents/algonquin_minnow/MixSIAR/shiner_model_comparison_table.csv")
