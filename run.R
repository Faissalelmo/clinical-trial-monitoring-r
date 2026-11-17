# Centralized statistical monitoring of the CDISC pilot study.
# Usage: Rscript run.R      (writes output/)

suppressPackageStartupMessages({
  library(pharmaversesdtm)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(jsonlite)
})
for (f in list.files("R", full.names = TRUE)) source(f)
dir.create("output", showWarnings = FALSE)

dm <- pharmaversesdtm::dm
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds
ae <- pharmaversesdtm::ae
lb <- pharmaversesdtm::lb

# 1. Analysis datasets ------------------------------------------------------
adsl <- build_adsl(dm, ex, ds)
adae <- build_adae(ae, adsl)
patients <- patient_table(adsl, adae)
write_csv(adsl, "output/adsl.csv")
write_csv(adae, "output/adae.csv")

# 2. Key risk indicators per site ------------------------------------------
counts <- site_counts(adsl, adae, lb)
kris <- assess_sites(counts)
teae_resampled <- resample_sites(patients)
# The AE rate is judged by patient resampling (see R/resample.R); the two
# proportions keep their exact binomial limits.
kris <- kris %>%
  left_join(select(teae_resampled, SITEID, rs_flag = flag, percentile), by = "SITEID") %>%
  mutate(flag = if_else(kri == "TEAE rate per patient-year", rs_flag, flag),
         z = if_else(kri == "TEAE rate per patient-year",
                     qnorm(pmin(pmax(percentile, 0.0005), 0.9995)), z)) %>%
  select(-rs_flag, -percentile)
write_csv(counts, "output/site_counts.csv")
write_csv(kris, "output/site_kris.csv")
write_csv(teae_resampled, "output/site_teae_resampled.csv")

# 3. Figures ----------------------------------------------------------------
flag_colours <- c(ok = "#7a8794", watch = "#e0a100", low = "#c0392b", high = "#c0392b",
                  `no data` = "#cccccc")
lt <- scale_linetype_manual(values = c(watch = "dashed", alert = "dotted"),
                            labels = c(alert = "99.8% limit", watch = "95% limit"))

# AE rate: limits from resampling n patients out of the whole study.
n_max <- max(teae_resampled$patients)
curve_rs <- bind_rows(lapply(2:(n_max + 5), function(n) {
  ref <- reference_rates(patients, site = "", n = n, draws = 4000, seed = n)
  q <- quantile(ref, c(0.001, 0.025, 0.975, 0.999), names = FALSE)
  tibble(denom = n, level = c("alert", "watch"), lo = q[1:2], hi = q[4:3])
}))
overall_rate <- sum(patients$teae) / sum(patients$years)
p <- ggplot() +
  geom_line(data = curve_rs, aes(denom, lo, linetype = level), colour = "#34495e") +
  geom_line(data = curve_rs, aes(denom, hi, linetype = level), colour = "#34495e") +
  geom_hline(yintercept = overall_rate, colour = "#34495e") +
  geom_point(data = teae_resampled, aes(patients, value, colour = flag), size = 3) +
  geom_text(data = teae_resampled, aes(patients, value, label = SITEID), size = 3, vjust = -1) +
  scale_colour_manual(values = flag_colours) + lt +
  coord_cartesian(ylim = c(0, 60)) +
  labs(title = "TEAE rate per patient-year by site",
       subtitle = "Limits from resampling patients (4,000 draws per site size)",
       x = "Patients treated at the site", y = "TEAEs per patient-year",
       colour = "Site", linetype = NULL) +
  theme_minimal(base_size = 11)
ggsave("output/funnel_teae.png", p, width = 8, height = 5, dpi = 150, bg = "white")

funnel_plot <- function(k, xlab, ylab, file) {
  d <- filter(kris, kri == k)
  curve <- funnel_curve(d$overall[1], max(d$denom), "proportion", phi = d$phi[1])
  p <- ggplot() +
    geom_line(data = curve, aes(denom, lo, linetype = level), colour = "#34495e") +
    geom_line(data = curve, aes(denom, hi, linetype = level), colour = "#34495e") +
    geom_hline(yintercept = d$overall[1], colour = "#34495e") +
    geom_point(data = d, aes(denom, value, colour = flag), size = 3) +
    geom_text(data = d, aes(denom, value, label = SITEID), size = 3, vjust = -1) +
    scale_colour_manual(values = flag_colours) + lt +
    scale_y_continuous(labels = scales::percent) +
    labs(title = k, x = xlab, y = ylab, colour = "Site", linetype = NULL) +
    theme_minimal(base_size = 11)
  ggsave(file, p, width = 8, height = 5, dpi = 150, bg = "white")
}
funnel_plot("Early discontinuation", "Patients treated", "Share discontinued",
            "output/funnel_discontinuation.png")
funnel_plot("Missing lab results", "Lab results expected", "Share missing",
            "output/funnel_missing_labs.png")

heat <- ggplot(kris, aes(kri, SITEID, fill = z)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(flag %in% c("low", "high", "watch"), flag, "")), size = 3) +
  scale_fill_gradient2(low = "#2e86c1", mid = "white", high = "#c0392b", midpoint = 0,
                       name = "Standardised\ndifference") +
  labs(title = "Site risk overview", x = NULL, y = "Site") +
  theme_minimal(base_size = 11)
ggsave("output/site_heatmap.png", heat, width = 7, height = 6, dpi = 150, bg = "white")

# 4. Does it catch a site that under-reports? ------------------------------
studies <- lapply(METHODS, function(m) detection_study(adsl, adae, lb, method = m))
names(studies) <- METHODS
write_csv(bind_rows(lapply(studies, function(s) mutate(s$runs, method = s$method))),
          "output/underreporting_runs.csv")

det <- bind_rows(lapply(studies, function(s) {
  s$summary %>% transmute(method = s$method, share, Alert = alert, `Alert or watch` = alert_or_watch)
})) %>%
  pivot_longer(c(Alert, `Alert or watch`), names_to = "level", values_to = "rate") %>%
  mutate(method = factor(method, METHODS,
                         c("Poisson (exact)", "Overdispersion-adjusted", "Patient resampling")))
p <- ggplot(det, aes(share, rate, colour = method)) +
  geom_line() + geom_point(size = 2) +
  facet_wrap(~level) +
  scale_x_continuous(labels = scales::percent, breaks = c(0.25, 0.5, 0.75)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(title = "Planted AE under-reporting: share of runs where the site is flagged",
       x = "Share of the site's TEAEs removed", y = "Runs flagged", colour = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom")
ggsave("output/underreporting_detection.png", p, width = 8, height = 4.5, dpi = 150, bg = "white")

# 5. Report -----------------------------------------------------------------
report <- list(
  data = list(patients_treated = sum(adsl$SAFFL == "Y"),
              sites = n_distinct(adsl$SITEID),
              adverse_events = nrow(adae),
              treatment_emergent = sum(adae$TRTEMFL == "Y", na.rm = TRUE),
              patient_years = round(sum(patients$years), 1),
              lab_results = sum(counts$lab_results)),
  study_rates = kris %>% distinct(kri, overall) %>% mutate(overall = round(overall, 4)),
  flags = kris %>% filter(flag != "ok") %>%
    transmute(site = SITEID, kri, flag, value = round(value, 4), study = round(overall, 4)),
  methods = lapply(studies, function(s) list(
    false_alarms_on_clean_data = s$false_alarms,
    watch_on_clean_data = s$watch_on_clean,
    sites_planted = s$sites,
    detection = s$summary
  )),
  overdispersion_phi = round(kris$phi[kris$kri == "TEAE rate per patient-year"][1], 2),
  draws_per_site_and_share = 20
)
write_json(report, "output/report.json", auto_unbox = TRUE, pretty = TRUE)
cat(toJSON(report, auto_unbox = TRUE, pretty = TRUE), "\n")
