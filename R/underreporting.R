# Does the method catch a site that under-reports adverse events?
#
# The pilot data carry no known misconduct, so the test is built: remove a
# share of one site's treatment-emergent AEs, re-run the assessment, and see
# whether that site is flagged "low". Repeated for every site large enough to
# matter, several shares and several random draws. The unmodified data give
# the false-alarm count.

suppressPackageStartupMessages(library(dplyr))

METHODS <- c("exact", "adjusted", "resample")

#' Drop `share` of one site's TEAEs at random.
plant_underreporting <- function(adae, site, share, seed) {
  set.seed(seed)
  target <- which(adae$SITEID == site & adae$TRTEMFL == "Y")
  drop <- target[runif(length(target)) < share]
  if (length(drop)) adae[-drop, ] else adae
}

#' Every site's AE-rate flag under one method.
teae_flags <- function(adsl, adae, lb, method = METHODS) {
  method <- match.arg(method)
  if (method == "resample") {
    r <- resample_sites(patient_table(adsl, adae))
    return(setNames(r$flag, r$SITEID))
  }
  k <- assess_sites(site_counts(adsl, adae, lb), adjust = method == "adjusted") %>%
    filter(kri == "TEAE rate per patient-year")
  setNames(k$flag, k$SITEID)
}

#' Detection table for one method. Only sites the method calls normal on the
#' untouched data are planted: a site already flagged would count as caught
#' without the test doing anything.
detection_study <- function(adsl, adae, lb, method = METHODS,
                            shares = c(0.25, 0.5, 0.75), min_patients = 10, draws = 20) {
  method <- match.arg(method)
  patients <- patient_table(adsl, adae)
  size <- table(patients$SITEID)
  base <- teae_flags(adsl, adae, lb, method)
  sites <- names(size)[size >= min_patients]
  sites <- sites[base[sites] %in% c("ok", "high")]

  # For resampling, a site's reference is drawn from the other sites only,
  # so it does not move when that site's own AEs are removed: build it once.
  refs <- if (method == "resample") {
    setNames(lapply(sites, function(s) {
      reference_rates(patients, s, size[[s]], seed = match(s, sort(unique(patients$SITEID))))
    }), sites)
  }

  runs <- expand.grid(site = sites, share = shares, draw = seq_len(draws),
                      stringsAsFactors = FALSE)
  runs$flag <- vapply(seq_len(nrow(runs)), function(i) {
    r <- runs[i, ]
    planted <- plant_underreporting(adae, r$site, r$share, seed = 1000 * i)
    if (method == "resample") {
      mine <- patient_table(adsl, planted) %>% filter(SITEID == r$site)
      resample_flag(sum(mine$teae) / sum(mine$years), refs[[r$site]])
    } else {
      teae_flags(adsl, planted, lb, method)[[r$site]]
    }
  }, character(1))

  list(
    method = method,
    false_alarms = sum(base %in% c("low", "high")),
    watch_on_clean = sum(base == "watch"),
    sites = sites,
    runs = runs,
    summary = runs %>%
      group_by(share) %>%
      summarise(runs = n(),
                alert = round(mean(flag == "low"), 3),
                alert_or_watch = round(mean(flag %in% c("low", "watch")), 3),
                .groups = "drop"),
    by_site = runs %>%
      group_by(site, share) %>%
      summarise(alert = mean(flag == "low"), .groups = "drop")
  )
}
