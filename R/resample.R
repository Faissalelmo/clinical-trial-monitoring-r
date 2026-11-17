# Patient-resampling reference for the AE rate.
#
# Adverse events cluster within patients: one patient with twelve AEs and
# another with none is ordinary, so the Poisson assumption (every event
# independent) makes the funnel limits far too narrow. Widening them by the
# spread between sites (overdispersion) fixes the false alarms but also
# absorbs the very site differences monitoring is meant to find.
#
# Here a site of n patients is compared with many "sites" of n patients drawn
# at random from the other sites' patients, each carrying their own exposure
# and AE count. That keeps the clustering inside patients and nothing else.

suppressPackageStartupMessages(library(dplyr))

#' One row per treated patient: site, years on treatment, TEAE count.
patient_table <- function(adsl, adae) {
  teae <- adae %>% filter(TRTEMFL == "Y") %>% count(USUBJID, name = "teae")
  adsl %>%
    filter(SAFFL == "Y", !is.na(TRTDURD), TRTDURD > 0) %>%
    transmute(USUBJID, SITEID, years = TRTDURD / 365.25) %>%
    left_join(teae, by = "USUBJID") %>%
    mutate(teae = coalesce(teae, 0L))
}

#' Reference distribution of the AE rate for a site of n patients, drawn
#' from the patients of every other site.
reference_rates <- function(patients, site, n, draws = 4000, seed = 1) {
  pool <- patients[patients$SITEID != site, ]
  set.seed(seed)
  idx <- matrix(sample.int(nrow(pool), n * draws, replace = TRUE), nrow = draws)
  rowSums(matrix(pool$teae[idx], nrow = draws)) /
    rowSums(matrix(pool$years[idx], nrow = draws))
}

#' Flag a site's AE rate against its resampled reference.
resample_flag <- function(observed, reference) {
  q <- quantile(reference, c(0.001, 0.025, 0.975, 0.999), names = FALSE)
  case_when(
    observed < q[1] ~ "low",
    observed > q[4] ~ "high",
    observed < q[2] | observed > q[3] ~ "watch",
    TRUE ~ "ok"
  )
}

#' Judge every site's AE rate by patient resampling.
resample_sites <- function(patients, draws = 4000) {
  sites <- sort(unique(patients$SITEID))
  bind_rows(lapply(seq_along(sites), function(i) {
    s <- sites[i]
    mine <- patients[patients$SITEID == s, ]
    ref <- reference_rates(patients, s, nrow(mine), draws, seed = i)
    obs <- sum(mine$teae) / sum(mine$years)
    tibble(SITEID = s, patients = nrow(mine), value = obs,
           ref_median = median(ref),
           ref_lo = quantile(ref, 0.001, names = FALSE),
           ref_hi = quantile(ref, 0.999, names = FALSE),
           percentile = mean(ref <= obs),
           flag = resample_flag(obs, ref))
  }))
}
