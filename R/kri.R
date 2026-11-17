# Site-level key risk indicators (KRIs) and funnel-plot control limits for
# centralized statistical monitoring (ICH E6(R2) 5.0, risk-based monitoring).
#
# Each KRI is a count against a denominator. A site is judged against the
# study-wide rate with exact control limits that widen as the site gets
# smaller, so a 3-patient site is not flagged for noise and a 50-patient site
# is not excused for a real difference.

suppressPackageStartupMessages(library(dplyr))

LIMITS <- c(watch = 0.95, alert = 0.998)

#' One row per site: exposure, TEAEs, discontinuations, lab results.
site_counts <- function(adsl, adae, lb) {
  saf <- filter(adsl, SAFFL == "Y")
  exposure <- saf %>%
    group_by(SITEID) %>%
    summarise(subjects = n(),
              patient_years = sum(TRTDURD, na.rm = TRUE) / 365.25,
              discontinued = sum(EOSSTT == "DISCONTINUED"),
              .groups = "drop")
  teae <- adae %>%
    filter(TRTEMFL == "Y") %>%
    count(SITEID, name = "teae")
  labs <- lb %>%
    semi_join(saf, by = "USUBJID") %>%
    mutate(SITEID = substr(USUBJID, 4, 6)) %>%
    group_by(SITEID) %>%
    summarise(lab_results = n(), lab_missing = sum(is.na(LBSTRESN)), .groups = "drop")

  exposure %>%
    left_join(teae, by = "SITEID") %>%
    left_join(labs, by = "SITEID") %>%
    mutate(across(c(teae, lab_results, lab_missing), ~ coalesce(.x, 0L)))
}

#' Exact Poisson limits for an event count given expected events.
poisson_limits <- function(expected, level) {
  a <- (1 - level) / 2
  cbind(lo = qpois(a, expected), hi = qpois(1 - a, expected))
}

#' Exact binomial limits for a count given n trials and the overall rate.
binomial_limits <- function(n, p, level) {
  a <- (1 - level) / 2
  cbind(lo = qbinom(a, n, p), hi = qbinom(1 - a, n, p))
}

#' Variance of one unit of the denominator under the null model.
unit_variance <- function(overall, type) {
  if (type == "rate") overall else overall * (1 - overall)
}

#' Overdispersion factor (Spiegelhalter 2005): the mean squared z-score of the
#' sites, winsorised at the 10th/90th percentiles so one extreme site cannot
#' widen the limits that should catch it. Never below 1.
overdispersion <- function(z) {
  z <- z[is.finite(z)]
  if (length(z) < 5) return(1)
  q <- quantile(z, c(0.1, 0.9), names = FALSE)
  max(1, mean(pmin(pmax(z, q[1]), q[2])^2))
}

#' Judge every site on one KRI.
#' type "rate": count per exposure (Poisson); "proportion": count of n (binomial).
#' adjust = FALSE uses exact Poisson/binomial limits; TRUE inflates the
#' variance by the overdispersion factor (normal approximation).
#' Flags: "low"/"high" outside the 99.8% limit, "watch" outside 95%, else "ok".
judge <- function(count, denom, type = c("rate", "proportion"), kri = "", adjust = TRUE) {
  type <- match.arg(type)
  overall <- sum(count) / sum(denom)
  expected <- overall * denom
  value <- ifelse(denom > 0, count / denom, NA_real_)
  z <- ifelse(denom > 0, (value - overall) / sqrt(unit_variance(overall, type) / denom), NA_real_)
  phi <- if (adjust) overdispersion(z) else 1

  # With no overdispersion the exact limits stand; the normal approximation
  # is only used to carry the inflation.
  if (phi > 1) {
    za <- z / sqrt(phi)
    k <- qnorm(1 - (1 - LIMITS) / 2)
    flag <- case_when(
      denom == 0 ~ "no data",
      za < -k[["alert"]] ~ "low",
      za > k[["alert"]] ~ "high",
      abs(za) > k[["watch"]] ~ "watch",
      TRUE ~ "ok"
    )
  } else {
    lims <- lapply(LIMITS, function(level) {
      if (type == "rate") poisson_limits(expected, level) else binomial_limits(denom, overall, level)
    })
    flag <- case_when(
      denom == 0 ~ "no data",
      count < lims$alert[, "lo"] ~ "low",
      count > lims$alert[, "hi"] ~ "high",
      count < lims$watch[, "lo"] | count > lims$watch[, "hi"] ~ "watch",
      TRUE ~ "ok"
    )
  }
  tibble(kri = kri, type = type, count = count, denom = denom, value = value,
         overall = overall, expected = expected, phi = phi,
         z = z / sqrt(phi), flag = flag)
}

#' All three KRIs, long format: one row per site and indicator.
assess_sites <- function(counts, adjust = TRUE) {
  one <- function(count, denom, type, kri) {
    bind_cols(SITEID = counts$SITEID, judge(count, denom, type, kri, adjust))
  }
  bind_rows(
    one(counts$teae, counts$patient_years, "rate", "TEAE rate per patient-year"),
    one(counts$discontinued, counts$subjects, "proportion", "Early discontinuation"),
    one(counts$lab_missing, counts$lab_results, "proportion", "Missing lab results")
  )
}

#' Funnel curve for plotting: the limits as a function of the denominator,
#' drawn the way the flags are decided (overdispersion-adjusted).
funnel_curve <- function(overall, max_denom, type = c("rate", "proportion"), phi = 1) {
  type <- match.arg(type)
  d <- seq(max_denom / 400, max_denom * 1.05, length.out = 400)
  sd <- sqrt(phi * unit_variance(overall, type) / d)
  bind_rows(lapply(names(LIMITS), function(nm) {
    k <- qnorm(1 - (1 - LIMITS[[nm]]) / 2)
    lo <- pmax(0, overall - k * sd)
    hi <- overall + k * sd
    if (type == "proportion") hi <- pmin(1, hi)
    tibble(denom = d, level = nm, lo = lo, hi = hi)
  }))
}
