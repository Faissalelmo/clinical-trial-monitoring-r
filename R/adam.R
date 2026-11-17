# ADaM-style analysis datasets built with {admiral} from the CDISC pilot SDTM
# data shipped in {pharmaversesdtm}.

suppressPackageStartupMessages({
  library(admiral)
  library(dplyr)
  library(stringr)
})

#' Subject-level dataset (ADSL): treatment dates, duration, safety flag,
#' end-of-study status and site.
build_adsl <- function(dm, ex, ds) {
  ex_ext <- ex %>%
    derive_vars_dtm(dtc = EXSTDTC, new_vars_prefix = "EXST") %>%
    derive_vars_dtm(dtc = EXENDTC, new_vars_prefix = "EXEN", time_imputation = "last")

  # A record counts as a dose when drug was given, or when it is placebo.
  is_dose <- expr(EXDOSE > 0 | (EXDOSE == 0 & str_detect(EXTRT, "PLACEBO")))

  eos <- ds %>%
    filter(DSCAT == "DISPOSITION EVENT", DSDECOD != "SCREEN FAILURE") %>%
    transmute(STUDYID, USUBJID,
              EOSSTT = if_else(DSDECOD == "COMPLETED", "COMPLETED", "DISCONTINUED"),
              DCSREAS = if_else(DSDECOD == "COMPLETED", NA_character_, DSDECOD))

  dm %>%
    filter(ARM != "Screen Failure") %>%
    derive_vars_merged(
      dataset_add = ex_ext,
      filter_add = !!is_dose & !is.na(EXSTDTM),
      new_vars = exprs(TRTSDTM = EXSTDTM),
      order = exprs(EXSTDTM, EXSEQ),
      mode = "first",
      by_vars = exprs(STUDYID, USUBJID)
    ) %>%
    derive_vars_merged(
      dataset_add = ex_ext,
      filter_add = !!is_dose & !is.na(EXENDTM),
      new_vars = exprs(TRTEDTM = EXENDTM),
      order = exprs(EXENDTM, EXSEQ),
      mode = "last",
      by_vars = exprs(STUDYID, USUBJID)
    ) %>%
    derive_vars_dtm_to_dt(source_vars = exprs(TRTSDTM, TRTEDTM)) %>%
    derive_var_trtdurd() %>%
    mutate(
      SAFFL = if_else(!is.na(TRTSDT), "Y", "N"),
      TRT01A = ARM
    ) %>%
    derive_vars_merged(dataset_add = eos, by_vars = exprs(STUDYID, USUBJID)) %>%
    mutate(EOSSTT = coalesce(EOSSTT, "ONGOING")) %>%
    select(STUDYID, USUBJID, SITEID, AGE, SEX, RACE, TRT01A, SAFFL,
           TRTSDT, TRTEDT, TRTDURD, EOSSTT, DCSREAS)
}

#' Adverse-event dataset (ADAE) with analysis dates and the
#' treatment-emergent flag (onset on or after first dose, up to 30 days
#' after the last).
build_adae <- function(ae, adsl) {
  ae %>%
    derive_vars_merged(
      dataset_add = adsl,
      new_vars = exprs(SITEID, TRT01A, SAFFL, TRTSDT, TRTEDT),
      by_vars = exprs(STUDYID, USUBJID)
    ) %>%
    filter(SAFFL == "Y") %>%
    derive_vars_dt(dtc = AESTDTC, new_vars_prefix = "AST",
                   highest_imputation = "M", min_dates = exprs(TRTSDT)) %>%
    derive_vars_dt(dtc = AEENDTC, new_vars_prefix = "AEN",
                   highest_imputation = "M", date_imputation = "last") %>%
    derive_var_trtemfl(
      start_date = ASTDT, end_date = AENDT,
      trt_start_date = TRTSDT, trt_end_date = TRTEDT,
      end_window = 30
    ) %>%
    select(STUDYID, USUBJID, SITEID, TRT01A, TRTSDT, TRTEDT, AESEQ, AEDECOD,
           AEBODSYS, AESEV, AESER, ASTDT, AENDT, TRTEMFL)
}
