adsl <- build_adsl(pharmaversesdtm::dm, pharmaversesdtm::ex, pharmaversesdtm::ds)
adae <- build_adae(pharmaversesdtm::ae, adsl)
patients <- patient_table(adsl, adae)

test_that("the patient table carries every treated patient's TEAEs", {
  expect_equal(sum(patients$teae),
               sum(adae$TRTEMFL == "Y" & adae$USUBJID %in% patients$USUBJID, na.rm = TRUE))
  expect_true(all(patients$years > 0))
})

test_that("a site's reference never draws its own patients", {
  # Give site 701 an absurd AE count: the reference must not move.
  inflated <- patients
  inflated$teae[inflated$SITEID == "701"] <- 1000L
  a <- reference_rates(patients, "701", 20, draws = 500, seed = 9)
  b <- reference_rates(inflated, "701", 20, draws = 500, seed = 9)
  expect_equal(a, b)
})

test_that("the reference is centred on the other sites' rate", {
  others <- patients[patients$SITEID != "710", ]
  ref <- reference_rates(patients, "710", 30, draws = 4000, seed = 2)
  expect_equal(median(ref), sum(others$teae) / sum(others$years), tolerance = 0.1)
})

test_that("a large site hiding 90% of its AEs is flagged low", {
  planted <- plant_underreporting(adae, "701", share = 0.9, seed = 4)
  res <- resample_sites(patient_table(adsl, planted), draws = 2000)
  expect_equal(res$flag[res$SITEID == "701"], "low")
})

test_that("on the untouched data at most one site is flagged", {
  res <- resample_sites(patients, draws = 2000)
  expect_lte(sum(res$flag %in% c("low", "high")), 1)
})
