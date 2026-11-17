dm <- pharmaversesdtm::dm
ex <- pharmaversesdtm::ex
ds <- pharmaversesdtm::ds
ae <- pharmaversesdtm::ae
adsl <- build_adsl(dm, ex, ds)
adae <- build_adae(ae, adsl)

test_that("ADSL has one row per randomised patient and no screen failures", {
  expect_equal(nrow(adsl), sum(dm$ARM != "Screen Failure"))
  expect_equal(anyDuplicated(adsl$USUBJID), 0L)
  expect_false("Screen Failure" %in% adsl$TRT01A)
})

test_that("treatment dates are ordered and duration counts both ends", {
  treated <- dplyr::filter(adsl, SAFFL == "Y", !is.na(TRTEDT))
  expect_true(all(treated$TRTSDT <= treated$TRTEDT))
  expect_equal(treated$TRTDURD, as.numeric(treated$TRTEDT - treated$TRTSDT) + 1)
})

test_that("every patient has an end-of-study status", {
  expect_true(all(adsl$EOSSTT %in% c("COMPLETED", "DISCONTINUED", "ONGOING")))
  expect_true(all(is.na(adsl$DCSREAS[adsl$EOSSTT == "COMPLETED"])))
})

test_that("an AE is treatment-emergent only inside the treatment window + 30 days", {
  te <- dplyr::filter(adae, TRTEMFL == "Y")
  expect_true(all(te$ASTDT >= te$TRTSDT))
  expect_true(all(te$ASTDT <= te$TRTEDT + 30, na.rm = TRUE))
  before <- dplyr::filter(adae, !is.na(ASTDT), ASTDT < TRTSDT)
  expect_true(all(is.na(before$TRTEMFL)))
})
