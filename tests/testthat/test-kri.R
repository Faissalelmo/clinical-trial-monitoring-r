test_that("a site at the study rate is ok and one with no events is low", {
  # Five sites at 10 a year; a sixth with one year of exposure reports none
  # where about 10 are expected (Poisson P(0) is about 1 in 18,000).
  count <- c(100, 100, 100, 100, 100, 0)
  denom <- c(10, 10, 10, 10, 10, 1)
  res <- judge(count, denom, "rate", adjust = FALSE)
  expect_equal(res$flag[1:5], rep("ok", 5))
  expect_equal(res$flag[6], "low")
})

test_that("a site with no denominator is not judged", {
  res <- judge(c(5, 5, 0), c(10, 10, 0), "proportion", adjust = FALSE)
  expect_equal(res$flag[3], "no data")
})

test_that("small sites get wider limits than large ones", {
  # Same 50% shortfall: 5 of 10 expected is noise, 50 of 100 is not.
  res <- judge(c(5, 50, 1000), c(1, 10, 100), "rate", adjust = FALSE)
  expect_equal(res$flag[1], "ok")
  expect_equal(res$flag[2], "low")
})

test_that("overdispersion is never below 1 and ignores one extreme site", {
  expect_equal(overdispersion(rep(0.1, 10)), 1)
  z <- c(rnorm(19, sd = 1.5), 40)
  set.seed(3)
  expect_lt(overdispersion(z), 10)
})

test_that("planting under-reporting removes only that site's treatment-emergent AEs", {
  adae <- data.frame(SITEID = c("A", "A", "A", "B"), TRTEMFL = c("Y", "Y", NA, "Y"))
  out <- plant_underreporting(adae, "A", share = 1, seed = 1)
  expect_equal(nrow(out), 2)
  expect_equal(sort(out$SITEID), c("A", "B"))
  expect_true(is.na(out$TRTEMFL[out$SITEID == "A"]))
})
