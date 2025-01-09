test_that("calcSampleSizeArcsine works with default parameters", {
  result <- calcSampleSizeArcsine(S0 = 0.90, S1 = 0.96)
  expect_type(result, "double")
  expect_gt(result, 0) # Sample size must be greater than 0
})

test_that("calcSampleSizeArcsine handles invalid S0 and S1", {
  expect_error(calcSampleSizeArcsine(S0 = -0.1, S1 = 0.96), "S0 must be strictly between 0 and 1")
  expect_error(calcSampleSizeArcsine(S0 = 0.90, S1 = 1.1), "S1 must be strictly between 0 and 1")
})


test_that("calcSampleSizeArcsine warns and returns NA when S1 <= S0", {
  expect_warning(result <- calcSampleSizeArcsine(S0 = 0.90, S1 = 0.90),
                 "effect size is non-positive or zero")
  expect_true(is.na(result))
})

test_that("calcSampleSizeArcsine handles extreme values", {
  result <- calcSampleSizeArcsine(S0 = 0.99, S1 = 0.999, timePoint = 12, accrual = 12, followup = 12)
  expect_gt(result, 0) # Should still calculate a positive sample size
})


