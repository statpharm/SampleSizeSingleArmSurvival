## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(SampleSizeSingleArmSurvival)
required_sample <- calcSampleSizeArcsine(S0 = 0.90, S1 = 0.96)
required_sample

## -----------------------------------------------------------------------------
custom_sample <- calcSampleSizeArcsine(
  S0 = 0.80,
  S1 = 0.85,
  accrual = 36,
  followup = 24,
  timePoint = 18
)
custom_sample

## -----------------------------------------------------------------------------
result <- calcSampleSizeArcsine(S0 = 0.90, S1 = 0.90)

