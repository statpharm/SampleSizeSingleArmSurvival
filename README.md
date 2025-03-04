# SampleSizeSingleArmSurvival
Provides methods to calculate sample size for single-arm survival      studies using the arcsine transformation, incorporating uniform accrual      and exponential survival assumptions. Includes functionality for      detailed numerical integration and simulation.




## Introduction

The `SampleSizeSingleArmSurvival` package provides a method for calculating sample sizes for single-arm survival studies using the arcsine transformation. This approach accounts for uniform accrual and exponential survival assumptions, leveraging the work of Nagashima et al. (2021).

This vignette demonstrates how to use the `calcSampleSizeArcsine()` function to calculate the required sample size based on specific survival probabilities, time points, and study designs.

## Installation

To install the package, you can use:

```r
# From source
devtools::install_github("statpharm/SampleSizeSingleArmSurvival")
```

## Usage

### Basic Example

The most basic use case calculates the sample size required when you specify the null and alternative survival probabilities:

```r
library(SampleSizeSingleArmSurvival)
required_sample <- calcSampleSizeArcsine(S0 = 0.90, S1 = 0.96)
required_sample
#> [1] 107
```

Here, `S0` represents the survival probability under the null hypothesis, and `S1` is the survival probability under the alternative hypothesis.

### Custom Study Parameters

You can adjust the study parameters to reflect different study designs, such as longer accrual periods or different time points of interest:

```r
custom_sample <- calcSampleSizeArcsine(
  S0 = 0.80,
  S1 = 0.85,
  accrual = 36,
  followup = 24,
  timePoint = 18
)
custom_sample
#> [1] 356
```

This example calculates the required sample size for a study with an accrual period of 36 months, an additional follow-up period of 12 months, and a time point of interest at 24 months.

### Handling Warnings

The function will provide warnings if parameters result in non-positive or zero effect sizes, or if the time point of interest exceeds the study duration:

```r
result <- calcSampleSizeArcsine(S0 = 0.90, S1 = 0.90)
#> Warning in calcSampleSizeArcsine(S0 = 0.9, S1 = 0.9): S1 <= S0: effect size is
#> non-positive or zero. Returning NA.
```

Warnings like these are meant to guide you in adjusting your parameters.

## Details of the Method

The method uses the arcsine transformation of survival probabilities to stabilize the variance and simplify comparisons between groups. Specifically:

1. The survival probabilities are converted to a hazard rate assuming exponential survival (`S(t) = exp(-lambda * t)`).
2. Variance at the time point of interest is calculated using the Greenwood formula.
3. The required sample size is solved iteratively using the `uniroot` function.

The arcsine transformation is particularly useful for stabilizing variance when survival probabilities are near 0 or 1.

## References

Nagashima, K., & colleagues (2021). Sample size calculations for single-arm survival studies using transformations of the Kaplan–Meier estimator. *Pharmaceutical Statistics*, 20(2), 228–241. [https://onlinelibrary.wiley.com/doi/10.1002/pst.2090](https://onlinelibrary.wiley.com/doi/10.1002/pst.2090)

## Conclusion

The `SampleSizeSingleArmSurvival` package offers a robust tool for calculating sample sizes for single-arm survival studies, accommodating various study designs and assumptions.

