
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reasonabletools

<!-- badges: start -->

<!-- badges: end -->

This package provides the user with convenience functions for cleaning
and summarizing data for use in Reasonable Potential Analyses (RPAs) and
Water Quality-based Effluent Limitation (WQBEL) development. Examples of
tasks one can do with reasonable tools include constructing composite
parameter measurement (e.g., total PCBs) from congener water samples,
finding maximum observed concentrations from partially censored
datasets, and projecting maximum effluent concentrations (MECs). Methods
for projecting MECs and computing effluent limitations are based on
EPA’s Technical Support Document for Water Quality-based Toxics
Control (TSD).

## Installation

You can install the current version of reasonabletools from
[github](https://github.com/mattreusswig/reasonabletools) using
devtools.

``` r
# If devetools is not installed, install it.
install.packages("devtools")

# Install reasonabletools from the github repository.
devtools::install_github("mattreusswig/reasonabletools")
```

## Usage

Functions in reasonabletools accept vectors as arguments–typically, a
vector of numeric data representing concentration point measurements or
detection limits, and a vector of character data indicating the
detected/non-detect/detected-but-not-quantified status of each numeric
value.

``` r
library(reasonabletools)

## Create an example dataset
zinc <- data.frame(qual   = c(rep("", 10), rep("<", 10)),
                   result = 1:20,
                   stringsAsFactors = FALSE)
zinc
#>    qual result
#> 1            1
#> 2            2
#> 3            3
#> 4            4
#> 5            5
#> 6            6
#> 7            7
#> 8            8
#> 9            9
#> 10          10
#> 11    <     11
#> 12    <     12
#> 13    <     13
#> 14    <     14
#> 15    <     15
#> 16    <     16
#> 17    <     17
#> 18    <     18
#> 19    <     19
#> 20    <     20

## Find observed MEC
find_mec(qual = zinc$qual, result = zinc$result)
#>   qual result
#> 1          10
```

### Non-Detects

By default, non-detect values and detected-but-not-quantified (DNQ or
j-flag) values are both treated as non-detect and are indicated by
either a “\<”, “nd”, or “ND” character string. Any string in the
non-detect qualifier vector which does not match these indicators is
treated as a detected value. However, the non-detect indicators can be
changed using the “nd” argument. In function output, “\<” is always used
to indicate a censored value and an empty string "" is used to indicate
a detected value.

``` r
## Create an example dataset with new custom non-detect indicator
zinc <- data.frame(qual   = c(rep("My Indicator", 20)),
                   result = 1:20,
                   stringsAsFactors = FALSE)
zinc
#>            qual result
#> 1  My Indicator      1
#> 2  My Indicator      2
#> 3  My Indicator      3
#> 4  My Indicator      4
#> 5  My Indicator      5
#> 6  My Indicator      6
#> 7  My Indicator      7
#> 8  My Indicator      8
#> 9  My Indicator      9
#> 10 My Indicator     10
#> 11 My Indicator     11
#> 12 My Indicator     12
#> 13 My Indicator     13
#> 14 My Indicator     14
#> 15 My Indicator     15
#> 16 My Indicator     16
#> 17 My Indicator     17
#> 18 My Indicator     18
#> 19 My Indicator     19
#> 20 My Indicator     20

## Find observed MEC
find_mec(qual = zinc$qual, result = zinc$result, 
         nd = c("My Indicator"))
#>   qual result
#> 1    <      1
```

## Assumptions for CVs and Distributions

### Minimum Information Content for CV

Consistent with the TSD, reasonabletools requires a minimum level of
information content in the dataset in order to estimate the Coefficient
of Variation (CV) for a pollutant. There must be at least 10
observations and fewer than 80 percent of the dataset must be non-detect
to calculate a CV value. If both those conditions are not met, a CV
value of 0.6 is assumed and output.

### Treatment of Non-Detects

Consistent with the TSD, reasonabletools uses substitution methods when
calculating mean and standard deviation values when non-detects are
present. When non-detects are present, result values corresponding to
non-detect qual values are multiplied by 0.5 prior to computing the mean
and standard deviation.

### Distributional Assumptions

Consistent with the TSD, reasonabletools assumes a lognormal
distribution is applicable to the data when performing effluent
limitation computations.
