
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reasonabletools

<!-- badges: start -->

<!-- badges: end -->

This package provides the user with convenience functions for cleaning
and summarizing data for use in Reasonable Potential Analyses (RPAs) and
Water Quality-based Effluent Limitation (WQBEL) development. Examples of
tasks one can do with reasonabletools include constructing composite
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

## Assumptions for CVs and Distributions

### Minimum Information Content for CV

Consistent with the TSD, reasonabletools requires a minimum level of
information content in the dataset in order to estimate the Coefficient
of Variation (CV) for a pollutant. There must be at least 10
observations and fewer than 80 percent of the dataset must be non-detect
to calculate a CV value. If either condition is not met, a CV value of
0.6 is assumed and output.

### Non-Detects Substitution

Consistent with the TSD, reasonabletools uses substitution methods when
calculating mean and standard deviation values when non-detects are
present. When non-detects are present, result values corresponding to
non-detect qual values are multiplied by 0.5 prior to computing the mean
and standard deviation.

### Distributional Assumptions

Consistent with the TSD, reasonabletools assumes a lognormal
distribution is applicable to the data when performing effluent
limitation computations.

## Usage

Functions in reasonabletools accept vectors as arguments–typically, a
vector of numeric data representing concentration point measurements or
detection limits, and a vector of character data indicating the
detected/non-detect/detected-but-not-quantified status of each numeric
value.

Combining the input vectors in a dataframe or a tibble is not necessary,
but can be useful when developing summaries for multiple pollutants at
once (see discussion below).

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

### Non-Detects, Detected Values, and DNQs

Functions in reasonabletools can distinguish between detected and
non-detect values. Non-detect values are indicated in the qualifier
vector with pre-specified character string. By default, this is a “\<”,
“nd”, or “ND” character string. Any string in the non-detect qualifier
vector which does not match these indicators is treated as a detected
value. However, the non-detect flags can be changed using the “nd”
argument (e.g., you could change it to “My ND Indicator”–see example
below). In function output, “\<” is always used to indicate a censored
value and an empty string "" is used to indicate a detected value.

``` r
## Create an example dataset with new custom non-detect indicator
zinc <- data.frame(qual   = c(rep("My ND Indicator", 20)),
                   result = 1:20,
                   stringsAsFactors = FALSE)
zinc
#>               qual result
#> 1  My ND Indicator      1
#> 2  My ND Indicator      2
#> 3  My ND Indicator      3
#> 4  My ND Indicator      4
#> 5  My ND Indicator      5
#> 6  My ND Indicator      6
#> 7  My ND Indicator      7
#> 8  My ND Indicator      8
#> 9  My ND Indicator      9
#> 10 My ND Indicator     10
#> 11 My ND Indicator     11
#> 12 My ND Indicator     12
#> 13 My ND Indicator     13
#> 14 My ND Indicator     14
#> 15 My ND Indicator     15
#> 16 My ND Indicator     16
#> 17 My ND Indicator     17
#> 18 My ND Indicator     18
#> 19 My ND Indicator     19
#> 20 My ND Indicator     20

## Find observed MEC
find_mec(qual = zinc$qual, result = zinc$result, 
         nd = c("My ND Indicator"))
#>   qual result
#> 1    <      1
```

When cleaning data for summary, recall that R is a case-sensitive
language and capitalization matters in non-detect flags. The default
censored flags in reasonabletools are “\<”, “nd”, and “ND”. This means
an observation in your data marked “nd” will be treated as a non-detect,
but an observation marked “nD” (notice the change in capitalization on
the second letter) will be treated as a detected value.

All elements in the results vector should have a numeric value–do not
include NA values. For non-detect values, the method detection limit
should be entered in the results column. NA values in the results column
are ignored and will result in incorrect results.

Functions in reasonabletools do not have specialized behavior for an
observation which is detected-but-not-quantified (DNQ or j-flag). The
user must decide when preprocessing the data whether these values should
be treated as non-detect values or as detected values. The observation
should be marked as non-detect or detected in the qualifier vector, and
the reporting limit or estimated concentration placed in result vector.

## Working with Multiple Pollutants

Real-life datasets used in NPDES permit development will typically
include many pollutants of concern which need to be summarized. This is
a split-apply-combine task which is cumbersome to address in Excel and
other spreadsheet-based programs for anything but simple summary
functions (e.g., computing averages or maximums without conditionality).

In base R (or using dplyr/tidyverse methods) this task becomes
straightforward. Here are some strategies from base R and dplyr for
summarizing effluent data with reasonabletools.

We will use the following toy dataset in these examples:

``` r
set.seed(10010)
data <- data.frame(pollutant = sort(rep(c("zinc", "arsenic", "copper",
                                          "lead", "nickel", "pcbs"), 5)),
                   qualifier = sample(c("<", ""), size = 30, 
                                      replace = TRUE, prob = c(0.8, 0.2)),
                   result = sample(seq(0.1, 0.5, 0.1), 30, replace = TRUE),
                   sampling.date = rep(seq.Date(from = as.Date("2012-03-10"), 
                                                to = as.Date("2012-03-10") + 365 * 5,
                                                by = 365), 5),
                   stringsAsFactors = FALSE)

data
#>    pollutant qualifier result sampling.date
#> 1    arsenic         <    0.1    2012-03-10
#> 2    arsenic              0.3    2013-03-10
#> 3    arsenic         <    0.4    2014-03-10
#> 4    arsenic         <    0.1    2015-03-10
#> 5    arsenic         <    0.3    2016-03-09
#> 6     copper         <    0.5    2017-03-09
#> 7     copper         <    0.4    2012-03-10
#> 8     copper         <    0.5    2013-03-10
#> 9     copper         <    0.5    2014-03-10
#> 10    copper         <    0.2    2015-03-10
#> 11      lead         <    0.1    2016-03-09
#> 12      lead         <    0.1    2017-03-09
#> 13      lead         <    0.2    2012-03-10
#> 14      lead         <    0.4    2013-03-10
#> 15      lead         <    0.4    2014-03-10
#> 16    nickel              0.5    2015-03-10
#> 17    nickel              0.2    2016-03-09
#> 18    nickel         <    0.4    2017-03-09
#> 19    nickel         <    0.1    2012-03-10
#> 20    nickel         <    0.3    2013-03-10
#> 21      pcbs         <    0.2    2014-03-10
#> 22      pcbs              0.4    2015-03-10
#> 23      pcbs         <    0.5    2016-03-09
#> 24      pcbs         <    0.2    2017-03-09
#> 25      pcbs         <    0.4    2012-03-10
#> 26      zinc              0.4    2013-03-10
#> 27      zinc         <    0.3    2014-03-10
#> 28      zinc         <    0.5    2015-03-10
#> 29      zinc              0.2    2016-03-09
#> 30      zinc         <    0.3    2017-03-09
```

### Base R

Use the split family of functions to create a list of dataframes grouped
by pollutant name, use lapply to apply the function to each element of
the list, and do.call + rbind to recombine the output.

``` r

data_split <- split.data.frame(x = data, f = data$pollutant)

## Note we need to use an anonymous function in lapply for functions with
## multiple non-default input arguments, like find_mec.
data_apply <- lapply(X = data_split, 
                     FUN = function(x) find_mec(x$qualifier, x$result))

combined_summary <- do.call(rbind, data_apply)

combined_summary
#>         qual result
#> arsenic         0.3
#> copper     <    0.2
#> lead       <    0.1
#> nickel          0.5
#> pcbs            0.4
#> zinc            0.4
```

Or, we can combine the apply and combine steps into a single line by
nesting the functions. Unfortunately, we can’t do this in fewer than two
steps.

``` r
data_split <- split.data.frame(x = data, f = data$pollutant)

combined_summary <- do.call(
       rbind, 
       lapply(X = data_split, 
              FUN = function(x) find_mec(x$qualifier, x$result))
  )

combined_summary
#>         qual result
#> arsenic         0.3
#> copper     <    0.2
#> lead       <    0.1
#> nickel          0.5
#> pcbs            0.4
#> zinc            0.4
```

### dplyr (v1.0 or later)

Using dplyr is often a more intuitive way to implement
split-apply-combine for many people. Using dplyr, we will use
dplyr::group\_by and wrap the reasonabletools summary function in
dplyr::summarise and it will automatically output all columns grouped by
pollutant. (Note that old versions of dplyr::summarise (pre-1.0) are not
set up to output two columns from a single summary function–see the next
section for pre-1.0 dplyr strategies.)

``` r
library(dplyr)

# data %>% 
#   group_by(pollutant) %>% 
#   summarise(find_mec(qualifier, result))
```

### dplyr (Pre-v1.0)

Older versions of dplyr were not set up to handle multiple column output
in the summarise functions. There were two main workarounds for this–one
built into reasonabletools and one using dplyr itself.

A number of data summarising functions in reasonabletools that output
two column values, like find\_mec or project\_mec, include an option for
automatically concatenating the two output vectors into a single
character vector. The “simple\_output” argument will produce a single
character vector when set to TRUE, and will produce two vectors when
left at its default FALSE value.

``` r
library(dplyr)

## simple_output = TRUE turns output into single character vector.
data %>% 
  group_by(pollutant) %>% 
  summarise(MEC = find_mec(qualifier, result, simple_output = TRUE))
#> # A tibble: 6 x 2
#>   pollutant MEC  
#>   <chr>     <chr>
#> 1 arsenic   0.3  
#> 2 copper    <0.2 
#> 3 lead      <0.1 
#> 4 nickel    0.5  
#> 5 pcbs      0.4  
#> 6 zinc      0.4
```

Note this setting is also useful when constructing summary tables for
inclusion in permit documents where one wants to present the qualifier
and value in the same table cell.

If one wants both output vectors retained in the output, use the
dplyr::do function as a wrapper around the desired reasonabletools
function. Use a dot prior to the $ operator and column name (e.g.,
.$qualifier) to indicate to dplyr::do it should look for the inputs in
the piped dataframe.

``` r

data %>% 
  group_by(pollutant) %>% 
  do(find_mec(.$qualifier, .$result)) 
#> # A tibble: 6 x 3
#> # Groups:   pollutant [6]
#>   pollutant qual  result
#>   <chr>     <chr>  <dbl>
#> 1 arsenic   ""       0.3
#> 2 copper    "<"      0.2
#> 3 lead      "<"      0.1
#> 4 nickel    ""       0.5
#> 5 pcbs      ""       0.4
#> 6 zinc      ""       0.4
```
