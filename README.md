
<!-- README.md is generated from README.Rmd. Please edit that file -->

# reasonabletools

<!-- badges: start -->

<!-- badges: end -->

This package provides the user with convenience functions for cleaning
and summarising data for use in Reasonable Potential Analyses (RPAs) and
Water Quality-based Effluent Limitation (WQBEL) development. Examples of
tasks one can do with reasonable tools include constructing composite
parameter measurment (e.g., total PCBs) from congener water samples,
finding maximum observed concentrations from partially censored
datasets, and projecting maximum effluent concentrations (MECs). Methods
for projecting MECs and computing effluent limitations are based on
EPA’s Technical Support Document for Water Quality-based Toxics
Control (TSD).

## Installation

You can install the current version of reasonabletools from
[github](https://github.com/mattreusswig/convertUnits) usng devtools.

``` r
# If devetools is not installed, install it.
install.packages("devtools")

# Install reasonabletools from the github repository.
devtools::install_github("mattreusswig/reasonabletools")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(reasonabletools)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
# summary(cars)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub\!
