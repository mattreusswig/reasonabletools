#' @title Find the adjusted coeffecient of variation
#' @description Calculates adjusted coeficient of variation (CV) according to methods described in EPA's Technical Support Document for Water Quality-based Toxics Control.
#'
#' @param qual  A character vector containing non-detect indicator strings, e.g., "<" or "ND". The strings used to indicate censored status can be edited in the "nd" argument.
#' @param result A numeric vector of concentration measurements.
#' @param nd A list indicating all the censoring flags included in the dataset. Defaults to "<", "nd", and "ND".
#' @param nd_adjustment Adjustment factor for non-dectect values. Non-detect values (as indicated in qual vector) are multiplied by nd_adjustment factor; i.e. result * nd_adjustment. Typically, method detection limits or reporting limits are used as result values for non-detects.
#'
#' @return A numeric coefficient of variation (CV) value
#' @export
#'
#' @examples 
#' # CV for all detected values 
#' cen_result <- rep("", 10)
#' result     <- c(1:10) 
#' cv_adj(cen_result, result)
#' 
#' # CV for all non-detected values
#' cen_result <- rep("<", 10)
#' cv_adj(cen_result, result)
#' 
#' # CV for fewer than 10 measurements
#' cen_result <- rep("", 5)
#' result <-     c(1:5)
#' cv_adj(cen_result, result)
#' 
#' # Change the default substitution value
#' cen_result <- c(rep("<", 5), rep("", 15))
#' result     <- c(101:120)
#' cv_adj(cen_result, result)   # Use default 0.5 multipler
#' cv_adj(cen_result, result, nd_adjustment = 1.0)  # Use 1.0 multiplier (equivalent to using MDL)
#' cv_adj(cen_result, result, nd_adjustment = 0)  # Use 0.0 multiplier (equivalent to zero substitution)
cv_adj <- function(qual, result, 
                   nd = c("<", "nd", "ND"), nd_adjustment = 0.5) {
  
  
  ## Error message is given if the the 'result' argument is not numeric.
  result_class <- class(result)
  if (!is.numeric(result)) stop("result arugment must be numeric. The result
                               argument entered is a ", result_class, " type.",
                                        call. = FALSE)
  

  ## Error message is given if the the 'result' and 'qual' arguments are of
  ## unequal lengths.
  if (length(result) != length(qual)) stop("The qual and result arguments must be equal in length. Qual is ", length(qual), " elements long, and result is ", length(result), ".", call. = FALSE)
  
  
  ## Main body of function
  
  # Helper function which creates an anti-union operator (opposite of %in%)
  "%not_in%" <- function(x, y) { !(x %in% y) }
  
  ## Index which observations are detects and which are nondetects
  detects <- which(qual %not_in% nd)             ## gives index of detects
  nondetects <- which(qual %in% nd)              ## gives index of nondetects
  
  
  ## The CV is conditional on data quality (i.e., n and proportion of censored
  ## values)
  n <- length(result)                           ## gives number of observations
  censored_prop <- length(nondetects) / n       ## gives fraction of non-detects
  
  
  ## Compute CV conditionally. If low data quality, CV == 0.6
  if ( (n < 10) | (censored_prop >= 0.8) ) { cv <- 0.6 
  
  
  ## If acceptable data quality, CV is the ratio of the sd to mean.
  ## Adjust censored values to the desired substitution values using
  ## the nd_adjustment argument.
  } else {
    
    ## Issue a warning when NA's are present which is unexpected behavior.
    if (anyNA(result)) {
      na_index <- which(is.na(result))
      na_index <- paste(na_index, collapse = ", ")
      warning(paste0("Your result argument contains NA values at index loacation(s) ", na_index, " which have been removed in this calculation. Your coefficient of variation estimate may contain errors."))
    } 
    
    ## Compute CV
    std_deviation <- sd( c(result[detects], 
                           nd_adjustment * result[nondetects]), 
                         na.rm = TRUE)
    
    average       <- mean( c(result[detects], 
                             nd_adjustment * result[nondetects]),
                           na.rm = TRUE)
    
    cv <- std_deviation / average  
    
  }
  
  return(cv)        ## Exit
  
}


#' @title Find the Maximum Observed Effluent Concentration (MEC)
#' @description Find the MEC (no projection) from the observed dataset using methods described in EPA's Technical Support Document for Water Quality-based Toxics Control.
#'
#' @param qual  A character vector containing non-detect indicator strings, e.g., "<" or "ND". The strings used to indicate censored status can be edited in the "nd" argument.
#' @param result A numeric vector of concentration measurements.
#' @param nd A list indicating all the censoring flags included in the dataset. Defaults to "<", "nd", and "ND". 
#' @param simple_output Logical value. If TRUE, the output columns are concatenated into a single character string (e.g., "<0.2") which can be useful for constructing summary tables. 
#'
#' @return Dataframe with a qualifier column (character) and a MEC column (numeric).
#' @export
#'
#' @examples
#' # Find observed MEC
#' cen_result <- c(rep("", 10), rep("<", 10))
#' result     <- 1:20
#' find_mec(cen_result, result)
#' 
#' cen_result <- rep("<", 20)
#' find_mec(cen_result, result)
#' 
#' cen_result <- rep("", 20)
#' find_mec(cen_result, result)
#' 
#' # Demonstrate simplified output
#' find_mec(cen_result, result, simple_output = TRUE)
#' 
#' # Define a set of custom non-detect flags
#' cen_result <- c(rep("non-detect", 5), rep("<", 10), rep("mdl", 5))
#' find_mec(cen_result, result, nd = c("non-detect", "<", "mdl"))
#' 
find_mec <- function(qual, result, 
                     nd = c("<", "nd", "ND"), simple_output = FALSE) {
  
  ## Error message is given if the the 'result' argument is not numeric.
  result_class <- class(result)
  
  if (!is.null(result)) {
    if (!is.numeric(result)) stop("The result arugment must be numeric. The result argument entered is a ", result_class, " type.",
                                          call. = FALSE)
  }
  
  ## Error message is given if the the 'result' and 'qual' arguments are of
  ## unequal lengths.
  if (length(result) != length(qual)) stop("The qual and result arguments must be equal in length. Qual is ", length(qual), " elements long, and result is ", length(result), ".", call. = FALSE)
  
  ## consider adding a helper function that looks for character strings not
  ## included in the provided "nd" list within the qual vector, and warn
  ## the user they are present.
  
  
  ## Main body of function
  
  # Helper function which creates an anti-union operator (opposite of %in%)
  "%not_in%" <- function(x, y) { !(x %in% y) }
  
  # Create indices for detected and non-detected results
  detects <- which(qual %not_in% nd)    ## gives index of detects
  nondetects <- which(qual %in% nd)     ## gives index of nondetects
  
  # How many observations are there
  n <- length(result)
  n_nd <- length(nondetects)
  
  all_nd <- if (n_nd == n) TRUE else FALSE
  
  if (all_nd == TRUE) {
    
    mec <- data.frame(qual = "<", result = min(result, na.rm = TRUE), 
                      stringsAsFactors = FALSE)

  }  
  
  if (all_nd == FALSE) {
    
    mec <- data.frame(qual = "", result = max(result[detects], 
                                              na.rm = TRUE), 
                      stringsAsFactors = FALSE)
    
  }
  
  if (isTRUE(simple_output)) mec <- paste0(mec[[1]], mec[[2]])
  
  return(mec)
  
}


#' @title Find the projected Maximum Observed Effluent Concentration (MEC)
#' @description Find the MEC projected from a lognormal distribution using methods described in EPA's Technical Support Document for Water Quality-based Toxics Control.
#'
#' @param qual  A character vector containing non-detect indicator strings, e.g., "<" or "ND". The strings used to indicate censored status can be edited in the "nd" argument.
#' @param result A numeric vector of concentration measurements.
#' @param nd A list indicating all the censoring flags included in the dataset. Defaults to "<", "nd", and "ND". 
#' @param percentile Location on lognormal distribution to estimate the projected MEC.
#' @param conf_level Confidence level of projected estimate.
#' @param nd_adjustment Adjustment factor for non-dectect values. Non-detect values (as indicated in qual vector) are multiplied by nd_adjustment factor; i.e. result * nd_adjustment. Typically, method detection limits or reporting limits are used as result values for non-detects.
#' @param simple_output Logical value. If TRUE, the output columns are concatenated into a single character string (e.g., "<0.2") which can be useful for constructing summary tables.
#'
#' @return  Dataframe with a qualifier column (character) and a MEC column (numeric).
#' @export
#'
#' @examples
#' # Find observed MEC
#' cen_result <- c(rep("", 10), rep("<", 10))
#' result     <- 1:20
#' project_mec(cen_result, result)
#' 
#' # Demonstrate simplified output
#' cen_result <- rep("<", 20)
#' project_mec(cen_result, result, simple_output = TRUE)
#' 
#' # Define a set of custom non-detect flags
#' cen_result <- c(rep("non-detect", 5), rep("<", 10), rep("mdl", 5))
#' project_mec(cen_result, result, nd = c("non-detect", "<", "mdl"))
#' 
#' # Change the substitution multiplier used for non-detect values
#' cen_result <- c(rep("<", 5), rep("", 15))
#' result     <- c(101:120)
#' project_mec(cen_result, result)   # Use default 0.5 multipler
#' project_mec(cen_result, result, nd_adjustment = 1.0)  # Use 1.0 multiplier (equivalent to using MDL)
#' project_mec(cen_result, result, nd_adjustment = 0)  # Use 0.0 multiplier (equivalent to zero substitution)
#' 
project_mec <- function(qual, result, 
                        nd = c("<", "nd", "ND"),
                        percentile = 0.95, conf_level = 0.99,
                        nd_adjustment = 0.5, simple_output = FALSE) {
  
  
  ## Error message is given if the the 'result' argument is not numeric.
  result_class <- class(result)
  
  if (!is.numeric(result)) stop("The result arugment must be numeric. The result argument entered is a ", result_class, " type.",
                                        call. = FALSE)
  
  
  ## Error message is given if the qual and result args do not have the same
  ## number of elements.
  n_qual = length(qual)
  n_rslt = length(result)
  
  if (n_qual != n_rslt) stop("The qual and result vectors must have an equal number of elements. The qual argument has ", n_qual, " elements and the result argument has ", n_rslt, " elements.")
  
  
  ## Find the indices of all the nondetect and detect observations
  detects    <- which(!(qual %in% nd))  ## gives index of detects
  nondetects <- which(qual %in% nd)  ## gives index of nondetects
  
  
  ## If all observations are non-dectect, then return the 
  ## minimum detection limit.
  all_nd <- if (length(nondetects) == length(qual)) TRUE else FALSE
  
  if (all_nd == TRUE) {
    
    mec = data.frame(qual = "<", 
                     result = min(result, na.rm = TRUE),
                     stringsAsFactors = FALSE)
  }   
  
  
  ## If at least one detected value then calculate a projected 
  ## MEC based on n and CV.
  if (all_nd == FALSE) {
    ## Compute the sample max effluent concentration.
    df <- find_mec(qual = qual, result = result, 
                   nd = nd, simple_output = FALSE)
    
    max_obs_result <- df[[2]]
    
    
    ## Calculate the total number of observations, and the location--as a
    ## percentile (p_n) of the sample mec.
    n_obs <- length(result)
    p_n <- (1 - conf_level)^(1/n_obs)
    
    
    ## Uses cv_adj() from this script to compute coeff of variation. 
    cv <- cv_adj(qual = qual, result = result, 
                 nd = nd, nd_adjustment = nd_adjustment)
    
    
    ## Compute sigma for a lognormal distribution, according to the TSD
    sigma <- sqrt(log(1 + cv^2))  ## log() is the natural log (i.e. ln) in R
    
    
    ## Compute the mec multiplier, according to the TSD method
    c_target   <- exp(qnorm(percentile) * sigma - 0.5 * sigma^2)
    c_observed <- exp(qnorm(p_n) * sigma - 0.5 * sigma^2)
    
    
    ## Compute the MEC
    mec <- data.frame(qual = "", 
                      result = max_obs_result * (c_target / c_observed), 
                      stringsAsFactors = FALSE)
  }
  
  ## If simple_output == TRUE, then convert the 2 column mec dataframe to
  ## a single character value
  if (isTRUE(simple_output)) mec <- paste0(mec[[1]], mec[[2]])
  
  return(mec)                                            ## EXIT
  
}



#' @title Compute Wasteload Allocation
#' @description Compute the wasteload allocation (WLA) used for effluent limit calculation.
#'
#' @param criteria Limiting water quality criterion. Must be in same units as background argument.
#' @param background Background pollutant concentration. Must be in same units as criteria argument.
#' @param Qrsw Upstream limiting/design receiving water flowrate. Must be in same units as Qeff argument. Flow arguments may be entered in ratio form.
#' @param Qeff Effluent limiting/design flowrate. Must be in same units as Qrsw argument.  Flow arguments may be entered in ratio form.
#'
#' @return WLA as numeric value in same units as criteria and background.
#' @export
#'
#' @examples
#' # WLA for pollutant with 2 ug/L acute criteria and upstream receiving water concentration 
#' # of 0.1 ug/L. The critical flows are 3 MGD (1Q10) and 0.5 MGD (max daily flow).
#' calc_WLA(2, 0.1, 3, 0.5) 
#' 
#' # When using dilution credits, put Qrsw and Qeff in terms of the dilution ratio (D).
#' D = 7   # Assume a jurisdiction that uses D = (Qrsw + Qeff) / Qeff
#' Qeff = 1 # Equal to 1 since its the denominator of the ratio, or you can use the critical flow
#' Qrsw = D - 1 # Same as the expression (D * Qeff) - Qeff
#' calc_WLA(2, 0.1, Qrsw = D - 1, Qeff = 1)
calc_WLA <- function(criteria, background, Qrsw, Qeff) {
  
  ## For dilution ratios, set Qeff = 1 and Qrsw to appropriate value
  
  Q <- Qrsw + Qeff
  
  WLA <- (Q * criteria - Qrsw * background) / Qeff
  
  return(WLA)                                           ## EXIT
   
}



#' @title Compute MDL
#' @description Compute Maximum Daily Effluent Limitation (MDL) according to methods in EPA's Technical Support Document for Water Quality-based Toxics Control.
#'
#' @param WLAa Numeric. Acute wasteload allocation (WLA). All WLA units should be identical.
#' @param WLAc Numeric. Chronic wasteload allocation (WLA). All WLA units should be identical.
#' @param WLAhh Numeric. Human health wasteload allocation (WLA). All WLA units should be identical.
#' @param cv Numeric. Coefficient of variation (CV) of effluent data. See cv_adj function.
#' @param n_samples Numeric. Number of sample observations.
#' @param prob_LTA Numeric (fraction). Allowable exceedance probability of the WLA used to estimate long-term average (LTA).
#' @param percentile_MDL Numeric (fraction). Lognormal distribution location for MDL. 
#' @param percentile_AML Numeric (fraction). Lognormal distribution location for AML.
#'
#' @return Numceric value in same units at WLA.
#' @export
#'
#' @examples
#' 
calc_MDL <- function(WLAa, WLAc, WLAhh, cv, 
                    n_samples = 4, prob_LTA = 0.99, 
                    percentile_MDL = 0.99, percentile_AML = 0.95) {
  
  ##TK add error function if arguments are non-numeric
  ##TK add error if arguments are not all of same length
  ##TK add error for negative values
  
  ## Step 1: Find the WLA for each criteria
  ## Given in function arguments
  
  ## Step 2: Compute a LTA for each 
  ### Compute the log variance (sigma_2) of the effuent data
  sigma_2  <- log(1 + cv^2)
  sigma4_2 <- log(1 + 0.25 * cv^2)
  
  ### Compute the acute and chronic multipliers
  acute_mult <- exp(0.5 * sigma_2 - qnorm(prob_LTA) * sqrt(sigma_2))
  chron_mult <- exp(0.5 * sigma4_2 - qnorm(prob_LTA) * sqrt(sigma4_2))
  
  ### Compute the LTAs
  LTAa <- WLAa * acute_mult
  LTAc <- WLAc * chron_mult
  
  ## Step 3: Find controlling LTA
  LTA <- min(LTAa, LTAc, na.rm = TRUE)
  
  ## Step 4: Compute the MDL for aquatic life
  sigman_2 <- log(1 + (1 / n_samples) * cv^2)
  MDL_mult <- exp(qnorm(percentile_MDL) * sqrt(sigma_2) - 0.5 * sigma_2)
  AML_mult <- exp(qnorm(percentile_AML) * sqrt(sigman_2) - 0.5 * sigma_2)
  MDLaq    <- LTA * MDL_mult
  
  ## Step 5: Compute the MDL for human health
  AMLhh <- WLAhh
  MDL_to_AML_ratio <- MDL_mult / AML_mult
  MDLhh <- AMLhh * MDL_to_AML_ratio
  
  ## Step 6: Select the controlling (minimum) MDL
  MDL <- min(MDLaq, MDLhh, na.rm = TRUE)
  
  
  return(MDL)                        ## EXIT
}


#' @title Compute AML
#' @description Compute average monthly effluent limitation (AML) according to methods in EPA's Technical Support Document for Water Quality-based Toxics Control.
#'
#' @param WLAa Numeric. Acute wasteload allocation (WLA). All WLA units should be identical.
#' @param WLAc Numeric. Chronic wasteload allocation (WLA). All WLA units should be identical.
#' @param WLAhh Numeric. Human health wasteload allocation (WLA). All WLA units should be identical.
#' @param cv Numeric. Coefficient of variation (CV) of effluent data. See cv_adj function.
#' @param n_samples Numeric. Number of sample observations.
#' @param prob_LTA Numeric (fraction). Allowable exceedance probability of the WLA used to estimate long-term average (LTA). 
#' @param percentile_AML Numeric (fraction). Lognormal distribution location for AML.
#'
#' @return
#' @export
#'
#' @examples
calc_AML <- function(WLAa, WLAc, WLAhh, cv, 
                    n_samples = 4, prob_LTA = 0.99, percentile_AML = 0.95) {
  
  ##TK add error function if arguments are non-numeric
  ##TK add error if arguments are not all of same length
  ##TK add error for negative values
  
  ## Step 1: Find the WLA for each criteria
  ## Given in function arguments
  
  ## Step 2: Compute a LTA for each 
  ### Compute the log variance (sigma_2) of the effuent data
  sigma_2  <- log(1 + cv^2)
  sigma4_2 <- log(1 + 0.25 * cv^2)
  
  ### Compute the acute and chronic multipliers
  M_acute <- exp(0.5 * sigma_2 - qnorm(prob_LTA) * sqrt(sigma_2))
  M_chron <- exp(0.5 * sigma4_2 - qnorm(prob_LTA) * sqrt(sigma4_2))
  
  ### Compute the LTAs
  LTAa <- WLAa * M_acute
  LTAc <- WLAc * M_chron
  
  ## Step 3: Find controlling LTA
  LTA <- min(LTAa, LTAc, na.rm = TRUE)
  
  ## Step 4: Compute the aquatic life AML
  ### The log variance is a function of the expected number of samples taken
  ### each month. Typically assume 4 for most toxics unless sampling is more
  ### frequent (e.g., daily = 30).
  sigman_2 <- log(1 + (1 / n_samples) * cv^2)
  AMLaq <- LTA * exp(qnorm(percentile_AML) * sqrt(sigman_2) - 0.5 * sigma_2)
  
  ## Step 5: Compute the human health AML
  AMLhh <- WLAhh
  
  ## Step 6: Return the controlling (minimum) AML
  AML <- min(AMLaq, AMLhh, na.rm = TRUE)
  
  return(AML)                        ## EXIT
}
