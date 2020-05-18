
#' @title Sum detected and non-detect concentrations.
#' @description Sum up congener sample concentrations to create a composite parameter value. Primarily intended as a helper function in fuse_samples.
#'
#' @param qual  A character vector containing non-detect indicator strings, e.g., "<" or "ND". The strings used to indicate censored status can be edited in the "nd" argument.
#' @param result A numeric vector of concentration measurements.
#' @param nd A list indicating all the censoring flags included in the dataset. Defaults to "<", "nd", and "ND".
#'
#' @return dataframe with two columns, qual (character) and result (numeric)
#'
#' @examples 
#' ## Example 1 - Some non-detects
#' result <- runif(5, 1, 10)
#' cen_result <- sample(c("<", ""), length(result), 
#'                      replace = TRUE, prob = c(0.2, 0.8))
#' cen_sum(cen_result, result)
#' 
#' ## Example 2 - All non-detects
#' result <- runif(5, 1, 10)
#' cen_result <- sample(c("<", ""), length(result), 
#'                      replace = TRUE, prob = c(1, 0))
#' cen_sum(cen_result, result)
#' 
cen_sum <- function(qual, result, nd = c("<", "nd", "ND")) {
  
  ## A function for adding up censored and non-cen result values
  ## When there are detected values, the MEC is assumed to be the sum
  ## of all detected values (non-detects are treated as zeroes).
  
  ## Error Checking Section ----------------
  ## Error message is given if the the 'result' argument is not numeric.
  result_class <- class(result)
  if (!is.numeric(result)) stop("result arugment must be numeric. The result
                               argument entered is a ", result_class, " type.",
                                call. = FALSE)
  
  ## Error message is given if the the 'result' and 'qual' arguments are of
  ## unequal lengths.
  if (length(result) != length(qual)) stop("The qual and result arguments must be equal in length. Qual is ", length(qual), " elements long, and result is ", length(result), ".", call. = FALSE)
  
  
  ## Main Body Section ----------------------
  detects <- which(!qual %in% nd)  ## gives index of detects
  nondetects <- which(qual %in% nd)  ## gives index of nondetects
  
  ## When all observations are non-detect, the minimum detection level
  ## is returned as a non-detect for the composite.
  all_nd <- if (length(nondetects) == length(qual)) TRUE else FALSE
  
  if (all_nd == TRUE) {
    
    value <- data.frame(qual = "<", result = min(result, na.rm = TRUE),
                        stringsAsFactors = FALSE) 
    
    return(value)
    
  } else {
    
    value <- data.frame(qual = "", result = sum(result[detects], na.rm = TRUE),
                        stringsAsFactors = FALSE)
    
    return(value) 
    
  }
  
} 


#' @title Combine values for a composite parameter
#' @description Calculate composite parameter concentrations using congener concentrations grouped by a sampling date vector. An example would application would be summing PCB congeners collected on a specific sampling date to produce a total PCBs concentration.
#'
#' @param date_grp A date vector to group the dataset.
#' @param qual  A character vector containing non-detect indicator strings, e.g., "<" or "ND". The strings used to indicate censored status can be edited in the "nd" argument.
#' @param result A numeric vector of concentration measurements.
#' @param nd A list indicating all the censoring flags included in the dataset. Defaults to "<", "nd", and "ND".
#'
#' @return A dataframe containing a column of sampling dates used to group data,
#' a qualifier column (character), and a MEC result column (numeric).
#' @export
#'
#' @examples
#' df <- data.frame(dates = rep(seq.Date(from = as.Date("1982-03-10"),
#'                                       to = as.Date("1982-03-15"),
#'                                       by = 1), 5),
#'                  congeners = sort(rep(LETTERS[1:6], 5)),
#'                  qualifier = sample(c("<", ""), size = 30, 
#'                                     replace = TRUE, prob = c(0.8, 0.2)),
#'                  result = sample(seq(0.1, 0.5, 0.1), 30, replace = TRUE))
#' 
#' fuse_samples(df$dates, df$qualifier, df$result)
#' 
#' 
fuse_samples <- function(date_grp, qual, result, nd = c("<", "nd", "ND")) {
  
  ## Error Handling ---------------------------------------------------
  if (!is.numeric(result)) stop("result arugment must be numeric. The result argument entered is a ", result_class, " type.",
                                        call. = FALSE)
  
  ## Main Function Body -----------------------------------------------
  ## Combine the selected vectors into a dataframe with named variables
  df <- data.frame(date = date_grp, qual = qual, result = result)
  

  ## Split and reaggregate the subsetted dataframe by date
  df_agg <- do.call(rbind, 
                    by(df, df[["date"]], 
                       function(x){data.frame(date = unique(x[["date"]]),
                                              cen_sum(qual = x[["qual"]], 
                                                      result = x[["result"]],
                                                      nd = nd))}))

  
  ## remove row names left over from split and by functions
  row.names(df_agg) <- NULL

  return(df_agg)          ## Exit
  
}



#' @title Merge together observed concentrations and detection limits 
#' @description Merge a results column with a detection limits column by overwriting the censored results values with the corresponding detection limit. 
#'
#' @param qual  A character vector containing non-detect indicator strings, e.g., "<" or "ND". The strings used to indicate censored status can be edited in the "nd" argument.
#' @param result A numeric vector of concentration measurements.
#' @param limit A numeric vecotr of method detection limit/reporting limit values.
#' @param nd A list indicating all the censoring flags included in the dataset. Defaults to "<", "nd", and "ND".
#'
#' @return A numeric vector. 
#' @export
#'
#' @examples
#' result <- runif(5, 1, 10)
#' cen_result <- sample(c("<", ""), length(result), replace = TRUE)
#' method_limit <- rep(0.1, length(result))
#' result <- mrg_cen_cols(cen_result, result, method_limit)
#' print(data.frame(cen_result, result, stringsAsFactors = FALSE))
mrg_cen_cols <- function(qual, result, limit,  
                         nd = c("<", "nd", "ND")){  
  
  ## This function returns a vector that should be copied over the results col. 
  ## It overwrites anything in the result column with a non-detect.
  ## TODO: Consider reworking function to override the overwrite behavior for
  ## situations where a portion of the dataset stores mdl's in results column.
  
  
  ## Error Checking Section ----------------
  ## Check that qual, result, and mdl/rl (limit) columns are the same length.
  n_qual   <- length(qual)
  n_result <- length(result)
  n_limit  <- length(limit)
  n_max    <- max(c(n_qual, n_result, n_limit))

  if (n_qual != n_max | n_result != n_max | n_limit != n_max) stop("The number of elements in qual, result, and limit arguments must be equal in length. Argument qual is ", n_qual, " elements long, result is ", n_result, " elements long, and limit is ", n_limit, " elements long.", call. = FALSE)
    
  ## Error message is given if the the 'result' argument is not numeric.
  result_class <- class(result)
  if (!is.numeric(result)) stop("result arugment must be numeric. The result
                               argument entered is a ", result_class, " type.",
                                call. = FALSE)
  
  
  ## Main Body Section ----------------------
  ## Helper function which creates an anti-union operator (opposite of %in%)
  "%not_in%" <- function(x, y) { !(x %in% y) }
  
  ## Index which observations are detects and which are nondetects
  detects <- which(qual %not_in% nd)             ## gives index of detects
  nondetects <- which(qual %in% nd)              ## gives index of nondetects  

  ## Replace non-detect result elements with corresponding limit values.
  result[nondetects] <- limit[nondetects]
  
  return(result)
}
