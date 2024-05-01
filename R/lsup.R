library(httr)
library(readr)
library(janitor)
library(dplyr)

#' Calculate Maximum Likelihood Estimate (MLE) for Bernoulli Distribution
#'
#' logLikBernoulli function calculates the parameter p that maximizes the log-likelihood for a Bernoulli distribution.
#'
#' @param data A vector of binary outcomes (0s and 1s).
#' @return The parameter p that maximizes the log-likelihood.
#' @examples
#' data <- c(1, 0, 0, 0, 1, 1, 1)
#' logLikBernoulli(data)
#' @export
logLikBernoulli <- function(data) {
  ll_hood <- function(p) {
    sum(data * log(p) + (1 - data) * log(1 - p))
  }
  probs <- seq(0, 1, by = 0.001)
  ll_hood_val <- sapply(probs, ll_hood)
  probs[which.max(ll_hood_val)]
}


#'
#' survCurv function calculates and plots a survival curve S(t) based on the provided status and time vectors.
#'
#' @param status A numerical vector representing censoring status (0 for censored, 1 for event).
#' @param time A numerical vector representing survival times.
#' @return Returns Kaplan-Meir survival curve plot.
#' @examples
#' status <- c(1, 1, 0, 1, 0, 1, 0)
#' time <- c(10, 15, 20, 25, 30, 35, 40)
#' survCurv(status, time)
#' @export
survCurv <- function(status, time) {
  data = data.frame(time, status)
  sorted_data = data[order(data$time, decreasing = F), ]

  n <- nrow(sorted_data)
  surv_prob <- rep(1, n)
  events <- sorted_data$status == 0

  for (i in 1:n) {
    if (events[i]) {
      surv_prob[i:n] <- surv_prob[i:n] * (1 - sum(events[1:(i-1)]) / (n - i + 1))
    }
  }
  plot(sorted_data$time, surv_prob, type = "s", ylim = c(0, 1),
       main = "Kaplan-Meier Survival Curve",
       xlab = "Time", ylab = "Survival Probability", col = "blue")
}


#' Reverse Centering and Scaling
#'
#' This function takes a vector that has been scaled and reverses the centering/scaling.
#'
#' @param x A numeric vector that has been scaled using 'scale' function.
#' @return The unscaled vector.
#' @examples
#' scaled_vec <- scale(c(1, 2, 3, 4, 5))
#' unscaled_vec <- unscale(scaled_vec)
#' unscaled_vec
#' @export
unscale <- function(x) {
  unscaled <- (x * attr(x, "scaled:scale")) + attr(x, "scaled:center")
  return(unscaled)
}

#' PC Approximation
#'
#' This function returns an approximation to the data `x` based on a specified number of principal components (`npc`).
#' The approximation is rescaled and centered to match the original data.
#'
#' @param x A numeric matrix or data frame representing the data.
#' @param npc Number of principal components to use for approximation.
#' @return An approximation of the data based on `npc` principal components.
#' @examples
#' # Create sample data matrix
#' set.seed(123)
#' data_matrix <- matrix(rnorm(100), nrow = 10)
#' # Approximate data using 3 principal components
#' approx_data <- pcApprox(data_matrix, 3)
#' approx_data
#' @export
pcApprox <- function(x, npc) {
  if (!is.numeric(x) && !is.matrix(x) && !is.data.frame(x)) {
    stop("Input 'x' must be a numeric matrix or data frame.")
  }
  if (npc >= ncol(x)) {
    stop("# of principal components is less than the number of variables in the data.")
  }

  pca_res <- prcomp(x, center = TRUE, scale. = TRUE)
  pc_app <- pca_res$x[, 1:npc]

  app_x <- pc_app %*% t(pca_res$rotation[, 1:npc])
  app_x_scale <- scale(app_x, center = TRUE, scale = apply(x, 2, sd))

  return(app_x_scale)
}


#' Standardize Variable Names in Tibble Data
#'
#' This function standardizes the variable names in a tibble data frame to a specified 'snake' case or 'small_camel' case using janitor and dplyr packages.
#'
#' @param data A tibble data frame.
#' @param case_ The case style to convert the variable names to (default is "snake").
#' @return A tibble data frame with standardized variable names.
#' @import dplyr
#' @import janitor
#' @import tibble
#' @examples
#' data <- tibble(`First Name` = c("Lucas", "Sean"), `Last Name` = c("McKay", "Jordan"))
#' # snake case
#' standardized_data <- standardizeNames(data, case_ = "snake")
#' standardized_data
#' # small_camel
#' standardized_data <- standardizeNames(data, case_ = "small_camel")
#' standardized_data
#'
#' @export
standardizeNames <- function(data, case_ = "snake") {

  if (!is.data.frame(data) || !inherits(data, "tbl_df")) {
    stop("Input 'data' must be a tibble data frame.")
  }

  if (case_ != "snake" & case_ != "small_camel"){
    stop("Selected case should be either 'small_camel' or 'snake'.")
  }

  clean_names <- make_clean_names(colnames(data), case = case_)
  renamed_data <- rename_with(data, .cols = everything(), .fn = ~ clean_names[match(.x, colnames(data))])

  return(renamed_data)

}

#' Minimum Sample Size Calculation for T-Test
#'
#' This function calculates the minimum sample size needed for a t-test with power 80% and significance level 0.05 (alpha).
#'
#' @param x1 Numeric vector representing the first sample.
#' @param x2 Numeric vector representing the second sample (optional).
#' @param power Desired statistical power (default is 0.8).
#' @param sig_level Significance level (default is 0.05).
#' @return Minimum sample size needed for the t-test.
#' @import pwr
#' @import fishmethods
#' @examples
#' # Calculate minimum sample size for one-sample t-test
#' min_n_1 <- minimumN(x1 = rnorm(30))
#' min_n_1
#' # Calculate minimum sample size for two-sample t-test
#' min_n_2 <- minimumN(x1 = rnorm(30), x2 = rnorm(30))
#' min_n_2
#' @export
minimumN <- function(x1, x2 = NULL, power = 0.8, sig_level = 0.05) {
  z_alpha2 <- qnorm(sig_level/2)
  z_beta <- qnorm(power, lower.tail = FALSE)
  if (is.null(x2)) {
    result <- pwr::pwr.t.test(d = abs(mean(x1) - 0) / sd(x1), sig.level = sig_level, power = power, type = "one.sample")
    min_n <- z_alpha2*z_beta * (sd(x1)/mean(x1))^2
  }
  else {
    xbars <- c(mean(x1), mean(x2))
    sds <-c(var(x1), var(x2))
    n <- c(length(x1), length(x2))

    var_x1x2 <- combinevar(xbars,sds,n)

    min_n <- 2*(z_alpha2 + z_beta)^2 * (sqrt(var_x1x2[2])/(mean(x1) - mean(x2)))^2
  }

  return(min_n)
}

#' Download Redcap Report
#'
#' Using the provided API token, URL, and report ID, this function queries a Redcap
#' server to download a specific report and returns the contents as a tibble.
#'
#' @param redcapTokenName The name of the environment variable containing the Redcap API token.
#' @param redcapUrl The URL of the Redcap server.
#' @param redcapReportId The ID of the Redcap report to download.
#'
#' @return A tibble containing the downloaded Redcap report data.
#' @import httr
#' @import readr
#'
#' @export
downloadRedcapReport <- function(redcapTokenName, redcapUrl, redcapReportId) {

  token <-redcapTokenName
  url <- redcapUrl
  formData <- list("token"=token,
                   content='report',
                   format='csv',
                   report_id=redcapReportId,
                   csvDelimiter='',
                   rawOrLabel='raw',
                   rawOrLabelHeaders='raw',
                   exportCheckboxLabel='false',
                   returnFormat='csv'
  )

  response <- POST(url, body = formData, encode = "form")

  if (http_status(response)$category != "Success") {
    stop("Error fetching data from Redcap.")
  }

  data <- read_csv(content(response, as = "text"))

  return(data)
}
