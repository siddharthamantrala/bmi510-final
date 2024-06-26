% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lsup.R
\name{standardizeNames}
\alias{standardizeNames}
\title{Standardize Variable Names in Tibble Data}
\usage{
standardizeNames(data, case_ = "snake")
}
\arguments{
\item{data}{A tibble data frame.}

\item{case_}{The case style to convert the variable names to (default is "snake").}
}
\value{
A tibble data frame with standardized variable names.
}
\description{
This function standardizes the variable names in a tibble data frame to a specified 'snake' case or 'small_camel' case using janitor and dplyr packages.
}
\examples{
data <- tibble(`First Name` = c("Lucas", "Sean"), `Last Name` = c("McKay", "Jordan"))
# snake case
standardized_data <- standardizeNames(data, case_ = "snake")
standardized_data
# small_camel
standardized_data <- standardizeNames(data, case_ = "small_camel")
standardized_data

}
