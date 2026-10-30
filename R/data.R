#' Extended Log-Quadratic coefficients
#'
#' This dataset contains the set of coefficients used for extending the
#' original Log-Quadratic Model from under-5 up to age 50.
#' @format A data frame (or tibble) with 7 variables:
#' \describe{
#'   \item{sex}{sex (female, male, or total)}
#'   \item{xi}{age index i}
#'   \item{x}{lower bound of age group (in days)}
#'   \item{ax}{ax coefficient}
#'   \item{bx}{bx coefficient}
#'   \item{cx}{cx coefficient}
#'   \item{vx}{vx coefficient; equal to 0 above age 5 years}
#' }
"extlogquad_coeffs"

#' Fertility Schedules Coale-Trussel
#'
#' This dataset contains fertility schedules derived by the Coale-Trussel
#' model.
#' @format A data frame with 4 variables:
#' \describe{
#'   \item{id}{schedule id, there are 795 schedules so id goes form 1 to 795}
#'   \item{tx_subset}{flag used to indicate whether this id should be used for estimating tx (tx_subset = 1)}
#'   \item{x}{lower bound of age group (in years)}
#'   \item{fx}{age-specific fertility rate at age x}
#' }
"fert_scheds_ct"


#' Fertility Schedules United Nations
#'
#' Fertility schedules derived from the UN World Population Prospects 2024 revision.
#' This dataset contains fertility schedules for every single country in the 2024 WPP from
#' 1950 to 2023. Each country-year correspond to one schedule.
#' @format A data frame with 7 variables:
#' \describe{
#'   \item{country}{country name}
#'   \item{country_code}{3-digit country code, from United Nations}
#'   \item{year}{year}
#'   \item{id}{schedule id, each country-year corresponds to a single schedule and there are 17464 schedules so id goes form 1 to 17464}
#'   \item{tx_subset}{flag used to indicate whether this id should be used for estimating tx (tx_subset = 1)}
#'   \item{x}{lower bound of age group (in years)}
#'   \item{fx}{age-specific fertility rate at age x}
#' }
"fert_scheds_un"
