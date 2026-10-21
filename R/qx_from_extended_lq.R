#' Estimate probability of dying at age a from the extended log-quadratic coefficients
#'
#' This function estimates the probability of dying at age x (qx) using the
#' extended (up to age 50) log-quadratic coefficients.
#'
#' @param sex sex ('male', 'female' or 'total')
#' @param q5 desired level of under-5 mortality rate
#' @param k value of k parameter in the log-quadratic model (only applied to ages below 5)
#'
#' @return A \code{data.frame} containing the 67 predicted values of probabilities
#' of dying between age 0 and x
#'
#' \describe{
#'   \item{x}{age group x, lower bound (in years)}
#'   \item{x_up}{upper bound of age group x (in years)}
#'   \item{x_d}{age group x, lower bound (in days)}
#'   \item{x_d_up}{upper bound of age group x (in days)}
#'   \item{qx}{Probability of dying from birth (age 0) to the end of group x (x_up)}
#'
#' }
#' @export
#'
#'

qx_from_extended_lq =
  function( sex = 'total', q5, k ){
    this_sex = sex
    mortality_lq = subset( extlogquad_coeffs, sex == this_sex )
    x_d_low = mortality_lq$x
    x_d_up   = mortality_lq$x + c( diff( mortality_lq$x ), 365.25 )
    x_low  = mortality_lq$x / 365.25
    x_up  = x_d_up / 365.25
    ax = mortality_lq$ax
    bx = mortality_lq$bx
    cx = mortality_lq$cx
    vx = mortality_lq$vx

    qx = exp( ax + bx * log(q5) + cx* (log(q5))^2 + vx*k )

    out_lq = data.frame( x = x_low, x_up, x_d = x_d_low, x_d_up, qx )

    return( out_lq )

  }
