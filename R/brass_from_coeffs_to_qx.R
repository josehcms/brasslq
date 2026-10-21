#' Estimate Brass original qx from Brass coefficients
#'
#' The goal of this function is to estimate the probability of dying at age a
#' for children born from mothers at each age group and the time to which each of
#' these qa's apply. given their reported proportions of deceased children 5PDX,
#' and adjusting the reference period to which each of these probabilities apply.
#'
#' It uses as main inputs PDX (reported proportions of deceased children) and
#' mean parities at ages 15-19 (P1), 20-24 (P2), and 25-29 (P3).
#'
#'
#' @param PDx5 proportion of deceased children by mother's age group
#' @param P1 Mean number of children ever born to women in ages 15-19
#' @param P2 Mean number of children ever born to women in ages 20-24
#' @param P3 Mean number of children ever born to women in ages 25-29
#' @param tsurv Reference time of survey/census
#' @param fertility_schedules a data frame with three columns 'id', which identifies
#' the id of a fertility schedule, 'x' for the lower bound of that age group,
#' and 'fx' for the fertility rate of age group x
#' @param q5 desired level of under-5 mortality rate
#' @param k value of k parameter in the log-quadratic model (only applied to ages below 5)
#' @param r stable population growth rate
#' @param Rd linear rate of mortality decline
#'
#' @return A \code{data.frame} containing the key regression coefficients used
#' in the Brass method to convert the proportion of deceased children by mother's
#' age group 5PDx into probabilities of dying at age a (qa) in a given time
#' tx years before the reference period.
#'
#' The regressions are set as 1) qa/5PDx = a0 + P1/P2 \* a1 + P2/P3 \* a2,
#'  2) tx = b0 + P1/P2 \* b1 + P2/P3 \* b2
#'
#'
#' \describe{
#'   \item{moage_group}{Age group of mother; in five-year groups}
#'   \item{a}{age corresponding to qa_hat estimate}
#'   \item{a0}{intercept coefficient for qa/5PDx regression on P1/P2 and P2/P3}
#'   \item{a1}{coefficient for P1/P2 for qa/5PDx regression on P1/P2 and P2/P3}
#'   \item{a2}{coefficient for P2/P3 for qa/5PDx regression on P1/P2 and P2/P3}
#'   \item{b0}{intercept coefficient for 5tx regression on P1/P2 and P2/P3}
#'   \item{b1}{coefficient for P1/P2 for 5tx regression on P1/P2 and P2/P3}
#'   \item{b2}{coefficient for P2/P3 for 5tx regression on P1/P2 and P2/P3}
#'   \item{PDx_5}{Proportion of deceased children by mother's age group (input)}
#'   \item{qa_hat}{Estimated probability of dying between birth and age a}
#'   \item{ref_date}{time at which the respective qa_hat value applies}
#' }
#' @export
#'

brass_from_coeffs_to_qx =
  function( PDx5,
            P1, P2, P3,
            tsurv,
            fertility_schedules,
            k, q5,
            r, Rd ){


    brass_coeffs = brass_estimate_coeffs( fertility_schedules = fertility_schedules,
                                          q5 = q5,
                                          k = k,
                                          r = r,
                                          Rd = Rd )

    if( nrow( brass_coeffs ) != length( PDx5 ) ){
      warning('Check 5PDx age groups!')
    }

    brass_coeffs$PDx_5 = PDx5

    p1p2 = P1 / P2
    p2p3 = P2 / P3

    brass_coeffs$qa_hat   = brass_coeffs$PDx_5  * ( brass_coeffs$a0  + brass_coeffs$a1 * p1p2 + brass_coeffs$a2 * p2p3 )
    brass_coeffs$ref_date = tsurv - ( brass_coeffs$b0  + brass_coeffs$b1 * p1p2 + brass_coeffs$b2 * p2p3 )

    return( brass_coeffs )
  }
