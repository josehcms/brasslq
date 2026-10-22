#' Run iterative Brass-Log-Quadratic approach
#'
#' The goal of this function is to iterate between the Brass Method and the
#' Log-Quadratic method to obtain a stable estimate for the under-5 mortality rate.
#' The mortality schedule to derive the coefficients (which replaces the model
#' life tables) is derived from the Log-Quadratic Method using q5 and k. Hence,
#' we need to define an initial mortality curve using the extended Log-Quadratic
#' version and these two parameters. We use a q5 initial value of 0.15 and we iterate
#' through the Brass method, getting the q5_hat output from Brass and pluggin-in
#' again into the Log-Quadratic method until we get a q5_hat value that is stable
#' (reaches some user-defined tolerance). With this approach, the only user required
#' mortality input is that one derived from the factor k.
#'
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
#' @param k value of k parameter in the log-quadratic model (only applied to ages below 5)
#' @param r stable population growth rate
#' @param Rd linear rate of mortality decline
#' @param q5_init initial q5 for starting iterative Brass-Log-Quadratic approach
#' @param tol tolerance level for q5 convergence (default = 10^-4) in iteration
#'
#' @return A \code{list} with two data frames: one containing the final iterative Brass results
#' (qa_hat, coefficients), and a second one containing the final dataset corresponding
#' to the age-specific mortality curve below age 5 for the q(5) estimate of the
#' iterative result.
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
#'
#' }
#' @export
#'

brass_logquad =
  function( PDx5,
            P1, P2, P3, tsurv,
            fertility_schedules,
            k = 0,
            r = 0.02, Rd = 0.03,
            q5_init = 0.150,
            tol = 10^-4 ){

    brass_result =
      brass_iterate( PDx5,
                     P1, P2, P3,
                     tsurv,
                     fertility_schedules,
                     k, q5 = q5_init,
                     r, Rd,
                     tol )

    brass_logquad_result =
      data.frame( ref_date = brass_result[ 4, ]$ref_date,
                  qx_from_extended_lq( q5 = brass_result[ 4, ]$qa_hat, k = k ) )

    brass_list =
      list(
        brasslq_iterative = brass_result,
        u5mr_schedule = brass_logquad_result[ 1:22, ]
      )

    return( brass_list )
  }
