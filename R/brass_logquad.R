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

    # iterate and get final results
    brass_result =
      brass_iterate( PDx5,
                     P1, P2, P3,
                     tsurv,
                     fertility_schedules,
                     k, q5 = q5_init,
                     r, Rd,
                     tol )

    # Use q(a) to retrieve U5MR curves using logquad5q0 package
    # We use q(2), q(3), q(5) only! and the selected k
    brass_logquad_result =
      do.call(
        rbind,
        lapply( c( 2, 3, 5 ),
                function( this_a ){
                  tmp_brass = brass_result[ brass_result$a == this_a, ]
                  tmp_logquad_input = logquad5q0::format_data(
                    rate      = tmp_brass$qa_hat,
                    lower_age = 0,
                    upper_age = 365.25*this_a,
                    type      = "qx",
                    sex       = "total" )
                  tmp_logquad_res = logquad5q0::lagrange5q0( data = tmp_logquad_input, k = k )$predictions

                  brasslq_out = data.frame( ref_date = tmp_brass$ref_date,
                                            tmp_logquad_res )
                  return( brasslq_out )
                })
      )

    brass_list =
      list(
        brasslq_iterative = brass_result,
        u5mr_schedule = brass_logquad_result
      )

    return( brass_list )
  }
