#' Generate ${_5}{PD}{_x}$ and ${_5}{t}{_x}$
#'
#' The goal of this function is to estimate the independent and dependent
#' variables to be used in the Brass regression model for each pair of
#' mortality and fertility schedules.This function estimates the proportion of
#' children dead by mother's 5-year age groups, ${_5}{PD}{_x}$,
#' and the time shift, ${_5}{t}{_x}$, associated with the linear period
#' mortality decline assumption. Also, it gives the mean parities
#' of age groups 15-19 (P1), 20-24 (P2), and 25-29 (P3).
#'
#' @param ages_qx vector of children's ages matching the qx vector
#' @param qx vector of children's probabilities of dying before age x
#' @param ages_fx vector of ages matching the fx vector
#' @param fx vector of age-specific fertility rates fx
#' @param ages_LxF vector of ages matching the LxF vector
#' @param LxF vector of life table person-years lived between ages x and x+1 for females
#' @param r stable population growth rate
#' @param Rd linear rate of mortality decline
#'
#' @return A \code{data.frame} containing the key regression variables of the Brass method
#' (qa, 5PDx, tx5 and P1, P2, P3) corresponding to the mortality (qx and LxF) and
#' fertility (fx) inputs given.
#'
#' \describe{
#'   \item{moage_group}{Age group of mother; in five-year groups}
#'   \item{a}{Upper age of children to which the probability of dying qa applies}
#'   \item{qa}{Probability of dying between birth (age 0) and age a}
#'   \item{5PDx}{Expected proportion of children deceased among mother's of that age group given the mortality and fertility inputs}
#'   \item{tx5}{Years before reference period to which the period qa value prevails}
#'   \item{P1} Mean number of children ever born to women in ages 15-19
#'   \item{P2} Mean number of children ever born to women in ages 20-24
#'   \item{P3} Mean number of children ever born to women in ages 25-29
#'
#' }
#' @export

generate_PDx5_tx5 = function( ages_qx, qx,
                              ages_fx, fx,
                              ages_LxF, LxF,
                              r,
                              Rd
){
  # set age ranges for analysis
  moage_range = 15:49 # mother's age range for output
  childbearing_age_start = 12 # reproductive age start
  childbearing_age_end = 50   # reproductive age ends

  # get vector lenghts m
  m_qx = length( qx )
  m_LxF = length( LxF )

  # age intervals of children's life table
  n_qx = c( diff( ages_qx ), 1 )

  # stop if lengths of age vectors and of contents do not match
  stopifnot( length( LxF ) == length( ages_LxF ) &
               length( fx ) == length( ages_fx ) &
               length( qx ) == length( ages_qx ) )

  # 1) Estimate PD(x,t) using granularity of logquad data output
  ## 1.A) create q(x,t) matrix applying linear rate of mortality decline
  ## t = 0: time of survey, t = 20: 20 years before survey
  n_t = childbearing_age_end
  qx_decline_trend = 1 + Rd * ( 0:( n_t - 1 ) )
  qxt = qx %*% t( qx_decline_trend )

  ## 1.B) Convert q(x,t) to L(x,t)
  Lxt = apply( qxt, 2, function( qx ){ qx_to_Lx( ages_qx, qx ) } )

  ## 1.C) Evaluate PD(x) = integral( f(a) * q(x-a) da ) / integral( f(a) da ) for each time t
  ##      We will use the granularity of Lx, but we need to approximate f(a) using linear interp
  reproductive_span = childbearing_age_start:childbearing_age_end
  # proportion of children dying for women at age x
  PDxt = matrix( NA_real_, nrow = length( reproductive_span ), ncol = n_t )
  # children ever born to women age x
  CEBx = rep( NA_real_, length = length( reproductive_span ) )
  for( i in 1:length( reproductive_span ) ){

    y = reproductive_span[i]

    # ages in which we need to evaluate fertility
    ages_fa = ( y - ages_qx )
    select_ages = which( ages_fa >= childbearing_age_start )

    # select ages for child survival probs
    ages_qa_start = ages_qx[ select_ages ]
    ages_qa_end = ages_qx[ select_ages + 1 ]
    n_qa = ages_qa_end - ages_qa_start

    # age start and age end of fertility interval
    n_fa = n_qa
    ages_fa_end = ages_fa[ select_ages ]
    ages_fa_start = ages_fa_end - n_fa


    # probability of surviving from birth to age group [a,a+n], l0 =1
    Pat = Lxt[ select_ages, ] / ( 1 * n_qa )

    # get granular fertility values
    fa = approx( x = ages_fx + 0.5,
                 y = fx,
                 xout = ages_fa_start + n_fa/2,
                 rule = 2 )$y

    # calculate number of children born within each age interval (fertility weights)
    w_fa = fa * n_fa

    # calculate number of children ever born until age y
    CEBx[i] = sum( w_fa )

    # proportion of children dying up to age y
    PDxt[i,] = 1 - ( w_fa %*% Pat ) / CEBx[i]

  }

  # 2) Estimate t(x)
  ## We need to convert L(x,t) into survival ratios to estimate a cohort PD(x) = PDc(x)
  ## For that the time intervals should match the age intervals, i.e., 1-year age groups

  ## 2.A) # Agregate age groups to 1-yr intervals for matrix L(x,t)
  index_ages_0_1yr  = which( ages_qx < 1 )
  index_ages_1_2yr  = which( ages_qx >= 1 & ages_qx < 2 )
  index_ages_2_50yr = which( ages_qx >= 2 )

  Lxt_1yr =   rbind( apply( Lxt[ index_ages_0_1yr, ], 2, sum ),
                     apply( Lxt[ index_ages_1_2yr, ], 2, sum ),
                     Lxt[ index_ages_2_50yr, ] )

  ## 2.B) calculate survival ratios S(x,t) = nLx+n / nLx, and S(0,t) = nL0 / n*l0
  Sxt = matrix( NA_real_, nrow = nrow( Lxt_1yr ), ncol = n_t )
  Sxt[ 1, ] = Lxt_1yr[ 1, ] / 1
  Sxt[ 2:nrow( Lxt_1yr ), ] = Lxt_1yr[ 2:nrow( Lxt_1yr ), ] / Lxt_1yr[ 1:(nrow( Lxt_1yr ) - 1), ]

  ## 2.C) calculate cohort survival probabilities Pcx for children born x years before survey
  # under a mortality decline context
  Pcx = rep( NA_real_, childbearing_age_end )
  Pcx[1] = Sxt[ 1, 1 ]
  for( j in 2:childbearing_age_end ){
    anti_diag_indexes = cbind( 1:j, j:1 )
    Pcx[ j ] = prod( Sxt[anti_diag_indexes] )
  }
  ages_Pcx = 1:childbearing_age_end - 1

  ## 2.D) Calculate proportion of children dying for the cohorts of children alive in time 0
  ## assuming linear mortality decline
  PDxc = rep( NA_real_, length = length( reproductive_span ) )
  for( i in 1:length( reproductive_span ) ){

    y = reproductive_span[i]

    # ages in which we need to evaluate fertility
    ages_fa = ( y -  ages_Pcx )
    select_ages = which( ages_fa >= childbearing_age_start )

    # select ages for child survival probs
    ages_Pca_start = ages_Pcx[ select_ages ]
    ages_Pca_end = ages_Pcx[ select_ages + 1 ]
    n_Pca = ages_Pca_end - ages_Pca_start

    # age start and age end of fertility interval
    n_fa = n_Pca
    ages_fa_end = ages_fa[ select_ages ]
    ages_fa_start = ages_fa_end - n_fa


    # probability of surviving from birth to age group [a,a+n], l0 =1
    Pca = Pcx[ select_ages ]

    # get granular fertility values
    fa = approx( x = ages_fx + 0.5,
                 y = fx,
                 xout = ages_fa_start + n_fa/2,
                 rule = 2 )$y

    # calculate number of children born within each age interval (fertility weights)
    w_fa = fa * n_fa

    # calculate number of children ever born until age y
    #CEBy[i] = sum( w_fa ) # we can use this from the previous loop

    # proportion of children dying up to age y
    PDxc[i] = 1 - ( w_fa %*% Pca ) / sum( w_fa )

  }

  ## 2.E) Estimate t(x) for each mother's age - time to which measure of qa will apply
  times = ( 1:n_t - 1 ) + 0.5
  tx = rep( NA, length(moage_range))
  for( i in 1:length(moage_range) ){
    y = moage_range[i]
    j = which( reproductive_span == y )
    tx[ i ] =
      approx( xout = PDxc[ j ],
              x = PDxt[ j, ],
              y = times,
              rule = 2 )$y
  }

  # 3) Calculate 5PDx and 5tx, for mother's age groups [x,x+5]

  ## 3.A) calculate weights 1ca assuming stable population and using women's LxF

  cxF =
    LxF[ which( ages_LxF %in% moage_range ) ] *
    exp( - r * ( ages_LxF[ which( ages_LxF %in% moage_range ) ] + 0.5 ) )
  cxF = cxF / sum( cxF )

  ## 3.B) Calculate 5PDx using cxF weights 5PDx = sum( PDx * cxF )
  PDx = PDxt[,1]
  x5 = unique( moage_range - moage_range %% 5 )
  age_groups = paste0( x5, '-', x5 + 4 )
  PDx5 = sapply( x5,
                 function( x ){
                   x_range = x:(x+4)
                   PDx5 =
                     weighted.mean( x = PDx[ which( reproductive_span %in% x_range ) ],
                                    w = cxF[ which( moage_range %in% x_range ) ] )
                   return( PDx5 )
                 })
  ## 3.C) Calculate 5tx using cxF weights 5tx = sum( tx * cxF )
  tx5 = sapply( x5,
                function( x ){
                  x_range = x:(x+4)
                  tx5 =
                    weighted.mean( x = tx[ which( moage_range %in% x_range )  ],
                                   w = cxF[ which( moage_range %in% x_range ) ] )
                  return( tx5 )
                })

  # 4) Get P1 (P20, cumulating 15-19), P2 (P25) and P3 (P30) from fertility schedules
  #CEBx
  P1 = weighted.mean( x = CEBx[ which( reproductive_span %in% 15:19 ) ],
                      w = cxF[ which( moage_range %in% 15:19 ) ] )
  P2 = weighted.mean( x = CEBx[ which( reproductive_span %in% 20:24 )],
                      w = cxF[ which( moage_range %in% 20:24 ) ] )
  P3 = weighted.mean( x = CEBx[ which( reproductive_span %in% 25:29 ) ],
                      w = cxF[ which( moage_range %in% 25:29 ) ] )

  # 5) get q(a) corresponding to each mother's age group
  map_moage_to_qx =
    data.frame( x = 0,
                n = c( 1, 2, 3, 5, 10, 15, 20 ),
                moage = seq( 15, 45, 5 ) )

  qa = qx[ which( ( ages_qx + n_qx ) %in% map_moage_to_qx$n ) ]
  names( qa ) = map_moage_to_qx$n

  data.frame( moage_group = age_groups,
              a = map_moage_to_qx$n,
              qa, PDx5, tx5, P1, P2, P3 )

}



#' Estimate Brass Coefficients
#'
#' The goal of this function is to estimate the Brass Method's coefficients
#' for converting proportions of deceased children by mother's age into
#' probabilities of dying, and for adjusting the reference period to which each
#' of these probabilities apply.
#'
#' It uses as main inputs q5 (under-5 mortality level) and k - for deriving the
#' reference life table to be used; and the input fertility schedules.
#'
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
#' (qa, 5PDx, tx5 and P1, P2, P3) corresponding to the mortality (qx and LxF) and
#' fertility (fx) inputs given.
#'
#' \describe{
#'   \item{moage_group}{Age group of mother; in five-year groups}
#'   \item{a0}{intercept coefficient for qa/5PDx regression on P1/P2 and P2/P3}
#'   \item{a1}{coefficient for P1/P2 for qa/5PDx regression on P1/P2 and P2/P3}
#'   \item{a2}{coefficient for P2/P3 for qa/5PDx regression on P1/P2 and P2/P3}
#'   \item{b0}{intercept coefficient for 5tx regression on P1/P2 and P2/P3}
#'   \item{b1}{coefficient for P1/P2 for 5tx regression on P1/P2 and P2/P3}
#'   \item{b2}{coefficient for P2/P3 for 5tx regression on P1/P2 and P2/P3}
#'
#' }
#' @export
#'
brass_estimate_coeffs =
  function( fertility_schedules,
            q5,
            k,
            r,
            Rd
  ){

    # map mother's age groups to qa
    map_moage_to_a = c( 1, 2, 3, 5, 10, 15, 20 )
    names( map_moage_to_a ) = paste0( seq( 15, 45, 5 ), '-', seq( 15, 45, 5 )+4 )
    # life table for both sexes for children - derived from logquad
    mlt_total = qx_from_extended_lq( sex = 'total', q5 = q5, k = k )
    ages_qx = mlt_total$x
    qx = mlt_total$qx
    # life table for women - derived from logquad
    mlt_females = qx_from_extended_lq( sex = 'female', q5 = q5, k = k )
    ages_LxF = mlt_females$x
    LxF = qx_to_Lx( x = ages_LxF, qx = mlt_females$qx )

    # fertility schedules' ids
    fert_ids = unique( fertility_schedules$id )
    pb = txtProgressBar( min = 1, max = length(fert_ids), initial = 1 )

    regression_inputs =
      do.call(
        rbind,
        lapply( fert_ids,
                function( fert_id ){
                  # progress bar
                  setTxtProgressBar( pb, fert_id )
                  this_fert_sched = fertility_schedules[ which( fertility_schedules$id == fert_id ), ]
                  fx = this_fert_sched$fx
                  ages_fx = this_fert_sched$x

                  PDx5_tx5_df =
                    data.frame( fert_id = fert_id,
                                generate_PDx5_tx5( ages_qx, qx,
                                                   ages_fx, fx,
                                                   ages_LxF, LxF,
                                                   r,
                                                   Rd ) )

                  return( PDx5_tx5_df )
                }
        )
      )

    # estimate coefficients for each mother's age group
    age_groups = unique( regression_inputs$moage_group )
    regression_inputs$k = regression_inputs$qa / regression_inputs$PDx5
    regression_inputs$P1P2 = regression_inputs$P1 / regression_inputs$P2
    regression_inputs$P2P3 = regression_inputs$P2 / regression_inputs$P3

    brass_coeffs =
      do.call(
        rbind,
        lapply( age_groups,
                function( this_ag ){
                  regression_inputs_ag = regression_inputs[ which( regression_inputs$moage_group == this_ag  ), ]

                  reg_ag_kx = lm( k ~ 1 + P1P2 + P2P3, data = regression_inputs_ag )
                  reg_ag_tx = lm( tx5 ~ 1 + P1P2 + P2P3, data = regression_inputs_ag )

                  brass_coeffs_k = coefficients( reg_ag_kx )
                  names( brass_coeffs_k ) = paste0( 'a', 0:2 )
                  brass_coeffs_t = coefficients( reg_ag_tx )
                  names( brass_coeffs_t ) = paste0( 'b', 0:2 )

                  data.frame( moage_group = this_ag,
                              a = map_moage_to_a[ this_ag ],
                              t( c( brass_coeffs_k, brass_coeffs_t ) ) )
                } )
      )

    return( brass_coeffs )

  }
