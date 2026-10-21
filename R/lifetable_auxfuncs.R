#' Convert qx into lx
#'
#' This function converts probabilities of dying at age x qx into number of
#' survivors at age x lx of a life table.
#'
#' @param x vector of ages matching the qx vector
#' @param qx vector of probabilities of dying before age x
#'
#' @return \code{lx} number of individuals surviving up to age x in life table
#' @export

qx_to_lx = function( x, qx ){
  stopifnot( length( x ) == length( qx ) )
  n_qx = length( qx )
  lx = numeric( n_qx )
  lx[1] = 1
  lx[2:n_qx] = 1-qx[1:(n_qx-1)]
  lx
}

#' Convert lx into Lx
#'
#' This function converts the number of survivors at age x lx
#' into the person-years lived between ages x and x+n nLx of a life table.
#' This function assumes constant force of mortality within each age group.
#'
#' @param x vector of ages matching the lx vector
#' @param lx vector of number of individuals surviving up to age x in life table
#'
#' @return \code{Lx} life table person-years lived between ages x and x+n
#' @export
#'
lx_to_Lx = function( x, lx ){
  stopifnot( length( x ) == length( lx ) )
  n_lx = length( lx )
  n = diff(x)
  Lx = rep( NA, n_lx )
  Lx[1:(n_lx-1)] =
    ( n * ( lx[ 2:n_lx ] - lx[ 1:(n_lx-1) ] ) ) /
    ( log( lx[ 2:n_lx ] / lx[ 1:(n_lx-1) ] ) )

  # for last, lx * n/2, assuming n = 1
  Lx[n_lx] = lx[n_lx] * 0.5
  Lx
}

#' Convert qx into Lx
#'
#' This function converts probabilities of dying at age x qx into number of
#' person-years lived between ages x and x+n nLx of a life table.
#'
#' @param x vector of ages matching the qx vector
#' @param qx vector of probabilities of dying before age x
#'
#' @return \code{Lx} person-years lived between ages x and x+n Lx of a life table
#' @export
#'
qx_to_Lx = function( x, qx ){
  lx_to_Lx( x, qx_to_lx( x, qx ) )
}

#' Convert Lx into Sx
#'
#' This function converts person-years lived between ages x and x+n nLx
#' of a life table into the probabilities of surviving from age group
#' x,x+n into age group x+n,x+2n
#'
#' @param x vector of ages matching the Lx vector
#' @param Lx vector of person-years lived between ages x and x+n nLx
#' @param l0 life table radix, set to 1
#'
#' @return \code{Sx} probability of surviving from one age group to the next age group
#' @export
#'
Lx_to_Sx = function( x, Lx, l0 = 1 ){
  n = diff(x)
  n_Lx = length( Lx )
  Sx = rep( NA, n_Lx )
  Sx[ 1 ] = Lx[ 1 ] / ( l0 * n[1] )
  Sx[ 2:n_Lx ] = Lx[ 2:n_Lx ] / Lx[ 1:(n_Lx-1) ]
  Sx
}
