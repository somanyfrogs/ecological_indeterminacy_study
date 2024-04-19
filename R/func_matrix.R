#' @title func_matrix.R
#' @description An R file describing functions of linear algebra for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221219 by K.Kawatsu.
#'      Last update: 20231227 by K.Kawatsu.

#' Calculate matrix power.
#'
#' @param M A real/complex matrix.
#' @param n A real, power number.
#' @param method A strings/numeric specifying the method for calculation.
#' @export
matPow <- function(M, n, method = "direct") {
    if(any(str_starts(method, c(1, "D", "d")))) {
        tmp <- diag(nrow(M))
        if(n > 0) for(i in 1:n) tmp <- tmp %*% M
    } else if(any(str_starts(method, c(2, "I", "i")))) {
        eig <- eigen(M)
        P <- eig$vectors
        tmp <- P %*% diag(eig$values^n) %*% solve(P)
    } else {
        tmp <- NA
    }

    if(!is.complex(M)) tmp <- Re(tmp)
    return(tmp)
}

#' Calculate k-th order approximation of Neumann series.
#'
#' @param M A real or complex matrix.
#' @param k An integer number specifying the approximation order.
matSeries <- function(M, k) {
    tmp <- diag(nrow(M))

    if(k > 0) for(i in 1:k) tmp <- tmp + matPow(M, i, method = 1)
    return(tmp)
}

#' Take off-diagonal elements.
#'
#' @param M A real/complex square matrix.
#' @param pos A strings/integer specifying elements to be accessed.
#' @export
offdiag <- function(M, pos = "all") {
    case <- case_when(any(str_starts(pos, c(1, "A", "a"))) ~ row(M) != col(M), any(str_starts(pos, c(2, "U", "u"))) ~ upper.tri(M), any(str_starts(pos, c(3, "L", "l"))) ~ lower.tri(M))
    return(M[case])
}

#' Calculate matrix trace
#'
#' @param M A real/complex matrix.
tr <- function(M) M |> diag() |> sum()

