#' @title func_main.R
#' @description An R file describing functions of main analysis for the Paper enetitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221219 by K.Kawatsu.
#'      Last update: 20240322 by K.Kawatsu.

## Load dependent packages
library(doParallel)
library(foreach)
library(ggraph)
library(patchwork)
library(tidygraph)
library(tidyverse)
library(viridis)

setwd("set the path of working directory")
source("R/func_matrix.R")
source("R/func_manage.R")

#' Function to calcluate the variables of sensitivity essence
#'
#' @param s An integer, system size.
#' @param c A double, connectance (link density).
#' @param sigma A double, the variance of off-diagonal entries of interaction matrix.
#' @param mu A double, the mean of off-diagonal entries of interaction matrix.
#' @param rho A double, the correlation of off-diagonal entries of interaction matrix.
calc_vars <- function(s, c, sigma, mu, rho, tbl = TRUE) {
    Va <- c * (sigma^2 + (1 - c) * mu^2)
    Pa <- (rho * sigma^2 + (1 - c) * mu^2) / (sigma^2 + (1 - c) * mu^2)

    if(is.infinite(s)) {
        d <- if_else(mu <= 0.0, -(1 + Pa) * sqrt(s * Va), -Inf)
    } else {
        d <- -(1 + max((1 + Pa) * sqrt(s * Va) - mu * c, mu * (s - 1) * c))
    }

    psi <- d + mu * c * (s - 1)

    if(tbl) {
        return(tibble(Va = Va, Pa = Pa, d = d, psi = psi))
    } else {
        return(c(Va = Va, Pa = Pa, d = d, psi = psi))
    }
}

#' Function to calculate the theoretical indeteRminacy criterion for finite-size networks
#'
#' @inheritParams calc_vars
criterion_fin <- function(s, c, sigma, mu, rho) {
    ## Specify the model variables
    vars <- calc_vars(s, c, sigma, mu, rho)
    Va <- vars$Va
    Pa <- vars$Pa
    d <- vars$d
    psi <- vars$psi

    ## Calculate theoretical maximum of the absolute eigenvalues
    Rmaj <- (1 + abs(Pa)) / abs(d - mu * c) * sqrt(s * Va)
    Rout <- abs(mu * c / psi) * sqrt(s * (s - 1) * Va) / abs(d - mu * c)

    if(mu * Pa < 0) {
        Rout <- Rout * sqrt((1 + Pa) / (2 * pi) + psi / (mu * c) * Pa)
    } else if(mu * Pa == 0) {
        Rout <- Rout * sqrt(2 / pi)
    } else {
        Rout <- Rout * (sqrt((1 + Pa) / (2 * pi)) + sqrt(-psi / (mu * c) * Pa))
    }

    return(max(Rmaj, Rout))
}

#' Function to calculate the theoretical indeteRminacy criterion for infinite-size networks
#'
#' @inheritParams calc_vars
criterion_inf <- function(c, sigma, mu, rho) {
    ## Specify the model variables
    vars <- calc_vars(Inf, c, sigma, mu, rho)
    Va <- vars$Va
    Pa <- vars$Pa

    ## Calculate theoreticadl maximum of the absolute eigenvalues
    if(mu <= 0.0) {
        Rmaj <- (1 + abs(Pa)) / (1 + Pa)
        Rout <- ifelse(Pa == 0.0, 0.0, sqrt(abs(Pa)) / (1 + Pa))
    } else {
        Rmaj <- 0.0
        Rout <- sqrt(Va)

        if(Pa < 0.0) {
            Rout <- Rout * sqrt((1 + Pa) / (2 * pi) - Pa / (mu * c))
        } else if(Pa > 0.0) {
            Rout <- Rout * (sqrt((1 + Pa) / (2 * pi)) + sqrt(Pa / (mu * c)))
        } else {
            Rout <- Rout * sqrt(2 / pi)
        }
    }

    return(max(Rmaj, Rout))
}

#' Cost function to deteRmine the criterion parameter (finite-size networks)
#'
#' @inheritParams criterion_fin
#' @param Rmax A double, criterion value to be calculated by 'criterion_fin'
func_cost <- function(s, c, sigma, mu, rho, Rmax) (Rmax - criterion_fin(s, c, sigma, mu, rho))^2

#' Define function obatine the eigenvalue distribution for sensitivity essence
#'
#' @inheritParams calc_vars
#' @inheritParams criterion_fin
#' @param replicate An integer, specify the replication number.
#' @param seed An integer, set the random seed number.
func_eig <- function(s, c, sigma, mu, rho, replicate, seed) {
    ## Reinstate system seed after simulation
    if(!is.null(seed)) {
        sysSeed <- .GlobalEnv$.Random.seed

        on.exit({
            if(!is.null(sysSeed)) {
                .GlobalEnv$.Random.seed <- sysSeed
            } else {
                rm(".Random.seed", envir = .GlobalEnv)
            }
        })

        set.seed(seed)
    }

    ## Specify the model variables
    vars <- calc_vars(s, c, sigma, mu, rho)
    Va <- vars$Va; Pa <- vars$Pa; d <- vars$d; psi <- vars$psi

    foreach(iter = 1:replicate, .combine = bind_rows) %do% {
        ## Make matrices A, E and Phi
        A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
        E <- A - c * mu - (vars$d - c * mu) * diag(s)
        Phi <- (E - (mu * c / vars$psi) * matrix(1, s, s) %*% E) / (vars$d - mu * c)

        eig <- eigen(Phi, only.values = TRUE)$values
        tibble(iter = iter, size = s, connectance = c, sigma = sigma, mu = mu, rho = rho, Va = vars$Va, Pa = vars$Pa, real = Re(eig), imag = Im(eig))
    }
}

#' Define function to calculate elliptic law
#'
#' @inheritParams criterion_fin
func_law <- function(s, c, sigma, mu, rho) {
    theta <- seq(0, 2 * pi, length.out = 250)
    vars <- calc_vars(s, c, sigma, mu, rho)

    tibble(size = s, connectance = c, sigma = sigma, mu = mu, rho = rho, Va = vars$Va, Pa = vars$Pa,
           real = (1 + vars$Pa) * sqrt(s * vars$Va) * cos(theta) / abs(vars$d - mu * c),
           imag = (1 - vars$Pa) * sqrt(s * vars$Va) * sin(theta) / abs(vars$d - mu * c))
}

#' Define function for the 1st mechanism of outlier eigenvalues
#'
#' @inheritParams criterion_fin
func_out <- function(s, c, sigma, mu, rho) {
    vars <- calc_vars(s, c, sigma, mu, rho)
    tmp <- sqrt(abs(-mu * c / vars$psi * s * (s - 1) * vars$Pa * vars$Va)) / abs(vars$d - mu * c)
    case_when(sign(mu * vars$Pa / vars$psi) > 0 ~ 0 + 1i, sign(mu * vars$Pa / vars$psi) < 0 ~ 1 + 0i, TRUE ~ 0+0i) * tmp
}

#' Make adjacency matrix
#'
#' @param s An integer specifying the network size.
#' @param c A double specifying the connectance probability.
get_adj_mat <- function(s, c) {
    adj <- matrix(0, s, s)
    adj[upper.tri(adj)] <- rbinom(s * (s - 1) / 2, 1, c)
    adj <- adj |> (\(M) M + t(M) + diag(s))()
    return(adj)
}

#' Make interaction matrix
#'
#' @inheritParams calc_vars
#' @param d A double specifying the diagonal entries.
get_int_mat <- function(s, sigma, mu, rho, d) {
    A <- matrix(0, s, s)
    A[upper.tri(A)] <- rnorm(s * (s - 1) / 2, mu, sigma)
    A <- A + apply(A, 1, \(.x) if_else(.x == 0, 0, mu + rho * (.x - mu) + sqrt(1 - rho^2) * rnorm(s, 0, sigma))) + d * diag(s)
    return(A)
}

#' Generalized Lotka-Volterra (GLV) model with press perturbation
#'
#' @param x A numeric vector setting community state.
#' @param params A list containing model parameters
#' A: Interaction matrix.
#' r: Intrinsic growth rate.
#' b: Press perturbation.
glv <- function(x, params) with(params, as.numeric((r + A %*% x + b) * x))

#' ODE solver with 4-th order Runge-Kutta approximation
#'
#' @param x_ini A numeric vector, set the initial states.
#' @param fun A functional object specifying the type of differential equation.
#' @param params A list containing required parameters in 'fun'.
#' @param t_end A double specifying the end time at the termination of calculation.
#' @param h A double, specifying integration time in RK-approximation.
#' @export
odeRK4 <- function(x_ini, fun, params, input = NULL, t_end, h = 0.01) {
    x <- x_ini
    t_seq <- seq(0, t_end, h)
    if(is.null(input)) input <- matrix(0, nrow = length(t_seq), ncol = length(x))

    foreach(i = icount(length(t_seq) - 1), .combine = rbind, .init = x_ini) %do% {
        k1 <- h * fun(x = x, params = c(params, u = list(input[i, ])))
        k2 <- h * fun(x = x + k1 / 2, params = c(params, u = list(input[i, ])))
        k3 <- h * fun(x = x + k2 / 2, params = c(params, u = list(input[i, ])))
        k4 <- h * fun(x = x + k3, params = c(params, u = list(input[i, ])))

        x <- x + (k1 + 2 * k2 + 2 * k3 + k4) / 6
    } |> as_tibbler(str_c("x", 1:length(x))) |> mutate(time = t_seq, .before = everything())
}

#' Function to calculate probability density function of absolute eigenvalues
#'
#' @inheritParams calc_vars
#' @param r A double the value of semi-major radius
pdf_abs_eigs <- function(r, s, c, sigma, mu, rho) {
    ## Specify the model variables, the semi-minor and semi-major axes
    vars <- calc_vars(s, c, sigma, mu, rho)
    Rmin <- (1 - abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - mu * c)
    Rmaj <- (1 + abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - mu * c)
    
    foreach(x = r, .combine = c) %do% {
        if(is.complex(x) | x < 0) {
            NA
        } else if(x == 0) {
            0
        } else if(x < Rmin) {
            2 * x / (Rmin * Rmaj)
        } else {
            2 * x / (Rmin * Rmaj) * (1 - 2 / pi * acos(sqrt((Rmaj^2 - x^2) / (Rmaj^2 - Rmin^2)) * Rmin / x))
        }
    }
}

#' Function to calcluate the relative volume of absolute eigenvalues exceeding the threshold |lambda| > 1
#'
#' @param s An integer, system size.
#' @param c A double, connectance (link density).
#' @param sigma A double, the variance of off-diagonal entries of interaction matrix.
#' @param mu A double, the mean of off-diagonal entries of interaction matrix.
#' @param rho A double, the correlation of off-diagonal pairs of interaction matrix.
vol_over_criterion <- function(s, c, sigma, mu, rho) {
    ## Specify the model variables, the semi-minor and semi-major axes
    vars <- calc_vars(s, c, sigma, mu, rho)
    Rmin <- (1 - abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - mu * c)
    Rmaj <- (1 + abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - mu * c)

    ## Expectation of the absolute eigenvalues within the ellipse
    K <- 1 - (Rmaj / Rmin)^2
    E1 <- 4 * Rmin / (3 * pi) * Carlson::elliptic_E(pi / 2, K)
    E2 <- 0

    ## Expectation of the absolute eigenvalues ranging within (1, Rmaj) if Rmaj >= 1
    if(Rmaj >= 1) {
        ## Specify the expectation variables
        R1 <- sqrt((1 - Rmin^2) / (Rmaj^2 - 1))
        R2 <- Rmin * sqrt((Rmaj^2 - 1) / (Rmaj^2 - Rmin^2))

        E2 <- 4 * Rmin / (3 * pi) * if_else(Rmin >= 1,
                                            Carlson::elliptic_E(pi / 2, K) - pi / (2 * Rmin^2 * Rmaj),
                                            Carlson::elliptic_E(pi / 2, K) - Carlson::elliptic_E(atan(R1), K) - (pi / 2 - acos(R2)) / (Rmin^2 * Rmaj))
    }

    ## Calculate the expectation of absolute eigenvalues
    lout <- abs(mu / (vars$psi * (vars$d - mu * c))) * c * sqrt(2 / pi * s * (s - 1) * (1 + vars$Pa) * vars$Va)

    if(mu * vars$Pa < 0) {
        tmp <- sqrt( mu * c / vars$psi * s * (s - 1) * vars$Pa * vars$Va) / abs(vars$d - mu * c) * (0 + 1i)
        lout <- lout / 2 + c(-tmp, tmp)
    } else if(mu * vars$Pa > 0) {
        tmp <- sqrt(-mu * c / vars$psi * s * (s - 1) * vars$Pa * vars$Va) / abs(vars$d - mu * c) * (1 + 0i)
        lout <- lout / 2 + c(-tmp, tmp)
    }

    lout <- abs(lout)
    if(any(lout > 1) & any(lout > Rmaj)) E2 <- (sum(lout) + (s - length(lout)) * E2) / s
    return(abs(E2 / E1))
}

