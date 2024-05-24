#' @title Figure_2.R
#' @description Makfing figure 2 for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221223 by K.Kawatsu.
#'      Last update: 20240524.

## Load R functions
source("R/functions.R")

## Simulation for Figure 2B-C: the relationship between theoretical value of criterion and approximation accuracy of netw interaction effects.
## In this simulation, the parameters sigma, mu and rho are randomly assigned from uniform distribution for each run.
set.seed(123)
s <- 250
c <- 0.5

system.time(sim_1 <- foreach(iter = 1:1000, .combine = bind_rows) %do% {
    ## Set the parameters
    sigma <- runif(1, 0.0, 1.0)
    mu <- runif(1, -0.5, 0.5)
    rho <- runif(1, -0.5, 0.5)

    ## Make matrices A, E, Di and Phi
    vars <- calc_vars(s, c, sigma, mu, rho)
    A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
    E <- A - mu * c - (vars$d - mu * c) * diag(s)
    Di <- ((vars$d + mu * c * (s - 1)) * diag(s) - mu * c * matrix(1, s, s)) / (vars$psi * (vars$d - mu * c))
    Phi <- Di %*% E

    Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()
    S_orig <- offdiag(solve(A))
    S_appr <-offdiag(matSeries(-Phi, 5) %*% Di)
    tibble(sigma = sigma, mu = mu, rho = rho, Va = vars$Va, Pa = vars$Pa, Rmax = criterion_fin(s, c, sigma, mu, rho), Robs = Robs, Rho = cor(S_orig, S_appr))
})

sim_1 |> write_csvr("output/sim_1.csv")

## Simulation for Figure 2D
## In this simulation, two of the mean, variance and correlation coefficients of interaction strengths are randomly assigned,
##      Then the remaining one is determined to satisfy the given criterion
set.seed(123)
lim <- 1000
s <- 250
c <- 0.5

system.time(sim_2 <- foreach(Rmax = c(0.74, 0.99, 1.24), .combine = bind_rows) %do% {
    ## Make the possible parameter set
    params <- tibble()

    repeat {
        mu <- runif(1, -0.5, 0.5)
        rho <- runif(1, -0.5, 0.5)
        tmp <- optimize(func_cost, interval = c(0.0, 0.25), s = s, c = c, mu = mu, rho = rho, Rmax = Rmax)

        if(tmp$objective <= 1e-9) params <- bind_rows(params, tibble(sigma = tmp$minimum, mu = mu, rho = rho))
        if(nrow(params) >= lim) break;
    }

    foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
        ## Set the parameters
        sigma <- tbl$sigma
        mu <- tbl$mu
        rho <- tbl$rho

        ## Make matrices A, E, Di and Phi
        vars <- calc_vars(s, c, sigma, mu, rho)
        A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
        E <- A - mu * c - (vars$d - mu * c) * diag(s)
        Di <- ((vars$d + mu * c * (s - 1)) * diag(s) - mu * c * matrix(1, s, s)) / (vars$psi * (vars$d - mu * c))
        Phi <- Di %*% E

        Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()
        S_orig <- offdiag(solve(A))

        foreach(k = 1:10, .combine = bind_rows) %do% {
            S_appr <- offdiag(matSeries(-Phi, k) %*% Di)
            tibble(sigma = sigma, mu = mu, rho = rho, Va = vars$Va, Pa = vars$Pa, Rmax = Rmax, Robs = Robs, order = k, Rho = cor(S_orig, S_appr))
        }
    }
})

sim_2 |> write_csvr("output/sim_2.csv")

## Make figure 2A-C
tbl1 <- read_csvr("output/sim_1.csv")
tbl2 <- read_csvr("output/sim_2.csv")
pal <- ggsci::pal_npg()(2)

gp1 <- tbl1 |> mutate(Rho = abs(Rho)) |> ggplot(aes(x = Rmax, y = Robs)) +
    geom_point(shape = 16, alpha = 0.25) + geom_abline(intercept = 0, slope = 1, color = pal[1], linewidth = 0.25) +
    xlab(expression(paste("Theoretical ", gamma[max]))) +
    ylab(expression(paste("Sampled ", max((group("|", lambda(Phi), "|")))))) + labs(tag = "A") + theme_st()

gp2 <- tbl1 |> mutate(Rho = abs(Rho)) |> ggplot(aes(x = Rmax, y = Rho)) +
    geom_point(shape = 16, alpha = 0.25) + geom_vline(xintercept = 1, color = pal[1], linewidth = 0.25) +
    xlab(expression(paste("Theoretical ", gamma[max]))) +
    ylab(expression(paste("Prediction skill ", italic(rho[k])))) + labs(tag = "B") + theme_st()

gp3 <- tbl2 |> mutate(Rho = abs(Rho), name = sprintf('gamma[max]==%.2f', Rmax)) |>
    ggplot(aes(x = factor(order), y = Rho)) + facet_wrap(. ~ name, labeller = label_parsed) +
    ggbeeswarm::geom_quasirandom(method = "pseudorandom", alpha = 0.25, size = 0.5) +
    stat_summary(geom = "point", fun = "mean", color = pal[2], fill = pal[2], size = 0.5, shape = 23) +
    scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "")) + scale_y_continuous(limits = c(0, 1), breaks = c(0, 1, 0.5)) +
    xlab(expression(paste("Approximation order ", italic(k)))) + ylab(expression(paste("Prediction skill ", italic(rho[k])))) + labs(tag = "C") + theme_st(lunit = 2)

((gp1 | gp2) / gp3) |> ggsaver("fig02", width = 8.7, height = 9, ext = "pdf")
