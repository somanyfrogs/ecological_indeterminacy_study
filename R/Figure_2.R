#' @title Figure_2.R
#' @description Makfing figure 2 for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221223 by K.Kawatsu.
#'      Last update: 20240408.

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

## Figure 1A: Demonstration of eigenvalue distributions of sensitivity essence Phi.
## Set the common parameters
s <- 250
c <- 1.0
sigma <- 0.1
seed <- 123

## Case 1: mu = -0.1, rho = 0.5
mu <- -0.1
rho <- 0.5

tbl_eig1 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 1, seed = seed) |> mutate(eig_type = "ellipse")
tbl_law1 <- func_law(s = s, c = c, sigma = sigma, mu = mu, rho = rho)

## Case 1: mu = -0.1, rho = 0.0
mu <- -0.1
rho <- 0.0

tbl_eig2 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 1, seed = seed) |> mutate(eig_type = "ellipse")
tbl_law2 <- func_law(s = s, c = c, sigma = sigma, mu = mu, rho = rho)

## Case 3: mu = -0.1, rho = -0.5
mu <- -0.1
rho <- -0.5

tbl_eig3 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 1, seed = seed) |> mutate(eig_type = "ellipse")
tbl_law3 <- func_law(s = s, c = c, sigma = sigma, mu = mu, rho = rho)

## Case 4: mu = 0.1, rho = 0.5
mu <- 0.1
rho <- 0.5
lout <- Re(func_out(s, c, sigma, mu, rho))

tbl_eig4 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 1, seed = seed) |> mutate(eig_type = c(rep("outlier", 2), rep("ellipse", s - 2)))
tbl_law4 <- func_law(s = s, c = c, sigma = sigma, mu = mu, rho = rho)
tbl_out4 <- tbl_eig4 |> reframe(across(1:7, mean), real = sum(real) / 2 + c(-lout, lout), imag = rep(0, 2))

## Case 5: mu = 0.1, rho = 0.0
mu <- 0.1
rho <- 0.0
lout <- 0

tbl_eig5 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 1, seed = 5) |> mutate(eig_type = c(rep("outlier", 1), rep("ellipse", s - 1)))
tbl_law5 <- func_law(s = s, c = c, sigma = sigma, mu = mu, rho = rho)
tbl_out5 <- tbl_eig5 |> reframe(across(1:7, mean), real = sum(real) + lout, imag = 0)

## Case 6: mu = 0.1, rho = -0.5
mu <- 0.1
rho <- -0.5
lout <- Im(func_out(s, c, sigma, mu, rho))

tbl_eig6 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 1, seed = seed) |> mutate(eig_type = c(rep("outlier", 2), rep("ellipse", s - 2)))
tbl_law6 <- func_law(s = s, c = c, sigma = sigma, mu = mu, rho = rho)
tbl_out6 <- tbl_eig6 |> reframe(across(1:7, mean), real = sum(real) / 2, imag = c(-lout, lout))

tbl_eig <- bind_rows(tbl_eig1, tbl_eig2, tbl_eig3, tbl_eig4, tbl_eig5, tbl_eig6) |> mutate(type = sprintf('mu==%.2f', sign(mu) * abs(mu))) |> mutate(type = factor(type, levels = unique(type)))
tbl_law <- bind_rows(tbl_law1, tbl_law2, tbl_law3, tbl_law4, tbl_law5, tbl_law6) |> mutate(type = sprintf('mu==%.2f', sign(mu) * abs(mu))) |> mutate(type = factor(type, levels = unique(type)))
tbl_out <- bind_rows(tbl_out4, tbl_out5, tbl_out6) |> mutate(type = sprintf('mu==%.2f', sign(mu) * abs(mu))) |> mutate(type = factor(type, levels = unique(type)))

gp1 <- tbl_eig |> mutate(rho = factor(rho, levels = c("0.5", "0", "-0.5"))) |>
    ggplot(aes(x = real, y = imag, color = rho)) + facet_wrap(. ~ type, labeller = label_parsed, scales = "free") +
    geom_point(pch = 21, aes(size = eig_type, fill = rho), alpha = 0.5) +
    geom_point(data = mutate(tbl_out, rho = factor(rho, levels = c("0.5", "0", "-0.5"))), pch = 4) +
    geom_path(data = mutate(tbl_law, rho = factor(rho, levels = c("0.5", "0", "-0.5"))), linewidth = 0.25) +
    ggsci::scale_color_npg(guide = "none") + ggsci::scale_fill_npg(label = c("0.5", "0.0", "-0.5")) +
    scale_size_discrete(range = c(0.1, 1.0), guide = "none") + labs(tag = "A", fill = expression(rho)) +
    xlab(expression(paste("Real of ", lambda(Phi)))) + ylab(expression(paste("Imaginary of ", lambda(Phi)))) +
    theme_st(lunit = 2, just = c(1, 0.01), pos = c(1, 0.01)) + 
    theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1, margin = margin(1, 1, 1, 1)), legend.title = element_text(margin = margin(1, 1, 1, 1)), legend.title.position = "left")

## Make figure 2B-D
tbl1 <- read_csvr("output/sim_1.csv")
tbl2 <- read_csvr("output/sim_2.csv")
pal <- ggsci::pal_npg()(2)

gp2 <- tbl1 |> mutate(Rho = abs(Rho)) |> ggplot(aes(x = Rmax, y = Robs)) +
    geom_point(shape = 16, alpha = 0.25) + geom_abline(intercept = 0, slope = 1, color = pal[1], linewidth = 0.25) +
    xlab(expression(paste("Theoretical ", gamma[max]))) +
    ylab(expression(paste("Sampled ", max((group("|", lambda(Phi), "|")))))) + labs(tag = "B") + theme_st()

gp3 <- tbl1 |> mutate(Rho = abs(Rho)) |> ggplot(aes(x = Rmax, y = Rho)) +
    geom_point(shape = 16, alpha = 0.25) + geom_vline(xintercept = 1, color = pal[1], linewidth = 0.25) +
    xlab(expression(paste("Theoretical ", gamma[max]))) +
    ylab(expression(paste("Prediction skill ", italic(rho[k])))) + labs(tag = "C") + theme_st()

gp4 <- tbl2 |> mutate(Rho = abs(Rho), name = sprintf('gamma[max]==%.2f', Rmax)) |>
    ggplot(aes(x = factor(order), y = Rho)) + facet_wrap(. ~ name, labeller = label_parsed) +
    ggbeeswarm::geom_quasirandom(method = "pseudorandom", alpha = 0.25, size = 0.5) +
    stat_summary(geom = "point", fun = "mean", color = pal[2], fill = pal[2], size = 0.5, shape = 23) +
    scale_x_discrete(labels = c(1, "", 3, "", 5, "", 7, "", 9, "")) + scale_y_continuous(limits = c(0, 1), breaks = c(0, 1, 0.5)) +
    xlab(expression(paste("Approximation order ", italic(k)))) + ylab(expression(paste("Prediction skill ", italic(rho[k])))) + labs(tag = "D") + theme_st(lunit = 2)

(gp1 / (gp2 | gp3) / gp4) |> ggsaver("fig02", width = 8.7, height = 13.5, ext = "pdf")
