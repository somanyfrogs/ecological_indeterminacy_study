#' @title Figure_5.R
#' @description Making figure 5 for the Paper entitled:
#'      "Unravaling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20230307 by K.Kawatsu.
#'      Last udpate: 20240409.

## Load R functions
source("R/functions.R")

#' Preparation of matrices A, E, Di and Phi
func_tmp <- function(s, c, sigma, mu, rho, seed = NULL) {
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

    ## Make matrices A, E, Di and Phi
    vars <- calc_vars(s, c, sigma, mu, rho)
    A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
    E <- A - mu * c - (vars$d - mu * c) * diag(s)
    Di <- ((vars$d + mu * (s - 1) * c) * diag(s) - mu * c * matrix(1, s, s)) / (vars$psi * (vars$psi - mu * c))
    Phi <- (E - (mu * c / vars$psi) * matrix(1, s, s) %*% E) / (vars$d - mu * c)

    return(list(A = A, E = E, Di = Di, Phi = Phi))
}

## Figure 1A: Demonstration of press perturbation dynamics close to different eeigenvectors
## Example: An IGP module with ecological indeterminacy (Case 1 for Figure 1A)
## Set the parameters
s <- 3
c <- 1.0
A <- rbind(c(0.00, -0.90, -0.20),
           c(0.60,  0.00, -0.90),
           c(0.10,  0.60,  0.00))
sigma <- sd(offdiag(A))
mu <- mean(offdiag(A))
rho <- cor(offdiag(A, 2), offdiag(t(A), 2))

## Make matrices A, E, Di, and Phi
vars <- calc_vars(s, c, sigma, mu, rho)
A <- A + vars$d * diag(s)
E <- A - mu * c - (vars$d - mu * c) * diag(s)
Di <- ((vars$d + mu * (s - 1) * c) * diag(s) - mu * c * matrix(1, s, s)) / (vars$psi * (vars$psi - mu * c))
Phi <- (E - (mu * c / vars$psi) * matrix(1, s, s) %*% E) / (vars$d - mu * c)

## Calculate eigenvalues and eigenvectors of Phi
eig <- eigen(Phi)
eval <- eig$values
evec <- eig$vectors

## Intrinsic growth rate to ensure the equal equilibrium density xe = 1/s
r <- -A |> rowMeans()

## Press perturbations close to the 1st (|lambda| = 1.0029) and 3rd (|lambda| = 0.0066) eigenvectors
set.seed(123)
delta <- runif(s, 0, 2 * pi) |> (\(theta) 2.5e-1 * complex(real = cos(theta), imaginary = sin(theta)))()
bs <- list(1e-3 * (evec[, 1] + delta), 1e-3 * (evec[, 3] + delta))

## Make true perturbation dynamics
foreach(i = 1:length(bs), .combine = bind_rows) %do% {
    odeRK4(x_ini = rep(1 / s, s), fun = glv, params = list(A = A, r = r, b = bs[[i]]), t_end = 1000, h = 0.01) |>
        mutate(across(everything(), abs), run = i, .before = everything())
} |> write_csvr("output/fig05A.csv")

## Make figure 5A
## Convert 3D data to 2D data with viewing transformation
theta <- 135; phi <- 45
pmat <- persp(z = diag(2), theta = theta, phi = phi)
orig <- tibble(x1 = 1/s, x2 = 1/s, x3 = 1/s) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat))

## 2D transformation of perturbation dynamics
tbl1 <- read_csvr("output/fig05A.csv") |> filter(run == 1) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
tbl2 <- read_csvr("output/fig05A.csv") |> filter(run == 2) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
tbl <- bind_rows(tbl1, tbl2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

## 2D transformation of eigenvectors
vec1 <- tibble(x1 = 1/s + c(0, 1e-3 * evec[1, 1]), x2 = 1/s + c(0, 1e-3 * evec[2, 1]), x3 = 1/s + c(0, 1e-3 * evec[3, 1])) |>
    mutate(across(everything(), Re)) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
vec2 <- tibble(x1 = 1/s + c(0, 1e-3 * evec[1, 3]), x2 = 1/s + c(0, 1e-3 * evec[2, 3]), x3 = 1/s + c(0, 1e-3 * evec[3, 3])) |>
    mutate(across(everything(), Re)) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
vec_tbl <- bind_rows(vec1, vec2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

## 2D transformation of perturbation vectors
bs1 <- tibble(x1 = 1/s + c(0, bs[[1]][1]) / 3, x2 = 1/s + c(0, bs[[1]][2]) / 3, x3 = 1/s + c(0, bs[[1]][3]) / 3) |> mutate(across(everything(), Re)) |>
    with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
bs2 <- tibble(x1 = 1/s + c(0, bs[[2]][1]) / 3, x2 = 1/s + c(0, bs[[2]][2]) / 3, x3 = 1/s + c(0, bs[[2]][3]) / 3) |> mutate(across(everything(), Re)) |>
    with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
bs_tbl <- bind_rows(bs1, bs2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

## 2D transformation of after perturbation equilibirum calculated with the true inverse interaction matrix A^-1
S_orig <- solve(A)
out_t1 <- abs(matrix(1/s, nrow = s) - S_orig %*% bs[[1]]) |> t() |> as_tibbler(str_c("x", 1:s)) |> with(trans3d(x = x1, y = x2, z = x3, pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
out_t2 <- abs(matrix(1/s, nrow = s) - S_orig %*% bs[[2]]) |> t() |> as_tibbler(str_c("x", 1:s)) |> with(trans3d(x = x1, y = x2, z = x3, pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
out_t <- bind_rows(out_t1, out_t2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

## 2D transformation of after perturbation outcomes calculated with the approximated inverse interaction matrix for the order k = 2, 3, 5
out_a1 <- foreach(k = c(2, 3, 5), .combine = rbind) %do% {abs(matrix(1/s, nrow = s) - (matSeries(-Phi, k) %*% Di) %*% bs[[1]]) |> t() |> as_tibbler(str_c("x", 1:s)) |>
    with(trans3d(x = x1, y = x2, z = x3, pmat)) |> (\(.lst) tibble(run = 1, order = k, x = .lst$x, y = .lst$y))()}
out_a2 <- foreach(k = c(2, 3, 5), .combine = rbind) %do% {abs(matrix(1/s, nrow = s) - (matSeries(-Phi, k) %*% Di) %*% bs[[2]]) |> t() |> as_tibbler(str_c("x", 1:s)) |>
    with(trans3d(x = x1, y = x2, z = x3, pmat)) |> (\(.lst) tibble(run = 2, order = k, x = .lst$x, y = .lst$y))()}
out_a <- bind_rows(out_a1, out_a2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

gp1 <- tbl |> ggplot(aes(x = x, y = y)) + facet_wrap(. ~ run, labeller = label_parsed) +
    annotate("segment", x = rep(0, 3), y = rep(0, 3) + 8e-5, xend = c(-4e-4, 0, 4e-4), yend = c(-3e-4, 4.5e-4, -3e-4), linewidth = 0.25, color = "grey50") +
    annotate("text", x = -3.6e-4, y = -2.0e-4, label = expression(italic(x)[2]), size = 2) + annotate("text", x = 0.5e-4, y = 4.3e-4, label = expression(italic(x)[3]), size = 2) + annotate("text", x = 3.6e-4, y = -2.0e-4, label = expression(italic(x)[1]), size = 2) +
    geom_path(data = vec_tbl, color = "orange") + geom_path(data = bs_tbl, linewidth = 0.5, arrow = arrow(length = unit(2, "mm"))) + geom_path() +
    annotate("point", x = 0, y = 0, shape = 21, color = "black", fill = "black") +
    geom_point(data = out_t, shape = 25, fill = "black") + geom_point(data = out_a, aes(fill = factor(order)), shape = 25) +
    ggsci::scale_fill_npg() + xlim(-4e-4, 4e-4) + ylim(-4e-4, 4.5e-4) + labs(tag = "A", fill = expression(italic(k))) + theme_st(lunit = 2) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.margin = margin(0, 0, 0, 0), panel.grid = element_blank())

## Figure 5B: Prediction skill of perturbation outcomes
## This demonstration illustrates the comparision of prediction skill between different complex networks with magnitude of Rmax > 1.
## Case 1: All eigenvalues follow the elliptic law.
## Case 2: Only outlier eigenvalues exceed Rmax > 1.

## Set the common parameters
s <- 250
seed <- 941
lim <- 250

## Case 1: All eigenvalues follow the elliptic law.
c <- 0.843
sigma <- 0.153
mu <- -0.0872
rho <- -0.530
vars <- calc_vars(s, c, sigma, mu, rho)

## Make theoretical and empirical PDF of absolute eigenvalues
Rmaj <- (1 + abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - mu * c)
dns1 <- tibble(case = 1, abs = seq(0, Rmaj, length.out = 250)) |> mutate(density = pdf_abs_eigs(abs, s, c, sigma, mu, rho))
eig1 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 100, seed = seed) |> transmute(case = 1, iter = iter, abs = sqrt(real^2 + imag^2))

## Make example matrices of A, E, Di, and Phi
set.seed(seed)

A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
E <- A - mu * c - (vars$d - mu * c) * diag(s)
Di <- ((vars$d + mu * (s - 1) * c) * diag(s) - mu * c * matrix(1, s, s)) / (vars$psi * (vars$d - mu * c))
Phi <- (E - (mu * c / vars$psi) * matrix(1, s, s) %*% E) / (vars$d - mu * c)

S_orig <- solve(A)

prd1 <- foreach(k = seq(2, 10, 2), .combine = bind_rows) %do% {
    S_appr <- matSeries(-Phi, k) %*% Di

    foreach(i = 1:lim, .combine = bind_rows) %do% {
        u <- rnorm(s) |> (\(v) sign(v) * sqrt(v^2 / sum(v^2)))() |> matrix(nrow = s)
        xo <- as.vector(S_orig %*% u)
        xa <- as.vector(S_appr %*% u)
        tibble(case = 1, order = k, run = i, orig = xo, appr = xa)
    }
}

## Case 2: Only outlier eigenvalues exceed the threshold
c <- 0.542
sigma <- 0.43
mu <- 0.0761
rho <- -0.90
vars <- calc_vars(s, c, sigma, mu, rho)

## Make theoretical and empirical PDF of absolute eigenvalues
Rmaj <- (1 + abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - mu * c)
dns2 <- tibble(case = 2, abs = seq(0, Rmaj, length.out = 250)) |> mutate(density = pdf_abs_eigs(abs, s, c, sigma, mu, rho))
eig2 <- func_eig(s = s, c = c, sigma = sigma, mu = mu, rho = rho, replicate = 100, seed = seed) |> transmute(case = 2, iter = iter, abs = sqrt(real^2 + imag^2))

## Make example matrices of A, E, Di, and Phi
set.seed(seed)
A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
E <- A - mu * c - (vars$d - mu * c) * diag(s)
Di <- ((vars$d + mu * (s - 1) * c) * diag(s) - mu * c * matrix(1, s, s)) / (vars$psi * (vars$d - mu * c))
Phi <- (E - (mu * c / vars$psi) * matrix(1, s, s) %*% E) / (vars$d - mu * c)

S_orig <- solve(A)

prd2 <- foreach(k = seq(2, 10, 2), .combine = bind_rows) %do% {
    S_appr <- matSeries(-Phi, k) %*% Di

    foreach(i = 1:lim, .combine = bind_rows) %do% {
        u <- rnorm(s) |> (\(v) sign(v) * sqrt(v^2 / sum(v^2)))() |> matrix(nrow = s)
        xo <- as.vector(S_orig %*% u)
        xa <- as.vector(S_appr %*% u)
        tibble(case = 2, order = k, run = i, orig = xo, appr = xa)
    }
}

dns <- bind_rows(dns1, dns2) |> mutate(case = str_c("Case: ", case))
eig <- bind_rows(eig1, eig2) |> mutate(case = str_c("Case: ", case))
prd <- bind_rows(prd1, prd2)

gp2 <- eig |> ggplot(aes(x = abs)) + facet_wrap(. ~ case, nrow = 2) + geom_histogram(aes(y = after_stat(density)), bins = 25, alpha = 0.5) +
    geom_line(data = dns, aes(y = density), linewidth = 0.5, color = ggsci::pal_npg()(4)[4]) +
    geom_vline(xintercept = 1, linewidth = 0.25, color = ggsci::pal_npg()(4)[1]) +
    scale_y_continuous(breaks = seq(0, 2, 1)) + xlab(expression(paste("Eigenvalues ", group("|", lambda(Phi), "|")))) +
    ylab("Density") + labs(tag = "B") + theme_st(lunit = 2)

gp3 <- prd |> group_by(case, order, run) |> summarize(MAE = mean(abs(orig - appr)), .groups = "drop") |> ggplot(aes(x = order, y = log(MAE), color = factor(case))) +
    ggbeeswarm::geom_quasirandom(method = "pseudorandom", shape = 16, size = 0.5, alpha = 0.25, fill = "grey25") +
    geom_smooth(method = "glm", linewidth = 0.5) + ggsci::scale_color_npg() +
    scale_x_continuous(breaks = seq(2, 10, 2)) + labs(color = "Case") + xlab(expression(paste("Approximation order ", italic(k)))) + ylab("Prediction error log(MAE)") + theme_st()

## Figure 5C: Relationship between prediction skill of perturbation outcomes and relative volume
## Simulations are conducted under different parameter set (c, sigma, mu, rho) to ensure the same magnitude of Rmax = 1.2, 1.3, and 1.4
set.seed(941)
s <- 250
lim <- 100

## Make parameter set satisfying the same theoretical criterion
params <- foreach(Rmax = c(1.2, 1.3, 1.4), .combine = bind_rows) %do% {
    tbl <- tibble()

    repeat {
        c <- runif(1, 0, 1)
        mu <- runif(1, -1, 1)
        rho <- runif(1, -1, 1)
        tmp <- optimize(func_cost, interval = c(0.0, 0.5), s = s, c = c, mu = mu, rho = rho, Rmax = Rmax)

        if(tmp$objective <= 1e-9) tbl <- bind_rows(tbl, tibble(Rmax = Rmax, s = s, c = c, sigma = tmp$minimum, mu = mu, rho = rho))
        if(nrow(tbl) >= lim) break
    }

    tbl |> mutate(run = 1:lim, seed = sample(32768, lim, replace = TRUE), .before = everything())
}

system.time(foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
    cat(sprintf("(Rmax, run) = (%.1f, %d): Start\n", tbl$Rmax, tbl$run))

    ## Make matrices A, E, Di and Phi
    tmp <- func_tmp(s = tbl$s, c = tbl$c, sigma = tbl$sigma, mu = tbl$mu, rho = tbl$rho, seed = tbl$seed)
    eval <- tmp$Phi |> eigen(only.values = TRUE) |> pullist(1) |> abs()

    ## Calculate theoretical and observed value of the relative volume of eigenspace over criterion
    obs <- sum(eval[eval > 1]) / sum(eval)
    prd <- vol_over_criterion(tbl$s, tbl$c, tbl$sigma, tbl$mu, tbl$rho)

    ## Calculate true and approximated inverse matrix
    S_orig <- solve(tmp$A)
    S_appr <- matSeries(tmp$Phi, 5) %*% tmp$Di

    Rho <- foreach(i = 1:250, .combine = c) %do% {
        ## Make a random perturbation vectors u with ||u|| = 1
        u <- rnorm(s) |> (\(v) sign(v) * sqrt(v^2 / sum(v^2)))() |> matrix(nrow = s)
        xo <- as.vector(S_orig %*% u)
        xa <- as.vector(S_appr %*% u)
        abs(cor(xo, xa))
    } |> mean()

    tibble(run = tbl$run, Rmax = tbl$Rmax, Robs = max(eval), size = tbl$s, connectance = tbl$c, sigma = tbl$sigma, mu = tbl$mu, obs = obs, pred = prd, Rho = Rho)
} |> write_csvr("output/sim_5.csv"))

## Make figure 5C
tbl <- read_csvr("output/sim_5.csv")
gp4 <- tbl |> ggplot(aes(x = obs, y = pred, color = factor(Rmax))) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_point(shape = 16, size = 1.0, alpha = 0.5) + geom_abline(slope = 1, intercept = 0, linetype = 2) + ggsci::scale_color_npg() +
    ylab(expression(paste("Sampled volume ", group("|", lambda(Phi), "|")>1))) +
    labs(tag = "C", color = expression(italic(R)[max])) + theme_st() + theme(axis.title.x = element_blank())

gp5 <- tbl |> ggplot(aes(x = obs, y = Rho, color = factor(Rmax))) + geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) +
    geom_point(shape = 16, size = 1.0, alpha = 0.5) + geom_smooth(linewidth = 0.25, method = "glm", method.args = list(family = Gamma)) + ggsci::scale_color_npg(guide = "none") +
    labs(color = expression(italic(R)[max])) + xlab(expression(paste("Theoretical volume ", Xi))) + ylab(expression(paste("Prediction skill ", rho[5]))) + theme_st()

gp_l <- gp1 / wrap_plots(gp2 | gp3)
gp_r <- (gp4 / gp5) + plot_layout(guide = "collect")
((gp_l | gp_r) + plot_layout(width = c(3, 2))) |> ggsaver("fig05", width = 17.8, height = 10, ext = "pdf")
