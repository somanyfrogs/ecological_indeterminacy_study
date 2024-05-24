#' @title Suppl.R
#' @descritption File for Supplementary Information for the Paper entitled:
#'      "Unravaling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20230518 by K.Kawatsu.
#'      Last update: 20240524.

## Load R functions
source("R/functions.R")

## Figure S1: Simulation analysis investigating the approximation accuracy of Stieltjes transformation
rslv <- function(E, z) solve(z * diag(ncol(E)) - E)

set.seed(123)
s <- 250

system.time(suppl_1 <- foreach(i = icount(100), .combine = rbind) %do% {
    c <- runif(1, 0, 1)
    sigma <- runif(1, 0, 1)
    rho <- runif(1, -1, 1)

    ## Make random matrix E
    E <- get_int_mat(s, sigma, 0.0, rho, 0.0) * get_adj_mat(s, c)
    Rmaj <- eigen(E, only.values = TRUE)$values |> abs() |> max()

    foreach(z = seq(0, 5, 0.1) + 1e-1, .combine = rbind) %do% {
        Ge <- rslv(E, Rmaj + z)
        tibble(z = z, tr = tr(Ge), total = sum(Ge))
    } |> mutate(run = i, size = s, connectance = c, sigma = sigma, rho = rho, Rmaj = Rmaj, .before = everything())
})

suppl_1 |> write_csvr("output/suppl_1.csv")

gp <- read_csvr("output/suppl_1.csv") |> mutate(diff = abs(total - tr) / total) |> ggplot(aes(x = z, y = diff)) + geom_point(shape = 16, size = 0.5, alpha = 0.25) +
    geom_smooth(method = "glm", method.args = list(family = Gamma), linewidth = 0.5, color = ggsci::pal_npg()(2)[1]) +
    xlab(expression(Delta*italic(z))) + ylab("Approximation error") + theme_st()
gp |> ggsaver("figS1", width = 8.7, height = 6, ext = "pdf")

## Figure S2: Demonstration of eigenvalu distributions of sensitivity essence Phi.
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

tbl_eig_a <- bind_rows(tbl_eig1, tbl_eig2, tbl_eig3) |> mutate(type = sprintf('(list(mu, rho))==(list(%.1f, %.1f))', sign(mu) * abs(mu), sign(rho) * abs(rho))) |> mutate(type = factor(type, levels = unique(type)))
tbl_eig_b <- bind_rows(tbl_eig4, tbl_eig5, tbl_eig6) |> mutate(type = sprintf('(list(mu, rho))==(list(%.1f, %.1f))', sign(mu) * abs(mu), sign(rho) * abs(rho))) |> mutate(type = factor(type, levels = unique(type)))

tbl_law_a <- bind_rows(tbl_law1, tbl_law2, tbl_law3) |> mutate(type = sprintf('(list(mu, rho))==(list(%.1f, %.1f))', sign(mu) * abs(mu), sign(rho) * abs(rho))) |> mutate(type = factor(type, levels = unique(type)))
tbl_law_b <- bind_rows(tbl_law4, tbl_law5, tbl_law6) |> mutate(type = sprintf('(list(mu, rho))==(list(%.1f, %.1f))', sign(mu) * abs(mu), sign(rho) * abs(rho))) |> mutate(type = factor(type, levels = unique(type)))

tbl_out_b <- bind_rows(tbl_out4, tbl_out5, tbl_out6) |> mutate(type = sprintf('(list(mu, rho))==(list(%.1f, %.1f))', sign(mu) * abs(mu), sign(rho) * abs(rho))) |> mutate(type = factor(type, levels = unique(type)))

gp1 <- tbl_eig_a |> ggplot(aes(x = real, y = imag)) + facet_wrap(. ~ type, nrow = 1, labeller = label_parsed) +
    geom_point(shape = 16, size = 0.25) +
    geom_path(data = tbl_law_a, linewidth = 0.25, color = ggsci::pal_npg()(3)[1]) +
    xlab(expression(paste("Real of ", lambda(Phi)))) + ylab(expression(paste("Imaginary of ", lambda(Phi)))) + theme_st(lunit = 2) 

gp2 <- tbl_eig_b |> ggplot(aes(x = real, y = imag)) + facet_wrap(. ~ type, nrow = 1, labeller = label_parsed) +
    geom_point(aes(size = eig_type), shape = 16) +
    geom_point(data = tbl_out_b, shape = 4, color = ggsci::pal_npg()(3)[3]) +
    geom_path(data = tbl_law_b, linewidth = 0.25, color = ggsci::pal_npg()(3)[1]) +
    scale_size_discrete(range = c(0.25, 1.5), guide = "none") +
    xlab(expression(paste("Real of ", lambda(Phi)))) + ylab(expression(paste("Imaginary of ", lambda(Phi)))) + theme_st(lunit = 2)

(gp1 / gp2 + plot_layout(axis_title = "collect")) |> ggsaver("figS2", width = 17.8, height = 10, ext = "pdf")

## Figure S3: Simulation analysis investigating the relationship between interaction mean and interaction correlation (mu & rho)
##      and the maximum absolute eigenvalue Rmax in finite-size networks (s = 250)
s <- 250
sigma <- 0.5

set.seed(123)
params <- expand_grid(c = c(0.1, 0.5, 0.9, 1.0), mu = seq(-1, 1, length.out = 21), rho = seq(-0.99, 0.99, length.out = 21)) |>
    mutate(run = 1:length(c), .before = mu)

system.time(suppl_2 <- foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
    ## Set the parameters
    vars <- calc_vars(s, tbl$c, sigma, tbl$mu, tbl$rho)
    cat(sprintf("%.2f%% start: (mu, rho) = (%.2f, %.2f)\n", 100 * tbl$run / nrow(params), tbl$mu, tbl$rho))

    foreach(i = icount(100), .combine = bind_rows) %do% {
        ## Make matrices A, E and Phi
        A <- get_int_mat(s, sigma, tbl$mu, tbl$rho, vars$d) * get_adj_mat(s, tbl$c)
        E <- A - tbl$mu * tbl$c - (vars$d - tbl$mu * tbl$c) * diag(s)
        Phi <- (E - tbl$mu * tbl$c / vars$psi * matrix(1, s, s) %*% E) / (vars$d - tbl$mu * tbl$c)
        Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()

        tibble(size = s, connectance = tbl$c, sigma = sigma, mu = tbl$mu, rho = tbl$rho, Rmax = criterion_fin(s, tbl$c, sigma, tbl$mu, tbl$rho), Robs = Robs)
    }
})

suppl_2 |> write_csvr("output/suppl_2.csv")

gp <- read_csvr("output/suppl_2.csv") |> group_by(size, connectance, sigma, mu, rho) |> summarize(across(c(Rmax, Robs), mean), .groups = "drop") |>
    mutate(connectance = sprintf('italic(C)==%.2f', connectance) |> (\(.x) factor(.x, levels = rev(unique(.x))))()) |>
    ggplot(aes(x = mu, y = rho, fill = log(Robs), z = log(Robs))) + facet_wrap(. ~ connectance, labeller = label_parsed) +
    geom_raster(interpolate = TRUE) + geom_contour(color = "grey60", linewidth = 0.1) + geom_contour(color = "grey40", linewidth = 0.25, breaks = log(1)) +
    scico::scale_fill_scico(palette = "vik", midpoint = 0, breaks = c(-1, 0, 1, 2)) + scale_x_continuous(expand = rep(1e-3, 2)) + scale_y_continuous(expand = rep(1e-3, 2)) +
    xlab(expression(paste("Mean interaction strength ", mu))) + ylab(expression(paste("Interaction correlation ", rho))) + labs(fill = expression(paste(log, group("|", lambda[1], "|")))) +
    theme_st() + theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1), panel.border = element_rect(color = "black", linewidth = 0.25), panel.spacing = unit(5, "mm"), strip.background = element_blank())

gp |> ggsaver("figS3", width = 8.7, height = 8.3, ext = "pdf")

## Figure S4: Simulation analysis investigating ecological indeterminacy of FR interaction systems
## Set the parameters
s <- 250
c <- 0.5
sigma <- 0.25
rho <- 0.0
mu <- 0.0
h <- 0.1
H <- h * (matrix(1, s, s) - diag(s))

## Type II functional response
zeta <- 1
set.seed(123)

dist1 <- foreach(i = 1:100, .combine = bind_rows) %do% {
    ## Make matrices A, D, E
    xe <- runif(s)
    A <- (get_int_mat(s, sigma, mu, rho, 0) * get_adj_mat(s, c)) |>
        (\(M) (zeta * (M %*% diag(xe^(zeta - 1)))) / (1 + H * (M %*% diag(xe^zeta)))^2)() |>
        (\(M) M - (max(Re(eigen(M, only.values = TRUE)$values + 1)) * diag(s)))()
    D <- matrix(mean(offdiag(A)), s, s) |> (\(M) M + diag(diag(A) - diag(M)))()
    E <- A - D

    ## calcualte eigenvalues of the sensitivity essence
    eig <- eigen(solve(D) %*% E, only.values = TRUE)$values
    tibble(zeta = zeta, run = i, Re = Re(eig), Im = Im(eig), abs = abs(eig))
}

## Type III functional response
zeta <- 2
set.seed(123)

dist2 <- foreach(i = 1:100, .combine = bind_rows) %do% {
    ## Make matrices A, D, E
    xe <- runif(s)
    A <- (get_int_mat(s, sigma, mu, rho, 0) * get_adj_mat(s, c)) |>
        (\(M) (zeta * (M %*% diag(xe^(zeta - 1)))) / (1 + H * (M %*% diag(xe^zeta)))^2)() |>
        (\(M) M - (max(Re(eigen(M, only.values = TRUE)$values + 1)) * diag(s)))()
    D <- matrix(mean(offdiag(A)), s, s) |> (\(M) M + diag(diag(A) - diag(M)))()
    E <- A - D

    ## calcualte eigenvalues of the sensitivity essence
    eig <- eigen(solve(D) %*% E, only.values = TRUE)$values
    tibble(zeta = zeta, run = i, Re = Re(eig), Im = Im(eig), abs = abs(eig))
}

## Make Figure S3A
gp1 <- bind_rows(dist1, dist2) |> mutate(zeta = sprintf('zeta==%1d', zeta)) |>
    ggplot(aes(x = abs)) + facet_wrap(. ~ zeta, nrow = 2, scales = "free_y", labeller = label_parsed) +
    geom_histogram(aes(y = after_stat(density)), bins = 25, alpha = 0.6) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, color = ggsci::pal_npg()(4)[1]) +
    scale_x_continuous(breaks = c(0.0, 0.5, 1.0)) + xlab(expression(paste("Absolute eigenvalues ", group("|", lambda(bold(Phi[B])), "|")))) + ylab("Density") + labs(tag = "A") + theme_st(lunit = 2)

## The relationship of maximum absolute eigenvalue and approximation accuracy of net interaction strengths
set.seed(123)
s <- 250
c <- 0.5
h <- 0.1
H <- h * (matrix(1, s, s) - diag(s))

system.time(suppl_3 <- foreach(zeta = c(rep(1, 500), rep(2, 500)), .combine = bind_rows) %do% {
    ## Set the parameters
    sigma <- runif(1, 0.0, 0.5)
    mu <- runif(1, -0.5, 0.5)
    rho <- runif(1, -0.5, 0.5)

    ## Make matrices A, D, E, Di, and Phi
    xe <- runif(s)
    A <- (get_int_mat(s, sigma, mu, rho, 0) * get_adj_mat(s, c)) |>
        (\(M) (zeta * (M %*% diag(xe^(zeta - 1)))) / (1 + H * (M %*% diag(xe^zeta)))^2)() |>
        (\(M) M - (max(Re(eigen(M, only.values = TRUE)$values + 1)) * diag(s)))()
    D <- matrix(mean(offdiag(A)), s, s) |> (\(M) M + diag(diag(A) - diag(M)))()
    E <- A - D
    Di <- solve(D)
    Phi <- Di %*% E

    Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()
    S_orig <- offdiag(solve(A))
    S_appr <- offdiag(matSeries(-Phi, 5) %*% Di)
    tibble(size = s, connectance = c, zeta = zeta, sigma = sigma, mu = mu, rho = rho, Robs = Robs, Rho = cor(S_orig, S_appr))
})

suppl_3 |> write_csvr("output/suppl_3.csv")
gp2 <- read_csvr("output/suppl_3.csv") |> ggplot(aes(x = Robs, y = abs(Rho), color = factor(zeta))) +
    geom_point(shape = 16, size = 0.5, alpha = 0.5) + geom_vline(xintercept = 1, linewidth = 0.25, color = ggsci::pal_npg()(4)[1]) +
    scale_color_brewer(palette = "Dark2") + scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) + xlab(expression(max((group("|", lambda(Phi[B]), "|"))))) +
    ylab(expression(paste("Prediction skill ", rho[italic(k)]))) + labs(tag = "B", color = expression(zeta)) + theme_st()

## The relationship between interaction mean/interaction correlation (mu & rho) and the largest eigenvalue in finite size networks with functional responses (s = 250)
s <- 250
sigma <- 0.5

h <- 0.1
H <- h * (matrix(1, s, s) - diag(s))

set.seed(123)
params <- expand_grid(zeta = c(1, 2), c = c(0.5, 1.0), mu = seq(-1, 1, length.out = 21), rho = seq(-0.99, 0.99, length.out = 21)) |>
    mutate(run = 1:length(c), .before = mu)

system.time(suppl_4 <- foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
    cat(sprintf("%.2f%% start: (mu, rho, zeta) = (%.2f, %.2f, %d)\n", 100 * tbl$run / nrow(params), tbl$mu, tbl$rho, tbl$zeta))

    x <- foreach(i = icount(100), .combine = bind_rows) %do% {
        ## Make matrices A, D, E
        xe <- runif(s)
        A <- (get_int_mat(s, sigma, tbl$mu, tbl$rho, 0) * get_adj_mat(s, tbl$c)) |>
            (\(M) (tbl$zeta * (M %*% diag(xe^(tbl$zeta - 1)))) / (1 + H * (M %*% diag(xe^tbl$zeta)))^2)() |>
            (\(M) M - (max(Re(eigen(M, only.values = TRUE)$values + 1))) * diag(s))()
        D <- matrix(mean(offdiag(A)), s, s) |> (\(M) M + diag(diag(A) - diag(M)))()
        E <- A - D
        Robs <- eigen(solve(D) %*% E, only.values = TRUE)$values |> abs() |> max()
        tibble(zeta = tbl$zeta, size = s, connectance = tbl$c, sigma = sigma, mu = tbl$mu, rho = tbl$rho, Robs = Robs)
    }
})

suppl_4 |> write_csvr("output/suppl_4.csv")

gp3 <- read_csvr("output/suppl_4.csv") |> group_by(zeta, size, connectance, sigma, mu, rho) |> summarize(Robs = mean(Robs), .groups = "drop") |>
    mutate(zeta = sprintf('zeta==%1d', zeta) |> (\(.x) factor(.x, levels = unique(.x)))(), connectance = sprintf('italic(C)==%.2f', connectance) |> (\(.x) factor(.x, levels = rev(unique(.x))))()) |>
    ggplot(aes(x = mu, y = rho, fill = log(Robs), z = log(Robs))) + facet_grid(zeta ~ connectance, labeller = label_parsed) +
    geom_raster(interpolate = TRUE) + geom_contour(color = "grey40", linewidth = 0.25, breaks = log(1)) +
    scico::scale_fill_scico(palette = "vik", midpoint = 0) + scale_x_continuous(expand = rep(1e-3, 2)) + scale_y_continuous(expand = rep(1e-3, 2)) +
    xlab(expression(paste("Mean interaction strength ", mu))) + ylab(expression(paste("Interaction correlation ", rho))) +
    labs(tag = "C", fill = expression(paste(log, group("|", lambda[1], "|")))) + theme_st(lunit = 2) + theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1), panel.border = element_rect(color = "black", linewidth = 0.25), panel.spacing = unit(5, "mm"))

(((gp1 / gp2 + plot_layout(height = c(2, 1))) | gp3) + plot_layout(width = c(2, 3))) |> ggsaver(name = "figS4", width = 17.8, height = 9, ext = "pdf")

## Figure S5: Simulation analysis investigating ecological indeterminacy of food webs
## A function to make a food web matrix
get_fwd_mat <- function(s, sigma, mu_u, mu_l, rho, d) {
    A <- matrix(0, s, s) + d * diag(s)
    A[upper.tri(A)] <- rnorm(s * (s - 1) / 2, mu_u, sigma)
    A <- A + apply(A, 1, \(x) if_else(x == 0, 0, mu_l + rho * (x - mu_u) + sqrt(1 - rho^2) * rnorm(s, 0, sigma)))
    return(A)
}

## The absolut eigenvalue distribution of sensitivity essence
## Set the parameters
s <- 250
c <- 0.5
sigma <- 0.1
rho <- 0.0

## Top-down regulated food webs
mu_u <- -0.3
mu_l <- 0.1

## Make theoretical PDF of absolute eigenvalues
vars <- calc_vars(s, c, sigma, (mu_u + mu_l) / 2, rho)
Rmaj <- (1 + abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - (mu_u + mu_l) / 2 * c)
dens1 <- tibble(case = 1, abs = seq(0, Rmaj, length.out = 250)) |> mutate(density = pdf_abs_eigs(abs, s, c, sigma, (mu_u + mu_l) / 2, rho))

## Make empirical PDF of absolute eigenvalues (food webs)
set.seed(123)

dist1 <- foreach(i = icount(100), .combine = bind_rows) %do% {
    ## Make matrices A, D, and E
    A <- get_fwd_mat(s, sigma, mu_u, mu_l, rho, 0) * get_adj_mat(s, c)
    diag(A) <- -(max(Re(eigen(A, only.values = TRUE)$values)) + 1)
    D <- get_fwd_mat(s, 0, c * mu_u, c * mu_l, 0, 0) + diag(diag(A))
    E <- A - D

    ## Calculate eigenavlues of the sensitivity essence
    eig <- eigen(solve(D) %*% E, only.values = TRUE)$values
    tibble(case = 1, run = i, Re = Re(eig), Im = Im(eig), abs = abs(eig))
}

## Bottom-up regulated food webs
mu_u <- -0.1
mu_l <- 0.3

## Make theoretical PDF of absolute eigenvalues
vars <- calc_vars(s, c, sigma, (mu_u + mu_l) / 2, rho)
Rmaj <- (1 + abs(vars$Pa)) * sqrt(s * vars$Va) / abs(vars$d - (mu_u + mu_l) / 2 * c)
dens2 <- tibble(case = 2, abs = seq(0, Rmaj, length.out = 250)) |> mutate(density = pdf_abs_eigs(abs, s, c, sigma, (mu_u + mu_l) / 2, rho))

## Make empirical PDF of absolute eigenvalues (food webs)
set.seed(123)

dist2 <- foreach(i = icount(100), .combine = bind_rows) %do% {
    ## Make matrices A, D, and E
    A <- get_fwd_mat(s, sigma, mu_u, mu_l, rho, 0) * get_adj_mat(s, c)
    diag(A) <- -(max(Re(eigen(A, only.values = TRUE)$values)) + 1)
    D <- get_fwd_mat(s, 0, c * mu_u, c * mu_l, 0, 0) + diag(diag(A))
    E <- A - D

    ## Calculate eigenavlues of the sensitivity essence
    eig <- eigen(solve(D) %*% E, only.values = TRUE)$values
    tibble(case = 2, run = i, Re = Re(eig), Im = Im(eig), abs = abs(eig))
}

## Make figure S5A
dist <- bind_rows(dist1, dist2) |> mutate(case = if_else(case == 1, "Top-down", "Bottom-up"))
dens <- bind_rows(dens1, dens2) |> mutate(case = if_else(case == 1, "Top-down", "Bottom-up"))

gp1 <- dist |> ggplot(aes(x = abs)) + facet_wrap(. ~ case, nrow = 1, scales = "free_y") +
    geom_histogram(aes(y = after_stat(density)), bins = 25, alpha = 0.6) +
    geom_line(data = dens, aes(y = density), linewidth = 0.25, color = ggsci::pal_npg()(4)[4]) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, color = ggsci::pal_npg()(4)[1]) +
    xlab(expression(paste("Absolute eigenvalues ", group("|", lambda(Phi[F]), "|")))) + ylab("Density") + labs(tag = "A") + theme_st()

## The relationship of maximum absolute eigenvalue and approximation accuracy of net interaction networks
set.seed(123)
s <- 250
c <- 0.5
rho <- 0.0

system.time(suppl_5 <- foreach(i = icount(1000), .combine = bind_rows) %do% {
    ## Set the parameters
    sigma <- runif(1, 0.0, 0.5)
    mu_u <- runif(1, -0.5, 0.0)
    mu_l <- runif(1,  0.0, 0.5)

    ## Make matrices A, D, E, Di and Phi
    A <- get_fwd_mat(s, sigma, mu_u, mu_l, rho, 0) * get_adj_mat(s, c)
    diag(A) <- -(max(Re(eigen(A, only.values = TRUE)$values)) + 1)
    D <- get_fwd_mat(s, 0, c * mu_u, c * mu_l, 0, 0) + diag(diag(A))
    E <- A - D
    Di <- solve(D)
    Phi <- Di %*% E

    Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()
    S_orig <- offdiag(solve(A))
    S_appr <- offdiag(matSeries(-Phi, 5) %*% Di)
    tibble(size = s, connectance = c, sigma = sigma, mu_u = mu_u, mu_l = mu_l, rho = rho, Robs = Robs, Rho = cor(S_orig, S_appr))
})

suppl_5 |> write_csvr("output/suppl_5.csv")
gp2 <- read_csvr("output/suppl_5.csv") |> ggplot(aes(x = Robs, y = Rho)) +
    geom_point(shape = 16, size = 0.5, alpha = 0.5) + geom_vline(xintercept = 1, linetype = 2, linewidth = 0.25, color = ggsci::pal_npg()(4)[1]) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), limits = c(0, 1)) + xlab(expression(max((group("|", lambda(bold(Phi[F])), "|"))))) +
    ylab(expression(paste("Prediction skill ", rho[italic(k)]))) + labs(tag = "B", color = expression(gamma)) + theme_st()

## The effect of mean interaction strengths of upper- and lower-triangular part
s <- 250
c <- 0.5
sigma <- 0.25
rho <- 0.0

set.seed(123)
params <- expand_grid(mu_u = seq(-1.05, -0.05, length.out = 21), mu_l = seq(0.05, 1.05, length.out = 21)) |>
    mutate(run = 1:length(mu_u), .before = mu_u)

system.time(suppl_6 <- foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
    cat(sprintf("%.2f%% start: (mu_u, mu_l) = (%.2f, %.2f)\n", 100 * tbl$run / nrow(params), tbl$mu_u, tbl$mu_l))

    foreach(i = icount(100), .combine = bind_rows) %do% {
        ## Make matrices A, D, E, and Phi
        A <- get_fwd_mat(s, sigma, tbl$mu_u, tbl$mu_l, rho, 0) * get_adj_mat(s, c)
        diag(A) <- -(max(Re(eigen(A, only.values = TRUE)$values)) + 1)
        D <- get_fwd_mat(s, 0, c * tbl$mu_u, c * tbl$mu_l, 0, 0) + diag(diag(A))
        E <- A - D
        Phi <- solve(D) %*% E

        Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()
        tibble(size = s, connectance = c, sigma = sigma, mu_u = tbl$mu_u, mu_l = tbl$mu_l, Robs = Robs)
    }
})

suppl_6 |> write_csvr("output/suppl_6.csv")

gp3 <- read_csvr("output/suppl_6.csv") |> group_by(size, connectance, sigma, mu_u, mu_l) |> summarize(Robs = mean(Robs), .groups = "drop") |>
    ggplot(aes(x = mu_u, y = mu_l, fill = log(Robs), z = log(Robs))) + geom_raster(interpolate = TRUE) + geom_contour(color = "grey40", linewidth = 0.25, breaks = log(1)) +
    geom_abline(slope = -1, intercept = 0, linetype = 2, linewidth = 0.25) +
    scico::scale_fill_scico(palette = "vik", midpoint = 0) + scale_x_continuous(expand = rep(1e-3, 2)) + scale_y_continuous(expand = rep(1e-3, 2)) +
    xlab(expression(mu[U])) + ylab(expression(mu[L])) + labs(tag = "C", fill = expression(paste(log, group("|", lambda[1], "|")))) + theme_st(lunit = 2) +
    theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1), panel.border = element_rect(color = "black", linewidth = 0.25), panel.spacing = unit(5, "mm"))

(gp1 / (gp2 + gp3 + plot_layout(width = c(6, 7)))) |> ggsaver(name = "figS5", width = 8.7, height = 8.5, ext = "pdf")

