#' @title Box_1.R
#' @description Making figures of Box 1 for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20240517 by K. Kawatsu based on Figure_1.R and Figure_5.R
#'      Last update: 20240522.

## Load R functions
source("R/functions.R")

## Panel A: Demonstration of direct and indirect effects in a small-network graph
## Example: An IGP (Intraguilde predation) module with two caess
## Press perturbation: addition of the top predator (x3) to the equilibrium
s <- 3
c <- 1.0
nodes <- tibble(id = str_c("italic(x)[", 1:3, "]"))
theta <- pi * seq(1 / 2 + 0.07, 3 / 2 - 0.12, length.out = 250)
arcs <- tibble(x = cos(theta) + 0.5, y = sin(theta) + 2)

## Case 1: Like food-chain model (The top predator consume the intermediate predator than the prey)
A <- rbind(c(0.00, -0.90, -0.20),
           c(0.60,  0.00, -0.90),
           c(0.10,  0.60,  0.00))
sigma <- sd(offdiag(A))
mu <- mean(offdiag(A))
rho <- cor(offdiag(A, 2), offdiag(t(A), 2))
d <- calc_vars(s, c, sigma, mu, rho)$d
A <- A + d * diag(s); r <- -A |> rowMeans()

## Eigenvalues of food-chain IGP (1.003, 1.003, 0.066)
D <- (d - mu * c) * diag(s) + mu * c
eigen(solve(D) %*% (A - D))$values |> abs()

edges <- tibble(from = c(2, 3, 3), to = c(1, 1, 2), weight = abs(offdiag(A, 2)), case = "food-chain IGP")
arcs <- arcs |> mutate(case = "food-chain IGP")

gp1 <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) |> ggraph() + facet_edges(~ case) +
    geom_edge_link(aes(width = weight), arrow = arrow(length = unit(1.5, "mm"), type = "closed"), start_cap = circle(2, "mm"), end_cap = circle(2, "mm")) +
    geom_path(data = arcs, aes(x= x, y = y), linetype = 2, linewidth = 0.7) + geom_node_point(size = 3, aes(color = id)) +
    geom_node_text(aes(label = id), size = 2, parse = TRUE) + annotate(geom = "point", x = cos(theta[250]) + 0.5, y = sin(theta[250]) + 2, shape = 16, size = 2) +
    ggsci::scale_color_npg(guide = "none") + scale_edge_width(range = c(0.20, 0.70), guide = "none") + scale_x_continuous(expand = rep(0.1, 2)) +
    scale_y_continuous(expand = rep(0.1, 2)) + labs(tag = "A") + theme_wt1() + theme(panel.background = element_blank(), panel.border = element_blank(), strip.background = element_blank(), strip.text = element_text(margin = margin(1, 1, 1, 1)))

## Dynamics after press perturbation
ts <- odeRK4(x_ini = rep(1 / 3, 3), fun = glv, params = list(A = A, r = r, b = c(0.0, 0.0, 1e-3)), t_end = 30, h = 0.001) |> mutate(across(starts_with("x"), \(.x) .x - 1 /3))
gp2 <- ts |> pivot_longer(!time, names_to = "species", values_to = "density") |> ggplot(aes(x = time, y = density, color = species)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) + geom_line(linewidth = 0.5) + scale_y_continuous(breaks = 0, limits = c(-4.5e-4, 7e-4)) +
    ggsci::scale_color_npg(guide = "none") + labs(color = NULL) + xlab("Time") + ylab("Density change") + theme_st() + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

## Case 2: Like exploitative competition (The top predator consume the prey than the intermediate predator)
A <- rbind(c(0.00, -0.90, -0.90),
           c(0.60,  0.00, -0.20),
           c(0.60,  0.10,  0.00))
sigma <- sd(offdiag(A))
mu <- mean(offdiag(A))
rho <- cor(offdiag(A, 2), offdiag(t(A), 2))
d2 <- calc_vars(s, c, sigma, mu, rho)$d
A <- A + d * diag(s); r <- -A |> rowMeans()

## Eigenvalues of competitive IGP (0.922, 0.922, 0.065)
D <- (d - mu * c) * diag(s) + mu * c
eigen(solve(D) %*% (A - D))$values |> abs()

edges <- tibble(from = c(2, 3, 3), to = c(1, 1, 2), weight = abs(offdiag(A, 2)), case = "competitive IGP")
arcs <- arcs |> mutate(case = "competitive IGP")

gp3 <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) |> ggraph() + facet_edges(~ case) +
    geom_edge_link(aes(width = weight), arrow = arrow(length = unit(1.5, "mm"), type = "closed"), start_cap = circle(2, "mm"), end_cap = circle(2, "mm")) +
    geom_path(data = arcs, aes(x= x, y = y), linetype = 2, linewidth = 0.2) + geom_node_point(size = 3, aes(color = id)) +
    geom_node_text(aes(label = id), size = 2, parse = TRUE) + annotate(geom = "point", x = cos(theta[250]) + 0.5, y = sin(theta[250]) + 2, shape = 16, size = 2) +
    ggsci::scale_color_npg(guide = "none") + scale_edge_width(range = c(0.20, 0.80), guide = "none") + scale_x_continuous(expand = rep(0.1, 2)) +
    scale_y_continuous(expand = rep(0.1, 2)) + theme_wt1() + theme(panel.border = element_blank(), panel.background = element_blank(), strip.background = element_blank(), strip.text = element_text(margin = margin(1, 1, 1, 1)))

## Dynamics after press perturbation
ts <- odeRK4(x_ini = rep(1 / 3, 3), fun = glv, params = list(A = A, r = r, b = c(0.0, 0.0, 1e-3)), t_end = 30, h = 0.001) |> mutate(across(starts_with("x"), \(.x) .x - 1 /3))
gp4 <- ts |> pivot_longer(!time, names_to = "species", values_to = "density") |> ggplot(aes(x = time, y = density, color = species)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) + geom_line(linewidth = 0.5) + scale_y_continuous(breaks = 0, limits = c(-4.5e-4, 7e-4)) +
    ggsci::scale_color_npg(guide = "none") + labs(color = NULL) + xlab("Time") + ylab("Density change") + theme_st() + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

## Penel B: Press perturbation dynamics close to different eigenvectors
## Example: A food-web-like IGP.
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
    odeRK4(x_ini = rep(1 / s, s), fun = glv, params = list(A = A, r = r, b  = bs[[i]]), t_end = 1000, h = 0.01) |>
        mutate(across(everything(), abs), run = i, .before = everything())
} |> write_csvr("output/boxB.csv")

## Make panel B
## Convert 3D data to 2D data with viewing transformation
theta <- 135; phi <- 45
pmat <- persp(z = diag(2), theta = theta, phi = phi)
orig <- tibble(x1 = 1/s, x2 = 1/s, x3 = 1/s) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat))

## 2D transformation of perturbation dynamics
tbl1 <- read_csvr("output/boxB.csv") |> filter(run == 1) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
tbl2 <- read_csvr("output/boxB.csv") |> filter(run == 2) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
tbl <- bind_rows(tbl1, tbl2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

## 2D transformation of eigenvectors
vec1 <- tibble(x1 = 1/s + c(0, 1e-3 * evec[1, 1]), x2 = 1/s + c(0, 1e-3 * evec[2, 1]), x3 = 1/s + c(0, 1e-3 * evec[3, 1])) |>
    mutate(across(everything(), Re)) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
vec2 <- tibble(x1 = 1/s + c(0, 1e-3 * evec[1, 3]), x2 = 1/s + c(0, 1e-3 * evec[2, 3]), x3 = 1/s + c(0, 1e-3 * evec[3, 3])) |>
    mutate(across(everything(), Re)) |> with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
vec_tbl <- bind_rows(vec1, vec2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

## 2D transformation of perturbation dynamics
bs1 <- tibble(x1 = 1/s + c(0, bs[[1]][1]) / 3, x2 = 1/s + c(0, bs[[1]][2]) / 3, x3 = 1/s + c(0, bs[[1]][3]) / 3) |> mutate(across(everything(), Re)) |>
    with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 1, x = .lst$x, y = .lst$y))()
bs2 <- tibble(x1 = 1/s + c(0, bs[[2]][1]) / 3, x2 = 1/s + c(0, bs[[2]][2]) / 3, x3 = 1/s + c(0, bs[[2]][3]) / 3) |> mutate(across(everything(), Re)) |>
    with(trans3d(x = x1, y = x2, z = x3, pmat = pmat)) |> (\(.lst) tibble(run = 2, x = .lst$x, y = .lst$y))()
bs_tbl <- bind_rows(bs1, bs2) |> mutate(run = sprintf('group("|", lambda[%d], "|")==%.3f', c(1, 3)[run], abs(eval[run + 1])), x = x - orig$x, y = y - orig$y)

# 2D transformation of after perturbation equilibirum calculated with the true inverse interaction matrix A^-1
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

gp5 <- tbl |> ggplot(aes(x = x, y = y)) + facet_wrap(. ~ run, nrow = 2, labeller = label_parsed) +
    annotate("segment", x = rep(0, 3), y = rep(0, 3) + 8e-5, xend = c(-4e-4, 0, 4e-4), yend = c(-3e-4, 4.5e-4, -3e-4), linewidth = 0.1, color = "grey50") +
    annotate("text", x = -3.6e-4, y = -2.0e-4, label = expression(italic(x)[2]), size = 2) + annotate("text", x = 0.5e-4, y = 4.3e-4, label = expression(italic(x)[3]), size = 2) + annotate("text", x = 3.6e-4, y = -2.0e-4, label = expression(italic(x)[1]), size = 2) +
    geom_path(data = vec_tbl, color = "orange") + geom_path(data = bs_tbl, linewidth = 0.5, arrow = arrow(length = unit(2, "mm"))) + geom_path(linewidth = 0.25) +
    annotate("point", x = 0, y = 0, shape = 21, color = "black", fill = "black") +
    geom_point(data = out_t, shape = 25, fill = "black") + geom_point(data = out_a, aes(fill = factor(order)), shape = 25) +
    ggsci::scale_fill_npg() + xlim(-4e-4, 4e-4) + ylim(-4e-4, 4.5e-4) + labs(tag = "B", fill = expression(italic(k))) + theme_st(lunit = 2) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), legend.margin = margin(0, 0, 0, 0), panel.grid = element_blank())

((gp1 / gp3 | gp2 / gp4 | gp5) + plot_layout(widths = c(1, 2, 2))) |> ggsaver("box", width = 11.4, height = 8, ext = "pdf")
