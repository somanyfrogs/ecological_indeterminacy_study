#' @title Figure_1.R
#' @description Making figure 1 for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20230104 by K.Kawatsu.
#'      Last update: 20240325.

## Load R functions
source("R/functions.R")

## Figure 1A: Demonstration of direct and indirect effects n a small-network graph
## Example: An IGP (Intraguild predation) module with two cases
## Press perturbation: addition of the top predator (x3) to the equilibrium
s <- 3
c <- 1.0
nodes <- tibble(id = str_c("italic(x)[", 1:3, "]"))
theta <- pi * seq(1 / 2 + 0.07, 3 / 2 - 0.12, length.out = 250)
arcs <- tibble(x = cos(theta) + 0.5, y = sin(theta) + 2)

## Case 1: Like food-chain model (The top predator coonsume the intermediate predator than the prey)
A <- rbind(c(0.00, -0.90, -0.20),
           c(0.60,  0.00, -0.90),
           c(0.10,  0.60,  0.00))
sigma <- sd(offdiag(A))
mu <- mean(offdiag(A))
rho <- cor(offdiag(A, 2), offdiag(t(A), 2))
d <- calc_vars(s, c, sigma, mu, rho)$d

edges <- tibble(from = c(2, 3, 3), to = c(1, 1, 2), weight = abs(offdiag(A, 2)), case = "Case 1")
arcs <- arcs |> mutate(case = "Case 1")

gp1 <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) |> ggraph() + facet_edges(~ case, strip.position = "left") +
    geom_edge_link(aes(width = weight), arrow = arrow(length = unit(1.5, "mm"), type = "closed"), start_cap = circle(2, "mm"), end_cap = circle(2, "mm")) +
    geom_path(data = arcs, aes(x= x, y = y), linetype = 2, linewidth = 0.7) + geom_node_point(size = 3, aes(color = id)) +
    geom_node_text(aes(label = id), size = 2, parse = TRUE) + annotate(geom = "point", x = cos(theta[250]) + 0.5, y = sin(theta[250]) + 2, shape = 16, size = 2) +
    ggsci::scale_color_npg(guide = "none") + scale_edge_width(range = c(0.20, 0.70), guide = "none") + scale_x_continuous(expand = rep(0.1, 2)) +
    scale_y_continuous(expand = rep(0.1, 2)) + labs(tag = "A") + theme_wt1() + theme(strip.text = element_text(margin = margin(1, 1, 1, 1)))

## Dynamics after press perturbation
A <- A + d * diag(s); r <- -A |> rowMeans()
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
d <- calc_vars(s, c, sigma, mu, rho)$d

edges <- tibble(from = c(2, 3, 3), to = c(1, 1, 2), weight = abs(offdiag(A, 2)), case = "Case 2")
arcs <- arcs |> mutate(case = "Case 2")

gp3 <- tbl_graph(nodes = nodes, edges = edges, directed = TRUE) |> ggraph() + facet_edges(~ case, strip.position = "left") +
    geom_edge_link(aes(width = weight), arrow = arrow(length = unit(1.5, "mm"), type = "closed"), start_cap = circle(2, "mm"), end_cap = circle(2, "mm")) +
    geom_path(data = arcs, aes(x= x, y = y), linetype = 2, linewidth = 0.2) + geom_node_point(size = 3, aes(color = id)) +
    geom_node_text(aes(label = id), size = 2, parse = TRUE) + annotate(geom = "point", x = cos(theta[250]) + 0.5, y = sin(theta[250]) + 2, shape = 16, size = 2) +
    ggsci::scale_color_npg(guide = "none") + scale_edge_width(range = c(0.20, 0.80), guide = "none") + scale_x_continuous(expand = rep(0.1, 2)) +
    scale_y_continuous(expand = rep(0.1, 2)) + theme_wt1() + theme(strip.text = element_text(margin = margin(1, 1, 1, 1)))

## Dynamics after press perturbation
A <- A + d * diag(s); r <- -A |> rowMeans()
ts <- odeRK4(x_ini = rep(1 / 3, 3), fun = glv, params = list(A = A, r = r, b = c(0.0, 0.0, 1e-3)), t_end = 30, h = 0.001) |> mutate(across(starts_with("x"), \(.x) .x - 1 /3))
gp4 <- ts |> pivot_longer(!time, names_to = "species", values_to = "density") |> ggplot(aes(x = time, y = density, color = species)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.25) + geom_line(linewidth = 0.5) + scale_y_continuous(breaks = 0, limits = c(-4.5e-4, 7e-4)) +
    ggsci::scale_color_npg(guide = "none") + labs(color = NULL) + xlab("Time") + ylab("Density change") + theme_st() + theme(axis.text.x = element_blank(), axis.ticks = element_blank())

## Figure 1B: Process of ecological indeterminacy analysis 1
## Decomposition of an interaction matrix A into a deterministic matrix D and a random matrix E.

## Set the parameters
set.seed(333)
s <- 25
c <- 0.5
sigma <- 0.2
mu <- -0.1
rho <- -0.75

## Make matrices A, D and E
d <- calc_vars(s, c, sigma, mu, rho)$d
A <- get_int_mat(s, sigma, mu, rho, d) * get_adj_mat(s, c)
D <- (d - mu * c) * diag(s) + mu * c
E <- A - D

tbl_A <- A |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "A")
tbl_D <- D |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "D")
tbl_E <- E |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "E")
tex <- tibble(type = c("A", "D", "E"), text = c("=", "+", NA))

gp5 <- bind_rows(tbl_A, tbl_D, tbl_E) |> ggplot() + facet_wrap(. ~ type) +
    geom_tile(aes(x = col, y = rev(row), fill = coeff)) + geom_tile(aes(x = col, y = rev(row), alpha = if_else(coeff == 0, "y", "n")), fill = grey(1)) +
    geom_text(data = tex, aes(label = text), x = 27.5, y = 13, size = 3) + scico::scale_fill_scico(palette = "vik", limits = c(-0.65, 0.65), breaks = c(-0.5, 0.0, 0.5), direction = -1) +
    scale_alpha_discrete(range = c(0, 1), guide = "none") + scale_x_continuous(expand = rep(0, 2)) + scale_y_continuous(expand = rep(0, 2)) + coord_cartesian(clip = "off") +
    labs(tag = "B") + theme_wt1() + theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1, margin = margin(0.0, 0, 0, 1)),
                                          panel.border = element_blank(), panel.grid = element_blank(), panel.spacing.x = unit(5, "mm"),
                                          strip.background = element_blank(), strip.text = element_text(size = 8, face = 2))

## Figure 1C: Process of ecological indeterminacy analysis 2
## Bulding sensitivity essence and Neumann series expansion

## Making matrices Di, Phi, True and approximated inverse interaction matrix
k <- 4
Di <- solve(D)
Phi <- Di %*% E
St <- solve(A)
Sa <- matSeries(-Phi, k) %*% Di

tbl_Phi <- Phi |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = 'bold(Phi)==bold(D)^-1*bold(E)')

gp6 <- tbl_Phi |> ggplot(aes(x = col, y = rev(row), fill = coeff)) + facet_wrap(. ~ type, labeller = label_parsed) + geom_tile() +
    scico::scale_fill_scico(palette = "vik", limits = c(-0.65, 0.65), direction = -1, guide = "none") +
    scale_x_continuous(expand = rep(1e-3, 2)) + scale_y_continuous(expand = rep(1e-3, 2)) + coord_cartesian(clip = "off") +
    labs(tag = "C") + theme_wt1() + theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1),
                                          panel.border = element_blank(), panel.grid = element_blank(), panel.spacing.x = unit(5, "mm"),
                                          strip.background = element_blank(), strip.text = element_text(size = 8, face = 2))

gp7 <- tibble(x = 0, y = 0, text = 'Neumann series\nexpansion') |> ggplot(aes(x = x, y = y)) + facet_wrap(. ~ text) +
    annotate(geom = "segment", x = -1, xend = 1, y = -0.5, yend = -0.5, arrow = arrow(length = unit(0.1, "inches"))) +
    annotate(geom = "label", x = 0, y = 0.4, label = expression(sum((-bold(Phi))^italic(i), italic(i)==0, italic(k))*bold(D)^-1%~~%bold(A)^-1), size = 2.3) +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + theme(plot.margin = unit(c(0, 0, 0, 0), "lines"), strip.text = element_text(size = 6))

gp8 <- tibble(true = offdiag(St), appr = offdiag(Sa), text = sprintf('italic(k)=="%d"', k)) |> ggplot(aes(x = true, y = appr)) + facet_wrap(. ~ text, labeller = label_parsed) +
    geom_point(shape = 16, size = 0.25) + geom_abline(intercept = 0, slope = 1, linewidth = 0.25, color = ggsci::pal_npg()(2)[1]) +
    xlab(expression(paste("True ", bold(A[italic(ij)]^-1)))) + ylab(expression(paste("Approximated ", bold(A[italic(ij)]^-1)))) +
    theme_st() + theme(axis.text = element_blank(), axis.ticks = element_blank(), strip.background = element_blank())

((((gp1 | gp2 | gp3 | gp4) + plot_layout(width = c(1, 2, 1, 2))) / gp5 / (gp6 + gp7 + gp8)) + plot_layout(height = c(2, 3, 3))) |> ggsaver("fig01a", width = 11.4, height = 12, ext = "pdf")
