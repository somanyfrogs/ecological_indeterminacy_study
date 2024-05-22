#' @title Figure_1.R
#' @description Making figure 1 for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20230104 by K.Kawatsu.
#'      Last update: 20240517.

## Load R functions
source("R/functions.R")

## Figure 1A: Proecess of ecological indeterminacy analysis 1
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

label <- c("'Interaction matrix'~bold(A)", "'Deterministic'~bold(D)", "'Random'~bold(E)")
tbl_A <- A |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "A")
tbl_D <- D |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "D")
tbl_E <- E |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "E")
tex <- tibble(type = c("A", "D", "E"), text = c("=", "+", NA)) |> mutate(type = factor(type, levels = unique(type), labels = label))

gp1 <- bind_rows(tbl_A, tbl_D, tbl_E) |> mutate(type = factor(type, levels = unique(type), labels = label)) |> ggplot() + facet_wrap(. ~ type, labeller = label_parsed) +
    geom_tile(aes(x = col, y = rev(row), fill = coeff)) + geom_tile(aes(x = col, y = rev(row), alpha = if_else(col == row, "y", "n")), fill = '#320000') +
    geom_tile(aes(x = col, y = rev(row), alpha = if_else(coeff == 0, "y", "n")), fill = grey(1)) +
    geom_text(data = tex, aes(label = text), x = 27.5, y = 13, size = 3) + scico::scale_fill_scico(palette = "vik", limits = c(-0.65, 0.65), breaks = c(-0.5, 0.0, 0.5), direction = -1) +
    scale_alpha_discrete(range = c(0, 1), guide = "none") + scale_x_continuous(expand = rep(0, 2)) + scale_y_continuous(expand = rep(0, 2)) + coord_cartesian(clip = "off") +
    labs(tag = "A") + theme_wt1() + theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1, margin = margin(0.0, 0, 0, 1)),
                                          panel.border = element_blank(), panel.grid = element_blank(), panel.spacing.x = unit(5, "mm"), strip.background = element_blank())

## Figure 1B: Process of ecological indeterminacy analysis 2
## Building sensitivity essence and Neumann series expansion

## Making matrices Di, Phi, True and approximated inverse interaction matrix
k <- 4
Di <- solve(D)
Phi <- Di %*% E
St <- solve(A)
Sa <- matSeries(-Phi, k) %*% Di

tbl_Phi <- Phi |> as_tibbler(1:s) |> mutate(row = 1:s) |> pivot_longer(!row, names_to = "col", values_to = "coeff") |> mutate(col = as.integer(col), type = "'Sensitivity\nessence'~Phi==bold(D)^-1*bold(E)")

gp2 <- tbl_Phi |> ggplot(aes(x = col, y = rev(row), fill = coeff)) + facet_wrap(. ~ type, labeller = label_parsed) + geom_tile() +
    scico::scale_fill_scico(palette = "vik", limits = c(-0.65, 0.65), direction = -1, guide = "none") +
    scale_x_continuous(expand = rep(1e-3, 2)) + scale_y_continuous(expand = rep(1e-3, 2)) + coord_cartesian(clip = "off") +
    labs(tag = "B") + theme_wt1() + theme(legend.margin = margin(0, 0, 0, 0), panel.border = element_blank(), panel.grid = element_blank(), panel.spacing.x = unit(5, "mm"), strip.background = element_blank())

gp3 <- tibble(x = 0, y = 0, text = 'Neumann series\nexpansion') |> ggplot(aes(x = x, y = y)) + facet_wrap(. ~ text) +
    annotate(geom = "segment", x = -1, xend = 1, y = -0.5, yend = -0.5, arrow = arrow(length = unit(0.1, "inches"))) +
    annotate(geom = "label", x = 0, y = 0.4, label = expression(sum((-bold(Phi))^italic(i), italic(i)==0, italic(k))*bold(D)^-1%~~%bold(A)^-1), size = 2.3) +
    annotate(geom = "text", x = 0, y = -0.1, label = expression(gamma[max]==0.756), size = 2.3) +
    xlim(-1, 1) + ylim(-1, 1) + theme_wt1() + theme(panel.background = element_blank(), panel.border = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "lines"), strip.text = element_text(size = 6), strip.background = element_blank())

gp4 <- tibble(true = offdiag(St), appr = offdiag(Sa), text = sprintf('italic(k)=="%d"', k)) |> ggplot(aes(x = true, y = appr)) + facet_wrap(. ~ text, labeller = label_parsed) +
    geom_point(shape = 16, size = 0.25) + geom_abline(intercept = 0, slope = 1, linewidth = 0.25, color = ggsci::pal_npg()(2)[1]) +
    xlab(expression(paste("True ", bold(A[italic(ij)]^-1)))) + ylab(expression(paste("Approximated ", bold(A[italic(ij)]^-1)))) +
    theme_st() + theme(axis.text = element_blank(), axis.ticks = element_blank(), strip.background = element_blank())

(gp1 / (gp2 + gp3 + gp4)) |> ggsaver("fig01a", width = 11.4, height = 8, ext = "pdf")
