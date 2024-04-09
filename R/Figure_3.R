#' @title Figure_3.R
#' @description Making figure 3 for the Paper entitled:
#'      "Unraveling emergenet network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221223 by K.Kawatsu.
#'      Last update: 20240408.

## Load R functions
source("R/functions.R")

## Simulation for Figure 3: investigating the relationship between ecological complexity (S & C) and the maximum absolute eigenvalue
## Set the parameters
set.seed(123)
params <- tibble(sim = 1:4, mu = c(-0.1, 0.0, 0.2, 0.1), sigma = 0.25, Pa = c(0.0, -0.25, 0.0, 0.25)) |>
    left_join(expand_grid(sim = 1:4, s = c(seq(10, 310, 10), seq(900, 1000, 50)), c = seq(0.4, 1.0, 0.3))) |> mutate(run = 1:length(sim), .before = mu)

tmp <- function(tbl, n_no, f_with) tbl |> mutate({{n_to}} := eval(parse(text = f_with)))

system.time(sim_3 <- foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
    ## Set the parameters
    s <- tbl$s
    c <- tbl$c
    sigma <- tbl$sigma
    mu <- tbl$mu
    rho <- tbl$Pa - (1 - tbl$Pa) * (1 - c) * (mu / sigma)^2

    vars <- calc_vars(s, c, sigma, mu, rho)
    cat(sprintf("%.2f%% start: (s, c, sigma, mu, Pa) = (%d, %.2f, %.2f, %.2f, %.2f)\n", 100 * tbl$run / nrow(params), s, c, sigma, mu, tbl$Pa))

    foreach(i = 1:100, .combine = bind_rows) %do% {
        ## Make matrices A, E and Phi
        A <- get_int_mat(s, sigma, mu, rho, vars$d) * get_adj_mat(s, c)
        E <- A - mu * c - (vars$d - mu * c) * diag(s)
        Phi <- (E - mu * c / vars$psi * matrix(1, s, s) %*% E) / (vars$d - mu * c)
        Robs <- eigen(Phi, only.values = TRUE)$values |> abs() |> max()

        tibble(sim = tbl$sim, size = s, connectance = c, sigma = sigma, mu = mu, rho = rho, Va = vars$Va, Pa = vars$Pa, Rmax = criterion_fin(s, c, sigma, mu, rho), Robs = Robs)
    }
})

sim_3 |> write_csvr("output/sim_3.csv")

## Make figure
tbl <- read_csvr("output/sim_3.csv") |> mutate(connectance = factor(connectance, levels = c("1", "0.7", "0.4"), labels = c("1.0", "0.7", "0.4")))
lgnd <- tbl |> ggplot(aes(x = size, y = Robs, color = connectance)) + facet_wrap(. ~ sim) +
    geom_point(shape = 16, size = 0.5, alpha = 0.25) + geom_line(aes(y = Rmax)) +
    ggsci::scale_color_npg() + labs(color = expression(italic(C))) + theme_st() +
    theme(legend.margin = margin(0, 0, 0, 0), legend.text = element_text(hjust = 1, margin = margin(1, 1, 1, 1)))
lgnd <- cowplot::get_legend(lgnd)

gp1 <- tbl |> filter(sim == 1) |> ggplot(aes(x = size, y = Robs, color = connectance)) + geom_hline(yintercept = 1, color = "grey50", linewidth = 0.25) +
    geom_point(shape = 16, size = 0.25, alpha = 0.1) + geom_line(aes(y = Rmax)) + ggsci::scale_color_npg(guide = "none") +
    ggbreak::scale_x_break(c(310, 900), space = 0.05) + scale_x_continuous(breaks = c(0, 100, 200, 900, 1000)) + scale_y_continuous(limits = c(0, 1.5), breaks = c(0.0, 0.5, 1.0, 1.5)) +
    labs(title = expression(paste((list(mu, italic(P)[bold(A)])), " = (-0.10, 0.00)"))) + theme_st() +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 6, hjust = 0.7, vjust = -2, margin = margin(0, 0, 0, 0)))

gp2 <- tbl |> filter(sim == 2) |> ggplot(aes(x = size, y = Robs, color = connectance)) + geom_hline(yintercept = 1, color = "grey50", linewidth = 0.25) +
    geom_point(shape = 16, size = 0.25, alpha = 0.1) + geom_line(aes(y = Rmax)) + ggsci::scale_color_npg(guide = "none") +
    ggbreak::scale_x_break(c(310, 900), space = 0.05) + scale_x_continuous(breaks = c(0, 100, 200, 900, 1000)) + scale_y_continuous(limits = c(0, 1.5), breaks = c(0.0, 0.5, 1.0, 1.5)) +
    labs(title = expression(paste((list(mu, italic(P)[bold(A)])), " = (0.00, -0.25)"))) + theme_st() +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 6, hjust = 0.7, vjust = -2, margin = margin(0, 0, 0, 0)))

gp3 <- tbl |> filter(sim == 3) |> ggplot(aes(x = size, y = Robs, color = connectance)) + geom_hline(yintercept = 1, color = "grey50", linewidth = 0.25) +
    geom_point(shape = 16, size = 0.25, alpha = 0.1) + geom_line(aes(y = Rmax)) + ggsci::scale_color_npg(guide = "none") +
    ggbreak::scale_x_break(c(310, 900), space = 0.05) + scale_x_continuous(breaks = c(0, 100, 200, 900, 1000)) + scale_y_continuous(limits = c(0, 1.5), breaks = c(0.0, 0.5, 1.0, 1.5)) +
    labs(title = expression(paste((list(mu, italic(P)[bold(A)])), " = (0.20, 0.00)"))) + theme_st() +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 6, hjust = 0.7, vjust = -2, margin = margin(0, 0, 0, 0)))

gp4 <- tbl |> filter(sim == 4) |> ggplot(aes(x = size, y = Robs, color = connectance)) + geom_hline(yintercept = 1, color = "grey50", linewidth = 0.25) +
    geom_point(shape = 16, size = 0.25, alpha = 0.1) + geom_line(aes(y = Rmax)) + ggsci::scale_color_npg(guide = "none") +
    ggbreak::scale_x_break(c(310, 900), space = 0.05) + scale_x_continuous(breaks = c(0, 100, 200, 900, 1000)) + scale_y_continuous(limits = c(0, 1.5), breaks = c(0.0, 0.5, 1.0, 1.5)) +
    labs(title = expression(paste((list(mu, italic(P)[bold(A)])), " = (0.10, 0.25)"))) + theme_st() +
    theme(axis.text.x.top = element_blank(), axis.ticks.x.top = element_blank(), axis.title = element_blank(), plot.title = element_text(size = 6, hjust = 0.7, vjust = -2, margin = margin(0, 0, 0, 0)))

gp <- print((gp1 | gp2) / (gp3 | gp4) + plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0))))
btm <- grid::textGrob(label = expression(paste("Community size ", italic(S))), y = 1.6, gp = grid::gpar(fontsize = unit(units = "char", 7)))
lft <- grid::textGrob(label = expression(paste("Theoretical ", gamma[max])), x = 1.6, rot = 90, gp = grid::gpar(fontsize = unit(units = "char", 7)))
gridExtra::grid.arrange(ggplotify::as.grob(gp), bottom = btm, left = lft, right = lgnd, padding = unit(-0.5, "mm")) |> ggsaver("fig03", width = 8.7, height = 9, ext = "pdf")

