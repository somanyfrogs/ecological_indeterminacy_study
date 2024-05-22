#' @title Figure_4.R
#' @description Making figure 4 for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221223 by K.Kawatsu.
#'      Last update: 20240522.

## Load R functions
source("R/functions.R")

## Simulation for Figure 4: Theoretical analysis investigating the relationship between interaction mean and interaction correlation (mu & rho)
##      and the maximum absolute eigenvalue in infinite networks
## Set the parameters
sigma <- 0.5
params <- expand_grid(c = c(0.1, 0.5, 0.9, 1.0), mu = seq(-1, 1, length.out = 81), rho = seq(-0.99, 0.99, length.out = 81))

system.time(sim_4 <- foreach(tbl = iter(params, by = "row"), .combine = bind_rows) %do% {
    vars <- calc_vars(Inf, tbl$c, sigma, tbl$mu, tbl$rho)
    Rmax <- criterion_inf(tbl$c, sigma, tbl$mu, tbl$rho)
    tibble(connectance = tbl$c, sigma = sigma, mu = tbl$mu, rho = tbl$rho, Va = vars$Va, Pa = vars$Pa, Rmax = Rmax)
})

sim_4 |> write_csvr("output/sim_4.csv")

## Make figure
tbl <- read_csvr("output/sim_4.csv") |> mutate(ex = if_else(Rmax > 1, "y", "n"))

gp <- tbl |> mutate(connectance = sprintf('italic(C)==%.2f', connectance) |> (\(.x) factor(.x, levels = rev(unique(.x))))()) |>
    ggplot(aes(x = mu, y = rho, fill = ex)) + facet_wrap(. ~ connectance, labeller = label_parsed) + geom_raster(interpolate = TRUE) +
    geom_vline(xintercept = 0, linetype = 2) + geom_hline(yintercept = 0, linetype = 2) + scale_x_continuous(expand = rep(1e-3, 2)) +
    scale_y_continuous(expand = rep(1e-3, 2)) + scale_fill_manual(values = c("grey90", "#c11000"), guide = "none") +
    xlab(expression(paste("Mean interaction strength ", mu))) + ylab(expression(paste("Interaction correlation ", rho))) +
    theme_st() + theme(panel.border = element_rect(color = "black", linewidth = 0.25), panel.spacing = unit(5, "mm"), strip.background = element_blank())

gp |> ggsaver("fig04", width = 8.7, height = 9.4, ext = "pdf")
