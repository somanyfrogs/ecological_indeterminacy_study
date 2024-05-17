#' @title func_manage.R
#' @description An R file describing functions of file/figure management for the Paper entitled:
#'      "Unraveling emergent network indeterminacy in complex ecosystems: a random matrix approach"
#'      Initially written on 20221219 by K.Kawatsu.
#'      Last update: 20231220 by K.Kawatsu.

#' Modified version of as_tibble
#'
#' @inheritParams tibble::as_tibble
#' @param coln A string vector setting column names.
as_tibbler <- function(x, coln = NULL) {
    suppressMessages(tbl <- as_tibble(x, .name_repair = "unique"))
    if(!is.null(coln)) tbl <- tbl |> setNames(coln)
    return(tbl)
}

#' Accessor function of list data
#'
#' @param lst A list.
#' @param name A strings for the name to be accessed.
pullist <- function(lst, name) lst[[name]]

#' Modified version of read_csv
#'
#' @param path A strings, setting the file path to be loaded from.
read_csvr <- function(path, ...) read_csv(path, progress = FALSE, show_col_types = FALSE, ...)

#' Modified version of write_csv
#'
#' @param x A tibble or data.frame to be writtne to.
#' @param path A strings, setting the file path to be written to.
write_csvr <- function(x, path) write_csv(x, path, progress = FALSE)

#' Modified version of ggsave
#'
#' @param gp ggplot object.
#' @param name A strings (file name).
#' @param width A numeric, figure width with cm scale.
#' @param height A numeric, figure height with cm scale.
#' @param ext A strings, specifying the file type (default = "eps").
#' @param rs An integer, specifying the value of fallback_resolution.
ggsaver <- function(gp, name, width, height, ext = "eps", rs = 900) {
    file <- str_c("fig/", name, ".", ext)

    if(ext == "eps") {
        ggsave(file, gp, device = cairo_ps, fallback_resolution = rs, family = "Helvetica", width = width, height = height, units = "cm")
    } else if(ext == "pdf") {
        ggsave(file, gp, device = cairo_pdf, fallback_resolution = rs, family = "Helvetica", width = width, height = height, units = "cm")
    } else {
        ggsave(file, gp, family = "Helvetica", width = width, height = height, units = "cm")
    }
}

#' Original ggplot theme1
#'
#' @param just A numeric vector of legend justification.
#' @param pos A numeric vector of legend position.
#' @param lunit A numeric specifying the line unit.
theme_st <- function(just = NULL, pos = NULL, lunit = 5.5) {
    mgn <- margin(0.5 * lunit, 0.5 * lunit, 0.5 * lunit, 0.5 * lunit)

    th <- theme_light() +
        theme(axis.text = element_text(size = 6),
              axis.title = element_text(size = 7),
              legend.background = element_blank(),
              legend.key.size = unit(0.25, "cm"),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 6, face = "bold"),
              panel.grid.minor = element_blank(),
              plot.tag = element_text(size = 10),
              plot.title = element_text(size = 7, face = "bold"),
              strip.text = element_text(size = 6, color = "black", margin = mgn))

    if(!is.null(just)) th <- th + theme(legend.justification = just)
    if(!is.null(pos)) th <- th + theme(legend.position = pos)
    return(th)
}

#' Original ggplot theme2
#'
#' @inheritParams theme_st
theme_wt0 <- function(just = NULL, pos = NULL) {
    th <- theme_st(just, pos) +
        theme(axis.text = element_blank(),
              axis.ticks = element_blank(),
              legend.title = element_blank(),
              panel.background = element_rect(color = "black", fill = "white"),
              panel.grid = element_blank())
    return(th)
}

#' Original ggplot theme3
#'
#' @inheritParams theme_st
theme_wt1 <- function(just = NULL, pos = NULL) theme_wt0(just, pos) + theme(axis.text = element_blank(), axis.title = element_blank())
