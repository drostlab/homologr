#' @title Plot pairwise dN, dS, or dN/dS distributions of orthologous genes
#' @description This function enables users to plot pairwise dN, dS, or dN/dS distributions of orthologous genes.
#' @param dnds_species_tbl a dnds table created with \code{\link{dnds_across_multiple_species}}.
#' @param type type of substitution rate quantification that shall be visualized on the y-axis. Options are:
#' \itemize{
#' \item \code{type = "dn"}
#' \item \code{type = "ds"}
#' \item \code{type = "dn/ds"}
#' \item \code{type = "dn+ds"}
#' }
#' @param order a character vector specifying the order of the species using their scientific names (which have to match the species names stored in the \code{subject_species}
#' column of the input \code{dnds_species_tbl}).
#' @param xlab label of x-axis.
#' @param ylab label of y-axis.
#' @param title title label.
#' @author Hajk-Georg Drost
#' @export

dnds_across_multiple_species_plot <- function(dnds_species_tbl, type, order, xlab = "", ylab = "", title = "") {

  if (!is.element(type, c("dn", "ds", "dn/ds", "dn+ds")))
    stop("Please specify a type that is supported by this function. E.g. dn, ds, dn+ds, or dn/ds.", call. = FALSE)

  subject_species <- dN <- dS <- dNdS <- dN_plus_dS <- ..density.. <- NULL

  if (type == "dn/ds") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = dNdS,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }

  if (type == "dn+ds") {
    blast_tbl <- dplyr::mutate(blast_tbl, dN_plus_dS = dN + dS)
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = dN_plus_dS,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }

  if (type == "dn") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = dN,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }

  if (type == "ds") {
    p <-
      ggplot2::ggplot(
        blast_tbl,
        ggplot2::aes(
          x = dS,
          y = factor(subject_species, levels = rev(order)),
          fill = subject_species,
          colour = subject_species,
          height = ..density..
        )
      )
  }
  p <- p +
    ggridges::geom_density_ridges(
      scale = 5,
      show.legend = FALSE,
      alpha = 0.4,
      size = 1.5,
      stat = "density", trim = TRUE
    ) + ggridges::theme_ridges() + ggsci::scale_colour_lancet() + ggsci::scale_fill_lancet() + ggplot2::labs(x = xlab,
                                                                                                             y = ylab) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 18)) +
    ggplot2::theme(
      axis.title  = ggplot2::element_text(size = 18, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 18, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 18, face = "bold"),
      panel.background = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(
        size = 18,
        colour = "black",
        face = "bold"
      )
    ) + ggplot2::ggtitle(title) + ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 8))

  return(p)
}
