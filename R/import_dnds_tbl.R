#' @title Import a dnds table generated with \code{\link{dnds}}
#' @description This function reads a file that stores a dnds table
#' generated with \code{\link{dnds}}.
#' @param file file path to dnds table with \code{sep} separated columns.
#' @param sep column separator used in the input \code{\link{dnds}} table.
#' @author Hajk-Georg Drost
#' @examples \dontrun{
#' # generate dnds table
#' dNdS_tbl <- dnds(query_file      = system.file('seqs/ortho_thal_cds.fasta', package = 'orthologr'),
#' subject_file    = system.file('seqs/ortho_lyra_cds.fasta', package = 'orthologr'),
#' ortho_detection = "DIAMOND-RBH",
#' aa_aln_type     = "pairwise",
#' aa_aln_tool     = "NW",
#' codon_aln_tool  = "pal2nal",
#' dnds_estimation = "Li",
#' cores           = 1 )
#'
#' # save dnds table as ';' column separated file
#' utils::write.table(
#' dNdS_tbl,
#' file.path(tempdir(), "dNdS_tbl.csv"), sep = ";",
#' col.names = TRUE,
#' row.names = FALSE,
#' quote     = FALSE )
#'
#' # import dNdS table into R session
#' dNdS_tbl_import <- import_dnds_tbl(file.path(tempdir(), "dNdS_tbl.csv"))
#' }
#' @export
import_dnds_tbl <- function(file, sep = ";") {

  if (!file.exists(file))
    stop("The file '",file,"' does not seem to exist. Please specify a valid path to your dnds table.", call. = FALSE)

  res <- tibble::as_tibble(readr::read_delim(
    file,
    col_names = TRUE,
    delim = sep,
    col_types = readr::cols("query_id" = readr::col_character(),
                            "subject_id"= readr::col_character(),
                            "dN" = readr::col_double(),
                            "dS"  = readr::col_double(),
                            "dNdS" = readr::col_double(),
                            "perc_identity" = readr::col_double(),
                            "alig_length" = readr::col_integer(),
                            "mismatches" = readr::col_integer(),
                            "gap_openings" = readr::col_integer(),
                            "q_start" = readr::col_integer(),
                            "q_end" = readr::col_integer(),
                            "s_start" = readr::col_integer(),
                            "s_end" = readr::col_integer(),
                            "evalue" = readr::col_double(),
                            "bit_score" = readr::col_double()
    ))
  )

  if (nrow(res) == 0)
    stop("The file '",file,"' seems to be empty. Please provide a non-empty dNdS file.", call. = FALSE)
  return(res)
}
