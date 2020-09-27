#' @title Import protein sequences from a \code{fasta} file
#' @description This function reads an organism specific proteome stored in a defined file format.
#' @param file a character string specifying the path to the file storing the proteome.
#' @param format a character string specifying the file format used to store the proteome, e.g. "fasta", "fastq".
#' @param ... additional arguments that are used by the \code{\link[Biostrings]{readAAStringSet}} function.
#' @author Hajk-Georg Drost
#' @details The \code{read.proteome} function takes a string specifying the path to the proteome file
#' of interest as first argument.
#'
#' It is possible to read in different proteome file standards such as \emph{fasta} or \emph{fastq}.
#'
#' Proteomes stored in fasta files can be downloaded from http://www.ebi.ac.uk/reference_proteomes.
#'
#' @examples \dontrun{
#' # reading a proteome stored in a fasta file
#' Ath_proteome <- import_proteome(system.file('seqs/ortho_thal_aa.fasta', package = 'homologr'),
#'                                format = "fasta")
#'
#' # look at results
#' Ath_proteome
#' }
#'
#' @return A data.table storing the gene id in the first column and the corresponding
#' sequence as string in the second column.
#' @export

import_proteome <- function(file, format, ...){

  if (!file.exists(file))
    stop("Please provide a valid path to the proteome fasta file. Your file '", file,"' does not seem to exist.", call. = FALSE)

  if(!is.element(format,c("fasta","fastq")))
    stop("Please choose a file format that is supported by this function.")

  geneids <- seqs <- NULL


  tryCatch({

    proteome <- Biostrings::readAAStringSet(filepath = file, format = format, ...)
    proteome_names <- as.vector(unlist(sapply(proteome@ranges@NAMES, function(x){return(strsplit(x, " ")[[1]][1])})))
    proteome.dt <- data.table::data.table(geneids = proteome_names,
                                          seqs = toupper(as.character(proteome)))

    data.table::setkey(proteome.dt, geneids)

  }, error = function(e){ stop("File ",file, " could not be read properly.", "\n",
                               "Please make sure that ",file," contains only amino acid sequences and is in ",format," format.")}
  )

  return(proteome.dt)
}
