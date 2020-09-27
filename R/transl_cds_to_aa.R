#' @title Translate CDS file to Amino Acids file
#' @description This function takes a \code{data.table} object as input and
#' translates the cds sequences stored as \code{seqs} column into the corresponding amino acid
#' sequence.
#' @param dt a \code{data.table} object storing cds sequences in a column named \code{seqs}.
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @import data.table
transl_cds_to_aa <- function(dt, delete_corrupt_cds = TRUE){

  if (!is.data.table(dt))
    stop("Your CDS file was not corretly transformed into a data.table object.", call. = FALSE)

  # define visible bindings for global variables
  seqs <- aa <- geneids <- NULL

  # When using data.tables within packaes, always make sure
  # 'data.table' is included in the DESCRIPTION file as 'Imports' AND
  # in the NAMESPACE file with 'imports(data.table)' -> @import when using roxygen2
  # http://stackoverflow.com/questions/10527072/using-data-table-package-inside-my-own-package
  # https://github.com/hadley/dplyr/issues/548

  # omit empty sequences
  dt <- dt[ , .SD[sapply(seqs, function(x) {
    return(!(is.na(x) || x == ""))
  })]]

  if (delete_corrupt_cds) {
    # omit sequences that are not multiples of 3
    dt <- dt[ , .SD[sapply(seqs, function(x) {
      return(nchar(x) %% 3 == 0)
    })]]
  }

  # omit sequences consisting of others than ACGT
  dt <- dt[ , .SD[sapply(seqs, is_dna_sequence)]]

  # translate cds to protein sequences
  tryCatch({
    dt[ , aa := transl(seqs), by = geneids]

  }, error = function(e) {
    stop(
      "The input coding sequences could not be translated properly to amino acid sequences.",
      "\n",
      " Please check whether ",
      file,
      " stores valid coding sequences.", call. = FALSE
    )
  })

  return(dt)
}
