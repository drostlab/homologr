is_codon_aln_tool <- function(tool = NULL){

  provided <-  c("pal2nal")

  if(is.null(tool)){
    message("The followig methods are provided: ", paste0(provided, collapse = ", "))
    return(FALSE)
  }

  return(is.element(tool, provided))
}
