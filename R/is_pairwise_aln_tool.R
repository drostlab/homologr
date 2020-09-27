is_pairwise_aln_tool <- function(tool = NULL){

  provided <- c("NW")

  if(is.null(tool)){
    message("The followig methods are provided: ", paste0(provided, collapse = ", "))
    return(FALSE)
  }

  return(is.element(tool, provided))
}
