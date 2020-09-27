is_dnds_estimator <- function(method = NULL){

  provided <- c("Li")

  if (is.null(method)){
    message("The followig methods are provided: ", paste0(provided, collapse = ", "))
    return(FALSE)
  }

  return(is.element(method, provided))
}
