is_ortho_detection_method <- function(method = NULL){

  provided <- c("DIAMOND-BH","DIAMOND-RBH")

  if(is.null(method)){
    message("The followig methods are provided: ", paste0(provided, collapse = ", "))
    return(FALSE)
  }

  return(is.element(method, provided))
}
