#' @title Internal function for dNdS computations
#' @description This function takes a pairwise alignment as input file and estimates the
#' dNdS ratio of the corresponding alignment. Nevertheless, this function is a helper function for
#' \code{\link{dNdS}}. For dNdS computations you should use the function: \code{\link{dNdS}}.
#' @param file a character string specifying the path to a codon alignment file
#' @param est.method a character string specifying the dNdS estimation method, e.g. "Comeron","Li" .
#' Note, that when using "Comeron" as dNdS estimation method, the program 'gestimator' is used to compute the
#' corresponding dNdS values from a given alignment. The program 'gestimator' can only read "fasta" files,
#' hence it is important to use format = "fasta" when choosing est.method = "Comeron".
#' @param format a character string specifying the file format in which the alignment is stored:
#' "mase", "clustal", "phylip", "fasta" , "msf"
#' @param quiet a logical value specifying whether the output of the coresponding interface shall be printed out.
#' @param subst_name a character string specifying the substitution name that shall be added to the internal folder path naming.
#' Default is \code{subst_name} = \code{NULL}.
#' @author Hajk-Georg Drost
#' @details This function takes a pairwise alignments file as input and estimates dNdS ratios
#' of the corresponding alignments using predefined dNdS estimation methods.
#'
#' The dNdS estimation methods available in this function are:
#'
#' \itemize{
#' \item "Comeron" : Comeron's method (1995)
#' }
#' @references
#' Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol., 36:96-99.
#'
#' Charif, D. and Lobry, J.R. (2007) SeqinR 1.0-2: a contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis.
#'
#' Thornton, K. (2003) libsequence: a C++ class library for evolutionary genetic analysis. Bioinformatics 19(17): 2325-2327
#'
#' Zhang Z, Li J, Zhao XQ, Wang J, Wong GK, Yu J: KaKs Calculator:
#' Calculating Ka and Ks through model selection and model averaging. Genomics Proteomics Bioinformatics 2006 , 4:259-263.
#'
#' https://code.google.com/p/kaks-calculator/wiki/KaKs_Calculator
#'
#' https://code.google.com/p/kaks-calculator/wiki/AXT
#'
#' @return A data.table storing the query_id, subject_id, dN, dS, and dNdS values or
#' a data.table storing the query_id, method, dN, dS, and dNdS values when using KaKs_Calculator.
#' If the dNdS value cannot be calculated NA is returned. This can happen because of constraints
#' of the used model. As each program throws different exception values we set all of them to NA instead.
#' @seealso \code{\link{dNdS}}, \code{\link{multi_aln}}, \code{\link{codon_aln}}, \code{\link{blast_best}},
#' \code{\link{blast_rec}}, \code{\link{read.cds}}
#' @export
substitutionrate <- function(file,
                             est.method,
                             format           = "fasta",
                             quiet            = FALSE,
                             subst_name       = NULL) {


  if (!is_dnds_estimator(est.method))
    stop("Please choose a dNdS estimation method that is supported by this function.",
         call. = FALSE)

  if (!is.element(format, c("mase", "clustal", "phylip", "fasta" , "msf")))
    stop("Please choose a format that is supported by seqinr::read.alignment.",
         call. = FALSE)

  # due to the discussion of no visible binding for global variable for
  # data.table objects see:
  # http://stackoverflow.com/questions/8096313/no-visible-binding-for-global-variable-note-in-r-cmd-check?lq=1
  query_id <- dN <- dS <- NULL

  # determine the file seperator of the current OS
  f_sep <- .Platform$file.sep

  if (!file.exists(file.path(tempdir(), "_calculation"))) {
    dir.create(file.path(tempdir(), "_calculation"))
  }

  if (est.method == "Li") {
    aln <- seqinr::read.alignment(file   = file,
                                  format = format)

    res_list <- seqinr::kaks(aln)

    res <-
      data.table::data.table(t(c(
        unlist(aln$nam),
        res_list$ka,
        res_list$ks,
        (res_list$ka / res_list$ks)
      )))

    # seqinr::kaks documentation:
    # When the alignment does not contain enough information (i.e we approach saturation),
    # the Ka and Ks values take the value 10. Negative values indicate that Ka and Ks can
    # not be computed.
    # we set both of them to NA

    data.table::setnames(
      x   = res,
      old = paste0("V", 1:5),
      new = c("query_id", "subject_id", "dN", "dS", "dNdS")
    )

    data.table::setkey(res, query_id)


    res[which(res[, dN < 0]), dNdS := NA]
    res[which(res[, dN < 0]), dN := NA]

    res[which(res[, dS < 0]), dNdS := NA]
    res[which(res[, dS < 0]), dS := NA]

    res[which(res[, dN > 9.9]), dNdS := NA]
    res[which(res[, dN > 9.9]), dN := NA]

    res[which(res[, dS > 9.9]), dNdS := NA]
    res[which(res[, dS > 9.9]), dS := NA]


    if (!quiet) {
      print("Substitutionrate successfully calculated using Li's method.")
    }

    return(res)
  }

  if (is.element(est.method, kaks_calc_methods)) {
    is_installed_kaks_calculator()

    operating_sys <- Sys.info()[1]

    if (operating_sys == "Darwin") {
      os <- "Mac"
      calc <- "KaKs_Calculator"
    }

    if (operating_sys == "Linux") {
      os <- "Linux"
      calc <- "KaKs_Calculator"
    }

    if (operating_sys == "Windows") {
      os <- "Windows"
      calc <- "KaKs_Calculator.exe"
    }


    # determine the file seperator of the current OS
    f_sep <- .Platform$file.sep
    fa2axt <-
      system.file(paste0("KaKs_Calc_parser", f_sep, "parseFastaIntoAXT.pl"),
                  package = "homologr")

    ### replace paths like this: path/to a folder/of_interest
    ### by: path/'to a folder'/of_interest
    file_name <- unlist(strsplit(file, f_sep))
    file_name <- file_name[length(file_name)]
    curr_wd <- unlist(strsplit(tempdir(), f_sep))
    wdir <- grepl(" ", curr_wd)

    if (any(wdir)) {
      curr_wd[wdir] <- stringr::str_replace(
        string     = curr_wd[wdir],
        replacement = paste0("'", curr_wd[wdir], "'"),
        pattern     = curr_wd[wdir]
      )
    }

    curr_wd <- paste0(curr_wd, collapse = f_sep)

    ###

    tryCatch({
      # converting to input file for KaKs_Calculator
      system(
        paste0(
          "perl ",
          fa2axt ,
          " ",
          file,
          " ",
          curr_wd,
          f_sep,
          "_calculation",
          f_sep,
          file_name
        )
      )


      # running KaKs_Calculator inside the homologr package
      #         KaKs_Calculator <- system.file(paste0("KaKs_Calculator1.2",f_sep,"bin",f_sep,os,f_sep,calc), package = "homologr")

      if (is.null(kaks_calc_path))
        KaKs_Calculator <- calc


      if (!is.null(kaks_calc_path))
        KaKs_Calculator <- paste0(kaks_calc_path, f_sep, calc)


      if (is.null(kaks_calc.params)) {
        system(
          paste0(
            KaKs_Calculator,
            " -i ",
            file.path(tempdir(),
                      "_calculation",
                      paste0(file_name, ".axt")),
            " -o ",
            file.path(
              tempdir(),
              "_calculation",
              paste0(file_name, ".axt.kaks")
            ),
            " -m ",
            est.method
          )
        )
      }


      if (!is.null(kaks_calc.params)) {
        system(
          paste0(
            KaKs_Calculator,
            " -i ",
            file.path(tempdir(),
                      "_calculation",
                      paste0(file_name, ".axt")),
            " -o ",
            file.path(
              tempdir(),
              "_calculation",
              paste0(file_name, ".axt.kaks")
            ),
            " ",
            kaks_calc.params
          )
        )
      }


    }, error = function(e) {
      stop("KaKs_Calculator couldn't run properly, please check your input files.")
    })

    tryCatch({
      kaks_tbl <-
        utils::read.csv(
          file.path(
            tempdir(),
            "_calculation",
            paste0(file_name, ".axt.kaks")
          ),
          sep = "\t",
          header = TRUE
        )
      kaks_tbl_res <- kaks_tbl[, 1:5]
      kaks_tbl_res <-
        data.frame(
          sapply(kaks_tbl_res[, 1], function(x)
            unlist(strsplit(
              as.character(x), "-"
            ))[1]),
          sapply(kaks_tbl_res[, 1], function(x)
            unlist(strsplit(
              as.character(x), "-"
            ))[2]) ,
          kaks_tbl_res[, c(3:5, 2)]
        )

      names(kaks_tbl_res) <-
        c("query_id",
          "subject_id",
          "dN",
          "dS",
          "dNdS",
          "method")
      kaks_tbl_res <- data.table::as.data.table(kaks_tbl_res)
      data.table::setkey(kaks_tbl_res, query_id)

      return(kaks_tbl_res)

    }, error = function(e) {
      stop(
        "Something went wront with KaKs_Calculator.",
        "\n",
        file.path(
          tempdir(),
          "_calculation",
          paste0(file_name, ".axt.kaks")
        ),
        " could not be read properly."
      )
    })
  }
}
