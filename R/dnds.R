#' @title Genome-wide pairwise inference of synonymous versus non-synonymous substitution rates (dnds)
#' @description This function takes the CDS files of two organisms of interest (query and subject)
#' and computes the dNdS estimation values for orthologous gene pairs between these organisms.
#' @param query a character string specifying the path to the CDS file of interest (query organism).
#' @param subject a character string specifying the path to the CDS file of interest (subject organism).
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format = "fasta"}, \code{format = "gbk"}. See
#'  \code{\link{import_proteome}} for more details.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed to detect orthologous genes. Options are:
#' \itemize{
#' \item \code{ortho_detection ="DIAMOND-BH"}: DIAMOND unidirectional best hit.
#' \item \code{ortho_detection = "DIAMOND-RBH"}: DIAMOND reciprocal/bidirectional best hit (Default).
#' }
#' @param delete_corrupt_cds a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.
#' @param store_locally a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.
#' @param sensitivity_mode specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
#' \itemize{
#'   \item \code{sensitivity_mode = "fast"} : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70% identity and short read alignment such as NGS reads.
#'   \item \code{sensitivity_mode = "mid-sensitive"} : fast alignments between the \code{fast} mode and the sensitive mode in sensitivity.
#'   \item \code{sensitivity_mode = "sensitive"} : fast alignments, but full sensitivity for hits >40% identity.
#'   \item \code{sensitivity_mode = "more-sensitive"} : more sensitive than the \code{sensitive} mode.
#'   \item \code{sensitivity_mode = "very-sensitive"} : sensitive alignment mode.
#'   \item \code{sensitivity_mode = "ultra-sensitive"} : most sensitive alignment mode (sensitivity as high as BLASTP).
#' }
#' @param out_format a character string specifying the format of the file in which the DIAMOND results shall be stored.
#' Available options are:
#'  \itemize{
#'  \item \code{out_format = "pair"} : Pairwise
#'  \item \code{out_format = "xml"} : XML
#'  \item \code{out_format = "csv"} : Comma-separated file
#'  }
#' @param hard_mask shall low complexity regions be hard masked with TANTAN? Default is \code{db_hard_mask = TRUE}.
#' @param diamond_exec_path a path to the DIAMOND executable or \code{conda/miniconda} folder.
#' @param add_makedb_options a character string specifying additional makedb options that shall be passed on to the diamond makedb command line call, e.g. \code{add_make_options = "--taxonnames"} (Default is \code{add_diamond_options = NULL}).
#' @param add_diamond_options a character string specifying additional diamond options that shall be passed on to the diamond command line call, e.g. \code{add_diamond_options = "--block-size 4.0 --compress 1 --no-self-hits"} (Default is \code{add_diamond_options = NULL}).
#' @param evalue Expectation value (E) threshold for saving hits (default: \code{evalue = 0.001}).
#' @param max_target_seqs maximum number of aligned sequences that shall be retained. Please be aware that \code{max_target_seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param cores number of cores for parallel DIAMOND searches.
#' @param aa_aln_tool a character string specifying the program that should be used e.g. "NW" (Needleman-Wunsch Global Alignments).
#' @param aa_aln_type a character string specifying the amino acid alignement type:
#' \itemize{
#' \item \code{aa_aln_type} = "pairwise" (default)
#' }
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used. Default is \code{codon_aln_tool = "pal2nal"}.
#' Right now only "pal2nal" can be selected as codon alignment tool.
#' @param dnds_estimation the dNdS estimation method that shall be used.
#' Options are:
#' \itemize{
#' \item \code{dnds_estimation = "Li"} (Default): Li's method
#' }
#' @param output_path a path to the location were the DIAMOND best hit output shall be stored. E.g. \code{output_path} = \code{getwd()}
#' to store it in the current working directory, or \code{output_path} = \code{file.path("put", "your", "path", "here")}.
#' @param quiet a logical value specifying whether the output of the corresponding alignment tool shall be printed out to the console.
#' Default is \code{quiet} = \code{FALSE}.
#' @param clean_folders a logical value specifying whether all internal folders storing the output of used programs
#' shall be removed. Default is \code{clean_folders} = \code{FALSE}.
#' @param print_citation a logical value indicating whether or not the citation message shall be printed.
#' @author Hajk-Georg Drost
#' @details
#' The dN/dS ratio quantifies the mode and strength of selection acting on a pair of orthologous genes.
#' This selection pressure can be quantified by comparing synonymous substitution rates (dS) that are assumed to be neutral
#' with nonsynonymous substitution rates (dN), which are exposed to selection as they
#' change the amino acid composition of a protein (Mugal et al., 2013 http://mbe.oxfordjournals.org/content/31/1/212).
#'
#' The \pkg{homologr} package provides the \code{\link{dnds}} function to perform dNdS estimation on pairs of orthologous genes.
#' This function takes the CDS files of two organisms of interest (\code{query} and \code{subject})
#' and computes the dN/dS estimation values for orthologous gene pairs between these organisms.
#'
#' The following pipeline resembles the dN/dS estimation process:
#'
#' \itemize{
#' \item 1) Orthology Inference: e.g. DIAMOND reciprocal best hit (DIAMOND-RBH)
#' \item 2) Pairwise sequence alignment: e.g. global pairwise alignments with Needleman-Wunsch
#' \item 3) Codon Alignment: e.g. pal2nal program
#' \item 4) dNdS estimation: e.g. Li's method
#' }
#' Note: it is assumed that when using \code{dnds()} all corresponding multiple sequence alignment programs you
#' want to use are already installed on your machine and are executable via either
#' the default execution \code{PATH} or you specifically define the location of the executable file
#' via the \code{aa_aln_path} or \code{blast_path} argument that can be passed to \code{dnds()}.
#'
#' The \code{dnds()} function can be used choosing the following options:
#'
#' \itemize{
#' \item \code{ortho_detection} :
#'    \itemize{
#'    \item \code{"DIAMOND-RBH"} (DIAMOND best reciprocal hit)
#'    \item \code{"DIAMOND-BH"} (DIAMOND best hit)
#'      }
#'
#' \item \code{aa_aln_type} :
#'  \itemize{
#'   \item \code{"pairwise"}
#'  }
#'
#' \item \code{aa_aln_tool} :
#'  \item \code{"NW"} (in case \code{aa_aln_type = "pairwise"})
#'
#' \item \code{codon_aln_tool} :
#' \itemize{
#'  \item \code{"pal2nal"}
#'  }
#' \item \code{dnds_estimation} :
#' \itemize{
#' \item "Li" : Li's method
#'
#' }
#'
#' }
#'
#' @references
#' seqinr: \url{http://seqinr.r-forge.r-project.org/}
#'
#' Li, W.-H., Wu, C.-I., Luo, C.-C. (1985) A new method for estimating synonymous and nonsynonymous rates of nucleotide substitution considering the relative likelihood of nucleotide and codon changes. Mol. Biol. Evol, 2:150-174
#'
#' Li, W.-H. (1993) Unbiased estimation of the rates of synonymous and nonsynonymous substitution. J. Mol. Evol., 36:96-99.
#'
#' Pamilo, P., Bianchi, N.O. (1993) Evolution of the Zfx and Zfy genes: Rates and interdependence between genes. Mol. Biol. Evol, 10:271-281
#' @return A tibble storing the dnds values of orthologous genes.
#' @import data.table
#' @examples \dontrun{
#'
#' # get a dnds table using:
#' # 1) reciprocal best hit for orthology inference (DIAMOND-RBH)
#' # 2) Needleman-Wunsch for pairwise amino acid alignments
#' # 3) pal2nal for codon alignments
#' # 4) Li for dNdS estimation
#' # 5) single core processing 'cores = 1'
#' result <- dnds(quer      = system.file('seqs/ortho_thal_cds.fasta', package = 'homologr'),
#'      subject    = system.file('seqs/ortho_lyra_cds.fasta', package = 'homologr'),
#'      ortho_detection = "DIAMOND-RBH",
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW",
#'      codon_aln_tool  = "pal2nal",
#'      dnds_estimation = "Li",
#'      cores      = 1 )
#' # look at output
#' result
#'
#' # specify the path to the DIAMOND executable if it is not in the system path
#' # using the 'diamond_exec_path' argument, e.g. diamond_exec_path = "/opt/miniconda3/bin/"
#' result <- dnds(query      = system.file('seqs/ortho_thal_cds.fasta', package = 'homologr'),
#'      subject    = system.file('seqs/ortho_lyra_cds.fasta', package = 'homologr'),
#'      diamond_exec_path = "/opt/miniconda3/bin/",
#'      ortho_detection = "DIAMOND-RBH",
#'      aa_aln_type     = "pairwise",
#'      aa_aln_tool     = "NW",
#'      codon_aln_tool  = "pal2nal",
#'      dnds_estimation = "Li",
#'      cores      = 1 )
#' # look at output
#' result
#' }
#' @seealso \code{\link{diamond_best_hits}}, \code{\link{diamond_reciprocal_best_hits}},
#' \code{\link{dnds_across_multiple_species}}, \code{\link{import_dnds_tbl}}
#' @export

dnds <- function(query,
                 subject,
                 format          = "fasta",
                 ortho_detection = "DIAMOND-RBH",
                 delete_corrupt_cds = FALSE,
                 store_locally   = FALSE,
                 sensitivity_mode = "ultra-sensitive",
                 out_format       = "csv",
                 evalue            = "1E-5",
                 max_target_seqs   = 5000,
                 cores           = 1,
                 hard_mask       = TRUE,
                 diamond_exec_path = NULL,
                 add_makedb_options = NULL,
                 add_diamond_options = NULL,
                 aa_aln_type     = "pairwise",
                 aa_aln_tool     = "NW",
                 codon_aln_tool  = "pal2nal",
                 dnds_estimation = "Li",
                 output_path = NULL,
                 quiet           = TRUE,
                 clean_folders   = FALSE,
                 print_citation = TRUE) {

  # determine the file separator of the current OS
  f_sep <- .Platform$file.sep

  message(
    "Starting orthology inference (",
    ortho_detection,
    ") and dN/dS estimation (",
    dnds_estimation,
    ") using the follwing parameters:"
  )
  message("query = '", basename(query), "'")
  message("subject = '", basename(subject), "'")
  message("e-value: ", evalue)
  message("aa_aln_type = '", aa_aln_type, "'")
  message("aa_aln_tool = '", aa_aln_tool, "'")
  message("cores = '", cores, "'")
  message("\n")

  if (!is_ortho_detection_method(ortho_detection))
    stop("Please choose a orthology detection method that is supported by this function.",
         call. = FALSE)

  if (!is.element(aa_aln_type, c("pairwise")))
    stop("Please choose the supported alignement type: 'pairwise'",
         call. = FALSE)


  if (store_locally) {
    if (file.exists("homologr_alignment_files")) {
      message(
        "The folder 'homologr_alignment_files' seems to exist already.",
        "Please make sure to delete this folder or store the previous results in a different location if you don't want files to be overwritten.",
        "The existing folder 'homologr_alignment_files' is used to store alignment files ..."
      )
    } else {
      message("Creating folder 'homologr_alignment_files' to store alignment files ...")
      dir.create("homologr_alignment_files")
      if (!file.exists(file.path(
        "homologr_alignment_files",
        "_pairwise_alignment_with_score"
      )))
        dir.create(file.path(
          "homologr_alignment_files",
          "_pairwise_alignment_with_score"
        ))
    }
  }

  if (!store_locally) {
    if (!file.exists(file.path(tempdir(), "_pairwise_alignment_with_score")))
      dir.create(file.path(tempdir(), "_pairwise_alignment_with_score"))
  }

  if (aa_aln_type == "pairwise") {
    if (!is_pairwise_aln_tool(aa_aln_tool))
      stop(
        "Please choose a pairwise alignment tool that is supported by this function." ,
        call. = FALSE
      )
  }

  if (!is_codon_aln_tool(codon_aln_tool))
    stop("Please choose a codon alignment tool that is supported by this function.",
         call. = FALSE)

  if (!is_dnds_estimator(dnds_estimation))
    stop("Please choose a dNdS estimation method that is supported by this function.",
         call. = FALSE)

  aa <- geneids <- NULL

  # run DIAMOND on each translated amino acid sequence against the related subject genome
  # to retrieve a hit table with pairs of geneids

  message("Step 1 - Running orthology inference using diamond ...")
  # use BLAST best hit as orthology inference method
  if (ortho_detection == "DIAMOND-BH") {

    qry_cds <- Biostrings::readDNAStringSet(query, format = "fasta")
    sbj_cds <- Biostrings::readDNAStringSet(subject, format = "fasta")

    qry_translated_path <- file.path(tempdir(), paste0("translated_", basename(query)))
    sbj_translated_path <- file.path(tempdir(), paste0("translated_", basename(subject)))

    Biostrings::writeXStringSet(Biostrings::translate(qry_cds),
                                filepath = qry_translated_path,
                                format = "fasta")
    Biostrings::writeXStringSet(Biostrings::translate(sbj_cds),
                                filepath = sbj_translated_path,
                                format = "fasta")

    # dnds() needs CDS files as input!
    diamond_output <-
      diamond_best_hits(
        query = qry_translated_path,
        subject = sbj_translated_path,
        is_subject_db = FALSE,
        format = format,
        sensitivity_mode = sensitivity_mode,
        out_format = out_format,
        evalue = evalue,
        max_target_seqs = max_target_seqs,
        cores = cores,
        hard_mask = hard_mask,
        diamond_exec_path = diamond_exec_path,
        add_makedb_options = add_makedb_options,
        add_diamond_options = add_diamond_options,
        output_path = output_path
      )

    data.table::setDT(diamond_output)
    data.table::setkeyv(diamond_output, c("query_id", "subject_id"))

    q_cds <- import_cds(file   = query,
                      format = format,
                      delete_corrupt_cds = delete_corrupt_cds)

    s_cds <- import_cds(file   = subject,
                      format = format,
                      delete_corrupt_cds = delete_corrupt_cds)

    filename_qry <-
      unlist(strsplit(
        query,
        f_sep,
        fixed = FALSE,
        perl = TRUE,
        useBytes = FALSE
      ))

    filename_qry <- filename_qry[length(filename_qry)]

    # input = paste0("query_", filename_qry, ".fasta")

    q_aa <-
      import_proteome(file = qry_translated_path,
                    format = "fasta")

    # translate coding sequences to amino acid sequences
    s_aa_tmp <-
      transl_cds_to_aa(s_cds, delete_corrupt_cds = delete_corrupt_cds)

    filename_subj <-
      unlist(strsplit(
        subject,
        f_sep,
        fixed = FALSE,
        perl = TRUE,
        useBytes = FALSE
      ))
    filename_subj <-
      filename_subj[length(filename_subj)]

    # seqinr::write.fasta(
    #   sequences = as.list(s_aa_tmp[ , aa]),
    #   names     = s_aa_tmp[ , geneids],
    #   nbchar    = 80,
    #   open      = "w",
    #   file.out  = file.path(tempdir(), "_diamond_db", filename_subj)
    # )

    s_aa <-
      import_proteome(file = sbj_translated_path,
                    format = "fasta")

  }

  # run DIAMOND best reciprocal hit as orthology inference method
  if (ortho_detection == "DIAMOND-RBH") {

    qry_cds <- Biostrings::readDNAStringSet(query, format = "fasta")
    sbj_cds <- Biostrings::readDNAStringSet(subject, format = "fasta")

    qry_translated_path <- file.path(tempdir(), paste0("translated_", basename(query)))
    sbj_translated_path <- file.path(tempdir(), paste0("translated_", basename(subject)))

    Biostrings::writeXStringSet(Biostrings::translate(qry_cds),
                                filepath = qry_translated_path,
                                format = "fasta")
    Biostrings::writeXStringSet(Biostrings::translate(sbj_cds),
                                filepath = sbj_translated_path,
                                format = "fasta")
    # dnds() needs CDS files as input!
    diamond_output <-
      diamond_reciprocal_best_hits(
        query = qry_translated_path,
        subject = sbj_translated_path,
        is_subject_db = FALSE,
        format = format,
        sensitivity_mode = sensitivity_mode,
        out_format = out_format,
        evalue = evalue,
        max_target_seqs = max_target_seqs,
        cores = cores,
        hard_mask = hard_mask,
        diamond_exec_path = diamond_exec_path,
        add_makedb_options = add_makedb_options,
        add_diamond_options = add_diamond_options,
        output_path = output_path
      )

    data.table::setDT(diamond_output)
    data.table::setkeyv(diamond_output, c("query_id", "subject_id"))

    q_cds <- import_cds(file   = query,
                      format = format,
                      delete_corrupt_cds = delete_corrupt_cds)

    s_cds <- import_cds(file   = subject,
                      format = format,
                      delete_corrupt_cds = delete_corrupt_cds)

    filename_qry <-
      unlist(strsplit(
        query,
        f_sep,
        fixed = FALSE,
        perl = TRUE,
        useBytes = FALSE
      ))

    filename_qry <- filename_qry[length(filename_qry)]

    # input_qry = paste0("query_", filename_qry, ".fasta")

    q_aa <-
      import_proteome(file = qry_translated_path,
                    format = "fasta")

    filename_subj <-
      unlist(strsplit(
        subject,
        f_sep,
        fixed = FALSE,
        perl = TRUE,
        useBytes = FALSE
      ))
    filename_subj <-
      filename_subj[length(filename_subj)]

    # input_subj = paste0("query_", filename_subj, ".fasta")

    s_aa <-
      import_proteome(file = sbj_translated_path,
                    format = "fasta")

  }

  data.table::setnames(
    q_cds,
    old = c("geneids", "seqs"),
    new = c("query_id", "query_cds")
  )
  data.table::setnames(
    s_cds,
    old = c("geneids", "seqs"),
    new = c("subject_id", "subject_cds")
  )
  data.table::setnames(q_aa,
                       old = c("geneids", "seqs"),
                       new = c("query_id", "query_aa"))
  data.table::setnames(
    s_aa,
    old = c("geneids", "seqs"),
    new = c("subject_id", "subject_aa")
  )

  # joining all tables to a final table containing: query_id, subject_id, query_aa_seq, subject_aa_seq, query_cds_seq, and subject_cds_seq
  query_tbl <-
    dplyr::inner_join(tibble::as_tibble(q_cds), tibble::as_tibble(q_aa), by = "query_id")
  subject_tbl <-
    dplyr::inner_join(tibble::as_tibble(s_cds), tibble::as_tibble(s_aa), by = "subject_id")
  joint_query_tbl <-
    dplyr::inner_join(tibble::as_tibble(diamond_output),
                      tibble::as_tibble(query_tbl),
                      by = "query_id")
  joint_subject_tbl <-
    dplyr::inner_join(tibble::as_tibble(diamond_output),
                      tibble::as_tibble(subject_tbl),
                      by = "subject_id")
  final_tbl <-
    dplyr::inner_join(
      tibble::as_tibble(joint_query_tbl),
      tibble::as_tibble(joint_subject_tbl),
      by = c("query_id", "subject_id")
    )

  if (nrow(final_tbl) == 0)
    stop("No orthologs could be found! Please check your input files!",
         call. = FALSE)

  message("\n")
  message("Orthology inference successfully completed.")
  message("\n")
  message("Step 2 - Running dnds estimation ...")
  dnds_output <-
    internal_dnds(
      complete_tbl    = data.table::as.data.table(final_tbl),
      aa_aln_type     = aa_aln_type,
      aa_aln_tool     = aa_aln_tool,
      codon_aln_tool  = codon_aln_tool,
      store_locally   = store_locally,
      dnds_estimation = dnds_estimation,
      quiet           = quiet,
      cores      = cores,
      clean_folders   = clean_folders
    )

  subject_id <- NULL
  hit.table_selected <- dplyr::select(diamond_output,-subject_id)

  res <-
    dplyr::inner_join(dnds_output, hit.table_selected, by = "query_id")

  message("\n")
  message("----------- dnds estimation successfully completed -----------")
  message("\n")

  if (print_citation) {
    message("\n")
    message("Citation:")
    message("Please cite the following paper when using homologr for your own research:")
    message(
      "B Buchfink and HG Drost. Methods in Molecular Biology (2021)."
    )
    message("\n")
  }

  # return the dNdS table for all query_ids and subject_ids
  return(res)
}
