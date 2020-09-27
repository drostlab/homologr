#' @title Genome-wide pairwise inference of synonymous versus non-synonymous substitution rates (dnds) across multiple subject organisms
#' @description This function allows you to compute dnds maps between a query organism
#' and a set subject organisms stored in the same folder. The corresponding dnds maps are then stored in an output folder.
#' @param query a character string specifying the path to the CDS file of the query organism.
#' @param subjects_folder a character string specifying the path to the folder where CDS files of the subject organisms are stored.
#' @param output_folder a character string specifying the path to the folder where output dnds maps should be stored.
#' @param format a character string specifying the file format of the sequence file, e.g. \code{format = "fasta"}, \code{format = "gbk"}. See
#'  \code{\link{import_proteome}} for more details.
#' @param ortho_detection a character string specifying the orthology inference method that shall be performed to detect orthologous genes. Options are:
#' \itemize{
#' \item \code{ortho_detection ="DIAMOND-BH"}: DIAMOND unidirectional best hit.
#' \item \code{ortho_detection = "DIAMOND-RBH"}: DIAMOND reciprocal/bidirectional best hit (default).
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
#' @param evalue a character string specifying the e-value for DIAMOND based orthology inference that is performed
#' in the process of dnds computations. Please use the scientific notation.
#' @param max_target_seqs maximum number of aligned sequences that shall be retained. Please be aware that \code{max_target_seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .
#' @param min_qry_coverage_hsp minimum \code{qcovhsp} (= query coverage of the HSP) of an orthologous hit (a value between 1 and 100).
#' @param min_qry_perc_identity minimum \code{perc_identity} (= percent sequence identity between query and selected HSP) of an orthologous hit (a value between 1 and 100).
#' @param aa_aln_type a character string specifying the amino acid alignement type: \code{aa_aln_type = "multiple"} or \code{aa_aln_type = "pairwise"}.
#' Default is \code{aa_aln_type = "pairwise"}.
#' @param aa_aln_tool a character string specifying the program that should be used e.g. "clustalw".
#' @param codon_aln_tool a character string specifying the codon alignment tool that shall be used.
#' Default is \code{codon_aln_tool = "pal2nal"}. Right now only "pal2nal" can be selected as codon alignment tool.
#' @param dnds_estimation a character string specifying the dnds estimation method, e.g. \code{dnds_estimation = "Li"} (default).
#' @param cores number of computing cores that shall be used to perform parallelized computations.
#' @param progress_bar should a progress bar be shown. Default is \code{progress_bar = TRUE}.
#' @param sep a file separator that is used to store maps as csv file.
#' @param diamond_hit_tables a file path to a folder where DIAMOND hit tables internally generated with \code{\link{diamond_best_hits}} or \code{\link{diamond_reciprocal_best_hits}} will be stored.
#' @param ... additional parameters that shall be passed to  \code{\link{dnds}}.
#' @details
#' Given a query organism and a set of subject organisms that are stored in the same folder,
#' this function crawls through all subject organism and as a first step computes the pairwise
#' dnds Maps between query and subject organism and as a second step stores the corresponding Map
#' in an output folder.
#' @author Hajk-Georg Drost
#' @seealso \code{\link{dnds}}, \code{\link{import_dnds_across_multiple_species}}, \code{\link{import_dnds_tbl}}
#' @examples
#' \dontrun{
#' # running dnds across several species
#' dnds_across_multiple_species(
#'    query      = system.file('seqs/ortho_thal_cds.fasta', package = 'homologr'),
#'    subjects_folder = system.file('seqs/map_gen_example', package = 'homologr'),
#'    aa_aln_type      = "pairwise",
#'    aa_aln_tool      = "NW",
#'    codon_aln_tool   = "pal2nal",
#'    dnds_estimation  = "Li",
#'    output_folder   = "homologr_dnds_maps",
#'    quiet           = TRUE,
#'    cores      = 1
#' )
#'
#' # running dnds across several species using DIAMOND executable
#' # path '/opt/miniconda3/bin/'
#' dnds_across_multiple_species(
#'    query      = system.file('seqs/ortho_thal_cds.fasta', package = 'homologr'),
#'    subjects_folder = system.file('seqs/map_gen_example', package = 'homologr'),
#'    diamond_exec_path = "/opt/miniconda3/bin/",
#'    aa_aln_type      = "pairwise",
#'    aa_aln_tool      = "NW",
#'    codon_aln_tool   = "pal2nal",
#'    dnds_estimation  = "Li",
#'    output_folder   = "homologr_dnds_maps",
#'    quiet           = TRUE,
#'    cores      = 1
#' )
#' }
#' @export

dnds_across_multiple_species <- function(query,
                               subjects_folder,
                               output_folder,
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
                               min_qry_coverage_hsp = 50,
                               min_qry_perc_identity = 10,
                               aa_aln_type      = "pairwise",
                               aa_aln_tool      = "NW",
                               codon_aln_tool   = "pal2nal",
                               dnds_estimation  = "Li",
                               progress_bar     = TRUE,
                               diamond_hit_tables = file.path(output_folder, "diamond_hit_tables"),
                               sep              = ";",
                               ... ){

  if (!file.exists(subjects_folder))
    stop("The specified folder path '", subjects_folder, "' does not seem to exist. Please provide a valid path to the folder storing all subject genomes.", call. = FALSE)
  # retrieve all subject files within a given folder

  subj.files <- list.files(subjects_folder)

  if (length(subj.files) == 0)
    stop("Your subjects_folder ", subjects_folder, " seems to be empty...", call. = FALSE)

  if (!dplyr::between(min_qry_coverage_hsp, 1, 100))
    stop("Please provide a valid min_qry_coverage_hsp value between 1 and 100.", call. = FALSE)

  if (!dplyr::between(min_qry_perc_identity, 1, 100))
    stop("Please provide a valid min_qry_perc_identity value between 1 and 100.", call. = FALSE)

  message("Starting pairwise genome comparisons (orthology inference and dnds estimation) between query species: ", basename(query), " and subject species: ", paste0(subj.files, collapse = ", "))

  # initialize progress bar
  if (progress_bar & (length(subj.files) > 1))
    pb <- utils::txtProgressBar(1, length(subj.files), style = 3)


  if (!file.exists(output_folder)){
    message("\n")
    message("Creating dnds output folder: ", output_folder)
    dir.create(output_folder)
  }

  if (!file.exists(diamond_hit_tables)){
    message("\n")
    message("Creating DIAMOND hits output folder: ", diamond_hit_tables)
    dir.create(diamond_hit_tables, recursive = TRUE)
  }

  qcovhsp <- perc_identity <- NULL
  message("\n")

  for (i in seq_len(length(subj.files))) {
    # compute pairwise dnds maps between query and all subject files
    OrgQuery_vs_OrgSubj <- dnds(
      query      = query,
      subject    = file.path(subjects_folder, subj.files[i]),
      format          = format,
      ortho_detection = ortho_detection,
      delete_corrupt_cds = delete_corrupt_cds,
      store_locally   = store_locally,
      sensitivity_mode = sensitivity_mode,
      out_format       = out_format,
      evalue            = evalue,
      max_target_seqs   = max_target_seqs,
      cores           = cores,
      hard_mask       = hard_mask,
      diamond_exec_path = diamond_exec_path,
      add_makedb_options = add_makedb_options,
      add_diamond_options = add_diamond_options,
      aa_aln_type      = aa_aln_type,
      aa_aln_tool      = aa_aln_tool,
      codon_aln_tool   = codon_aln_tool,
      dnds_estimation  = dnds_estimation,
      output_path = diamond_hit_tables,
      print_citation = FALSE,
      ...
    )

    message("Filtering for DIAMOND hits with min_qry_coverage_hsp >= ", min_qry_coverage_hsp, " and min_qry_perc_identity >= ", min_qry_perc_identity, " ...")
    OrgQuery_vs_OrgSubj <- dplyr::filter(OrgQuery_vs_OrgSubj, qcovhsp >= min_qry_coverage_hsp, perc_identity >= min_qry_perc_identity)

    dnds_output_i <- file.path(
      output_folder,
      paste0(
        "dnds_query=",
        basename(query),
        "_subject=",
        subj.files[i],
        "_ortho_method=",
        ortho_detection,
        "_evalue=",
        evalue,
        "_dnds_estimation=",
        dnds_estimation,
        "_min_qry_coverage_hsp=",
        min_qry_coverage_hsp,
        "_min_qry_perc_identity=",
        min_qry_perc_identity,
        ".csv"
      ))

    message("\n")
    message("Storing final dnds table in: ", dnds_output_i)

    readr::write_delim(
      OrgQuery_vs_OrgSubj,
      dnds_output_i,
      delim       = sep,
      col_names = TRUE
    )

    if (progress_bar)
      utils::setTxtProgressBar(pb, i)

  }


  message("\n")
  message("Citation:")
  message("Please cite the following paper when using homologr for your own research:")
  message(
    "B Buchfink and HG Drost. Methods in Molecular Biology (2021)."
  )
  message("\n")

  message("All maps are stored in ", output_folder, ".")

}




