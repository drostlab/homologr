% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dnds_across_multiple_species.R
\name{dnds_across_multiple_species}
\alias{dnds_across_multiple_species}
\title{Genome-wide pairwise inference of synonymous versus non-synonymous substitution rates (dnds) across multiple subject organisms}
\usage{
dnds_across_multiple_species(
  query,
  subjects_folder,
  output_folder,
  format = "fasta",
  ortho_detection = "DIAMOND-RBH",
  delete_corrupt_cds = FALSE,
  store_locally = FALSE,
  sensitivity_mode = "ultra-sensitive",
  out_format = "csv",
  evalue = "1E-5",
  max_target_seqs = 5000,
  cores = 1,
  hard_mask = TRUE,
  diamond_exec_path = NULL,
  add_makedb_options = NULL,
  add_diamond_options = NULL,
  min_qry_coverage_hsp = 50,
  min_qry_perc_identity = 10,
  aa_aln_type = "pairwise",
  aa_aln_tool = "NW",
  codon_aln_tool = "pal2nal",
  dnds_estimation = "Li",
  progress_bar = TRUE,
  diamond_hit_tables = file.path(output_folder, "diamond_hit_tables"),
  sep = ";",
  ...
)
}
\arguments{
\item{query}{a character string specifying the path to the CDS file of the query organism.}

\item{subjects_folder}{a character string specifying the path to the folder where CDS files of the subject organisms are stored.}

\item{output_folder}{a character string specifying the path to the folder where output dnds maps should be stored.}

\item{format}{a character string specifying the file format of the sequence file, e.g. \code{format = "fasta"}, \code{format = "gbk"}. See
\code{\link{import_proteome}} for more details.}

\item{ortho_detection}{a character string specifying the orthology inference method that shall be performed to detect orthologous genes. Options are:
\itemize{
\item \code{ortho_detection ="DIAMOND-BH"}: DIAMOND unidirectional best hit.
\item \code{ortho_detection = "DIAMOND-RBH"}: DIAMOND reciprocal/bidirectional best hit (default).
}}

\item{delete_corrupt_cds}{a logical value indicating whether sequences with corrupt base triplets should be removed from the input \code{file}. This is the case when the length of coding sequences cannot be divided by 3 and thus the coding sequence contains at least one corrupt base triplet.}

\item{store_locally}{a logical value indicating whether or not alignment files shall be stored locally rather than in \code{tempdir()}.}

\item{sensitivity_mode}{specify the level of alignment sensitivity. The higher the sensitivity level, the more deep homologs can be found, but at the cost of reduced computational speed.
\itemize{
  \item \code{sensitivity_mode = "fast"} : fastest alignment mode, but least sensitive (default). Designed for finding hits of >70% identity and short read alignment such as NGS reads.
  \item \code{sensitivity_mode = "mid-sensitive"} : fast alignments between the \code{fast} mode and the sensitive mode in sensitivity.
  \item \code{sensitivity_mode = "sensitive"} : fast alignments, but full sensitivity for hits >40% identity.
  \item \code{sensitivity_mode = "more-sensitive"} : more sensitive than the \code{sensitive} mode.
  \item \code{sensitivity_mode = "very-sensitive"} : sensitive alignment mode.
  \item \code{sensitivity_mode = "ultra-sensitive"} : most sensitive alignment mode (sensitivity as high as BLASTP).
}}

\item{out_format}{a character string specifying the format of the file in which the DIAMOND results shall be stored.
Available options are:
 \itemize{
 \item \code{out_format = "pair"} : Pairwise
 \item \code{out_format = "xml"} : XML
 \item \code{out_format = "csv"} : Comma-separated file
 }}

\item{evalue}{a character string specifying the e-value for DIAMOND based orthology inference that is performed
in the process of dnds computations. Please use the scientific notation.}

\item{max_target_seqs}{maximum number of aligned sequences that shall be retained. Please be aware that \code{max_target_seqs} selects best hits based on the database entry and not by the best e-value. See details here: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty833/5106166 .}

\item{cores}{number of computing cores that shall be used to perform parallelized computations.}

\item{hard_mask}{shall low complexity regions be hard masked with TANTAN? Default is \code{db_hard_mask = TRUE}.}

\item{diamond_exec_path}{a path to the DIAMOND executable or \code{conda/miniconda} folder.}

\item{add_makedb_options}{a character string specifying additional makedb options that shall be passed on to the diamond makedb command line call, e.g. \code{add_make_options = "--taxonnames"} (Default is \code{add_diamond_options = NULL}).}

\item{add_diamond_options}{a character string specifying additional diamond options that shall be passed on to the diamond command line call, e.g. \code{add_diamond_options = "--block-size 4.0 --compress 1 --no-self-hits"} (Default is \code{add_diamond_options = NULL}).}

\item{min_qry_coverage_hsp}{minimum \code{qcovhsp} (= query coverage of the HSP) of an orthologous hit (a value between 1 and 100).}

\item{min_qry_perc_identity}{minimum \code{perc_identity} (= percent sequence identity between query and selected HSP) of an orthologous hit (a value between 1 and 100).}

\item{aa_aln_type}{a character string specifying the amino acid alignement type: \code{aa_aln_type = "multiple"} or \code{aa_aln_type = "pairwise"}.
Default is \code{aa_aln_type = "pairwise"}.}

\item{aa_aln_tool}{a character string specifying the program that should be used e.g. "clustalw".}

\item{codon_aln_tool}{a character string specifying the codon alignment tool that shall be used.
Default is \code{codon_aln_tool = "pal2nal"}. Right now only "pal2nal" can be selected as codon alignment tool.}

\item{dnds_estimation}{a character string specifying the dnds estimation method, e.g. \code{dnds_estimation = "Li"} (default).}

\item{progress_bar}{should a progress bar be shown. Default is \code{progress_bar = TRUE}.}

\item{diamond_hit_tables}{a file path to a folder where DIAMOND hit tables internally generated with \code{\link{diamond_best_hits}} or \code{\link{diamond_reciprocal_best_hits}} will be stored.}

\item{sep}{a file separator that is used to store maps as csv file.}

\item{...}{additional parameters that shall be passed to  \code{\link{dnds}}.}
}
\description{
This function allows you to compute dnds maps between a query organism
and a set subject organisms stored in the same folder. The corresponding dnds maps are then stored in an output folder.
}
\details{
Given a query organism and a set of subject organisms that are stored in the same folder,
this function crawls through all subject organism and as a first step computes the pairwise
dnds Maps between query and subject organism and as a second step stores the corresponding Map
in an output folder.
}
\examples{
\dontrun{
# running dnds across several species
dnds_across_multiple_species(
   query      = system.file('seqs/ortho_thal_cds.fasta', package = 'homologr'),
   subjects_folder = system.file('seqs/map_gen_example', package = 'homologr'),
   aa_aln_type      = "pairwise",
   aa_aln_tool      = "NW",
   codon_aln_tool   = "pal2nal",
   dnds_estimation  = "Li",
   output_folder   = "homologr_dnds_maps",
   quiet           = TRUE,
   cores      = 1
)

# running dnds across several species using DIAMOND executable
# path '/opt/miniconda3/bin/'
dnds_across_multiple_species(
   query      = system.file('seqs/ortho_thal_cds.fasta', package = 'homologr'),
   subjects_folder = system.file('seqs/map_gen_example', package = 'homologr'),
   diamond_exec_path = "/opt/miniconda3/bin/",
   aa_aln_type      = "pairwise",
   aa_aln_tool      = "NW",
   codon_aln_tool   = "pal2nal",
   dnds_estimation  = "Li",
   output_folder   = "homologr_dnds_maps",
   quiet           = TRUE,
   cores      = 1
)
}
}
\seealso{
\code{\link{dnds}}, \code{\link{import_dnds_across_multiple_species}}, \code{\link{import_dnds_tbl}}
}
\author{
Hajk-Georg Drost
}
