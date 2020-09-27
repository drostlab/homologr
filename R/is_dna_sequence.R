is_dna_sequence <- function(seq){

  return((length(grep("[^(A|C|G|T|a|c|g|t)]", seq)) == 0))
}
