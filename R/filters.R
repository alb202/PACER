filter_unique_positions <- function(gr){
  gr <- sort.GenomicRanges(gr)
  return(sort.GenomicRanges(gr[!(GenomicRanges::duplicated.GenomicRanges(x = gr))]))
}

filter_by_gene <- function(gr, gene_list, invert=FALSE){

  ### If the gene list is empty, just return the full set of alignments
  if(length(gene_list)==1 & gene_list[1]==""){
    return(gr)
  }

  ### Check if the gene list are ensembl IDs or external IDs
  if(isTRUE(is_ensembl_gene_name(gr = gr, gene_list = gene_list))){
    # If the gene list is Ensembl, use the Ensembl ID for comparison
    matches <- mcols(gr)$ensembl_gene_id %in% gene_list
  } else {
    # If the gene list is not Ensembl, use the external ID for comparison
    matches <- mcols(gr)$external_gene_name %in% gene_list
  }

  ### If you want genes that are not in the list, invert the matching index
  if(isTRUE(invert))
    matches <- !matches
  return(gr[matches])
}

is_ensembl_gene_name <- function(gr, gene_list){
  ensembl_matches <- gene_list %in% mcols(gr)$ensembl_gene_id
  external_matches <- gene_list %in% mcols(gr)$external_gene_name
  if(sum(ensembl_matches)>=sum(external_matches)){
    return(TRUE)
  } else{
    return(FALSE)
  }
}

assign_5prime_to_a_length <- function(gr, primary_length){
  # Filter alignments by strand
  pos <- gr[strand(gr)=="+"]
  neg <- gr[strand(gr)=="-"]
  # Filter the primary length alignments
  pos_primary <- pos[width(pos)==primary_length]
  neg_primary <- neg[width(neg)==primary_length]
  # Filter the secondary length alignments
  pos_secondary <- pos[width(pos)!=primary_length]
  neg_secondary <- neg[width(neg)!=primary_length]
  # Filter secondary lengths by overlap
  pos_secondary_filtered <- pos_secondary[-subjectHits(findOverlaps(query = pos_primary, subject = pos_secondary, minoverlap = 1, type = "start", select = "all"))]
  neg_secondary_filtered <- neg_secondary[-subjectHits(findOverlaps(query = neg_primary, subject = neg_secondary, minoverlap = 1, type = "end", select = "all"))]
  # Combine the primary length alignments with the filtered secondary length alignments, sort and return
  return(sort.GenomicRanges(c(pos_primary, neg_primary, pos_secondary_filtered, neg_secondary_filtered)))
}
#
# assign_5prime_to_longer_slower <- function(gr){
#   chromosomes <- sort(unique.Vector(seqnames(gr)))
#   lengths <- sort(unique.Vector(width(gr)), decreasing = TRUE)
#   results <- GAlignments()
#   for(j in 1:length(chromosomes)){
#     #print(chromosomes[j])
#     # Filter alignments by strand
#     pos <- gr[strand(gr)=="+" & seqnames(gr)==chromosomes[j]]
#     neg <- gr[strand(gr)=="-" & seqnames(gr)==chromosomes[j]]
#     #results <- c(gr[width(gr)==lengths[1]])
#     for (i in 1:(length(lengths))){
#       #print(lengths[i])
#       #print(length(pos))
#       #print(length(neg))
#       # Filter the primary length alignments
#       pos_primary <- pos[width(pos)==lengths[i]]
#       neg_primary <- neg[width(neg)==lengths[i]]
#       #print(length(pos_primary))
#       #print(length(neg_primary))
#       # Remove the alignments that share a 5' end with a primary alignment
#       pos <- pos[!start(pos) %in% as.vector(start(pos_primary))]
#       neg <- neg[!end(neg) %in% as.vector(end(neg_primary))]
#       # Concatenate, sort and return results
#       results <- c(results, pos_primary, neg_primary)
#     }
#   }
#   #print(class(results))
#   return(sort.GenomicRanges(results))
# }
assign_5prime_to_longer <- function(gr){
  lengths <- sort(unique.Vector(width(gr)), decreasing = TRUE)
  results <- c(gr[width(gr)==lengths[1]])
  for (i in 1:(length(lengths)-1)){
    #higher <- alignments[qwidth(alignments) %in% lengths[1:i]]
    lower <- gr[width(gr) %in% lengths[(i+1)]]
    #positive <- results[strand(results)=="+"]
    #negative <- results[strand(results)=="-"]
    positive_results <- findOverlaps(query = lower,
                                     subject = results[strand(results)=="+"],
                                     type = "start",
                                     select = "first",
                                     ignore.strand=FALSE)
    negative_results <- findOverlaps(query = lower,
                                     subject = results[strand(results)=="-"],
                                     type = "end",
                                     select = "first",
                                     ignore.strand=FALSE)
    positive_results <- replace(x = positive_results, !is.na(positive_results), FALSE)
    positive_results <- replace(x = positive_results, is.na(positive_results), TRUE)
    negative_results <- replace(x = negative_results, !is.na(negative_results), FALSE)
    negative_results <- replace(x = negative_results, is.na(negative_results), TRUE)
    #full_results <- positive_results & negative_results
    results <- c(results, lower[positive_results & negative_results])
  }
  return(sort.GenomicRanges(results))
}

filter_alignments_by_size <- function(alignments, minimum=10, maximum=30){
  # results <- alignments[qwidth(alignments)>=minimum & qwidth(alignments)<=maximum]
  return(sort.GenomicRanges(alignments[qwidth(alignments)>=minimum & qwidth(alignments)<=maximum]))
}


filter_BAM_tags <- function(gr){
  # Get the index for alignments with no mismatches
  no_mismatches_index <- mcols(gr)$NM==0

  # Get the index for alignments with up to 2 mismatches
  two_mismatches_index <- mcols(gr)$NM<=2

  # Get the index for alignments with no mismatches in the first 22 bases
  MD_split <- strsplit(mcols(gr)$MD, split = "[A-Z]")
  setA <- ifelse(as.numeric(mcols(gr)$NM)<=0, TRUE, FALSE)
  setB <- ifelse(strand(gr)=="+"&unlist(lapply(X = MD_split, FUN = function(x) as.numeric(x[[1]][1])>=22)), TRUE, FALSE)
  setC <- ifelse(strand(gr)=="-"&unlist(lapply(X = MD_split, FUN = function(x) as.numeric(x[[length(x)]][1])>=22)), TRUE, FALSE)
  no_mismatches_in_seed_index <- setA | setB | setC
  return(list(no_mm=no_mismatches_index, two_mm=two_mismatches_index, no_mm_seed=no_mismatches_in_seed_index))
}
#
# filter_MD_tag <- function(strand, NM, MD){
#   if(NM==0)
#     return(TRUE)
#   MD_split <- strsplit(x = MD, split = "[A-Z]")
#   # return(filter_MD_tags3(strand = strand, MD = MD_split))
#   if(strand=="-")
#     MD_split <- rev(MD_split[[1]])
#   if((as.numeric(MD_split[[1]][1])>=as.numeric(22)))
#     return(TRUE)
#   else
#     return(FALSE)
  # MD <- strsplit(x = mcols(gr)$MD, split = "[A-Z]")
  # ((strand(gr)=="+" & unlist(lapply(X = MD, FUN = function(x) as.numeric(x[[1]])))>=22) |
  #    (strand(gr)=="-" & unlist(lapply(X = rev(MD), FUN = function(x) as.numeric(x[[1]])))>=22) | )
  #else()
#   #  return(FALSE)
# }


filter_by_metadata <- function(target, source, column){

  matches <- mcols(target)[,column] %in% mcols(source)[,column]
  results <- target[matches]
  return(sort.GenomicRanges(results))
}

filter_by_regions <- function(gr, regions, type=c("both", "sense", "antisense"), invert=FALSE){
  if (type=="both") {
    results <- subsetByOverlaps(query = gr, subject = regions, invert = invert, ignore.strand=TRUE)
  }
  if (type=="sense") {
    results <- subsetByOverlaps(query = gr, subject = regions, invert = invert, ignore.strand=FALSE)
  }
  if (type=="antisense") {
    strand(regions) <- invert_vector(as.character(strand(regions)))
    results <- subsetByOverlaps(query = gr, subject = regions, invert = invert, ignore.strand=FALSE)
  }
  return(sort.GenomicRanges(results))
}


filter_RNA_from_intervals <- function(gr){
  results <- subset(x = gr,
                    gene_biotype!="snoRNA" &
                      gene_biotype!="miRNA" &
                      gene_biotype!="rRNA" &
                      gene_biotype!="tRNA" &
                      gene_biotype!="snRNA")
  return(sort.GenomicRanges(results))
}

remove_overrepresented_sequences <- function(alignments, cutoff=0.001){
  counts <- rle(sort(as.character(mcols(alignments)$seq)))
  counts_df <- data.frame(values=counts$values, lengths=counts$lengths)
  overreppresented <- subset(counts_df, lengths>(cutoff*length(alignments)))
  '%nin%' <- Negate('%in%')
  results <- subset(alignments, seq %nin% overreppresented$values)
  return(sort.GenomicRanges(results))
}

subsample_gr <- function(gr, size){
  return(sort.GenomicRanges(gr[sample(x = 1:length(gr), size = size, replace = FALSE)]))
}

filter_ambiguous_bases <- function(seqs){
  return(!unlist(mclapply(X = as.character(seqs),
                         FUN = function(x) grepl(pattern = "R|Y|S|W|K|M|B|D|H|V|N|\\.",
                                                 x = x,
                                                 ignore.case = TRUE,
                                                 fixed = FALSE))))
}
# Test the 5' filter
# table(start(two_mm_5prime_filtered[width(two_mm_5prime_filtered)!=22 & strand(two_mm_5prime_filtered)=="+" &
# seqnames(two_mm_5prime_filtered)=="I"]) %in% start(two_mm_5prime_filtered[width(two_mm_5prime_filtered)==22 &
# strand(two_mm_5prime_filtered)=="+" & seqnames(two_mm_5prime_filtered)=="I"]))
