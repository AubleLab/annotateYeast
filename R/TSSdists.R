#' Find the distance to the nearest genomic feature
#'
#' For a given query set of genomic regions, and a given feature set of
#' regions, this function will return the distance for each query region to its
#' closest feature. It ignores strand and returns the distance as positive or
#' negative, depending on whether the feature is upstream or downstream
#'
#' This function is similar to the bioconductor distanceToNearest function, but
#' returns negative values for downstream distances instead of absolute values.
#' This allows you to assess the relative location.
#'
#' @param query A data.frame with following 3 columns: chr, start, end
#' @param features A data.frame with gene coordinates. Must contain:
#'                 chr, start, end, gene, geneNum, strand
#'
#' @return A data frame with distances to the nearest gene
#' @export
#' @examples
#' data("testData")
#' data("yeastGenes")
#' calcFeatureDist_aY(testData, yeastGenes)
calcFeatureDist_aY = function(query, features) {
  if (is(query, "GRangesList")) {
    # Recurse over each GRanges object
    x = lapply(query, calcFeatureDist, features)
    return(x)
  }
  
  # # splid the tables by chromosomes and make sure that mitochondrial
  # # chromosome is called the same
  # if ("chrmt" %in% query[, "chr"]){
  #   query[query[,"chr"]=="chrmt","chr"] = "chrM"
  # }
  #
  # if ("chrmt" %in% features[, "chr"]){
  #   features[features[,"chr"]=="chrmt","chr"] = "chrM"
  # }
  
  query = as.data.table(query)
  features = as.data.table(features)
  
  queryDTs = splitDataTable_aY(query, split_factor="chr")
  featureDTs = splitDataTable_aY(features, split_factor="chr")
  
  annotatedPeaks = mapply(queryDTs, featureDTs[names(queryDTs)],
                          FUN=DTNearest_aY)
  if (class(annotatedPeaks) == "list"){
    finaltable = do.call("rbind", annotatedPeaks)
  } else {
    for (i in names(queryDTs)){
      annotTable = DTNearest_aY(queryDTs[[i]], featureDTs[[i]])
      
      if (i == names(queryDTs)[1]){
        finaltable = annotTable
      } else {
        finaltable = rbind(finaltable, annotTable)
      }
    }
  }
  
}

#' Efficiently split a data.table by a column in the table
#'
#' @param DT Data.table to split
#' @param split_factor Column to split, which can be a character vector
#'        or an integer.
#' @return List of data.table objects, split by column
# @examples
# DT = data.table::data.table(letters, grp = rep(c("group1", "group2"), 13))
# splitDataTable(DT, "grp")
# splitDataTable(DT, 2)
splitDataTable_aY = function(DT, split_factor) {
  factor_order = unique(DT[, get(split_factor)])
  if (is.numeric(split_factor)) {
    split_factor = colnames(DT)[split_factor]
    message("Integer split_factor, changed to: ", split_factor)
  }
  l = lapply(split(seq_len(nrow(DT)), DT[, get(split_factor)]),
             function(x) DT[x])
  return(l[factor_order])
}

# Function uses data.table rolling join to identify the nearest features
# really quickly.
#
# @param DT1 A data.table object to be joined to a second data.table object.
# @param DT2 A second data.table object to join with DT1.
#
# @return A rolling joined data.table object.
DTNearest_aY = function(DT1, DT2) {
  if (is.null(DT1)) {
    return(NULL)
  }
  if (is.null(DT2)) {
    return(rep(NA, nrow(DT1)))
  }
  # get middle point in
  DT1[, mid:=start + round((end-start)/2)]
  # mid here is TSS: on + strand start, on - strand end
  DT2[, mid:=ifelse(strand == "+", start, end)]
  data.table::setorder(DT1, mid)
  data.table::setorder(DT2, mid)
  data.table::setattr(DT1, "sorted", "mid")
  data.table::setattr(DT2, "sorted", "mid")
  data.table::setkey(DT1, mid)
  data.table::setkey(DT2, mid)
  merged = DT2[DT1, roll="nearest"]
  merged[, distance:=ifelse(strand == "+", (mid-start), (end-mid))]

  merged = as.data.frame(merged)
  colnames(merged) = c("chrGene", "startGene", "endGene", "gene", "geneNum", "strand",
                       "peakCenter", "chr", "start", "end", "distTSStoPeakCenter")
  return(merged)
}
