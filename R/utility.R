#' Checks class of the list of variables. To be used in functions
#'
#' @param checkList list of object to check, e.g.
#' list(varname=c("data.frame", "numeric")).
#' Multiuple strings in the vector are treated as OR.
#' @return A warning if the wrong input class is provided.
#' @examples
#' x <- function(var1) {
#'     cl = list(var1=c("numeric","character"))
#'     .validateInputs(cl)
#'     return(var1^2)
#' }
.validateInputs_aY <- function(checkList) {
  nms = names(checkList)
  for(i in seq_along(checkList)){
    fail = FALSE
    clss = checkList[[i]]
    x = get(nms[i], envir=parent.frame(1))
    for(cls in clss){
      if (is(x, cls)) fail = append(fail, TRUE)
    }
    if(!any(fail))
      stop(paste0(nms[i], " must be a ", paste(clss, collapse=" or "),
                  ".  Got: ", class(x)))
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
