#' Units' Border
#'
#' Identify the border of geomorphic units. It takes as input a list of elements each representing a single unit.
#' Each element can consist of the raster cell numbers composing the geomorphic unit or the indices of the raster cells
#' in the attribure table (see \code{\link{attTbl}})
#'
#'
#' elements represent the
#'
#' @param group named list,
#' @param ngbList
#' @param silent
#'
#' @return
#' @export


obj.border <- function(group, ngbList, silent = FALSE){

  if(is.null(names(group))){stop("Groups should be named")}

  grp.bord <- list()
  N <- length(group)

  for(g in seq_along(group)){

    if(!silent){
      cat("\r", g, "/", N)
    }


    G <- group[[g]]

    if(length(G) == 1){

      grp.bord[[g]] <- G

    } else {

      nbs_G   <- unique(unlist(ngbList[ G ]))
      out_brd <- setdiff(nbs_G, G)

      grp.bord[[g]] <- setdiff(intersect(nbs_G, unlist(ngbList[ out_brd ])), out_brd)

    } #else

  } #for

  names(grp.bord) <- names(group)
  return(grp.bord)

}
