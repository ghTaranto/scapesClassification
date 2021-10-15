#' Geomorphic Units' Borders
#'
#' Takes as input a list of elements representing geomorphic units and identify
#' the border of each unit.
#'
#' @param group named list, each element represents a geomorphic unit. Each unit
#'   is identified by the indices pointing to the raster cells in the attribute
#'   table that compose that unit (see \code{\link{attTbl}}).
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{attTbl}}). This list has to be
#'   generated setting the argument of the function \code{\link{ngbList}}
#'   \code{index = TRUE} (see \code{\link{ngbList}}).
#' @param silent logic, progress bar is not printed on the console.
#'
#' @seealso [attTbl()], [ngbList()]
#'
#' @note Note that \code{group} has to be a named list whose names correspond to
#'   geomorphic unit-IDs.
#'
#' @return The function returns a named list. Each element of the list
#'   represents a geomorphic unit as identified by the list's names, while the
#'   values of each element of the list consist of indices. These indices point
#'   to the rows in the attribute table that correspond to the raster cells
#'   comprising the borders. The actual raster cells comprising borders can be
#'   extracted from the attribute table (e.g. \code{attTbl$Cell[indices]}) (see
#'   \code{\link{attTbl}}).
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


#' Geomorphic Units' Neighbors
#'
#' Identify all the units adjacent to geomorphic units of interest. For each
#' adjacent unit, individual borders are identified. Both function inputs and
#' outputs use indices pointing to raster cells in the attribute table (see
#' \code{\link{attTbl}} and \code{\link{ngbList}} with the argument \code{index
#' = TRUE}).
#'
#' @param grp.bord named list returned by the function \code{\link{obj.border}}.
#'   Each element of the list represents a geomorphic unit identified by the
#'   list's names, while the values of each element of the list consist of the
#'   indices pointing to the rows in the attribute table that correspond to the
#'   raster cells constituting the border.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{attTbl}}). This list has to be
#'   generated setting the argument of the function \code{\link{ngbList}}
#'   \code{index = TRUE}. (see \code{\link{ngbList}}).
#' @param only_grp character vector, if geomorphic unit-IDs are provided as a
#'   character vector, then these are the only ones considered. If this argument
#'   is set as \code{NULL}, then all the geomorphic units included in the
#'   argument \code{grp.bord} are considered.
#' @param silent logic, progress bar is not printed on the console.
#'
#' @return The function returns a named list. Each element of the list represent
#'   a geomorphic unit as identified by the list's names. Each element is
#'   constituted by a nested named list. The names of the nested lists are the
#'   IDs of all the adjacent geomorphic units. The values of each nested list
#'   consist of indices. These indices point to the rows in the attribute table
#'   (see \code{\link{attTbl}}) that correspond to the raster cells that mark
#'   the border between the geomorphic unit and the adjacent unit. The actual
#'   raster cells comprising borders can be extracted from the attribute table
#'   (e.g. \code{attTbl$Cell[indices]}) (see \code{\link{attTbl}}).
#'
#' @export


obj.nbs <- function(grp.bord, ngbList, only_grp = NULL, silent = FALSE){

  if(is.null(names(grp.bord))){stop("Groups should be named")}

  fct       <- rep(names(grp.bord), lengths(grp.bord))
  grp.unlst <- unlist(grp.bord)

  N <- length(grp.bord)

  grp.nbs <- list()
  if(!is.null(only_grp)){
    grp_list <- only_grp
  } else {
    grp_list <- names(grp.bord)
  }

  for(g in grp_list){ #for1

    if(!silent){
      cat("\r", which(g == names(grp.bord)), "/", N)
    }


    G        <- grp.bord[[g]]

    nbs_G    <- unique(unlist(ngbList[ G ]))

    shar_brd <- intersect(nbs_G, grp.unlst[fct != g])
    nbs_grp  <- unique(fct[grp.unlst %in% shar_brd])


    l <- list()
    for(i in nbs_grp){ #for2

      tmp <- shar_brd[shar_brd %in% grp.bord[[i]]]
      tmp <- unique(unlist(ngbList[ tmp ]))
      l[[i]] <- intersect(tmp, G)

    } #for2

    grp.nbs[[g]] <- l

  } #for1

  return(grp.nbs)

}

