#' Borders of raster objects
#'
#' Identify the borders of raster objects.
#'
#' @param group named list, each element represents a raster object composed of
#'   several raster cells. If there are NA values on the raster surface, raster
#'   cells must be identified by attribute table row indices (each corresponding
#'   to a raster cell) (see \code{\link{attTbl}}).
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}. The list of neighborhoods has to be computed setting
#'   the argument \code{rNumb=TRUE}.
#' @param silent logic, progress bar is not printed on the console.
#'
#' @seealso [attTbl()], [ngbList()], [obj.nbs()]
#'
#' @note \itemize{ \item Note that \code{group} has to be a **named** list whose
#'   names correspond to raster object IDs.
#'
#'   \item If there are NA values on the raster surface, raster cells must be
#'   identified by attribute table row indices (each corresponding to a raster
#'   cell). Row indices can be converted into raster cells using the \code{Cell}
#'   column of the attribute table (e.g. \code{attTbl$Cell[indices]}) (see
#'   \code{\link{attTbl}}).}
#'
#' @return The function returns a named list of object borders. List names
#'   identify the objects; list element values identify the raster cells
#'   comprising the borders.
#'
#' @export
#' @examples
#' # DUMMY DATA
#' ######################################################################################
#' # LOAD LIBRARIES
#' library(scapesClassification)
#' library(terra)
#'
#' # LOAD THE DUMMY RASTER
#' r <- list.files(system.file("extdata", package = "scapesClassification"),
#'                 pattern = "dummy_raster\\.tif", full.names = TRUE)
#' r <- terra::rast(r)
#'
#' # ADD NA-VALUE
#' r[11] <- NA
#'
#' # COMPUTE THE ATTRIBUTE TABLE
#' at <- attTbl(r, "dummy_var")
#'
#' # COMPUTE THE LIST OF NEIGBORHOODS
#' nbs <- ngbList(r, rNumb=TRUE, attTbl=at) # rnumb MUST be true to use obj.border
#'
#' ################################################################################
#' # COMPUTE RASTER OBJECTS
#' ################################################################################
#' at$cv <- anchor.seed(at, nbs, silent=TRUE, class = NULL, rNumb=TRUE,
#'                      cond.filter = "dummy_var > 1",
#'                      cond.seed   = "dummy_var==max(dummy_var)",
#'                      cond.growth = "dummy_var<dummy_var[]",
#'                      lag.growth  = 0)
#'
#' # Raster objects
#' RO <- split(1:NROW(at), at$cv)
#'
#' print(RO) # values are attribute table row indices
#'
#' ################################################################################
#' # COMPUTE BORDERS
#' ################################################################################
#' RO_bd <- obj.border(RO, nbs, silent = TRUE)
#'
#' RO_bd1 <- at$Cell[RO_bd[["1"]]] # Convert row numbers to cell numbers
#' RO_bd2 <- at$Cell[RO_bd[["2"]]] # Convert row numbers to cell numbers
#'
#' print(RO_bd)  # attribute table row indices
#' print(RO_bd1) # cell numbers
#' print(RO_bd2) # cell numbers
#'
#' ################################################################################
#' # PLOT BORDERS
#' ################################################################################
#' plot(cv.2.rast(r,at$cv), type="classes", col=c("#E6E600","#00A600"),
#'      main="Borders")
#' points(terra::xyFromCell(r, RO_bd1), pch=20, col="blue")
#' points(terra::xyFromCell(r, RO_bd2), pch=20, col="red")
#' text(xyFromCell(r, 11), "NA\nvalue")

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


#' Shared borders of raster objects
#'
#' Identify the shared borders of neighboring raster objects.
#'
#' @param grp.bord named list, the list of borders returned by the function
#'   \code{\link{obj.border}}.
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}. The list of neighborhoods has to be computed setting
#'   the argument \code{rNumb=TRUE}.
#' @param only_grp character vector. If \code{NULL}, all IDs in \code{grp.bord}
#'   are considered. If IDs are provided, then they are the only ones
#'   considered.
#' @param silent logic, progress bar is not printed on the console.
#'
#' @seealso [attTbl()], [ngbList()], [obj.border()]
#'
#' @note If there are NA values on the raster surface, raster cells are
#'   identified by attribute table row indices (each corresponding to a raster
#'   cell). Row indices can be converted into raster cells using the \code{Cell}
#'   column of the attribute table (e.g. \code{attTbl$Cell[indices]}) (see
#'   \code{\link{attTbl}}).
#'
#' @return The function returns a named list. Each element represents a raster
#'   object as identified by the list names and contains a nested named list.
#'   The names of the nested lists are the IDs of neighboring raster objects and
#'   their values identify the raster cells comprising the shared borders.
#'
#' @export
#' @examples
#' # DUMMY DATA
#' ######################################################################################
#' # LOAD LIBRARIES
#' library(scapesClassification)
#' library(terra)
#'
#' # LOAD THE DUMMY RASTER
#' r <- list.files(system.file("extdata", package = "scapesClassification"),
#'                 pattern = "dummy_raster\\.tif", full.names = TRUE)
#' r <- terra::rast(r)
#'
#' # ADD ONE NA VALUE
#' r[11] <- NA
#'
#' # COMPUTE THE ATTRIBUTE TABLE
#' at <- attTbl(r, "dummy_var")
#'
#' # COMPUTE THE LIST OF NEIGBORHOODS
#' nbs <- ngbList(r, rNumb=TRUE, attTbl=at) # rnumb MUST be true to use obj.nbs
#'
#' ################################################################################
#' # COMPUTE RASTER OBJECTS
#' ################################################################################
#' at$cv <- anchor.seed(at, nbs, silent=TRUE, class = NULL, rNumb=TRUE,
#'                      cond.filter = "dummy_var > 1",
#'                      cond.seed   = "dummy_var==max(dummy_var)",
#'                      cond.growth = "dummy_var<dummy_var[]",
#'                      lag.growth  = 0)
#'
#' RO <- split(1:NROW(at), at$cv)
#'
#' print(RO)
#'
#' ################################################################################
#' # COMPUTE BORDERS
#' ################################################################################
#' RO_bd <- obj.border(RO, nbs, silent = TRUE)
#'
#' ################################################################################
#' # COMPUTE SHARED BORDERS
#' ################################################################################
#' RO_sbd <- obj.nbs(RO_bd, nbs, silent = TRUE)
#'
#' # Convert row indices to cell numbers
#' RO_sbd1 <- RO_sbd[["1"]]
#' RO_sbd1 <- at$Cell[unlist(RO_sbd1)]
#'
#' RO_sbd2 <- RO_sbd[["2"]]
#' RO_sbd2 <- at$Cell[unlist(RO_sbd2)]
#'
#' # Borders in `RO_sbd` are identified by row indices
#' print(RO_sbd[["1"]]) # Row indices
#' print(RO_sbd1) # Cell numbers
#'
#' print(RO_sbd[["2"]])  # Row indices
#' print(RO_sbd2) # Cell numbers
#'
#' # Neighbor objects
#' names(RO_sbd[["1"]]) # RO1 has one neighbor, RO2
#' names(RO_sbd[["2"]]) # RO2 has one neighbor, RO1
#'
#' ################################################################################
#' # PLOT BORDERS
#' ################################################################################
#' plot(cv.2.rast(r,at$cv), type="classes", col=c("#E6E600","#00A600"),
#'      main="Shared borders")
#' points(terra::xyFromCell(r, RO_sbd1), pch=20, col="blue")
#' points(terra::xyFromCell(r, RO_sbd2), pch=20, col="red")
#' text(xyFromCell(r, 11), "NA\nvalue")

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

