#' Relative position index
#'
#' Compute the relative position index of raster objects.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param RO column name, the name of the column with the raster object IDs.
#' @param el column name, the name of column with the elevation values on
#'   which the relative position index is computed.
#' @param type character, defines if position index values are _standardized_
#'   (\code{"s"}) or _normalized_ (\code{"n"}).
#' @param plot logic, plot the results.
#' @param r a \code{SpatRaster} object, the raster used to compute the
#'   attribute table. Required only if \code{plot = TRUE}.
#'
#' @details Position index values are computed only for cells that belong to a
#'   raster object.
#'
#'   * _Standardized position index values_ (\code{type="s"}) are computed with
#'   the formula \code{( x - mean(x) ) / sd(x)};
#'
#'   * _Normalized position index values_ (\code{type="n"}) are computed with
#'   the formula \code{( x - min(x) ) / ( max(x) - min(x) )};
#'
#'   *  Variable \code{x} represents the elevation values
#'   of individual raster object.
#'
#' @return The function returns a vector with relative position index values.
#'   The vector has length equal to the number of rows of the attribute table.
#'   NA values are assigned to cells that do not belong to any raster object.
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
#' # COMPUTE THE ATTRIBUTE TABLE
#' at <- attTbl(r, "dummy_var")
#'
#' # COMPUTE THE LIST OF NEIGBORHOODS
#' nbs <- ngbList(r, rNumb=TRUE, attTbl=at) # rnumb MUST be true to use obj.nbs
#'
#' ################################################################################
#' # COMPUTE RASTER OBJECTS
#' ################################################################################
#' at$RO <- anchor.seed(at, nbs, silent=TRUE, class = NULL, rNumb=TRUE,
#'                      cond.filter = "dummy_var > 1",
#'                      cond.seed   = "dummy_var==max(dummy_var)",
#'                      cond.growth = "dummy_var<dummy_var[]",
#'                      lag.growth  = 0)
#'
#' ################################################################################
#' # STANDARDIZED RELATIVE POSITION INDEX
#' ################################################################################
#' relPI <- rel.pi(attTbl = at, RO = "RO", el = "dummy_var",
#'                 type = "s",
#'                 plot = TRUE, r = r)
#'
#' points(terra::xFromCell(r, at$Cell[which(at$RO==1)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==1)]) - 0.04,
#'        pch=20, col="yellow")
#' points(terra::xFromCell(r, at$Cell[which(at$RO==2)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==2)]) - 0.04,
#'        pch=20, col="darkgreen")
#' text(xyFromCell(r,at$Cell), as.character(round(relPI,2)))
#'
#' ################################################################################
#' # NORMALIZED RELATIVE POSITION INDEX
#' ################################################################################
#' relPI <- rel.pi(attTbl = at, RO = "RO", el = "dummy_var",
#'                 type = "n",
#'                 plot = TRUE, r = r)
#'
#' points(terra::xFromCell(r, at$Cell[which(at$RO==1)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==1)]) - 0.04,
#'        pch=20, col="yellow")
#' points(terra::xFromCell(r, at$Cell[which(at$RO==2)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==2)]) - 0.04,
#'        pch=20, col="darkgreen")
#' text(xyFromCell(r,at$Cell), as.character(round(relPI,2)))

rel.pi <- function(attTbl,
                   RO,
                   el,
                   type = "s",
                   plot = FALSE,
                   r=NULL){

  # TEST FOR COLUMN CELL IN attTbl
  if (!("Cell" %in% names(attTbl))){
    stop("attribute table mising 'Cell' column. Check ?attTbl")
  }

  # TEST ARGUMENT RO, el
  if( !all(c(RO,el) %in% names(attTbl)) ){
    stop("'RO', 'el' must be columns of 'attTbl'")
  }

  # TEST relPI TYPES
  if( !(type %in% c("s","n")) ){
    stop("type must be s (standardized) or n (normalized)")
  }

  # TEST FOR PLOT
  if(is.null(r) & plot){stop("r must be provided if plot = TRUE")}

  ROcell <- split(attTbl[["Cell"]], attTbl[[RO]])
  ROcell <- unlist(ROcell)

  if(type == "s"){ # standardized
    relPI <- unlist( lapply(split(attTbl[[el]], attTbl[[RO]]),
                            function(x) (x-mean(x))/stats::sd(x)) )
  }

  if(type == "n"){ # normalized
    relPI <- unlist( lapply(split(attTbl[[el]], attTbl[[RO]]),
                            function(x) (x-min(x))/(max(x)-min(x)) ))
  }

  attTbl$relPI <- as.numeric(NA)
  attTbl$relPI[match(ROcell, attTbl$Cell)] <- relPI

  if(plot){

    graphics::layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
    m <- c(1,1,1,3)

    r_RO  <- cv.2.rast(r = r, classVector = attTbl[[RO]])
    terra::plot(r_RO, type="classes", main="Raster objects", mar=m,
                plg=list(x=1, y=1, cex=0.9))

    r_rPI <- cv.2.rast(r = r, classVector = attTbl$relPI)

    # breaks
    if(type == "s"){

      brk <- seq(floor(min(relPI, na.rm = T)),
                 ceiling(max(relPI, na.rm = T)),
                 0.5)
      tt <- "Standardized relPI"
    }

    if(type == "n"){
      brk <- seq(0,1,0.1)

      tt <- "Normalized relPI"
    }

    terra::plot(r_rPI, type="interval", main=tt, mar=m, breaks=brk,
                plg=list(x=1, y=1, cex=0.9))
  }

  on.exit(graphics::layout(1))
  return(attTbl$relPI)

}
