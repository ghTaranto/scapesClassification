#' Standardized relative position index
#'
#' Compute the standardized relative position index of raster objects
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param RO_col column name, the name of the attribute table column with the
#'   raster object IDs.
#' @param el_col column name, the name of the attribute table column with the
#'   elevation values on which the relative position index is computed.
#' @param plot logic, plot the results.
#' @param SpatRaster a \code{SpatRaster} object, the raster used to compute the
#'   attribute table. Required only if \code{plot=TRUE}
#'
#' @details For every raster object mean and standard deviation for elevation
#'   values are computed. Then, the _standardized relative position index_ is
#'   computed by subtracting the mean and dividing by the standard deviation for
#'   each elevation value within raster objects.
#'
#' @return The function returns a vector with standardized relative position
#'   index values. The vector has length equal to the number of rows of the
#'   attribute table. NA values are assigned to cells that do not belong to any
#'   raster object.
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
#' relPI <- rel.pi(attTbl = at, RO_col = "RO", el_col = "dummy_var",
#'                 plot = TRUE, SpatRaster = r)
#'
#' points(terra::xFromCell(r, at$Cell[which(at$RO==1)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==1)]) - 0.04,
#'        pch=20, col="yellow")
#' points(terra::xFromCell(r, at$Cell[which(at$RO==2)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==2)]) - 0.04,
#'        pch=20, col="darkgreen")
#' text(xyFromCell(r,at$Cell), as.character(round(relPI,2)))

rel.pi <- function(attTbl,
                   RO_col,
                   el_col,
                   plot=FALSE,
                   SpatRaster=NULL){

  # TEST FOR COLUMN CELL IN attTbl
  if (!("Cell" %in% names(attTbl))){
    stop("attribute table mising 'Cell' column. Check ?attTbl")
  }

  # TEST ARGUMENT RO_col, regPI_col, el_col and locPI_col
  if( !all(c(RO_col,el_col) %in% names(attTbl)) ){
    stop("'RO_col', 'el_col' must be columns of 'attTbl'")
  }

  ROcell <- split(attTbl[["Cell"]], attTbl[[RO_col]])
  ROcell <- unlist(ROcell)

  relPI         <- numeric(length = NROW(attTbl))
  relPI[]       <- NA
  relPI[ROcell] <- unlist( lapply(split(attTbl[[el_col]], attTbl[[RO_col]]),
                                  function(x) (x-mean(x))/stats::sd(x)) )
  if(plot){

    graphics::layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
    m <- c(1,1,1,1)

    r_RO  <- cv.2.rast(r = SpatRaster, classVector = attTbl[[RO_col]])
    terra::plot(r_RO, type="classes", main="Raster objects", mar=m)

    r_rPI <- cv.2.rast(r = SpatRaster, classVector = relPI)
    terra::plot(r_rPI, type="interval", main="Standardized relPI",mar=m)
  }

  return(relPI)

}
