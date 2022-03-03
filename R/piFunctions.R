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
#' graphics::layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
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
#' graphics::layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
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


#' Position index segmentation
#'
#' Segment raster objects based on position index values.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}.
#' @param rNumb logic, the neighborhoods of the argument \code{ngbList} are
#'   identified by cell numbers (\code{rNumb=FALSE}) or by row numbers
#'   (\code{rNumb=TRUE}) (see \code{\link{ngbList}}). It is advised to use row
#'   numbers for large rasters.
#' @param RO column name, the name of the column with the raster object IDs.
#' @param mainPI column name, the name of the column with main position index
#'   values.
#' @param secPI column name, the name of the column with secondary position
#'   index values.
#' @param cut.mPI numeric, threshold of the main position index values below
#'   which cells are excluded from raster objects.
#' @param cut.sPI numeric, threshold of the secondary position index values
#'   below which cells are excluded from raster objects.
#' @param min.N numeric, the minimum number of cells a raster object has to have
#'   to be included in the function output.
#' @param plot logic, plot the results.
#' @param r a \code{SpatRaster} object, the raster used to compute the attribute
#'   table. Required only if \code{plot = TRUE}.
#'
#' @return The function returns a class vector with raster objects IDs. The
#'   vector has length equal to the number of rows of the attribute table. NA
#'   values are assigned to cells that do not belong to any raster object.
#'
#' @details Raster objects are segmented based on position index values. Two
#'   different position indices can be passed to the function (\code{mainPI} and
#'   \code{secPI}).
#'
#'   * Cells with Values below \code{mainPI} **or** below \code{mainPI} are
#'   excluded from raster objects;
#'
#'   * Input raster objects are assigned to the same class to flag cells that
#'   are part of a raster objects;
#'
#'   * Each non-continuous group of raster object cells will identify the
#'   an output raster objects.
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
#' nbs <- ngbList(r, attTbl=at)
#'
#' ################################################################################
#' # COMPUTE RASTER OBJECTS
#' ################################################################################
#' at$RO <- anchor.seed(at, nbs, silent=TRUE, class = NULL, rNumb=TRUE,
#'                      cond.filter = "dummy_var > 1",
#'                      cond.seed   = "dummy_var==max(dummy_var)",
#'                      cond.growth = "dummy_var<dummy_var[]",
#'                      lag.growth  = Inf)
#'
#' # One raster object
#' unique(at$RO)
#'
#' ################################################################################
#' # NORMALIZED RELATIVE POSITION INDEX
#' ################################################################################
#' at$relPI <- rel.pi(attTbl = at, RO = "RO", el = "dummy_var", type = "n")
#'
#' ################################################################################
#' # POSITION INDEX SEGMENTATION
#' ################################################################################
#' RO1 <- pi.sgm(at, nbs,
#'               RO = "RO",        # Raster objects
#'               mainPI = "relPI", # PI segmentation layer
#' #'               cut.mPI = 0,      # segment on relPI values <= 0
#'              plot = TRUE, r = r)
#'
#' graphics::layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
#' text(xyFromCell(r,at$Cell), as.character(at$relPI)) # visualize relPI values
#'
#' # Two raster object
#' unique(at$RO1)

pi.sgm <- function(attTbl,
                   ngbList,
                   rNumb = FALSE,
                   RO,
                   mainPI,
                   secPI = NULL,
                   cut.mPI = NULL,
                   cut.sPI = NULL,
                   min.N = NULL,
                   plot = FALSE,
                   r = NULL){


  # TEST ARGUMENT RO, mainPI
  if( !all(c(RO, mainPI) %in% names(attTbl)) ){
    stop("'RO', 'mainPI' must be columns of 'attTbl'")
  }

  if( !is.null(secPI) ){
    if(!secPI %in% names(attTbl)){
      stop("'secPI' must be a column of 'attTbl'")
    }
  }

  # TEST FOR COLUMN CELL IN attTbl
  if (!("Cell" %in% names(attTbl))){
    stop("attribute table mising 'Cell' column. Check ?attTbl")
  }

  # TEST FOR CORRESPONDENCE attTbl, ngbList
  if (length(ngbList) != nrow(attTbl)) {
    stop("ngbList and attTbl shoud have the same length/nrows")
  }

  # TEST FOR PLOT
  if(is.null(r) & plot){stop("r must be provided if plot = TRUE")}

  # CONVERT NBS FORM CELL IDS TO CELL INDECES
  if(!rNumb){
    fct     <- rep(seq_along(lengths(ngbList)), lengths(ngbList))
    ngbList <- match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- split(ngbList, fct)

    rm(fct, no_nas)
  }

  # SEGMENT ########################################################################
  RO1 <- attTbl[[RO]]

  if(!is.null(cut.mPI) | !is.null(cut.sPI)){

    noNA <- which(!is.na(attTbl[[RO]]))

    # Segmentation rules
    if(is.null(cut.mPI)){cut.mPI <- min(attTbl[[mainPI]], na.rm = TRUE)-0.1}
    if(is.null(cut.sPI)){

      if(is.null(secPI)){secPI <- mainPI}
      cut.sPI <- min(attTbl[[secPI]], na.rm = TRUE)-0.1

    }

    ind <- which( attTbl[[mainPI]][noNA]<cut.mPI | attTbl[[secPI]][noNA]<cut.sPI )
    RO1[ match(attTbl[noNA,][["Cell"]][ind], attTbl$Cell) ] <- NA

    # Compute new raster objects
    RO1 <- anchor.seed(data.frame(Cell=attTbl$Cell, RO1), ngbList,
                       rNumb = TRUE, class = NULL, silent = TRUE,
                       cond.filter = "!is.na(RO1)",
                       cond.seed   = "!is.na(RO1)",
                       cond.growth = "!is.na(RO1)")
  }

  # Remove small raster objects
  if(!is.null(min.N)){

    l <- lengths(split(RO1, RO1)) >= min.N
    l <- as.numeric(names(l)[l])

    RO1[!(RO1 %in% l)] <- NA
  }

  # Plot
  if(plot){
    graphics::layout(matrix(c(1, 2), nrow=1, byrow=TRUE))
    m <- c(1.2,1.2,1.2,3)

    r_RO  <- cv.2.rast(r = r, classVector = attTbl[[RO]])
    terra::plot(r_RO, type="classes", main="Raster objects - Input", mar=m)

    r_RO1 <- cv.2.rast(r = r, attTbl$Cell, classVector = RO1)
    terra::plot(r_RO1, type="classes", main="Raster objects - Output", mar=m)
  }

  on.exit(graphics::layout(1))
  return(RO1)

}
