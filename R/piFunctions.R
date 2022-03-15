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
#' @param r single or multi-layer raster of the class \code{SpatRaster} (see
#'   \code{help("rast", terra)}) used to compute the attribute table. Required
#'   only if \code{plot = TRUE}.
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
#' @seealso [attTbl()], [ngbList()], [pi.add()], [pi.sgm()]
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
#' # Convert class vector at$RO to raster and plot
#' r_RO  <- cv.2.rast(r = r, classVector = at$RO)
#' terra::plot(r_RO, type="classes", main="Raster objects",
#'             plg=list(x=1, y=1, cex=0.9))
#'
#' ################################################################################
#' # STANDARDIZED RELATIVE POSITION INDEX
#' ################################################################################
#' relPI <- rel.pi(attTbl = at, RO = "RO", el = "dummy_var",
#'                 type = "s",
#'                 plot = TRUE, r = r)
#'
#' # Annotate relPI
#' points(terra::xFromCell(r, at$Cell[which(at$RO==1)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==1)]) - 0.04,
#'        pch=20, col="yellow")
#' points(terra::xFromCell(r, at$Cell[which(at$RO==2)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==2)]) - 0.04,
#'        pch=20, col="darkgreen")
#' text(xyFromCell(r,at$Cell), as.character(round(relPI,2)))
#' legend(1.02, 0.4, legend=c("1", "2"), bty = "n", title="RO:", xpd=TRUE,
#' col=c("#E6E600", "#00A600"), pch=20, cex=0.9, pt.cex = 1.5)
#'
#' ################################################################################
#' # NORMALIZED RELATIVE POSITION INDEX
#' ################################################################################
#' # Compute normalized relative position index
#' relPI <- rel.pi(attTbl = at, RO = "RO", el = "dummy_var",
#'                 type = "n",
#'                 plot = TRUE, r = r)
#'
#' # Annotate relPI
#' points(terra::xFromCell(r, at$Cell[which(at$RO==1)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==1)]) - 0.04,
#'        pch=20, col="yellow")
#' points(terra::xFromCell(r, at$Cell[which(at$RO==2)]),
#'        terra::yFromCell(r, at$Cell[which(at$RO==2)]) - 0.04,
#'        pch=20, col="darkgreen")
#' text(xyFromCell(r,at$Cell), as.character(round(relPI,2)))
#' legend(1.02, 0.4, legend=c("1", "2"), bty = "n", title="RO:", xpd=TRUE,
#' col=c("#E6E600", "#00A600"), pch=20, cex=0.9, pt.cex = 1.5)

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

    terra::plot(r_rPI, type="interval", main=tt, breaks=brk,
                plg=list(x=terra::ext(r_rPI)[2], y=terra::ext(r_rPI)[4], cex=0.9))
  }

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
#' @param cut.mPI numeric, threshold of main position index values. Cells with
#'   values below the threshold are excluded from raster objects.
#' @param cut.sPI numeric, threshold of secondary position index values. Cells
#'   with values below the threshold are excluded from raster objects.
#' @param min.N numeric, the minimum number of cells a raster object has to have
#'   to be included in the function output.
#' @param plot logic, plot the results.
#' @param r single or multi-layer raster of the class \code{SpatRaster} (see
#'   \code{help("rast", terra)}) used to compute the attribute table. Required
#'   only if \code{plot = TRUE}.
#'
#' @return The function returns a class vector with raster objects IDs. The
#'   vector has length equal to the number of rows of the attribute table. NA
#'   values are assigned to cells that do not belong to any raster object.
#'
#' @details Raster objects are segmented based on position index values. Two
#'   different position indices can be passed to the function (\code{mainPI} and
#'   \code{secPI}).
#'
#'   * Input raster objects are assigned to the same class to flag cells that
#'   are part of raster objects;
#'
#'   * Cells with values below \code{mainPI} **OR** below \code{mainPI} are
#'   flagged as not being part of any raster object;
#'
#'   * Each non-continuous group of raster object cells will identify an output
#'   raster object.
#'
#'   * Only raster objects with at least as many cells as specified by the
#'   argument \code{min.N} are included in the function output.
#'
#'   * If both \code{mainPI} and \code{secPI} are equal to \code{NULL}, the
#'   function will exclusively filter raster objects based on their size
#'   (\code{min.N}).
#'
#' @seealso [attTbl()], [ngbList()], [rel.pi()], [pi.add()]
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
#' # One input raster object
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
#'               cut.mPI = 0,      # segment on relPI values <= 0
#'               plot = FALSE, r = r)
#'
#' ################################################################################
#' # PLOT
#' ################################################################################
#' # Convert class vectors to raster
#' r_RO  <- cv.2.rast(r = r, classVector = at$RO)
#' r_RO1 <- cv.2.rast(r = r, classVector = RO1)
#'
#' # Plot
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(4.5, 0.5, 2, 3.2)
#'
#' terra::plot(r_RO, type="classes", main="Raster objects - Input", mar=m,
#' plg=list(x=1, y=1, cex=0.9))
#'
#' terra::plot(r_RO1, type="classes", main="Raster objects - Output", mar=m,
#'             plg=list(x=1, y=1, cex=0.9))
#' text(xyFromCell(r,at$Cell), as.character(round(at$relPI,2))) # visualize relPI
#' text(0.01, 1, "Cut on relPI <= 0", adj=c(0,1), cex = 0.8)
#' par(oldpar)
#'
#' # Two output raster objects
#' unique(RO1)

pi.sgm <- function(attTbl,
                   ngbList,
                   rNumb = FALSE,
                   RO,
                   mainPI = NULL,
                   secPI = NULL,
                   cut.mPI = NULL,
                   cut.sPI = NULL,
                   min.N = NULL,
                   plot = FALSE,
                   r = NULL){


  if( !RO %in% names(attTbl) ){
    stop("'RO', 'mainPI' must be columns of 'attTbl'")
  }

  if( !is.null(mainPI) ){
    if(!mainPI %in% names(attTbl)){
      stop("'mainPI' must be a column of 'attTbl'")
    }
  }

  if( !is.null(secPI) ){
    if(!secPI %in% names(attTbl)){
      stop("'secPI' must be a column of 'attTbl'")
    }
  }

  if(all( sapply(list(mainPI, secPI, min.N), is.null) )){
    stop("'mainPI', 'secPI' and 'min.N' are all NULL")
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

    ind <- which( attTbl[[mainPI]][noNA]<=cut.mPI | attTbl[[secPI]][noNA]<=cut.sPI )
    RO1[ match(attTbl[noNA,][["Cell"]][ind], attTbl$Cell) ] <- NA

    if(all(is.na(RO1))){
      stop("All raster object cells below PI values")
    }

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

    r_RO1 <- cv.2.rast(r = r, attTbl$Cell, classVector = RO1)
    terra::plot(r_RO1, type="classes", main="Raster objects - Output",
                plg=list(x=terra::ext(r_RO1)[2], y=terra::ext(r_RO1)[4], cex=0.9))

  }

  return(RO1)

}

#' Position index addition
#'
#' Add new raster objects based on position index values.
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
#' @param add.mPI numeric, threshold of main position index values. Cells with
#'   values above the threshold are flagged as cells potentially being part of
#'   new raster objects.
#' @param add.sPI numeric, threshold of secondary position index values. Cells
#'   with values above the threshold flagged as cells potentially being part of
#'   new raster objects.
#' @param cond.filter character string, defines what cells have to be considered
#'   by the function the arguments. Test cell absolute conditions can be used
#'   (see \code{\link{conditions}}).
#' @param min.N numeric, the minimum number of cells a raster object has to have
#'   to be included in the function output.
#' @param plot logic, plot the results.
#' @param r single or multi-layer raster of the class \code{SpatRaster} (see
#'   \code{help("rast", terra)}) used to compute the attribute table. Required
#'   only if \code{plot = TRUE}.
#'
#' @return The function returns a class vector with raster objects IDs. The
#'   vector has length equal to the number of rows of the attribute table. NA
#'   values are assigned to cells that do not belong to any raster object.
#'
#' @details New raster objects are added based on position index values. Two
#'   different position indices can be passed to the function (\code{mainPI} and
#'   \code{secPI}).
#'
#'   * Input raster objects are assigned to the same class to flag cells that
#'   are part of raster objects;
#'
#'   * Cells with values above \code{mainPI} **OR** above \code{mainPI} are
#'   flagged as cells potentially being part of new raster objects;
#'
#'   * If not connected to any of the existing raster objects, the groups of
#'   cells above position index values are assigned to new raster objects.
#'
#'   * Only raster objects with at least as many cells as specified by the
#'   argument \code{min.N} are included in the function output.
#'
#'   * If both \code{mainPI} and \code{secPI} are equal to \code{NULL}, the
#'   function will exclusively filter raster objects based on their size
#'   (\code{min.N}).
#'
#' @note Raster objects are added only if they do not share any border with
#'   input raster objects.
#'
#' @seealso [attTbl()], [ngbList()], [rel.pi()], [pi.sgm()], [conditions()]
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
#' at$RO[at$dummy_var==8] <- 1
#' at$RO <- cond.4.nofn(at, nbs, classVector = at$RO, class=1, nbs_of = 1,
#'                      cond = "dummy_var < dummy_var[] & dummy_var > 2")
#'
#' # One raster object
#' unique(at$RO)
#'
#' ################################################################################
#' # POSITION INDEX
#' ################################################################################
#' at$PI <- (at$dummy_var - mean(at$dummy_var))/stats::sd(at$dummy_var)
#'
#' ################################################################################
#' # POSITION INDEX ADDITION
#' ################################################################################
#' RO1 <- pi.add(at, nbs,
#'               RO = "RO",     # Raster objects
#'               mainPI = "PI", # PI addition layer
#'               add.mPI = 1,   # add disjoint objects with PI values > 1
#'               plot = FALSE, r = r)
#'
#' ################################################################################
#' # PLOT
#' ################################################################################
#' # Convert class vectors to raster
#' r_RO  <- cv.2.rast(r = r, classVector = at$RO)
#' r_RO1 <- cv.2.rast(r = r, classVector = RO1)
#'
#' # Plot
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(4.5, 0.5, 2, 3.2)
#'
#' terra::plot(r_RO, type="classes", main="Raster objects - Input", mar=m,
#'             plg=list(x=1, y=1, cex=0.9))
#'
#' terra::plot(r_RO1, type="classes", main="Raster objects - Output", mar=m,
#'             plg=list(x=1, y=1, cex=0.9))
#' text(xyFromCell(r,at$Cell), as.character(round(at$PI,2)),
#' cex = 0.8) # visualize relPI
#' text(0.01, 1, "Add on PI >= 1", adj=c(0,0), cex = 0.8)
#' par(oldpar)
#'
#' # Two raster object
#' unique(RO1)

pi.add <- function(attTbl,
                   ngbList,
                   rNumb = FALSE,
                   RO,
                   mainPI = NULL,
                   secPI = NULL,
                   add.mPI = NULL,
                   add.sPI = NULL,
                   cond.filter = NULL,
                   min.N = NULL,
                   plot = FALSE,
                   r = NULL){


  # TEST ARGUMENT RO, mainPI
  if( !RO %in% names(attTbl) ){
    stop("'RO', 'mainPI' must be columns of 'attTbl'")
  }

  if( !is.null(mainPI) ){
    if(!mainPI %in% names(attTbl)){
      stop("'mainPI' must be a column of 'attTbl'")
    }
  }

  if( !is.null(secPI) ){
    if(!secPI %in% names(attTbl)){
      stop("'secPI' must be a column of 'attTbl'")
    }
  }

  if(all( sapply(list(mainPI, secPI, min.N), is.null) )){
    stop("'mainPI', 'secPI' and 'min.N' are all NULL")
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

  # ADD SEGMENTS ###################################################################
  RO1 <- attTbl[[RO]]

  if(!is.null(add.mPI) | !is.null(add.sPI)){

    if(is.null(add.mPI)){add.mPI <- max(attTbl[[mainPI]], na.rm = TRUE)+0.1}
    if(is.null(add.sPI)){

      if(is.null(secPI)){secPI <- mainPI}
      add.sPI <- max(attTbl[[secPI]], na.rm = TRUE)+0.1

    }

    RO1[!is.na(RO1)] <- 1
    RO1[(attTbl[[mainPI]]>=add.mPI | attTbl[[secPI]]>=add.sPI ) & is.na(RO1)] <- 0

    # REMOVE SEGMENT CONTINUOUS WITH EXISTING OBJECTS
    RO1 <- reclass.nbs(attTbl, ngbList, rNumb = rNumb, classVector = RO1,
                       nbs_of = 1, class = 0, reclass = NA, reclass_all = TRUE)

  }

  # APPLY COND.FILTER ##############################################################
  if(!is.null(cond.filter)){
    cf <- paste0("!is.na(RO1)&", cond.filter)
  } else{
    cf <- "!is.na(RO1)"
  }

  # RASTER OBJECTS #################################################################
  RO1 <- anchor.seed(cbind(attTbl, RO1), ngbList,
                     rNumb = TRUE, class = NULL, silent = TRUE,
                     cond.filter = cf,
                     cond.seed   = cf,
                     cond.growth = cf)

  # REMOVE SMALL RASTER OBJECTS
  if(!is.null(min.N)){

    l <- lengths(split(RO1, RO1)) >= min.N
    l <- as.numeric(names(l)[l])

    RO1[!(RO1 %in% l)] <- NA
  }

  if(plot){

    r_RO1 <- cv.2.rast(r = r, attTbl$Cell, classVector = RO1)
    terra::plot(r_RO1, type="classes", main="Raster objects - Output",
                plg=list(x=terra::ext(r_RO1)[2], y=terra::ext(r_RO1)[4], cex=0.9))
  }

  return(RO1)

}
