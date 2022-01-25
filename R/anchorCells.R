#' Anchor cells from spatial vector objects
#'
#' Returns a vector of raster cell numbers extracted at the locations of a
#' spatial object.
#'
#' @param rstack \code{Raster*} object.
#' @param dsn data source name (filename) or an `sf` or a `Spatial` object.
#' @param only_NAs logic, cell numbers extracted only for incomplete cases at
#'   the locations of a spatial object. Incomplete cases are cells having an
#'   NA-value in one or more layers of the \code{Raster*} object.
#' @param fill_NAs logic, cell numbers extracted at the locations of a spatial
#'   object _and_ at contiguous locations that are incomplete cases.
#' @param plot logic, plot anchor points.
#' @param saveRDS filename, if a file name is provided save the anchor cell
#'   vector as an RDS file.
#' @param writeRaster filename, if a raster name is provided save the anchor
#'   cell vector as a raster file.
#' @param overWrite logic, if RDS and raster names already exist, existing files
#'   are overwritten.
#'
#' @return Numeric vector of raster cell numbers.
#'
#' @details When the arguments \code{only_NA} and \code{fill_NAs} are FALSE the
#'   numeric output is equivalent to the output of the function
#'   \code{raster::extract()} over a raster whose values are cell numbers.
#'
#' @export
#' @examples
#'
#' # NOT RUN
#' \dontrun{
#' # LOAD LIBRARIES AND DATA
#' library(raster)
#' library(scapesClassification)
#'
#' # CELL NUMBERS OF A DUMMY RASTER (7X7)
#' r_cn <- raster(matrix(1:49, nrow = 7, byrow = TRUE))
#'
#' # SET SOME NA-VALUE
#' r_cn[c(9, 10, 11, 17, 18)] <- NA
#'
#' # BULD A DUMMY POLYGON
#' pol <- rbind(c(0,0.95), c(0.28,1), c(0.24, 0.72), c(0.05,0.72), c(0,0.95))
#' pol <- spPolygons(pol)
#'
#' # EMPTY RASTER TO PLOT anchor.svo RESULTS
#' r2 <- r_cn; r2[] <- NA
#'
#' # SET PLOT LAYOUT
#' par(mfrow=c(2,2), mar=c(1, 0, 2, 0))
#'
#' # PLOT 1
#' ################################################################################
#' # only_NAs = FALSE; fill_NAs = FALSE
#' ac     <- anchor.svo(r_cn, pol, only_NAs = FALSE, fill_NAs = FALSE)
#' r2[ac] <- 1
#'
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE,
#'      main = "only_NAs = FALSE; fill_NAs = FALSE", xlim = c(0,1))
#' text(r_cn); plot(pol, add = TRUE, lwd = 2.5, border = "red")
#' ################################################################################
#'
#' # PLOT 2
#' ################################################################################
#' # only_NAs = TRUE; fill_NAs = FALSE
#' ac     <- anchor.svo(r_cn, pol, only_NAs = TRUE, fill_NAs = FALSE)
#' r2[]   <- NA; r2[ac] <- 1
#'
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE,
#'      main = "only_NAs = TRUE; fill_NAs = FALSE")
#' text(r_cn); plot(pol, add = TRUE, lwd = 2.5, border = "red")
#' ################################################################################
#'
#' # PLOT 3
#' ################################################################################
#' # only_NAs = FALSE; fill_NAs = TRUE
#' ac     <- anchor.svo(r_cn, pol, only_NAs = FALSE, fill_NAs = TRUE)
#' r2[]   <- NA; r2[ac] <- 1
#'
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE,
#'      main = "only_NAs = FALSE; fill_NAs = TRUE")
#' text(r_cn); plot(pol, add = TRUE, lwd = 2.5, border = "red")
#' ################################################################################
#'
#' # PLOT 4
#' ################################################################################
#' # only_NAs = TRUE; fill_NAs = TRUE
#' ac     <- anchor.svo(r_cn, pol, only_NAs = TRUE, fill_NAs = TRUE)
#' r2[]   <- NA; r2[ac] <- 1
#'
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE,
#'      main = "only_NAs = TRUE; fill_NAs = TRUE")
#' text(r_cn); plot(pol, add = TRUE, lwd = 2.5, border = "red")
#' }

anchor.svo <- function(rstack,
                       dsn,
                       only_NAs = FALSE,
                       fill_NAs = FALSE,
                       plot = FALSE,
                       saveRDS = NULL,
                       writeRaster = NULL,
                       overWrite = FALSE) {
  if (!overWrite) {
    # RDS
    if (!is.null(saveRDS)) {
      if (file.exists(saveRDS))
        stop("RDS filename exists; use a different name")
    }
    # Raster
    if (!is.null(writeRaster)) {
      if (file.exists(writeRaster))
        stop("raster filename exists; use a different name")
    }

  }

  if(is.character(dsn)){

    p <- rgdal::readOGR(dsn = dsn, verbose = FALSE)
    p <- sp::spTransform(p, raster::crs(rstack))

  } else if(methods::is(dsn, "sf")|methods::is(dsn, "Spatial")){

    p <- dsn

  } else {

    stop("spatial_vector_name must be a data source name OR an 'sf' object OR a 'Spatial' object")

  }

  r2     <-
    raster::raster(
      ext = raster::extent(rstack),
      crs = raster::crs(rstack),
      res = raster::res(rstack)
    )
  r2[]   <- seq(raster::ncell(rstack))
  i_cell <- unique(unlist(raster::extract(r2, p)))

  if (is.null(i_cell))
    stop("no overlap between raster and shape files")

  # only_NA & fillNAs arguments
  v           <- as.data.frame(raster::values(rstack))
  v[["Cell"]] <- seq(raster::ncell(rstack))
  c_cc        <- v[stats::complete.cases(v), "Cell"]

  if (only_NAs) {
    i_cell <- i_cell[!(i_cell %in% c_cc)]
  }

  if (fill_NAs) {
    nbs <- nbg8(raster::nrow(rstack) , raster::ncol(rstack))

    continue <- T
    i_cell0  <- i_cell

    while (continue) {
      continue <- F

      na_nbs <- unlist(nbs[i_cell])
      i_cell <- unique(na_nbs[!(na_nbs %in% c_cc)])

      i_cell <- setdiff(i_cell, i_cell0)

      if (length(i_cell) > 0) {
        continue <- T
        i_cell0  <- unique(c(i_cell, i_cell0))

      }

    }

    i_cell <- i_cell0

  }

  # outputs
  r2[]       <- NA
  r2[i_cell] <- 1

  if (plot)
    raster::plot(r2)
  if (!is.null(saveRDS))
    saveRDS(i_cell, saveRDS)
  if (!is.null(writeRaster))
    raster::writeRaster(r2, writeRaster, overwrite = overWrite)

  return(i_cell)

}


#' Cell numbers to class vector
#'
#' Converts a vector of cell numbers into a class vector.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param rstack the \code{Raster*} object used to compute the
#'   \code{\link{attTbl}}.
#' @param anchor integer vector of raster cell numbers.
#' @param class numeric, the classification number to attribute to all cells
#'   that meet the function conditions.
#' @param classVector numeric vector, if provided, it defines the cells in the
#'   attribute table that have already been classified and that have to be
#'   ignored by the function (unless the argument \code{overwrite_class =
#'   TRUE}).
#' @param class2cell logic, attribute the classification number to the cells of
#'   the argument \code{anchor}. If there is a \code{classVector} input, the
#'   classification number is only assigned to \code{classVector} NA-cells.
#' @param class2nbs logic, attribute the classification number to cells adjacent
#'   to the ones of the argument \code{anchor}. If there is a \code{classVector}
#'   input, the classification number is only assigned to \code{classVector}
#'   NA-cells.
#' @param overwrite_class logic, if there is a \code{classVector} input,
#'   reclassify cells that were already classified and that meet the function
#'   conditions.
#' @param plot logic, plot the class vector output.
#' @param writeRaster filename, if a raster name is provided, save the class
#'   vector in a raster file.
#' @param overWrite logic, if the raster names already exist, the existing file
#'   is overwritten.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. If there is no \code{classVector} input, the function return
#'   a new class vector.
#'
#' @details Converts a vector of cell numbers into a class vector. If there is a
#'   \code{classVector} input, then the class vector is updated assigning a
#'   classification number to all cells that meet the function conditions.
#'
#' @seealso [conditions()], [anchor.svo()], [attTbl()]
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # LOAD LIBRARIES AND DATA
#' library(raster)
#' library(scapesClassification)
#'
#' # CELL NUMBERS OF A DUMMY RASTER (7X7)
#' r_cn <- raster(matrix(1:49, nrow = 7, byrow = TRUE))
#'
#' # COMPUTE ATTRIBUTE TABLE AND LIST OF NEIGHBORHOODS
#' at  <- attTbl(r_cn, "dummy_var")
#' nbs <- ngbList(r_cn)
#'
#' # SET PLOT LAYOUT
#' par(mfrow=c(2,2), mar=c(1, 0, 2, 0))
#'
#' # PLOT 1
#' ################################################################################
#' # class2cell   = TRUE & class2nbs = FALSE
#' # anchor cells = 1:7
#' cv <- anchor.cell(attTbl = at,
#'                   rstack = r_cn,
#'                   anchor = 1:7,
#'                   class  = 10,
#'                   class2cell = TRUE,
#'                   class2nbs  = FALSE)
#'
#' r_cv <- cv.2.rast(r = r_cn, index = at$Cell, classVector = cv)
#'
#' plot(r_cv, axes=FALSE, box=FALSE, legend = FALSE,
#'      colNA="#818792", col="#78b2c4", main = "class2cell=TRUE; class2nbs=FALSE")
#' text(r_cn)
#' ################################################################################
#'
#' # PLOT 2
#' ################################################################################
#' # class2cell   = FALSE & class2nbs = TRUE
#' # anchor cells = 1:7
#' cv <- anchor.cell(attTbl = at,
#'                   rstack = r_cn,
#'                   anchor = 1:7,
#'                   class  = 10,
#'                   class2cell = FALSE,
#'                   class2nbs  = TRUE)
#'
#' r_cv <- cv.2.rast(r = r_cn, index = at$Cell, classVector = cv)
#'
#' plot(r_cv, axes=FALSE, box=FALSE, legend = FALSE,
#'      colNA="#818792", col="#78b2c4", main = "class2cell=FALSE; class2nbs=TRUE")
#' text(r_cn)
#' ################################################################################
#'
#' # PLOT 3
#' ################################################################################
#' # class2cell   = TRUE & class2nbs = TRUE
#' # anchor cells = 1:7
#' cv <- anchor.cell(attTbl = at,
#'                   rstack = r_cn,
#'                   anchor = 1:7,
#'                   class  = 10,
#'                   class2cell = TRUE,
#'                   class2nbs  = TRUE)
#'
#' r_cv <- cv.2.rast(r = r_cn, index = at$Cell, classVector = cv)
#'
#' # PLOT RASTER AND CELL NUMBERS
#' plot(r_cv, axes=FALSE, box=FALSE, legend = FALSE,
#'      colNA="#818792", col="#78b2c4", main = "class2cell=TRUE; class2nbs=TRUE")
#' text(r_cn)
#' }

anchor.cell <-
  function(attTbl,
           rstack,
           anchor,
           class,
           classVector = NULL,
           class2cell = TRUE,
           class2nbs = TRUE,
           overwrite_class = FALSE,
           plot = FALSE,
           writeRaster = NULL,
           overWrite = FALSE) {
    ###
    if (!("Cell" %in% names(attTbl))){
      stop(
        "the attribute table should have one column named 'Cell' with cell numbers
                                       that indicate the position of each row in the original Raster object"
      )
    }

    if (!overWrite) {
      # Raster
      if (!is.null(writeRaster)) {
        if (file.exists(writeRaster))
          stop("raster filename exists; use a different name")
      }
    }


    if (is.null(classVector)) {
      classVector <- rep(as.integer(NA), NROW(attTbl))
    }


    ###
    if (class2cell) {
      if (!overwrite_class) {
        classVector[which(attTbl$Cell %in% anchor &
                            is.na(classVector))] <- class

      } else {
        classVector[which(attTbl$Cell %in% anchor)] <- class

      }

    }

    if (class2nbs) {
      nbs_anchor <-
        nbg8(n_row = raster::nrow(rstack), raster::ncol(rstack))
      nbs_anchor <- unlist(nbs_anchor[anchor])

      if (!overwrite_class) {
        classVector[which(attTbl$Cell %in% nbs_anchor &
                            is.na(classVector))] <- class

      } else {
        classVector[which(attTbl$Cell %in% nbs_anchor)] <- class

      }

      if(!class2cell){
        classVector[which(attTbl$Cell %in% anchor)] <- NA
      }

    }

    r2   <- rstack[[1]]
    r2[] <- NA
    r2[attTbl$Cell[!is.na(classVector)]] <-
      classVector[!is.na(classVector)]

    if (plot)
      raster::plot(r2)
    if (!is.null(writeRaster))
      raster::writeRaster(r2, writeRaster, overwrite = overWrite)

    return(classVector)

  }



