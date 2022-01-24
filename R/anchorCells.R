#' Anchor cells from spatial vector objects
#'
#' Returns a vector of raster cell numbers extracted at the locations of a
#' spatial vector object.
#'
#' @param rstack \code{Raster*} object.
#' @param spatial_vector_name data source name (filename), sf or Spatial object.
#' @param only_NAs logic, cell numbers extracted only at the locations of a
#'   spatial vector object that are not complete cases (i.e. have some missing
#'   values in some of the layers of the \code{Raster*} object).
#' @param fill_NAs logic, cell numbers extracted also at locations contiguous to
#'   those of the spatial vector object that are not complete cases (i.e. have
#'   some missing values in some of the layers of the \code{Raster*} object).
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
#'   output is equivalent to the output of the function \code{raster::extract()}
#'   over a raster whose values are cell numbers.
#'
#' @export
#'

anchor.svo <- function(rstack,
                       spatial_vector_name,
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

  if(is.character(spatial_vector_name)){

    p <- rgdal::readOGR(dsn = spatial_vector_name, verbose = FALSE)
    p <- sp::spTransform(p, raster::crs(rstack))

  } else if(methods::is(spatial_vector_name, "sf")|methods::is(spatial_vector_name, "Spatial")){

    p <- spatial_vector_name

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


#' Class vector from cell numbers
#'
#' Converts a vector of cell numbers into a class vector.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}} (see \code{\link{attTbl}}).
#' @param rstack The \code{Raster*} object used to compute the attribute table
#'   (see \code{\link{attTbl}}).
#' @param anchor integer vector, vector of raster cell numbers.
#' @param class numeric, the classification number to attribute to cells
#'   indicated by the argument \code{anchor}.
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified. If provided, the function do not
#'   attribute a class to cells that were already classified unless the argument
#'   \code{overwrite_class = TRUE}.
#' @param class2cell logic, attribute the classification number to the cells
#'   indicated by argument \code{anchor}.
#' @param class2nbs logic, attribute the classification number to the cells
#'   adjacent to the ones indicated by argument \code{anchor}.
#' @param overwrite_class logic, reclassify cells that had already been
#'   classified.
#' @param plot logic, plot \code{classVector}.
#' @param writeRaster filename, if a raster name is provided save the class
#'   vector in a raster file.
#' @param overWrite logic, if the raster names already exist, the existing file
#'   is overwritten.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. If no \code{classVector} was provided, the function return a
#'   new class vector.
#'
#' @details Converts a vector of cell numbers into a class vector. If a class
#'   vector is provided as an input (argument \code{classVector}), then this
#'   class vector is updated assigning a classification number to the cell
#'   indicated by the argument \code{anchor}. The classification vector can be
#'   indexed into the \code{rstack} using the \code{Cell} column of the
#'   attribute table (see \code{\link{attTbl}}).
#'
#' @seealso [anchor.seed()], [anchor.svo()], [attTbl()]
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # LOAD LIBRARIES AND DATA
#' library(raster)
#' library(scapesClassification)
#'
#' r <- list.files(system.file("extdata", package = "scapesClassification"),
#'                 pattern = "dummy_raster\\.tif", full.names = TRUE)
#' r <- raster(r)
#'
#' # SET RASTER DATA EQUAL TO CELL NUMBERS
#' r[]<- 1:49
#'
#' # COMPUTE ATTRIBUTE TABLE AND LIST OF NEIGHBORHOODS
#' at  <- attTbl(r, "dummy_var")
#' nbs <- ngbList(r)
#'
#' # CELLS FROM 1:7 IN THE `anchor` ARGUMENT
#' anch <-  1:7
#'
#'
#' # class2cell = TRUE & class2nbs = FALSE
#' cv <- anchor.cell(attTbl = at,
#'                   rstack = r,
#'                   anchor = anch,
#'                   class  = 10,
#'                   class2cell = TRUE,
#'                   class2nbs  = FALSE)
#'
#' r_cv <- cv.2.rast(r = r, index = at$Cell, classVector = cv)
#'
#' # PLOT RASTER AND CELL NUMBERS
#' plot(r_cv, axes=FALSE, legend = FALSE,
#'      colNA="#818792", col="#78b2c4", main = "class2cell=TRUE & class2nbs=FALSE")
#' text(r)
#'
#'
#' # class2cell = FALSE & class2nbs = TRUE
#' cv <- anchor.cell(attTbl = at,
#'                   rstack = r,
#'                   anchor = anch,
#'                   class  = 10,
#'                   class2cell = FALSE,
#'                   class2nbs  = TRUE)
#'
#' r_cv <- cv.2.rast(r = r, index = at$Cell, classVector = cv)
#'
#' # PLOT RASTER AND CELL NUMBERS
#' plot(r_cv, axes=FALSE, box=FALSE, legend = FALSE,
#'      colNA="#818792", col="#78b2c4", main = "class2cell=FALSE & class2nbs=TRUE")
#' text(r)
#'
#'
#' # class2cell = TRUE & class2nbs = TRUE
#' cv <- anchor.cell(attTbl = at,
#'                   rstack = r,
#'                   anchor = anch,
#'                   class  = 10,
#'                   class2cell = TRUE,
#'                   class2nbs  = TRUE)
#'
#' r_cv <- cv.2.rast(r = r, index = at$Cell, classVector = cv)
#'
#' # PLOT RASTER AND CELL NUMBERS
#' plot(r_cv, axes=FALSE, box=FALSE, legend = FALSE,
#'      colNA="#818792", col="#78b2c4", main = "class2cell=TRUE & class2nbs=TRUE")
#' text(r)
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



