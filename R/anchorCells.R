#' Set Anchor Cells from Spatial Vector Objects
#'
#' Returns a vector of cell numbers overlapping with the locations of spatial
#' vector data. These cells can be used as anchor cells in other
#' \code{scapesClassification} functions.
#'
#' @param rstack \code{RasterStack} or \code{RasterLayer} object.
#' @param spatial_vector_name OGR data source name.
#' @param only_NAs logic, return only cell numbers overlapping with the spatial
#'   vector data that have missing values in rstack.
#' @param fill_NAs logic, return cell numbers of cells adjacent to those
#'   overlapping with the spatial vector data that have missing values in
#'   rstack.
#' @param plot logic, plot anchor points.
#' @param saveRDS filename, if a file name is provided save the anchor cell
#'   vector as an RDS file.
#' @param writeRaster filename, if a raster name is provided save the anchor
#'   cell vector as a raster file.
#' @param overWrite logic, if RDS and raster names already exist, existing files
#'   are overwritten.
#'
#' @return anchor cell vector.
#'
#' @details When the arguments \code{only_NA} and \code{fillNAs} are FALSE the
#'   output is equivalent to the output of the function \code{raster::extract()}
#'   over a raster whose values are cell numbers. Cell numbers start with 1 in
#'   the upper-left corner and increase from left to right and from top to
#'   bottom.
#'
#' @export
#'


anchor.svo <- function(rstack,
                       spatial_vector_name,
                       only_NAs = FALSE,
                       fill_NAs = FALSE,
                       plot = TRUE,
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

  p <- rgdal::readOGR(dsn = spatial_vector_name, verbose = F)
  p <- sp::spTransform(p, raster::crs(rstack))

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


#' Classify Raster Cells Based on Anchor Cells
#'
#' Return a classification vector based on anchor cells.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param rstack \code{RasterStack} or \code{RasterLayer} object used to compute
#'   the attribute table (see \code{\link{attTbl}}).
#' @param anchor integer vector, contains the cell numbers to be considered as
#'   anchor points.
#' @param class numeric, the class to attribute to cells meeting contiguity
#'   conditions.
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified.
#' @param class2cell logic, attribute the new classification, defined by the
#'   argument \code{class}, to \code{anchor} cells.
#' @param class2nbs logic, attribute the new classification, defined by the
#'   argument \code{class}, to the cells adjacent to \code{anchor} cells.
#' @param overwrite_class logic, reclassify cells that have already been
#'   classified.
#' @param plot logic, plot \code{classVector}.
#' @param writeRaster filename, if a raster name is provided save
#'   \code{classVector} in a raster file.
#' @param overWrite logic, if the raster names already exist, the existing file
#'   is overwritten.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. If no \code{classVector} was provided, the function return a
#'   new classification vector.
#'
#' @details Classify cells of the \code{rstack} argument based on the anchor
#'   cells provided by \code{anchor}. Return a classification vector that can be
#'   indexed into \code{rstack} using the \code{Cell} column of the attribute
#'   table.
#'
#' @seealso [anchor.seed()], [anchor.svo()], [attTbl()]
#'
#' @export

anchor.cell <-
  function(attTbl,
           rstack,
           anchor,
           class,
           classVector = NULL,
           class2cell = TRUE,
           class2nbs = TRUE,
           overwrite_class = FALSE,
           plot = TRUE,
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



