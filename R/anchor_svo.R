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
#' par(mfrow=c(2,2), mar=c(3, 2, 2, 2))
#'
#' # PLOT 1, only_NAs = FALSE; fill_NAs = FALSE
#' ################################################################################
#' # ANCHOR.SVO
#' ac     <- anchor.svo(r_cn, pol, only_NAs = FALSE, fill_NAs = FALSE)
#'
#' # ADD ANCHOR CELL TO RASTER FOR PLOTTING
#' r2[ac] <- 1
#'
#' # PLOT
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE, asp = NA)
#' plot(pol, add = TRUE, lwd = 2.5, border = "red")
#'
#' # REFERENCE PLOT
#' text(r_cn)
#' title(adj = 0.0, line = 1,
#'       sub = paste0("only_NAs = FALSE; fill_NAs = FALSE\nac=",
#'       paste0(sort(ac), collapse = ",")))
#' legend("bottomleft", ncol = 1, bg = "white",
#'        legend = c("Anchor cell", "Polygon"), fill = c("#78b2c4", "red"))
#' ################################################################################
#'
#' # PLOT 2, only_NAs = TRUE; fill_NAs = FALSE
#' ################################################################################
#' # ANCHOR.SVO
#' ac     <- anchor.svo(r_cn, pol, only_NAs = TRUE, fill_NAs = FALSE)
#'
#' # ADD ANCHOR CELL TO RASTER FOR PLOTTING
#' r2[]   <- NA
#' r2[ac] <- 1
#'
#' # PLOT
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE, asp = NA)
#' plot(pol, add = TRUE, lwd = 2.5, border = "red")
#'
#' # REFERENCE PLOT
#' text(r_cn)
#' title(adj = 0.0, line = 1,
#'       sub = paste0("only_NAs = TRUE; fill_NAs = FALSE\nac=",
#'       paste0(sort(ac), collapse = ",")))
#' legend("bottomleft", ncol = 1, bg = "white",
#'        legend = c("Anchor cell", "Polygon"), fill = c("#78b2c4", "red"))
#' ################################################################################
#'
#' # PLOT 3, only_NAs = FALSE; fill_NAs = TRUE
#' ################################################################################
#' # ANCHOR.SVO
#' ac     <- anchor.svo(r_cn, pol, only_NAs = FALSE, fill_NAs = TRUE)
#'
#' # ADD ANCHOR CELL TO RASTER FOR PLOTTING
#' r2[]   <- NA
#' r2[ac] <- 1
#'
#' # PLOT
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE, asp = NA)
#' plot(pol, add = TRUE, lwd = 2.5, border = "red")
#'
#' # REFERENCE PLOT
#' text(r_cn)
#' title(adj = 0.0, line = 1,
#'       sub = paste0("only_NAs = FALSE; fill_NAs = TRUE\nac=",
#'       paste0(sort(ac), collapse = ",")))
#' legend("bottomleft", ncol = 1, bg = "white",
#'        legend = c("Anchor cell", "Polygon"), fill = c("#78b2c4", "red"))
#' ################################################################################
#'
#' # PLOT 4, only_NAs = TRUE; fill_NAs = TRUE
#' ################################################################################
#' # ANCHOR.SVO
#' ac     <- anchor.svo(r_cn, pol, only_NAs = TRUE, fill_NAs = TRUE)
#'
#' # ADD ANCHOR CELL TO RASTER FOR PLOTTING
#' r2[]   <- NA
#' r2[ac] <- 1
#'
#' # PLOT
#' plot(r2, col="#78b2c4", colNA= "grey", axes=FALSE, box=FALSE, legend = FALSE, asp = NA)
#' plot(pol, add = TRUE, lwd = 2.5, border = "red")
#'
#' # REFERENCE PLOT
#' text(r_cn)
#' title(adj = 0.0, line = 1,
#'       sub = paste0("only_NAs = TRUE; fill_NAs = TRUE\nac=",
#'       paste0(sort(ac), collapse = ",")))
#' legend("bottomleft", ncol = 1, bg = "white",
#'        legend = c("Anchor cell", "Polygon"), fill = c("#78b2c4", "red"))

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