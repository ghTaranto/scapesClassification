#' Attribute table
#'
#' Converts a single or a multi-layer raster into an attribute table
#' (\code{data.frame}).
#'
#' @param r single or multi-layer raster of the class \code{SpatRaster} (see
#'   \code{help("rast", terra)}).
#' @param var_names character vector, raster layers' names in the attribute
#'   table. If \code{NULL}, then the original layers' names are used.
#'
#' @encoding UTF-8
#'
#' @return data.frame
#'
#' @details Attribute tables come with a column named **\code{"Cell"}** which
#'   stores raster cell numbers and associate each row of the attribute table
#'   with a cell of the raster object. The remaining columns of the attribute
#'   table store the data contained in the raster layers. Note that only raster
#'   cells having no missing value in no layer (**complete cases**) are included
#'   in the attribute table.
#'
#' @note **Attribute table contains only complete cases**, i.e., raster cells
#'   having a value for every layer of the stack.
#'
#' @export
#' @examples
#' library(scapesClassification)
#' library(terra)
#'
#' ## CREATE A DUMMY RASTER ##
#' r <- terra::rast(matrix(c(NA,100,100,NA,100,100,0,0,0),
#'                         nrow = 3,
#'                         ncol = 3,
#'                         byrow = TRUE))
#'
#' ## RASTER CELL NUMBERS ##
#' rcn <- r; rcn[] <- 1:9
#'
#' ## PLOT DATA AND CELL NUMBERS ##
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(4, 1, 4, 1)
#'
#' plot(r, col="grey90", colNA="red3", mar=m, asp = NA, axes=FALSE, legend=FALSE)
#' text(r)
#' lines(r)
#' mtext(side=3, line=0.2, adj=0, cex=1.5, font=2, "Dummy_var")
#' legend("bottomright", ncol=1, bg="white", fill=c("red3"),
#'        legend = c("NA cells (1 and 4)"))
#'
#' plot(rcn, col="grey90", mar=m, asp=NA, axes=FALSE, legend=FALSE)
#' text(rcn)
#' lines(rcn)
#' mtext(side=3, line=0.2, adj=0, cex=1.5, font=2, "Cell numbers")
#' par(oldpar)
#'
#' ## VISUALIZE ATTRIBUTE TABLE ##
#'
#' at <- attTbl(r, var_names = c("dummy_var"))
#' at
#'
#' # Note that cells 1 and 4 have missing values and therefore are not included in the table
#' any(at$Cell %in% c(1,4))

attTbl <- function(r, var_names = NULL){

  dt <- data.frame(Cell=1:terra::ncell(r), terra::values(r))
  dt <- dt[stats::complete.cases(dt),]
  row.names(dt) <- NULL

  if( !is.null(var_names) ) {

    if( length(var_names) != length(names(dt)) -1 ){

      stop("var_names length should be equal to the number of layers in r")

    } else {

      names(dt)[names(dt) != "Cell"] <- var_names

    }

  }

  return(dt)

}


#' List of neighborhoods
#'
#' Computes the neighborhoods of the cells of a raster. Neighborhoods are not
#' computed for cells with missing values.
#'
#' @param r single or multi-layer raster of the class \code{SpatRaster} (see
#'   \code{help("rast", terra)}).
#' @param rNumb logic, the neighbors of a raster cell are identified by **cell
#'   numbers (\code{rNumb=FALSE})** or by **row numbers (\code{rNumb=TRUE})**.
#'   If true, the argument \code{attTbl} cannot be NULL.
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}} (see \code{\link{attTbl}}). It is required only if the
#'   argument \code{rNumb=TRUE}.
#'
#' @encoding UTF-8
#'
#' @return Named list of integer vectors.
#'
#' @details **Definition of neighborhood**
#'
#'   * A cell with coordinates \code{(x, y)} has 8 neighbors with coordinates:
#'   \code{(x±1, y)},  \code{(x, y±1)} and \code{(x±1, y±1)}. Cells on the edge
#'   of a raster have less than 8 neighbors.
#'
#'   **Neighborhoods (\code{rNumb=FALSE})**
#'
#'   * Neighbors are identified by their cell numbers if the argument
#'   \code{rNumb=FALSE}.
#'
#'   **Neighborhoods (\code{rNumb=TRUE})**
#'
#'   * Neighbors are identified by their positions in the attribute table (i.e.
#'   row numbers) if the argument \code{rNumb=TRUE};
#'
#'   * When the argument \code{rNumb = TRUE}, neighbors with missing values are
#'   omitted;
#'
#'   * \code{(scapes)Classifications} are faster when the list of neighborhoods
#'   uses row numbers.
#'
#'   **Neighborhood names**
#'
#'   The list of neighborhoods is named.
#'
#'   * When \code{rNumb = FALSE}, the element name identifies the raster cell to
#'   which the neighborhood refers. For instance, the element with name
#'   \code{"n"} stores the neighborhood of the raster cell \code{n}.
#'
#'   * When \code{rNumb = TRUE}, the element name identifies the row number to
#'   which the neighborhood refers. For instance, the element with name
#'   \code{"n"} stores the neighborhood of the raster cell located in the
#'   \code{nth} row of the attribute table (\code{attTbl$Cell[n]}).
#'
#' @seealso [ngb8()], [attTbl()]
#'
#' @note
#'   * There is always a correspondence between the indices of the attribute
#'   table (\code{\link{attTbl}}) and the indices of the list of neighborhoods:
#'   the 1st element of the list corresponds to the neighbors of the cell stored
#'   in the 1st row of the attribute table; the 2nd element corresponds to the
#'   2nd row; etc.
#'
#'   * There is a correspondence between the raster cell number and the indices
#'   of the list of neighborhoods only when no missing value is present in the
#'   raster.
#'
#' @export
#' @examples
#' library(scapesClassification)
#' library(terra)
#'
#' ## CREATE A DUMMY RASTER AND COMPUTE ATTRIBUTE TABLE ##
#' r <- terra::rast(matrix(c(NA,100,100,NA,100,100,0,0,0),
#'                         nrow = 3,
#'                         ncol = 3,
#'                         byrow = TRUE))
#'
#' at <- attTbl(r, var_names = c("dummy_var"))
#'
#' ## RASTER CELL NUMBERS ##
#' rcn <- r; rcn[] <- 1:9
#'
#' ## PLOT DATA AND CELL NUMBERS ##
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(4, 1, 4, 1)
#'
#' plot(r, col="grey90", colNA="red3", mar=m, asp=NA, axes=FALSE, legend=FALSE)
#' text(r)
#' lines(r)
#' mtext(side=3, line=0.2, adj=0, cex=1.5, font=2, "Dummy_var")
#' legend("bottomright", ncol = 1, bg = "white", fill = c("red3"),
#'        legend = c("NA cells (1 and 4)"))
#'
#' plot(rcn, col="grey90", mar=m, asp = NA, axes=FALSE, legend=FALSE)
#' text(rcn)
#' lines(rcn)
#' mtext(side=3, line=0.2, adj=0, cex=1.5, font=2, "Cell numbers")
#' par(oldpar)
#'
#' ## NEIGHBORHOODS - CELL NUMBERS ##
#'
#' # Cells 1 and 4 are omitted because they are NAs
#' nbs_CELL <- ngbList(r, rNumb = FALSE)
#' nbs_CELL
#'
#'
#' ## NEIGHBORHOODS - ROW NUMBERS ##
#'
#' # Cells 1 and 4 are omitted because they are NAs
#' nbs_ROW <- ngbList(r, rNumb = TRUE, attTbl = at)
#' nbs_ROW
#'
#' # Numbers in 'nbs_ROW' refer to row numbers
#' # (e.g. number 1 refers to the cell #2)
#' at$Cell[1]
#'
#' # (e.g. number 2 refers to the cell #3)
#' at$Cell[2]
#'
#' # (e.g. number 5 refers to the cell #7)
#' at$Cell[5]
#'
#'
#' ## CONSIDER THE NEIGHBORHOOD OF CELL #2 ##
#'
#' # Cell #2 corresponds to the 1st element of both 'nbs_CELL' and 'nbs_ROW'
#' # because raster cell 1 is an NA-cell
#' r[1]
#'
#' # Neighborhood cell #2 corresponds to cells:
#' nbs_CELL[1]
#'
#' # Neighborhood cell #2 corresponds to rows:
#' nbs_ROW[1]
#'
#' # Rows can be converted to cell numbers
#' at$Cell[ nbs_ROW[[1]] ]
#'
#' # Note that 'at$Cell[ nbs_ROW[[1]] ]' is not equal to 'nbs_CELL'
#' identical( at$Cell[ nbs_ROW[[1]] ] , nbs_CELL[[1]] )
#'
#' # This is because raster cells 1 and 4 (NA-cells) are omitted in 'nbs_ROW'
#' setdiff(nbs_CELL[[1]], at$Cell[ nbs_ROW[[1]] ])
#' r[c(1,4)]


ngbList <- function(r, rNumb = FALSE, attTbl = NULL){

  ind <- which(stats::complete.cases( as.data.frame(terra::values(r)) ))
  nbs <- ngb8(terra::nrow(r[[1]]), terra::ncol(r[[1]]))

  nbs <- nbs[ind]

  # CONVERT nbs FORM CELL IDS TO CELL INDECES
  if(rNumb){
    if(is.null(attTbl)){stop("When the argument rNumb = TRUE, an attribute table must be provided")}
    fct    <- rep(seq_along(lengths(nbs)), lengths(nbs))
    nbs    <- match(unlist(nbs), attTbl$Cell)
    no_nas <- !is.na(nbs)
    nbs    <- nbs[no_nas]
    fct    <- fct[no_nas]

    nbs    <- split(nbs, fct)
  }

  return(nbs)

}
