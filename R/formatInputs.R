#' Attribute table
#'
#' Converts a \code{Raster*} object into an attribute table (\code{data.frame}).
#'
#' @param rstack \code{Raster*} object.
#' @param var_names character vector, names of the \code{Raster*} object layers
#'   in the attribute table. If \code{NULL} layers' names are used.
#'
#' @return data.frame
#'
#' @details Attribute tables come with a column named **\code{"Cell"}** which
#'   stores raster cell numbers and associate each row of the attribute table
#'   with a cell of the raster object. Each of the remaining columns stores the
#'   values of a layer of the \code{Raster*} object. Note that only raster cells
#'   having no missing value in no layer (**complete cases**) are included in
#'   the attribute table.
#'
#' @note **Attribute table contains only complete cases**, i.e., raster cells
#'   having a value for every layer in the stack.
#'
#' @export
#' @examples
#' library(ggplot2)
#' library(reshape2)
#' library(raster)
#'
#' ## CREATE A DUMMY RASTER ##
#' r <- raster(
#'
#'   matrix(c(NA,100,100,NA,100,100,0,0,0),
#'          nrow = 3,
#'          ncol = 3,
#'          byrow = TRUE) )
#'
#'
#' ## VISUALIZE RASTER CELL NUMBERS ##
#' m <- matrix(1:9, nrow=3, ncol=3, byrow = TRUE)
#'
#' m <- t(m)[,nrow(m):1]
#' m <- melt(m, value.name = "cell")
#'
#' ggplot(m, aes(x=Var1, y=Var2)) +
#'   geom_tile(colour="gray90", lwd=1.5, show.legend = FALSE) +
#'   coord_fixed(ratio=1) +
#'   geom_text(aes(label=cell), color = "white", family=c("serif"), size=8) +
#'   theme_void()
#'
#'
#' ## VISUALIZE RASTER DUMMY LAYER ##
#' r_plot <- as.matrix(r)
#' r_plot <- t(r_plot)[,nrow(r_plot):1]
#' r_plot <- melt(r_plot, value.name = "dummy")
#'
#' r_plot$nas <- NA
#' r_plot$nas[which(is.na(r_plot$dummy))] <- 'NA'
#'
#'
#' ggplot(r_plot, aes(x=Var1, y=Var2)) +
#'   geom_tile(colour="gray90", lwd=1.5, show.legend = FALSE) +
#'   coord_fixed(ratio=1) +
#'   geom_text(aes(label=dummy), color = "white", family=c("serif"), size=8, na.rm=TRUE) +
#'   geom_text(aes(label=nas), color = "red", family=c("serif"), size=8, fontface='bold', na.rm=TRUE) +
#'   theme_void()
#'
#'
#' ## VISUALIZE ATTRIBUTE TABLE ##
#'
#' # Note that cells 1 and 4 have missing values and therefore are not included in the table
#' attTbl(r, var_names = c("dummy_var"))


attTbl <- function(rstack, var_names = NULL){

  dt <- data.frame(Cell=1:raster::ncell(rstack), raster::values(rstack))
  dt <- dt[stats::complete.cases(dt),]
  row.names(dt) <- NULL

  if( !is.null(var_names) ) {

    if( length(var_names) != length(names(dt)) -1 ){

      stop("var_names length should be equal to the number of layers in rstack")

    } else {

      names(dt)[names(dt) != "Cell"] <- var_names

    }

  }

  return(dt)

}


#' List of neighborhoods
#'
#' Computes the neighborhoods of the cells of a \code{Raster*} object.
#' Neighborhoods are not computed for cells with missing values.
#'
#' @param rstack \code{Raster*} object.
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
#'   * When the argument \code{rNumb = TRUE}, neighbors with one or more missing
#'   values are omitted;
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
#' @seealso [nbg8()], [attTbl()]
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
#'   Raster* object.
#'
#' @export
#' @examples
#' library(ggplot2)
#' library(reshape2)
#' library(raster)
#'
#' ## CREATE A DUMMY RASTER AND COMPUTE ITS ATTRIBUTE TABLE##
#' r <- raster(
#'
#'   matrix(c(NA,100,100,NA,100,100,0,0,0),
#'          nrow = 3,
#'          ncol = 3,
#'          byrow = TRUE) )
#'
#' at <- attTbl(r, var_names = c("dummy_var"))
#'
#'
#' ## VISUALIZE RASTER CELL NUMBERS ##
#' m <- matrix(1:9, nrow=3, ncol=3, byrow = TRUE)
#' m <- t(m)[,nrow(m):1]
#' m <- melt(m, value.name = "cell")
#'
#' ggplot(m, aes(x=Var1, y=Var2)) +
#'   geom_tile(colour="gray90", lwd=1.5, show.legend = FALSE) +
#'   coord_fixed(ratio=1) +
#'   geom_text(aes(label=cell), color = "white", family=c("serif"), size=8) +
#'   theme_void()
#'
#'
#' ## VISUALIZE RASTER DUMMY LAYER ##
#' r_plot <- as.matrix(r)
#' r_plot <- t(r_plot)[,nrow(r_plot):1]
#' r_plot <- melt(r_plot, value.name = "dummy")
#' r_plot$nas <- NA
#' r_plot$nas[which(is.na(r_plot$dummy))] <- 'NA'
#'
#' ggplot(r_plot, aes(x=Var1, y=Var2)) +
#'   geom_tile(colour="gray90", lwd=1.5, show.legend = FALSE) +
#'   coord_fixed(ratio=1) +
#'   geom_text(aes(label=dummy), color = "white", family=c("serif"), size=8, na.rm=TRUE) +
#'   geom_text(aes(label=nas), color = "red", family=c("serif"), size=8, fontface='bold', na.rm=TRUE) +
#'   theme_void()
#'
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


ngbList <- function(rstack, rNumb = FALSE, attTbl = NULL){

  ind <- which( stats::complete.cases( as.data.frame(raster::values(rstack)) ) )
  nbs <- nbg8(raster::nrow(rstack[[1]]), raster::ncol(rstack[[1]]))

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
