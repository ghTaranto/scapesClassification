#' Attribute table
#'
#' Converts a \code{Raster*} object into an attribute table (\code{data.frame}).
#'
#' @param rstack \code{Raster*} object.
#' @param var_names character vector, names of variables used to store the
#'   values of the layers of the \code{Raster*} object.
#'
#' @return data.frame
#'
#' @details Attribute tables come with a column named **\code{"Cell"}** which
#'   stores raster cell numbers and associate each row of the attribute table
#'   with a cell of the raster object. Each of the remaining columns stores the
#'   values of a layer of the \code{Raster*} object. Note that only raster cells
#'   having no missing value in any layer (**complete cases**) are included in
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
#' Excluding cells on the edge of a raster, a cell with coordinates \code{(x,
#' y)} has 8 neighbors with coordinates: \code{(x±1, y)},  \code{(x, y±1)} and
#' \code{(x±1, y±1)}. The function computes the neighborhoods of all the cells
#' of a \code{Raster*} object. Cell with missing values are omitted.
#'
#' @param rstack \code{Raster*} object.
#' @param rNumb logic, the neighbors of a raster cell are identified by cell
#'   numbers if \code{rNumb=FALSE} or by row numbers if \code{rNumb=TRUE}. If
#'   the argument is true the argument \code{attTbl} cannot be NULL.
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}} (see \code{\link{attTbl}}). It is required only if the
#'   argument \code{rNumb=TRUE}.
#'
#'
#' @return Named list of integer vectors.
#'
#' @details **Neighborhoods (\code{rNumb=FALSE})**
#'
#'   * Neighbors are identified by their cell numbers if the argument
#'   \code{rNumb=FALSE}.
#'
#'   **Neighborhoods (\code{rNumb=TRUE})**
#'
#'   * Neighbors are identified by their positions in the attribute table (i.e.
#'   row numbers) if the argument \code{rNumb=TRUE} and the argument
#'   \code{attTbl!=NULL} (see \code{\link{attTbl}});
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
#'   * When \code{rNumb = FALSE}, the element name identify the raster cell to
#'   which the neighborhood refers. For instance, the element with name
#'   \code{"6"} stores the neighborhood of the raster cell \code{6}.
#'
#'   * When \code{rNumb = TRUE}, the element name identify the row number to
#'   which the neighborhood refers. For instance, the element with name
#'   \code{"6"} stores the neighborhood of the raster cell located in the 6th
#'   row of the attribute table (see \code{\link{attTbl}}).\cr
#'
#' @seealso [nbg8()], [attTbl()]
#'
#' @note There is a correspondence between the indices of the attribute table
#'   (\code{\link{attTbl}}) and the indices of the list of neighborhoods. For
#'   instance, the first element of the list corresponds to the neighbors of the
#'   cell stored in the first row of the attribute table.
#'
#' @export
#' @examples
#' ## r1, raster with cell numbers as values
#' r1 <- raster::raster(matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE))
#'
#' ## r2, raster with missing values in cell numbers 1 and 2
#' r2 <- raster::raster(matrix(1, nrow = 3, ncol = 4, byrow = TRUE))
#' r2[c(1,2)] <- NA
#'
#'
#' ngbList(raster::stack(r1,r2))
#' matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)


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
