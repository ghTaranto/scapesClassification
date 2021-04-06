#' Raster Object to Data Frame
#'
#' Returns the attribute table (\code{data.frame}) used by other
#' \code{scapesClassification} functions. Each row corresponds to a cell on the
#' \code{Raster*} object.
#'
#' @param rstack \code{Raster*} object.
#' @param var_names character vector, re-name the variables of the
#'   \code{Raster*} object in the \code{data.frame}.
#'
#' @return data.frame
#'
#' @details Returns a \code{data.frame} with a number of columns equal to the
#'   number of layers in \code{rstack} plus one and with a maximum number of
#'   rows equal to the number of cells on \code{rstack}. Its first column, named
#'   \code{'Cell'}, refers to positions on the \code{Raster*} object with
#'   complete cases, i.e., raster cells having a value for every layer in the
#'   stack.
#'
#' @note The attribute table contains only complete cases, i.e., raster cells
#'   having a value for every layer in the stack.
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
#' attTbl(raster::stack(r1,r2))
#' attTbl(raster::stack(r1,r2), var_names = c("var_r1", "var_r2"))


attTbl <- function(rstack, var_names = NULL){

  dt <- data.frame(Cell=1:raster::ncell(rstack), raster::values(rstack))
  dt <- dt[stats::complete.cases(dt),]

  if( !is.null(var_names) ) {
    names(dt)[names(dt) %in% names(rstack)] <- var_names
  }

  return(dt)

}


#' Eight Neighbors for Complete Cases
#'
#' Return the 8-neighbors, as cell numbers (or index), of cells on a
#' \code{Raster*} object with complete cases, i.e., raster cells having a value
#' for every layer in the stack.
#'
#' @param rstack \code{Raster*} object.
#' @param index Return neighbor list as index position in the attribute table.
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}. It is required only if the argument \code{index =
#'   TRUE}
#'
#' @return List of integer vectors.
#'
#' @details The names of the list correspond to focal cell numbers
#'   (\code{as.character}), the integer vectors within the list contain the
#'   8-neighbors of each focal cell.
#'
#' @seealso [nbg8()], [attTbl()]
#'
#' @note There is a correspondence between the indices of the attribute table
#'   (\code{\link{attTbl}}) and the indices of the 8-neighbors list. For
#'   instance, the first element in the 8-neighbors list corresponds to the
#'   neighbors of the cell whose values are stored in the first row of the
#'   attribute table.
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


ngbList <- function(rstack, index = FALSE, attTbl = NULL){

  ind <- which( stats::complete.cases( as.data.frame(raster::values(rstack)) ) )
  nbs <- nbg8(raster::nrow(rstack[[1]]), raster::ncol(rstack[[1]]))

  nbs <- nbs[ind]

  # CONVERT nbs FORM CELL IDS TO CELL INDECES
  if(index){
    fct    <- rep(seq_along(lengths(nbs)), lengths(nbs))
    nbs    <- match(unlist(nbs), attTbl$Cell)
    no_nas <- !is.na(nbs)
    nbs    <- nbs[no_nas]
    fct    <- fct[no_nas]

    nbs    <- split(nbs, fct)
  }


  return(nbs)

}
