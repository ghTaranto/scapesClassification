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


#' Eight Neighbors for Complete Cases
#'
#' The function \code{ngbList} returns the list of 8-neighbors for all the cells
#' of \code{rstack} with complete cases (i.e., cells without any missing value).
#' Neighbors are identified by cell numbers or by row numbers in the attribute
#' table (see \code{\link{attTbl}}) if the argument \code{rNumb = TRUE}. \cr\cr
#' The list of 8-neighbors is named. When \code{rNumb = FALSE}, names identify
#' the cell number for which the neighborhood was computed. When \code{rNumb =
#' TRUE}, names refers to row indices. For instance, the name \code{"6"} can be
#' used to call the neighborhood of cell number 6 when \code{rNumb = FALSE}.
#' However, when \code{rNumb = TRUE}, the name \code{"6"} can be used to call
#' the neighborhood of the raster cell stored in the 6th row of the attribute
#' table \code{attTbl}.\cr\cr When the argument \code{rNumb = TRUE}, only
#' neighbors with complete cases are considered, i.e. if a neighbor position
#' correspond to a raster cell with one or more missing values it is not
#' considered. Therefore, if a cell has all neighboring cells with missing
#' values, the 8-neighbor vector of that cell will have zero-length.
#'
#' @param rstack \code{Raster*} object.
#' @param rNumb logic, the 8-neighbors of each cell are identified by their row
#'   number in the attribute table. For instance, if cell 3 is located in row 2
#'   it will be identified by the number 2.
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}. It is required only if the argument \code{rNumb =
#'   TRUE}
#'
#' @return Named list of integer vectors.
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


ngbList <- function(rstack, rNumb = FALSE, attTbl = NULL){

  ind <- which( stats::complete.cases( as.data.frame(raster::values(rstack)) ) )
  nbs <- nbg8(raster::nrow(rstack[[1]]), raster::ncol(rstack[[1]]))

  nbs <- nbs[ind]

  # CONVERT nbs FORM CELL IDS TO CELL INDECES
  if(rNumb){
    fct    <- rep(seq_along(lengths(nbs)), lengths(nbs))
    nbs    <- match(unlist(nbs), attTbl$Cell)
    no_nas <- !is.na(nbs)
    nbs    <- nbs[no_nas]
    fct    <- fct[no_nas]

    nbs    <- split(nbs, fct)
  }


  return(nbs)

}
