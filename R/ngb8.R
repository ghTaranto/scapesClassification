#' Eight neighbors
#'
#' Return the 8 neighbors, as cell numbers, of each cell on a raster.
#'
#' @param n_row Integer. The number of rows of a Raster or object.
#' @param n_col Integer. The number of columns of a Raster object.
#'
#' @details A cell with coordinates \code{(x, y)} has 8 neighbors with
#'   coordinates: \code{(x±1, y)},  \code{(x, y±1)} and \code{(x±1, y±1)}. Cells
#'   on the edge of a raster have less than 8 neighbors. The function identifies
#'   the neighbors of a cell as cell numbers.
#'
#' @return Named list, the \code{nth} element of the list corresponds to the 8
#'   adjacent cell numbers of the \code{nth} cell on the \code{Raster*} object.
#'
#' @seealso [ngbList()]
#'
#' @export
#' @examples
#' ## Matrix m mocking a raster of 3 rows and 4 columns
#' m <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
#' m
#'
#' nbs <- ngb8(3, 4)
#' nbs


ngb8 <- function(n_row, n_col){

  n_cel <- n_row * n_col
  m     <- matrix(1:n_cel, n_row, n_col, byrow = T)

  nbs_m <- array(dim = c(n_row, n_col, 8))

  dimnames(nbs_m)[[3]] <- c("TL", "T", "TR", "L", "R", "BL", "B", "BR")

  # TL
  nbs_m[, , 1] <- m - n_col - 1
  # T
  nbs_m[, , 2] <- m - n_col
  # TR
  nbs_m[, , 3] <- m - n_col + 1
  # L
  nbs_m[, , 4] <- m - 1
  # R
  nbs_m[, , 5] <- m + 1
  #BL
  nbs_m[, , 6] <- m  + n_col - 1
  #B
  nbs_m[, , 7] <- m + n_col
  # BR
  nbs_m[, , 8] <- m + n_col + 1

  # Top row
  nbs_m[1, , c("TL", "T", "TR")]      <- NA
  # Bottom row
  nbs_m[n_row, , c("BL", "B", "BR")]  <- NA
  # Left column
  nbs_m[, 1 , c("TL", "L", "BL")]     <- NA
  # Right column
  nbs_m[, n_col , c("TR", "R", "BR")] <- NA


  # Array to vector
  nbs_m <- as.integer(nbs_m[1:(n_cel*8)])

  # Vector split index
  ind <- rep(seq_along(1:n_cel), 8)

  # Remove NAS
  no_nas <- which(!is.na(nbs_m))
  nbs_m  <- nbs_m[no_nas]
  ind    <- ind[no_nas]

  # Vector to list
  nbs_m <- split(nbs_m, ind)

  # Reorder list by row
  rbyrow <- m[1:n_cel]
  names(nbs_m) <- as.character(rbyrow)

  nbs_m <- nbs_m[order(as.numeric(names(nbs_m)))]

  return(nbs_m)

}
