#' Eight Neighbors
#'
#' Return the 8 neighbors, as cell numbers, of each cell on a raster. Cell numbers
#' start with 1 in the upper-left corner and increase from left to right and
#' from top to bottom.
#'
#' @param n_row Integer. The number of rows of a Raster or RasterStack object.
#' @param n_col Integer. The number of columns of a Raster or RasterStack
#'   object.
#'
#' @return List of length equal to the number of cells on a raster. The nth
#'   element of the list corresponds to the 8 adjacent cell numbers of the nth
#'   cell on the RasterStack object.
#'
#' @seealso [ngbList()]
#'
#' @export
#' @examples
#' ## Matrix m mocking a raster of 3 rows and 4 columns
#' m <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
#' m
#'
#' nbs <- nbg8(3, 4)
#' nbs


nbg8 <- function(n_row, n_col){

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
