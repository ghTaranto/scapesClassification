#' Reclassify Neighbors
#'
#' If members of two classes are contiguous, one of them is reclassified.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{ngbList}}).
#' @param nbsIndex logic, \code{ngbList} contain the neighbors index position
#'   in the attribute table (see \code{\link{ngbList}}).
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified.
#' @param nbs_of numeric or numeric vector, if a cell classified as \code{class}
#'   is adjacent to a cell classified as indicated by the argument
#'   \code{nbs_of}, then it is reclassified as indicated by the argument
#'   \code{reclass}.
#' @param class numeric or numeric vector, indicates the classes for which
#'   adjacent neighbors have to be evaluated.
#' @param reclass numeric, the numeric class to attribute to the cells meeting
#'   contiguity conditions.
#' @param reclass_all logic, all cells adjacent to a reclassified cell
#'   classified as \code{class} are reclassified as indicated by the argument
#'   \code{reclass}.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function.
#'
#' @details The function evaluates if a cell of class \code{class} is adjacent
#'   to one of the classes of the argument \code{nbs_of} and, if \code{TRUE},
#'   reclassifies it as indicated by the argument \code{reclass}. If the argument
#'   \code{reclass_all} is \code{TRUE}, then all cells adjacent to a
#'   reclassified cell classified as \code{class} are also reclassified.
#'
#' @seealso [attTbl()], [ngbList()], [cond.reclass()] and contiguity functions
#'   in [conditions()]
#'
#' @export

reclass.nbs  <- function(attTbl,
                         ngbList,
                         nbsIndex = FALSE,
                         classVector,
                         nbs_of,
                         class,
                         reclass,
                         reclass_all = TRUE) {

  if(!nbsIndex){
  # CONVERT NBS FORM CELL IDS TO CELL INDECES
    fct     <- rep(seq_along(lengths(ngbList)), lengths(ngbList))
    ngbList <- match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- split(ngbList, fct)

    rm(fct, no_nas)
  }

  # INITIALIZE ALGORITHM
  ind   <- which(classVector %in% class)

  ### NO MIN BORDER CONDITION ############################################################
  for (c in ind) {
    cind <- ngbList[[c]]

    if (length(cind) == 0) {
      next
    }
    if (isTRUE(any(classVector[cind] %in% nbs_of))) {
      classVector[c] <- reclass

      continue <- TRUE

      new_cell_ind0 <- c(c, cind)
      new_cell_ind  <- c(c, cind)

      while (continue & reclass_all) {
        continue <- F
        i        <- which(classVector[new_cell_ind] %in% class)

        if (length(i) > 0) {
          continue <- T
          classVector[new_cell_ind[i]] <- reclass


          new_cell_ind  <-
            setdiff(unlist(ngbList[new_cell_ind[i]]), new_cell_ind0)
          new_cell_ind0 <- unique(c(new_cell_ind, new_cell_ind0))

        } #if

      }#while

    } #if

  } #for
  ########################################################################################
  return(classVector)

}


#' Classify All Unclassified Cells
#'
#' Classify all cells in \code{classVector} that have not yet been classified
#' based on contiguity conditions.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{ngbList}}).
#' @param nbsIndex logic, \code{ngbList} contain the neighbors index position
#'   in the attribute table (see \code{\link{ngbList}}).
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function.
#'
#' @details For each unclassified cell its 8-cell neighborhood is considered.
#'   Among neighbors, the class with the highest number of members is assigned
#'   to the unclassified cell. If two or more classes have the same number of
#'   members then one of these classes is assigned randomly to the unclassified
#'   cell. If none of the neighbors have been classified then no class is
#'   assigned to the unclassified cell.
#'
#' @seealso [attTbl()], [ngbList()] and contiguity functions in [conditions()]
#'
#' @export

classify.all  <- function(attTbl,
                          ngbList,
                          nbsIndex = FALSE,
                          classVector){

  if(!nbsIndex){
  # CONVERT NBS FORM CELL IDS TO CELL INDECES
    fct     <- rep(seq_along(lengths(ngbList)), lengths(ngbList))
    ngbList <- match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- split(ngbList, fct)

    rm(fct, no_nas)
  }


  # INITIALIZE ALGORITHM, CONSIDER ONLY UNCLASSIFIED CELLS WITH SOME NEIGHBOURING CELL WITH CLASS
  ind <- which(is.na(classVector))

  nbs_ind1 <- ngbList[ind]
  fct      <- rep(seq_along(lengths(nbs_ind1)), lengths(nbs_ind1))

  nbs_ind1_vect <- !is.na(classVector[unlist(nbs_ind1)])
  anynbs_class  <- split(nbs_ind1_vect, fct)
  anynbs_class  <- sapply(anynbs_class, any, USE.NAMES = F)

  ind <- ind[anynbs_class]
  continue <- TRUE

  while(continue){

    continue <- FALSE
    k = 1
    list_new_cell_ind <- list()

    for( c in ind ){

      nind <- ngbList[[ c ]]

      cv_nind <- classVector[nind]
      freq    <- table(classVector[nind])
      fq_max  <- which(freq == max(freq))[1]

      classVector[ c ] <- as.numeric(names(fq_max))

      i <- which(is.na(cv_nind))

      if(length(i) > 0){

        list_new_cell_ind[[k]] <- nind[i]
        k <- k + 1

      } #if

    } #for

    if(!is.null(unlist(list_new_cell_ind))){

      ind <- unique(unlist(list_new_cell_ind))
      ind <- ind[ which(is.na(classVector[ind])) ]

      if(length(ind) > 0){
        continue <- TRUE

      }

    }

  } #while
  ########################################################################################
  return(classVector)

}
