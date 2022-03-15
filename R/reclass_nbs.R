#' Reclassify neighbors
#'
#' Evaluate if members of two classes are contiguous and, if they are, one of
#' them is reclassified.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}.
#' @param rNumb logic, the neighborhoods of the argument \code{ngbList} are
#'   identified by cell numbers (\code{rNumb=FALSE}) or by row numbers
#'   (\code{rNumb=TRUE}) (see \code{\link{ngbList}}). It is advised to use row
#'   numbers for large rasters.
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified. See \code{\link{conditions}} for more
#'   information about class vectors.
#' @param nbs_of numeric or numeric vector, indicates the class(es) of focal and
#'   anchor cells.
#' @param class numeric or numeric vector, cells of classes \code{class}
#'   adjacent to cells belonging to one of the classes of \code{nbs_of} are
#'   reclassified as indicated by the argument \code{reclass}.
#' @param reclass numeric, the classification number to assign to all cells that
#'   meet the function conditions.
#' @param reclass_all logic, all cells of class \code{class} are also
#'   reclassified if they are connected to a reclassified cell.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. See \code{\link{conditions}} for more information about class
#'   vectors.
#'
#' @details \itemize{ \item The function evaluates if a cell of class
#'   \code{class} is adjacent to a cell of class \code{nbs_of} and, if it is, it
#'   is reclassifies as indicated by the argument \code{reclass}.
#'
#'   \item If the argument \code{reclass_all = TRUE}, all cells of class
#'   \code{class} are also reclassified if they are connected to a reclassified
#'   cell.}
#'
#' @seealso [attTbl()], [ngbList()], [cond.reclass()], [classify.all()]
#'
#' @export
#' @examples
#' library(scapesClassification)
#' library(terra)
#'
#' # LOAD THE DUMMY RASTER
#' r <- list.files(system.file("extdata", package = "scapesClassification"),
#'                 pattern = "dummy_raster\\.tif", full.names = TRUE)
#' r <- terra::rast(r)
#'
#' # COMPUTE THE ATTRIBUTE TABLE
#' at <- attTbl(r, "dummy_var")
#'
#' # COMPUTE THE LIST OF NEIGBORHOODS
#' nbs <- ngbList(r)
#'
#' ################################################################################
#' # RECLASS.NBS
#' ################################################################################
#'
#' # Compute an inital class vector with `cond.4.all`
#' cv <- cond.4.all(attTbl = at, cond = "dummy_var > 5", class = 1)
#'
#' # Update the class vector with a second class
#' cv <- cond.4.all(attTbl = at, cond = "dummy_var >= 2", class = 2,
#'                  classVector = cv)
#'
#'
#' # Reclassify cells of class 2 adjacent to cells of class 1
#'
#' # reclass_all = FALSE
#' rc1 <- reclass.nbs(attTbl = at, ngbList = nbs,
#'
#'                    # CLASS VECTOR `cv`
#'                    classVector = cv,
#'
#'                    # CELLS OF CLASS...
#'                    class = 2,
#'
#'                    # ...ADJACENT TO CELLS OF ANOTHER CLASS...
#'                    nbs_of = 1,
#'
#'                    # ...WILL BE RECLASSIFIED...
#'                    reclass = 3,
#'
#'                    # NO MORE RECLASSIFICATIONS
#'                    reclass_all = FALSE)
#'
#' # reclass_all = TRUE
#' rc2 <- reclass.nbs(attTbl = at, ngbList = nbs,
#'
#'                    # CLASS VECTOR `cv`
#'                    classVector = cv,
#'
#'                    # CELLS OF CLASS...
#'                    class = 2,
#'
#'                    # ...ADJACENT TO CELLS OF ANOTHER CLASS...
#'                    nbs_of = 1,
#'
#'                    # ...WILL BE RECLASSIFIED...
#'                    reclass = 3,
#'
#'                    # ...AND SO ALL CELLS OF CLASS 1 CONNECTED TO A RECLASSIFIED CELL
#'                    reclass_all = TRUE)
#'
#' # Convert class vectors to rasters
#' r_cv  <- cv.2.rast(r, at$Cell,classVector = cv, plot = FALSE)
#' r_rc1 <- cv.2.rast(r, at$Cell,classVector = rc1, plot = FALSE)
#' r_rc2 <- cv.2.rast(r, at$Cell,classVector = rc2, plot = FALSE)
#'
#' ################################################################################
#' # PLOTS
#' ################################################################################
#' oldpar <- par(mfrow = c(2,2))
#' m = c(0.1, 3.5, 3.2, 3.5)
#'
#'
#' # 1)
#' plot(r_cv, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar=m,
#'      colNA="#818792", col=c("#1088a0", "#78b2c4"))
#' text(r)
#' mtext(side=3, line=2, adj=0, cex=1, font=2, "COND.4.ALL")
#' mtext(side=3, line=1, adj=0, cex=0.9, "Step1: 'dummy_var > 5', class: 1")
#' mtext(side=3, line=0, adj=0, cex=0.9, "Step2: 'dummy_var > 3', class: 2")
#' legend("bottomright", ncol = 1, bg = "white", y.intersp= 1.2,
#'        legend = c("Class 1", "Class 2", "Unclassified cells"),
#'       fill = c("#1088a0", "#78b2c4", "#818792"))
#'
#' # 2)
#' plot(r_rc1, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar=m,
#'      colNA="#818792", col=c("#1088a0", "#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=2, adj=0, cex=1, font=2, "RECLASS.NBS")
#' mtext(side=3, line=1, adj=0, cex=0.9, "Reclass: class 2 adjacent to class 1")
#' mtext(side=3, line=0, adj=0, cex=0.9, "reclass_all = FALSE")
#' legend("bottomright", ncol = 1, bg = "white", y.intersp= 1.2,
#'        legend = c("Reclassified cells"), fill = c("#cfad89"))
#'
#' # 3)
#' plot(r_rc2, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar=m,
#'      colNA="#818792", col=c("#1088a0", "#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=2, adj=0, cex=1, font=2, "RECLASS.NBS")
#' mtext(side=3, line=1, adj=0, cex=0.9, "Reclass: class 2 adjacent to class 1")
#' mtext(side=3, line=0, adj=0, cex=0.9, "reclass_all = TRUE")
#' legend("bottomright", ncol = 1, bg = "white", y.intersp= 1.2,
#'        legend = c("Reclassified cells"), fill = c("#cfad89"))
#' par(oldpar)

reclass.nbs  <- function(attTbl,
                         ngbList,
                         rNumb = FALSE,
                         classVector,
                         nbs_of,
                         class,
                         reclass,
                         reclass_all = TRUE) {

  if(!rNumb){
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
