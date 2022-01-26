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
#'   that have already been classified.
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
#'   the function.
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
#'
#' library(raster)
#' library(scapesClassification)
#'
#' # LOAD THE DUMMY RASTER
#' r <- list.files(system.file("extdata", package = "scapesClassification"),
#'                 pattern = "dummy_raster\\.tif", full.names = TRUE)
#' r <- raster(r)
#'
#' # COMPUTE THE ATTRIBUTE TABLE
#' at <- attTbl(r, "dummy_var")
#'
#' # COMPUTE THE LIST OF NEIGBORHOODS
#' nbs <- ngbList(r)
#'
#' # COMPUTE A CLASS VECTOR
#' ################################################################################
#' # conditions: "dummy_var > 5"
#' # class: 1
#'
#' cv <- cond.4.all(attTbl = at, conditions = "dummy_var > 5", class = 1)
#'
#' # UPDATE THE CLASS VECTOR
#' ################################################################################
#' # conditions: "dummy_var > 3"
#' # class: 2
#'
#' cv <- cond.4.all(attTbl = at, conditions = "dummy_var >= 2", class = 2,
#'
#'                  classVector = cv)
#'
#'
#' # RECLASSIFY CELL OF CLASS 2 ADJACENT TO CELL OF CLASS 1
#' ################################################################################
#' # class: 2
#' # adjacent to class: 1
#' # new class: 3
#' # reclass_all = FALSE
#'
#' # RECLASSIFY NEIGHBORS
#' rc1 <- reclass.nbs(attTbl = at,
#'                    ngbList = nbs,
#'
#'                    # CLASS VECTOR COMPUTED WITH THE RULE "dummy_var > dummy_var{}"
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
#'
#' # RECLASSIFY ALL NEIGHBORS
#' rc2 <- reclass.nbs(attTbl = at,
#'                    ngbList = nbs,
#'
#'                    # CLASS VECTOR COMPUTED WITH THE RULE "dummy_var > dummy_var{}"
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
#' ################################################################################
#'
#' # CONVERT CLASS VECTORs INTO RASTERS
#' r_cv  <- cv.2.rast(r, at$Cell,classVector = cv, plot = FALSE)
#' r_rc1 <- cv.2.rast(r, at$Cell,classVector = rc1, plot = FALSE)
#' r_rc2 <- cv.2.rast(r, at$Cell,classVector = rc2, plot = FALSE)
#'
#' # SET PLOT LAYOUT
#' par(mfrow=c(2,2), mar=c(3, 2, 4, 2))
#'
#' # PLOT RESULTS
#'
#' # 1.
#' plot(r_cv, axes=FALSE, box=FALSE, legend = FALSE, asp = NA,
#'      colNA="#818792", col=c("#1088a0", "#78b2c4"))
#'
#' # REFERENCE PLOT 1
#' text(r)
#' title("STEPS 1&2. COND.4.ALL", adj = 0.0, line = 1,
#' sub="Step1. Rule: 'dummy_var > 5'; class: 1\nStep2. Rule: 'dummy_var > 3'; class: 2")
#'
#' legend("bottomright", ncol = 1, bg = "white",
#'        legend = c("Class 1", "Class 2", "Unclassified cells"),
#'        fill = c("#1088a0", "#78b2c4", "#818792"))
#'
#' # 2.
#' plot(r_rc1, axes=FALSE, box=FALSE, legend = FALSE, asp = NA,
#'      colNA="#818792", col=c("#1088a0", "#78b2c4", "#cfad89"))
#'
#' # REFERENCE PLOT 2
#' text(r)
#' title("STEP 3a. RECLASS.NBS", adj = 0.0, line = 1,
#' sub="reclass_all = FALSE ->\nreclass based on 'cell contiguity'")
#'
#' legend("bottomright", ncol = 1, bg = "white",
#'        legend = c("Class 1", "Class 2", "Reclassified cells",
#'                   "Unclassified cells"),
#'        fill = c("#1088a0", "#78b2c4", "#cfad89", "#818792"))
#'
#' # 3.
#' plot(r_rc2, axes=FALSE, box=FALSE, legend = FALSE, asp = NA,
#'      colNA="#818792", col=c("#1088a0", "#78b2c4", "#cfad89"))
#'
#' # REFERENCE PLOT 3
#' text(r)
#' title("STEP 3b. RECLASS.NBS", adj = 0.0, line = 1,
#' sub="reclass_all = TRUE ->\nreclass based on 'cell contiguity' and 'cell continuity'")
#' legend("bottomright", ncol = 1, bg = "white",
#'        legend = c("Class 1", "Reclassified cells",
#'                   "Unclassified cells"),
#'        fill = c("#1088a0", "#cfad89", "#818792"))

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


#' Classify All Unclassified Cells
#'
#' Classify all cells in \code{classVector} that have not yet been classified
#' based on contiguity conditions.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{ngbList}}).
#' @param rNumb logic, \code{ngbList} contain the neighbors index position
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
                          rNumb = FALSE,
                          classVector){

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
