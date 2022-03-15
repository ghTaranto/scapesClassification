#' Classify All Unclassified Cells
#'
#' Classify all cells in \code{classVector} that have not yet been classified
#' based on contiguity and continuity conditions.
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
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. See \code{\link{conditions}} for more details about class
#'   vectors.
#'
#' @details The neighborhood of unclassified cells is considered. Among
#'   neighbors, the class with the highest number of members is assigned to the
#'   unclassified cell. If two or more classes have the same number of members,
#'   then one of these classes is assigned randomly to the unclassified cell.
#'   The function considers _class continuity_, thus, even cells that at first
#'   were not contiguous to any class will be classified if continuous with at
#'   least one cell having a class (see \code{\link{conditions}}).
#'
#' @seealso [attTbl()], [ngbList()], [conditions()]
#'
#' @export
#' @examples
#' # DUMMY DATA
#' ############################################################################
#' # LOAD LIBRARIES
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
#' # CLASSIFY.ALL
#' ################################################################################
#' # compute example class vector
#' cv <- cond.4.all(attTbl = at, cond = "dummy_var <= 1", class = 1)
#' # update example calss vector
#' cv <- cond.4.all(attTbl = at, cond = "dummy_var <= 3", class = 2,
#'                  classVector = cv) # input previous class vector
#'
#' # classify all unclassified cells
#' ca <- classify.all(attTbl = at, ngbList = nbs, rNumb = TRUE, classVector = cv)
#'
#'
#' # Convert class vectors into rasters
#' r_cv <- cv.2.rast(r, at$Cell, classVector = cv)
#' r_ca <- cv.2.rast(r, at$Cell,classVector = ca)

#' ################################################################################
#' # PLOTS
#' ################################################################################
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(3, 1, 5, 1)
#'
#' # 1)
#' plot(r_cv, type="classes", axes=FALSE, legend=FALSE, asp=NA, mar=m,
#'      colNA="#818792", col=c("#78b2c4", "#cfc1af"))
#' text(r)
#' mtext(side=3, line=2, adj=0, cex=1, font=2, "COND.4.ALL")
#' mtext(side=3, line=1, adj=0, cex=0.9, "Step1: 'dummy_var<=1', Class: 1")
#' mtext(side=3, line=0, adj=0, cex=0.9, "Step2: 'dummy_var<=3', Class: 2")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfc1af", "#818792"),
#'        legend = c("Class 1", "Class 2", "Unclassified cells"))
#'
#' # 2)
#' plot(r_ca, type="classes", axes=FALSE, legend=FALSE, asp=NA,  mar=m,
#'      colNA="#818792", col=c("#78b2c4", "#cfc1af"))
#' text(r)
#' mtext(side=3, line=2, adj=0, cex=1, font=2, "CLASSIFY.ALL")
#' mtext(side=3, line=1, adj=0, cex=0.9, "Classify all unclassified cells")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfc1af", "#818792"),
#'        legend = c("Class 1", "Class 2"))
#' par(oldpar)

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
