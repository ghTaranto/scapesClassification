#' Test conditions for all cells
#'
#' Evaluate conditions for unclassified cells and classify them if conditions
#' are true.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#' \code{\link{attTbl}}.
#' @param cond character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{class}. If there is a
#'   \code{classVector} input, the classification number is only assigned to
#'   unclassified cells unless the argument \code{ovw_class = TRUE}. See
#'   \code{\link{conditions}} for more details.
#' @param classVector numeric vector, if provided, it defines the cells in the
#'   attribute table that have already been classified. See
#'   \code{\link{conditions}} for more information about class vectors.
#' @param class numeric, the classification number to assign to all cells that
#'   meet the function conditions.
#' @param ovw_class logic, if there is a \code{classVector} input,
#'   reclassify cells that were already classified and that meet the function
#'   conditions.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. If there is no \code{classVector} input, the function returns
#'   a new class vector. See \code{\link{conditions}} for more details about
#'   class vectors.
#'
#' @details \itemize{ \item The function evaluates the conditions of the
#'   argument \code{conditions} for all unclassified cells (i.e.,
#'   \code{classVector} NA-cells).
#'
#'   \item Cells that meet the function conditions are classified as indicted by
#'   the argument \code{class}.
#'
#'   \item _Absolute test cell conditions_ can be used (see
#'   \code{\link{conditions}}).}
#'
#' @seealso [conditions()], [attTbl()], [cond.4.nofn()], [cond.reclass()]
#'
#' @export
#' @examples
#' # DUMMY DATA
#' ################################################################################
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
#' # COND.4.ALL
#' ################################################################################
#' # compute new class vector
#' # conditions: "dummy_var == 1"
#' cv1 <- cond.4.all(attTbl = at, cond = "dummy_var <= 1", class = 1)
#'
#' unique(cv1) # one class (class 1)
#'
#' # update class vector `cv1`
#' # conditions: "dummy_var <= 3"
#' cv2   <- cond.4.all(attTbl = at, cond = "dummy_var <= 3", class = 2,
#'                     classVector = cv1) # input previous class vector
#'
#' unique(cv2) # two classes (class 1 and class 2)
#'
#' # convert class vector 2 raster
#' r_cv1 <- cv.2.rast(r, at$Cell, classVector = cv1)
#' r_cv2 <- cv.2.rast(r, at$Cell, classVector = cv2)
#'
#' ################################################################################
#' # PLOTS
#' ################################################################################
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(4.5, 0.5, 2, 3.2)
#'
#' # 1.
#' r_cv1[which(is.na(values(r_cv1)))] <- 10
#' plot(r_cv1, type="classes", mar=m, col=c("#78b2c4","#818792"), axes=FALSE,
#'      plg=list(x=1, y=1, cex=.80, title="Classes",legend=c("1", "NA")))
#' text(r); lines(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "1. COND.4.ALL")
#' mtext(side=3, line=0, adj=0, cex=0.9, "New class vector")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var <= 1'")
#' mtext(side=1, line=1, cex=0.9, adj=0, "Class: 1")
#'
#' # 2.
#' r_cv2[which(is.na(values(r_cv2)))] <- 10
#' plot(r_cv2, type="classes", mar=m, col=c("#78b2c4","#cfad89","#818792"), axes=FALSE,
#'      plg=list(x=1, y=1, cex=.80, title="Classes",legend=c("1", "2", "NA")))
#' text(r); lines(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "2. COND.4.ALL")
#' mtext(side=3, line=0, adj=0, cex=0.9, "Update class vector (class 1 not overwritten)")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var <= 3'")
#' mtext(side=1, line=1, cex=0.9, adj=0, "Class: 2")
#' par(oldpar)

cond.4.all <- function(attTbl,
                       cond,
                       classVector = NULL,
                       class,
                       ovw_class = FALSE) {

  if (is.null(classVector)) {
    classVector <- rep(as.integer(NA), NROW(attTbl))
  }

  # PARSE CONDITION
  cond <- cond.parse(names(attTbl), cond)

  ctype <- names(cond[[2]])[lengths(cond[[2]])>0]
  if(!all("v_ab" == ctype)){
    stop("Only absolute test cell conditions allowed for `cond.4.all` (see ?conditions)")
  }

  cond <- gsub("\\bl_ab\\b", "attTbl", cond[[1]], perl = TRUE)
  if(!ovw_class){
    cond <- paste(cond, "& is.na(classVector)")
  }
  cond <- parse(text = cond)

  # EVALUATE CONDITION
  classVector[which(eval(parse(text = cond)))] <- class

  return(classVector)
}

