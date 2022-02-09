#' Test conditions for all cells
#'
#' Evaluate conditions for unclassified cells and classify them if conditions
#' are true.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param conditions character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{class}. If there is a
#'   \code{classVector} input, the classification number is only assigned to
#'   \code{classVector} NA-cells unless the argument \code{overwrite_class =
#'   TRUE}. See \code{\link{conditions}} for more details.
#' @param classVector numeric vector, if provided, it defines the cells in the
#'   attribute table that have already been classified. See
#'   \code{\link{conditions}} for more information about class vectors.
#' @param class numeric, the classification number to assign to all cells that
#'   meet the function conditions.
#' @param overwrite_class logic, if there is a \code{classVector} input,
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
#' ################################################################################
#' # LOAD DATA
#' ################################################################################
#' library(terra)
#' library(scapesClassification)
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
#' cv1 <- cond.4.all(attTbl = at, conditions = "dummy_var == 1", class = 1)
#'
#' unique(cv1) # one class (class 1)
#'
#' # update class vector `cv1`
#' # conditions: "dummy_var <= 3"
#' cv2   <- cond.4.all(attTbl = at, conditions = "dummy_var <= 3", class = 2,
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
#' par(mar = c(4, 0.5, 4, 0.5), mfrow=c(1,2))
#'
#' # 1)
#' plot(r_cv1, type="classes", axes=FALSE, legend=FALSE, asp=NA, col="#78b2c4", colNA="#818792")
#' text(r)
#' title("COND.4.ALL", adj = 0.0, line = 1,
#'       sub = "New class vector
#'       rule: 'dummy_var == 1'; class: 1")
#' legend("bottomleft", bg = "white", fill = c("#78b2c4", "#818792"),
#'        legend = c("Class 1", "Unclassified cells"))
#'
#' # 2)
#' plot(r_cv2, type="classes", axes=FALSE, legend=FALSE, asp=NA, col=c("#78b2c4","#cfad89"),
#'      colNA="#818792")
#' text(r)
#' title("COND.4.ALL", adj = 0.0, line = 1,
#'       sub = "Update class vector (class 1 is not overwritten)
#'       rule: 'dummy_var <= 3'; class: 2")
#' legend("bottomleft", bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Class 1", "Class 2", "Unclassified cells"))

cond.4.all <- function(attTbl,
                       conditions,
                       classVector = NULL,
                       class,
                       overwrite_class = FALSE) {
  if (is.null(classVector)) {
    classVector <- rep(as.integer(NA), NROW(attTbl))
  }

  if (!overwrite_class) {
    conditions <- paste("(", conditions, ")", "& is.na(classVector)")
  }

  v_ab <-
    names(attTbl)[stringr::str_detect(conditions, paste0("\\b", names(attTbl), "\\b"))]
  for (v in v_ab) {
    conditions <-
      stringr::str_replace_all(conditions, paste0("\\b", v, "\\b"), paste0("attTbl$", v))
  }

  classVector[which(eval(parse(text = conditions)))] <- class

  return(classVector)

}
