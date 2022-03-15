#' Identify local maxima or minima
#'
#' Identify local maxima or minima on a raster surface.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}.
#' @param rNumb logic, the neighborhoods of the argument \code{ngbList} are
#'   identified by cell numbers (\code{rNumb=FALSE}) or by row numbers
#'   (\code{rNumb=TRUE}) (see \code{\link{ngbList}}). It is advised to use row
#'   numbers for large rasters.
#' @param p_col character, the column of the attribute table over which maxima
#'   or minima are searched.
#' @param p_fun character, if 'max' the function searches for local maxima; if
#'   'min' the function searches for local minima.
#' @param p_edge logic, if false local maxima or minima are not searched on edge
#'   cells. Edge cells are considered cells on the edge of the raster and cell
#' neighboring NA-cells.
#'
#' @details \itemize{ \item A cell constitutes a _local maximum_ if its
#'   elevation value is larger than the values of all the cells in its
#'   neighborhood (see \code{\link{ngbList}}).
#'
#'   \item A cell constitutes a _local minimum_ if its elevation value is
#'   smaller than the values of all the cells in its neighborhood (see
#'   \code{\link{ngbList}}).}
#'
#' @return A \code{classVector} with peak cells identified by the numeric class
#'   \code{1}. See \code{\link{conditions}} for more details about class
#'   vectors.
#'
#' @seealso [conditions()], [attTbl()], [ngbList()]
#'
#' @export
#'
#' @examples
#' # DUMMY DATA
#' ################################################################################
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
#' ################################################################################
#'
#' # PEAK.CELL
#' ################################################################################
#' # p_edge = FALSE
#' pc_a <- peak.cell(attTbl = at, ngbList = nbs, rNumb = FALSE,
#'                   p_col = "dummy_var", p_fun = "max", p_edge = FALSE)
#'
#' # p_edge = TRUE
#' pc_b <- peak.cell(attTbl = at, ngbList = nbs, rNumb = FALSE,
#'                   p_col = "dummy_var", p_fun = "max", p_edge = TRUE)
#'
#' # CONVERT THE CLASS VECTORS INTO RASTERS
#' r_pca <- cv.2.rast(r, at$Cell, classVector = pc_a, plot = FALSE)
#' r_pcb <- cv.2.rast(r, at$Cell, classVector = pc_b, plot = FALSE)
#' ################################################################################
#'
#' #PLOTS
#' ###############################################################################
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(4, 1, 4, 1)
#'
#' # PLOT 1 - p_edge = FALSE
#' plot(r_pca, axes=FALSE, legend=FALSE, asp=NA, mar=m,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "PEAK.CELL")
#' mtext(side=3, line=0, adj=0, cex=0.9, "p_edge = FALSE")
#' legend("bottomright", bg = "white",
#'        legend = c("Peak cell", "Unclassified cells"),
#'        fill = c("#cfad89", "#818792"))
#'
#' # PLOT 2 - p_edge = TRUE
#' plot(r_pcb, axes=FALSE, legend=FALSE, asp=NA, mar=m,
#'     colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "PEAK.CELL")
#' mtext(side=3, line=0, adj=0, cex=0.9, "p_edge = TRUE")
#' legend("bottomright", bg = "white",
#'        legend = c("Peak cell", "Unclassified cells"),
#'        fill = c("#cfad89", "#818792"))
#' par(oldpar)


peak.cell <- function(attTbl,
                      ngbList,
                      rNumb = FALSE,
                      p_col,
                      p_fun = "max",
                      p_edge = FALSE){


  # TEST ARGUMENT p_col
  if( !(p_col %in% names(attTbl)) ){
    stop("'p_col' must be a column of 'attTbl'")
  }

  # TEST ARGUMENT p_fun
  if( !(p_fun %in% c("min", "max")) ){
    stop("'p_fun' can be either 'max' or 'min'")
  }

  # TEST CELL VECTOR
  tc <- attTbl[[p_col]]

  # TEST EDGE (CELLS WITH < 8 NEIGHBORS)
  lnbs <- lengths(ngbList)

  if(!p_edge){
    ind_noEdge <- which(lnbs == max(lnbs))
    ind_Edge   <- which(lnbs != max(lnbs))

    tc      <- tc[ind_noEdge]
    ngbList <- ngbList[ind_noEdge]

    lnbs <- lengths(ngbList)
  }

  ## CONVERT NBS FORM CELL IDS TO CELL INDECES
  fct  <- base::rep(seq_along(lnbs), lnbs)

  if(!rNumb){
    ngbList <- base::match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- base::split(ngbList, fct)
  }

  # NEIGHBORS VECTORS
  nb <- attTbl[[p_col]][unlist(ngbList)]
  nb <- base::split(nb, fct)

  # FIND LOCAL PEAKS
  if(p_fun == "max"){
    lp <- mapply(function(x,y) all(x>y), tc, nb)
  } else {
    lp <- mapply(function(x,y) all(x<y), tc, nb)
  }

  # TEST EDGE (CELLS WITH < 8 NEIGHBORS)
  if(!p_edge){
    lp0 <- logical(length = NROW(attTbl))
    lp0[ind_noEdge] <- lp
    lp0[!lp0] <- NA
    lp <- as.numeric(lp0)

  } else{
    lp[!lp] <- NA
    lp <- as.numeric(lp)
  }

  return(lp)
}
