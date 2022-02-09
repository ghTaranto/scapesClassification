#' Test conditions and reclassify
#'
#' Evaluate conditions for cells of a class and reclassify them if conditions
#' are true.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}. Only necessary if using a _neighborhood condition_
#'   (see \code{\link{conditions}}).
#' @param rNumb logic, the neighborhoods of the argument \code{ngbList} are
#'   identified by cell numbers (\code{rNumb=FALSE}) or by row numbers
#'   (\code{rNumb=TRUE}) (see \code{\link{ngbList}}). It is advised to use row
#'   numbers for large rasters.
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified. See \code{\link{conditions}} for more
#'   information about class vectors.
#' @param conditions character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{reclass}. See
#'   \code{\link{conditions}} for more details.
#' @param class numeric or numeric vector, indicates the class(es) for which
#'   conditions have to be evaluated.
#' @param reclass numeric, the classification number to assign to all cells that
#'   meet the function conditions.
#' @param fn_perc numeric value between 0 and 1. If a _neighborhood condition_
#'   is considered, test cells are classified if the conditions are true for at
#'   least as many evaluations as the ones specified by the argument
#'   \code{fn_perc} (see \code{\link{conditions}}).
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. See \code{\link{conditions}} for more information about class
#'   vectors.
#'
#' @details \itemize{ \item The function evaluates the conditions of the
#'   argument \code{conditions} for all cells in the classes of the
#'   argument \code{class}.
#'
#'   \item Cells that meet the function conditions are classified as indicted by
#'   the argument \code{reclass}.
#'
#'   \item Absolute test cell and neighborhood conditions can be used. The
#'   condition string can only include one neighborhood condition (\code{'{}'})
#'   (see \code{\link{conditions}}).}
#'
#' @seealso [conditions()], [attTbl()], [ngbList()]
#'
#' @export
#' @examples
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
#'
#
#' ################################################################################
#' # RECLASS.NBS
#' ################################################################################
#'
#' # Compute an example class vector
#' cv <- cond.4.all(attTbl = at, conditions = "dummy_var > 1", class = 1)
#'
#' # Reclassify cells
#' cr <- cond.reclass(attTbl = at, ngbList = nbs,
#'
#'                    # CLASS VECTOR COMPUTED WITH THE RULE "dummy_var > 1"
#'                    classVector = cv,
#'
#'                    # CELLS TO RECLASSIFY HAVE THIS CLASS
#'                    class = 1,
#'
#'                    # ABSOLUTE CONDITION
#'                    conditions = "dummy_var > 5",
#'
#'                    # NEW CLASSIFICATION NUMBER
#'                    reclass = 2)
#'
#' # Convert class vectors to rasters
#' r_cr <- cv.2.rast(r, at$Cell,classVector = cr, plot = FALSE)
#'
#'
#' ################################################################################
#' # PLOTS
#' ################################################################################
#' plot(r_cr, type="classes", axes=FALSE, legend = FALSE, asp=NA,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' title("COND.RECLASS", adj = 0.0, line = 1,
#'  sub =
#'  "Rule: 'dummy_var > 1'; Function: cond.4.all; Class: 1
#' Rule: 'dummy_var > 5'; Function: cond.reclass; Class: 2")
#' legend("bottomright", ncol = 1, bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Class 1","Class 2 (reclass)","Unclassified cells"))

cond.reclass <- function(attTbl,
                         ngbList = NULL,
                         rNumb = FALSE,
                         classVector,
                         conditions,
                         class,
                         reclass,
                         fn_perc = 1) {

  ### PARSE CONDITIONS
  v_ab <-
    names(attTbl)[stringr::str_detect(conditions, paste0("\\b", names(attTbl), "(?!\\[|\\{)", "\\b"))]
  v_fn <-
    names(attTbl)[stringr::str_detect(conditions, paste0(names(attTbl), "\\{\\}"))]

  for (v in v_ab) {
    conditions <-
      stringr::str_replace_all(conditions, paste0("\\b", v, "(?!\\[|\\{)", "\\b"), paste0("l_ab$", v))
  }
  for (v in v_fn) {
    conditions <-
      stringr::str_replace_all(conditions, paste0(v, "\\{\\}"), paste0("l_fn$", v))
  }

  cond_parsed <- parse(text = conditions)

  if (length(v_fn) != 0) {
    fn = TRUE
  } else{
    fn = FALSE
  }


  if(!fn){

    c <- which(classVector %in% class)

    l_ab <- lapply(as.list(v_ab), function(x)
      attTbl[[x]][c])
    names(l_ab) <- v_ab

    c <- c[eval(cond_parsed)]

    classVector[c] <- reclass

    return(classVector)

  }

  # CONVERT NBS FORM CELL IDS TO CELL INDECES
  if(!rNumb){
    fct     <- rep(seq_along(lengths(ngbList)), lengths(ngbList))
    ngbList <- match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- split(ngbList, fct)

    rm(fct, no_nas)
  }

  ###INITIALIZE ALGORITHM #########################################################################
  new_cell_id       <- which(classVector %in% class)

  l_ab <- list()
  l_fn <- list()

  ### RUN ALGORITHM #################################################################### while ####
  for (c in new_cell_id) {
    #

    n_ind    <- ngbList[[c]]

    if (length(n_ind) == 0 & fn) {
      next
    }


    ### TEST FOR CONDITIONS #################### while//for//if//relative_condition ####
    if (fn) {
      # IF CONSIDERING FOCAL NEIGHBORHOOD
      l_fn <- lapply(as.list(v_fn), function(x)
        attTbl[[x]][n_ind])
      names(l_fn) <- v_fn
    }

    l_ab <- lapply(as.list(v_ab), function(x)
      attTbl[[x]][c])
    names(l_ab) <- v_ab


    ##################################################### while//for//if//relative_condition ####

    ### TEST FOR CELLS MEETING CONDITIONS ########################### while//for//conditions ####
    ev_cond <- eval(cond_parsed)
    ev_cond <- ev_cond[!is.na(ev_cond)]

    if(length(ev_cond) == 0){next}

    if (fn & length(ev_cond) > 0) {
      rc <- sum(ev_cond) / length(n_ind) >= fn_perc
    } else {
      rc <- ev_cond
    }

    if (rc & length(ev_cond) > 0) {
      classVector[c] <- reclass
    }

  } #FOR ENDS

  return(classVector)

}
