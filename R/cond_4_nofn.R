#' Test conditions for neighbors and neighbors of neighbors
#'
#' Evaluate conditions for cells neighboring specific classes and classify them
#' if conditions are true.
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
#' @param class numeric, the classification number to assign to all cells that
#'   meet the function conditions.
#' @param nbs_of numeric or numeric vector, indicates the class(es) of focal and
#'   anchor cells. Conditions are only evaluated at positions adjacent to anchor
#'   and focal cells. If the classification number assigned with the argument
#'   \code{class} is also included in the argument \code{nbs_of}, the function
#'   takes into account _class continuity_ (see \code{\link{conditions}}).
#' @param cond character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{class}. The classification
#'   number is only assigned to unclassified cells unless the argument
#'   \code{ovw_class = TRUE}. See \code{\link{conditions}} for more details.
#' @param min.bord numeric value between 0 and 1. A test cell is classified if
#'   conditions are true and if among its bordering cells a percentage equal
#'   or greater than \code{min.bord} belong to one of the classes of
#'   \code{nbs_of}. Percentages are computed counting only valid neighbors
#'   (i.e., neighbors with complete cases).
#' @param max.iter integer, the maximum number of iterations.
#' @param peval numeric value between 0 and 1. If _absolute or relative
#'   neighborhood conditions_ are considered, test cells are classified if the
#'   number of positive evaluations is equal or greater than the percentage
#'   specified by the argument \code{peval} (see \code{\link{conditions}}).
#' @param directional logic, absolute or relative neighborhood conditions are
#'   tested using the _directional neighborhood_ (see \code{\link{conditions}}).
#' @param ovw_class logic, reclassify cells that were already classified and
#'   that meet the function conditions.
#' @param hgrowth logic, if true the classes in \code{nbs_of} are treated as
#'   discrete raster objects and the argument \code{class} is ignored.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. See \code{\link{conditions}} for more details about class
#'   vectors.
#'
#' @details \itemize{ \item The function evaluates the conditions of the
#'   argument \code{cond} for all unclassified cells in the neighborhood of
#'   focal and anchor cells (specified by the argument \code{nbs_of}).
#'   Unclassified cells are NA-cells in \code{classVector}.
#'
#'   \item Cells that meet the function conditions are classified as indicted by
#'   the argument \code{class}.
#'
#'   \item _Class continuity_ is considered if the classification number
#'   assigned with the argument \code{class} is also included in the argument
#'   \code{nbs_of}. This means that, at each iteration, newly classified cells
#'   become focal cells and conditions are tested in their neighborhood.
#'
#'   \item All types of conditions can be used. The condition string can only
#'   include one neighborhood condition (\code{'{}'}) (see
#'   \code{\link{conditions}}).}
#'
#'   **Homogeneous growth (\code{hgrowth})**
#'
#'   If the argument \code{hgrowth} is true the classes in \code{nbs_of} are
#'   treated as discrete raster objects and the argument \code{class} is
#'   ignored. Iterations proceed as follow:
#'
#'   * cells contiguous to the first element of \code{nbs_of} are evaluated
#'   against the classification rules and, when evaluations are true, cells are
#'   assigned to that element;
#'
#'   * the same process is repeated for cells contiguous to the second element
#'   of \code{nbs_of}, then for cells contiguous to the third element and so on
#'   until the last element of \code{nbs_of};
#'
#'   * once cells contiguous to the last element of \code{nbs_of} are evaluated
#'   the iteration is complete;
#'
#'   * cells classified in one iteration become focal cells in the next
#'   iteration;
#'
#'   * a new iteration starts as long as new cells were classified in the
#'   previous iteration and if the iteration number < \code{max.iter}.
#'
#' @seealso [conditions()], [attTbl()], [ngbList()]
#'
#' @export
#' @examples
#' # DUMMY DATA
#' ######################################################################################
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
#' # SET A DUMMY FOCAL CELL (CELL #25)
#' at$cv[at$Cell == 25] <- 0
#'
#' # SET FIGURE MARGINS
#' m <- c(2, 8, 2.5, 8)
#'
#' ######################################################################################
#' # ABSOLUTE TEST CELL CONDITION - NO CLASS CONTINUITY
#' ######################################################################################
#'
#' # conditions: "dummy_var >= 3"
#' cv1 <- cond.4.nofn(attTbl = at, ngbList = nbs,
#'
#'                    # CLASS VECTOR - INPUT
#'                    classVector = at$cv,
#'
#'                    # CLASSIFICATION NUMBER
#'                    class = 1,
#'
#'                    # FOCAL CELL CLASS
#'                    nbs_of = 0,
#'
#'                    # ABSOLUTE TEST CELL CONDITION
#'                    cond = "dummy_var >= 3")
#'
#' # CONVERT THE CLASS VECTOR INTO A RASTER
#' r_cv1 <- cv.2.rast(r, at$Cell,classVector = cv1, plot = FALSE)
#'
#' # PLOT
#' plot(r_cv1, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "CONDITION: ABSOLUTE TEST CELL")
#' mtext(side=3, line=0, adj=0, cex=1, "Class continuity: NO")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var >= 3'")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Focal cell", "Classified cells", "Unclassified cells"))
#'
#' ######################################################################################
#' # ABSOLUTE TEST CELL CONDITION - WITH CLASS CONTINUITY
#' ######################################################################################
#'
#' # conditions: "dummy_var >= 3"
#' cv2 <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = at$cv,
#'
#'                   # CLASSIFICATION NUMBER
#'                    class = 1,
#'
#'                    nbs_of = c(0,  # FOCAL CELL CLASS
#'                               1), # CLASSIFICATION NUMBER
#'
#'                    # ABSOLUTE CONDITION
#'                    cond = "dummy_var >= 3")
#'
#' # CONVERT THE CLASS VECTOR INTO A RASTER
#' r_cv2 <- cv.2.rast(r, at$Cell,classVector = cv2, plot = FALSE)
#'
#' # PLOT
#' plot(r_cv2, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "CONDITION: ABSOLUTE TEST CELL")
#' mtext(side=3, line=0, adj=0, cex=1, "Class continuity: YES")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var >= 3'")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Focal cell", "Classified cells", "Unclassified cells"))
#'
#' ######################################################################################
#' # ABSOLUTE NEIGHBORHOOD CONDITION
#' ######################################################################################
#'
#' # conditions: "dummy_var{} >= 3"
#' cv3 <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = at$cv, nbs_of = c(0,1), class = 1,
#'
#'                    # ABSOLUTE NEIGHBORHOOD CONDITION
#'                    cond = "dummy_var{} >= 3",
#'
#'                    # RULE HAS TO BE TRUE FOR 100% OF THE EVALUATIONS
#'                    peval = 1)
#'
#' # CONVERT THE CLASS VECTOR INTO A RASTER
#' r_cv3 <- cv.2.rast(r, at$Cell,classVector = cv3, plot = FALSE)
#'
#' #PLOT
#' plot(r_cv3, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "CONDITION: ABSOLUTE NEIGHBORHOOD")
#' mtext(side=3, line=0, adj=0, cex=1, "Class continuity: YES")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var{ } >= 3'")
#' mtext(side=1, line=0, cex=0.9, adj=1, "('{ }' cell neighborhood)")
#' mtext(side=1, line=1, cex=0.9, adj=0, "Fn_perc: 1 (100%)")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Focal cell", "Classified cells", "Unclassified cells"))
#'
#' ######################################################################################
#' # RELATIVE NEIGHBORHOOD CONDITION
#' ######################################################################################
#'
#' # conditions: "dummy_var > dummy_var{}"
#' cv4 <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = at$cv, nbs_of = c(0,1), class = 1,
#'
#'                    # RELATIVE NEIGHBORHOOD CONDITION
#'                    cond = "dummy_var > dummy_var{}",
#'
#'                    # RULE HAS TO BE TRUE FOR AT LEAST 60% OF THE EVALUATIONS
#'                    peval = 0.6)
#'
#'
#' # CONVERT THE CLASS VECTOR INTO A RASTER
#' r_cv4 <- cv.2.rast(r, at$Cell, classVector = cv4, plot = FALSE)
#'
#' #PLOT
#' plot(r_cv4, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "CONDITION: RELATIVE NEIGHBORHOOD")
#' mtext(side=3, line=0, adj=0, cex=1, "Class continuity: YES")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var > dummy_var{ }'")
#' mtext(side=1, line=0, cex=0.9, adj=1, "('{ }' cell neighborhood)")
#' mtext(side=1, line=1, cex=0.9, adj=0, "Fn_perc: 0.6 (60%)")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Focal cell", "Classified cells", "Unclassified cells"))
#'
#' ######################################################################################
#' # RELATIVE FOCAL CELL CONDITION
#' ######################################################################################
#'
#' # conditions: "dummy_var > dummy_var[]"
#' cv5 <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = at$cv, nbs_of = c(0,1), class = 1,
#'
#'                    # RELATIVE FOCAL CELL CONDITION
#'                    cond = "dummy_var > dummy_var[]")
#'
#'
#' # CONVERT THE CLASS VECTOR INTO A RASTER
#' r_cv5 <- cv.2.rast(r, at$Cell,classVector = cv5, plot = FALSE)
#'
#' #PLOT
#' plot(r_cv5, type="classes", axes=FALSE, legend = FALSE, asp = NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "CONDITION: RELATIVE FOCAL CELL")
#' mtext(side=3, line=0, adj=0, cex=1, "Class continuity: YES")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Rule: 'dummy_var > dummy_var[ ]'")
#' mtext(side=1, line=0, cex=0.9, adj=1, "('[ ]' focal cell)")
#' legend("bottomright", bg = "white", fill = c("#78b2c4", "#cfad89", "#818792"),
#'        legend = c("Focal cell", "Classified cells", "Unclassified cells"))
#'
#' ######################################################################################
#' # HOMOGENEOUS GROWTH
#' ######################################################################################
#'
#' # Dummy raster objects 1 and 2
#' ro <- as.numeric(rep(NA, NROW(at)))
#' ro[which(at$dummy_var == 10)] <- 1
#' ro[which(at$dummy_var == 8)] <- 2
#'
#' # Not homogeneous growth
#' nhg <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = ro,
#'                    nbs_of = 1, class = 1, # GROWTH ROBJ 1
#'                    cond = "dummy_var <= dummy_var[] & dummy_var != 1")
#'
#' nhg <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = nhg, # UPDATE nhg
#'                    nbs_of = 2, class = 2, # GROWTH ROBJ 2
#'                    cond = "dummy_var <= dummy_var[] & dummy_var != 1")
#'
#'
#' # Homogeneous growth
#' hg <- cond.4.nofn(attTbl = at, ngbList = nbs, classVector = ro,
#'                   nbs_of = c(1, 2), class = NULL,
#'                   cond = "dummy_var <= dummy_var[] & dummy_var != 1",
#'                   hgrowth = TRUE) # HOMOGENEOUS GROWTH
#'
#' # Convert class vectors into rasters
#' r_nhg <- cv.2.rast(r, at$Cell,classVector = nhg, plot = FALSE)
#' r_hg  <- cv.2.rast(r, at$Cell,classVector = hg, plot = FALSE)
#'
#' # Plots
#' oldpar <- par(mfrow = c(1,2))
#' m <- c(3, 1, 5, 1)
#'
#' # Original raster objects (for plotting)
#' r_nhg[at$dummy_var == 10] <- 3
#' r_nhg[at$dummy_var == 8]  <- 4
#'
#' r_hg[at$dummy_var == 10] <- 3
#' r_hg[at$dummy_var == 8]  <- 4
#' #t
#' # 1)
#' plot(r_nhg, type="classes", axes=FALSE, legend=FALSE, asp=NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfc1af", "#1088a0", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "RASTER OBJECTS GROWTH")
#' mtext(side=3, line=0, adj=0, cex=0.9, "Not homogeneous (hgrowth = FALSE)")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Growth rule:")
#' mtext(side=1, line=1, cex=0.9, adj=0, "'dummy_var<=dummy_var[ ] & dummy_var!=1''")
#' legend("topleft", bg = "white", y.intersp= 1.3,
#'        fill = c("#1088a0", "#cfc1af", "#78b2c4", "#cfc1af", "#818792"),
#'        legend = c("RO1", "RO2", "RO1 - growth", "RO2 - growth", "Unclassified cells"))
#' # 2)
#' plot(r_hg, type="classes", axes=FALSE, legend=FALSE, asp=NA, mar = m,
#'      colNA="#818792", col=c("#78b2c4", "#cfc1af", "#1088a0", "#cfad89"))
#' text(r)
#' mtext(side=3, line=1, adj=0, cex=1, font=2, "RASTER OBJECTS GROWTH")
#' mtext(side=3, line=0, adj=0, cex=0.9, "Homogeneous (hgrowth = TRUE)")
#' mtext(side=1, line=0, cex=0.9, adj=0, "Growth rule:")
#' mtext(side=1, line=1, cex=0.9, adj=0, "'dummy_var<=dummy_var[ ] & dummy_var!=1''")
#' legend("topleft", bg = "white", y.intersp= 1.3,
#'        fill = c("#1088a0", "#cfc1af", "#78b2c4", "#cfc1af", "#818792"),
#'        legend = c("RO1", "RO2", "RO1 - growth", "RO2 - growth", "Unclassified cells"))
#' par(oldpar)

cond.4.nofn <- function(attTbl,
                        ngbList,
                        rNumb = FALSE,
                        classVector,
                        class,
                        nbs_of,
                        cond,
                        min.bord = NULL,
                        max.iter = +Inf,
                        peval = 1,
                        directional = FALSE,
                        ovw_class = FALSE,
                        hgrowth = FALSE) {

  # TEST FOR COLUMN CELL IN attTbl
  if (!("Cell" %in% names(attTbl))){
    stop("attribute table mising 'Cell' column. Check ?attTbl")
  }

  # TEST FOR CORRESPONDENCE attTbl, ngbList
  if (length(ngbList) != nrow(attTbl)) {
    stop("ngbList and attTbl shoud have the same length/nrows")
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

  # HANDLE CONDITION STRING
  cond        <- cond.parse(names(attTbl), cond)
  cond_parsed <- cond[[1]]

  ## CONDITIONS TYPE CONTROLS
  v_ab  <- cond[[2]][["v_ab"]]
  v_fc  <- cond[[2]][["v_fc"]]
  v_n   <- cond[[2]][["v_n"]]
  v_nAB <- cond[[2]][["v_nAB"]]

  if(length(v_ab) > 0)  {fa = TRUE} else {fa = FALSE}
  if(length(v_fc) > 0)  {fc = TRUE} else {fc = FALSE}
  if(length(v_n) > 0)   {tn = TRUE} else {tn = FALSE}
  if(length(v_nAB) > 0) {tnAB = TRUE} else {tnAB = FALSE}

  ## OVERWRITE CLASSES
  if(!ovw_class | is.na(ovw_class)) {
    flt <- parse(text = "is.na(classVector[n_ind])")
  } else {
    flt <-
      parse(text = "(classVector[n_ind] != class | is.na(classVector[n_ind]))")
  }

  ###INITIALIZE ALGORITHM #########################################################################
  if(is.null(min.bord)){tb <- F}
  if(!is.null(min.bord)){
    if (min.bord > 1|min.bord < 0){
      stop("min.bord have to be a value between 0 and 1")}
    tb <- T
  }

  if(peval>1|peval<0){stop("peval have to be a value between 0 and 1")}

  # HOMOGENEOUS GROWTH
  if(hgrowth){

    new_cell_id_list       <- list()
    classification_t0_list <- list()

    for(n in seq_along(nbs_of)){

      new_cell_id_list[[n]] <- which(classVector %in% nbs_of[n])
      classification_t0_list[[n]] <- new_cell_id_list[[n]]

    }

    # FOCAL CELLS BY RASTER OBJECT
    new_cell_id       <- new_cell_id_list[[1]]
    classification_t0 <- classification_t0_list[[1]]
    class   <- nbs_of[1]
    obj_ind <- 1

    # NORMAL BEHAVIOUR
  } else {

    # FOCAL CELLS
    new_cell_id       <- which(classVector %in% nbs_of)
    classification_t0 <- new_cell_id

  }

  # check if 'class' in 'nbs_of' (CLASS CONTINUITY)
  if(class %in% nbs_of){
    class_continuity <- TRUE
  } else {
    class_continuity <- FALSE
  }

  ### RUN ALGORITHM #################################################################### while ####
  continue <- TRUE
  itr <- 0

  while (continue & itr < max.iter) {
    continue <- FALSE

    k = 1
    list_new_cell_ind <- list()

    for (c in new_cell_id) {

      n_ind    <- ngbList[[c]]
      n_indAll <- n_ind
      if (length(n_ind) == 0) {
        next
      }

      n_ind <- n_ind[eval(flt)]
      if (length(n_ind) == 0) {
        next
      }


      ### TEST FOR CONDITIONS #################### while//for//if//relative_condition ####
      fct <- seq_along(n_ind)

      # CONSIDERING TEST NEIGHBORHOOD
      if (tn) {

        if (directional) {
          fn_ind <- lapply(ngbList[n_ind], intersect, y = c(c, n_indAll))
        } else {
          fn_ind <- ngbList[n_ind]
        }

        fn_ind0 <- which(lengths(fn_ind) != 0)

        if (length(fn_ind0) == 0) {
          next
        }

        n_ind  <- n_ind[fn_ind0]
        fn_ind <- fn_ind[fn_ind0]

        fct <-
          rep(seq_along(lengths(fn_ind)), lengths(fn_ind)) # NUMBER OF FOCAL NEIGHBORS FOR EACH N_IND

        l_n <-
          lapply(v_n, function(x)
            attTbl[[x]][unlist(fn_ind)])
        names(l_n) <- v_n
      }

      if (tnAB) {
        # IF CONSIDERING ABSOLUTE CONDITION FOCAL NEIGHBORHOOD
        fn_ind <- mapply(c, ngbList[n_ind], n_ind, SIMPLIFY=FALSE)

        if (directional) {
          fn_ind <- lapply(fn_ind, intersect, y = c(n_ind, n_indAll))
        }

        fn_ind0 <- which(lengths(fn_ind) != 0)

        if (length(fn_ind0) == 0) {
          next
        }

        n_ind  <- n_ind[fn_ind0]
        fn_ind <- fn_ind[fn_ind0]

        fct <-
          rep(seq_along(lengths(fn_ind)), lengths(fn_ind)) # NUMBER OF FOCAL NEIGHBORS FOR EACH N_IND

        l_nAB <-
          lapply(v_nAB, function(x)
            attTbl[[x]][unlist(fn_ind)])
        names(l_nAB) <- v_nAB
      }

      # FOCAL CELL CONDITION
      if (fc) {
        l_fc <-
          lapply(v_fc, function(x)
            rep(attTbl[[x]][c], length(fct)))
        names(l_fc) <- v_fc
      }

      # ABSOLUTE CONDITION
      if (fa) {
        l_ab <-
          lapply(v_ab, function(x)
            attTbl[[x]][n_ind][fct])
        names(l_ab) <- v_ab
      }

      ##################################################### while//for//if//relative_condition ####

      ### TEST FOR CELLS MEETING CONDITIONS ########################### while//for//conditions ####
      ev_cond <- eval(cond_parsed)

      if (tn|tnAB) {

        ev_cond <-
          sapply(split(ev_cond, fct), function(x)
            sum(x) / length(x), USE.NAMES = F)
        i       <- which(ev_cond >= peval)
      } else {
        i <- which(ev_cond)
      }

      if (length(i) == 0) {
        next
      }
      ################################################################# while//for//conditions ####


      ### TEST FOR MIN BORDER CONDITION ########################## while//for//if//test_border ####
      n_ind <- n_ind[i]

      if (tb) {

        if(hgrowth){
          nbs_itr <- nbs_of[obj_ind]
        } else {
          nbs_itr <- nbs_of
        }

        test_min_border <- rep(FALSE, length(i))
        for (mb in seq_along(i)) {
          nbg_index <- ngbList[[n_ind[mb]]]
          test_min_border[mb] <-
            sum(classVector[nbg_index] %in% nbs_itr) / 8 >= min.bord

        }

        i <- i[test_min_border]

        if (length(i) == 0) {
          next
        }

        n_ind <- n_ind[i]
      }
      ############################################################ while//for//if//test_border ####

      ### ASSIGN CELLS TO NEW CLASS ####################################### while//for//assign ####
      classVector[n_ind] <- class

      list_new_cell_ind[[k]]  <- n_ind
      k <- k + 1

    } #FOR ENDS

    new_cell_id <-
      setdiff(unlist(list_new_cell_ind), classification_t0)

    ### TEST IF NEW CELLS CHANGED CLASS ############################### while//if//new_cell_id ####
    if (length(new_cell_id) != 0 & class_continuity & !hgrowth) {
      classification_t0 <- c(new_cell_id , classification_t0)
      continue          <- TRUE
      itr <- itr + 1
    }

    # IF HOMOGENEOUS GROWTH
    if(hgrowth){

      if(length(new_cell_id) != 0){
        new_cell_id_list[[obj_ind]]       <- new_cell_id
        classification_t0_list[[obj_ind]] <- classification_t0

        # UPDATE RASTER OBJECT INDEX
        obj_ind <- obj_ind + 1
        if(obj_ind > length(nbs_of)){
          obj_ind <- 1
          itr <- itr + 1
        }

        # MOVE TO THE NEXT RASTER OBJECT
        new_cell_id       <- new_cell_id_list[[obj_ind]]
        classification_t0 <- classification_t0_list[[obj_ind]]

        class <- nbs_of[obj_ind]
        continue <- TRUE
      }

      if(length(new_cell_id) == 0){

        # REMOVE COMPLETE RASTER OBJECT
        new_cell_id_list <- new_cell_id_list[-obj_ind]
        classification_t0_list <- classification_t0_list[-obj_ind]
        nbs_of <- nbs_of[-obj_ind]

        if(length(nbs_of) == 0){next}

        if(obj_ind > length(nbs_of)){
          obj_ind <- 1
          itr <- itr + 1
        }

        # MOVE TO THE NEXT RASTER OBJECT
        new_cell_id       <- new_cell_id_list[[obj_ind]]
        classification_t0 <- classification_t0_list[[obj_ind]]

        class <- nbs_of[obj_ind]
        continue <- TRUE

      }

    }

    ################################################################### while//if//new_cell_id ####

  } #WHILE ENDS

  return(classVector)

}
