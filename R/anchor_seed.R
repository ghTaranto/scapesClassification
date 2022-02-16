#' Identify seed cells
#'
#' Returns a vector of cell numbers at the locations of seed cells.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, the list of neighborhoods returned by the function
#'   \code{\link{ngbList}}.
#' @param rNumb logic, the neighborhoods of the argument \code{ngbList} are
#'   identified by cell numbers (\code{rNumb=FALSE}) or by row numbers
#'   (\code{rNumb=TRUE}) (see \code{\link{ngbList}}). It is advised to use row
#'   numbers for large rasters.
#' @param class numeric, the classification number to assign to all cells that
#'   meet the function conditions. If \code{NULL}, a new class numbers is
#'   assigned every time a new seed cell is identified. Growth buffers have the
#'   same classification number as the seed cell to which they refer.
#' @param cond.filter character string, the conditions to define for what cells
#'   the arguments \code{cond.seed}, \code{cond.growth} and \code{cond.isol}
#'   have to be evaluated. It can be \code{NULL}. Absolute conditions can be
#'   used (see \code{\link{conditions}}).
#' @param cond.seed character string, the conditions to identify seed cells.
#'   Absolute conditions can be used (see \code{\link{conditions}}), including
#'   \code{"variable_x == max(variable_x)"} or \code{"variable_x ==
#'   min(variable_x)"} if seed cell correspond to local maxima or minima. It
#'   cannot be \code{NULL}.
#' @param sort.seed character, sort seeds based on column values. If
#'   \code{"max"}, seeds are evaluated from the maximum to the minimum. If
#'   \code{"min"}, seeds are evaluated from the minimum to the maximum.
#' @param sort.col character, the column name in the \code{attTbl} on which the
#'   \code{sort.seed} is based on.
#' @param cond.growth character string, the conditions to define a growth buffer
#'   around seed cells. It can be \code{NULL}. Absolute and focal cell
#'   conditions can be used (see \code{\link{conditions}}).
#' @param lag.growth numeric, it defines the lag on which focal cell conditions
#'   in \code{cond.growth} are evaluated.
#' @param cond.isol character string, the conditions to define an isolation
#'   buffer around seed cells. It can be \code{NULL}. Absolute and focal cell
#'   conditions can be used (see \code{\link{conditions}}).
#' @param lag.isol numeric, it defines the lag on which focal cell conditions in
#'   \code{cond.isol} are evaluated.
#' @param saveRDS filename, if a file name is provided save the class vector as
#'   an RDS file.
#' @param overWrite logic, if the RDS names already exist, existing files are
#'   overwritten.
#' @param isolationBuff logic, return the isolation buffer (class = -1).
#' @param silent logic, progress is not printed on the console.
#'
#' @return Classification vector.
#'
#' @details This function implements an algorithm to identify seed cells, growth
#'   buffers and isolation buffers.
#'
#'   \cr**Condition arguments**
#'
#'   The function takes as inputs four sets of conditions with
#'   \code{cond.growth} and \code{cond.isol} taking into account class
#'   contiguity and continuity (see \code{\link{conditions}}):
#'
#'   1. **\code{cond.filter}**, the conditions to define what cells have to be
#'   evaluated by the function.
#'
#'   2. **\code{cond.seed}**, the conditions to identify, at each iteration, the
#'   seed cell. The seed cell is the cell around which growth and isolation
#'   conditions are applied.
#'
#'   3. **\code{cond.growth}**, the conditions to define a buffer around the
#'   seed cell.
#'
#'   4. **\code{cond.isol}**, the conditions to isolate one seed cell (and its
#'   growth buffer) from another.
#'
#'   \cr**Iterations**
#'
#'   * The argument \code{cond.filter} defines the set of cells to be considered
#'   by the function.
#'
#'   * A seed cell is identified based on \code{cond.seed} and receives a
#'   classification number as specified by the argument \code{class}. If
#'   \code{class=NULL}, then a new class is assigned to every new seed cell.
#'
#'   * Cells continuous and continuous to the seed cell meeting the conditions
#'   specified by \code{cond.growth} are assigned to the same class of the seed
#'   cell (growth buffer).
#'
#'   * Cells continuous and continuous to the seed cell (and its growth buffer)
#'   meeting the conditions specified by \code{cond.isol} are assigned to the
#'   isolation buffer (\code{class = -1}).
#'
#'   * A new seed cell is identified and a new iteration starts. Seed, growth
#'   and isolation cells identified in previous iteration are ignored in
#'   successive iterations.
#'
#'   * The function stops when it cannot identify any new seed cell.
#'
#'   \cr**Node structure**
#'
#'   * The evaluation of \code{cond.growth} and \code{cond.isol} follows a node
#'   structure. Let us consider \code{cond.growth} as an example.
#'
#'   * At \code{level_0} there is a single node, a cell previously classified as
#'   \code{seed}. Every unclassified cell adjacent to the node in \code{level_0}
#'   that evaluates positively to \code{cond.growth} becomes a \code{level_1}
#'   node. The \code{level_0} cell becomes the \code{referenceNode} for
#'   \code{level_1}
#'
#'   * Then, the algorithm considers the first node at \code{level_1}.
#'   Unclassified cells adjacent to the first node evaluating positively to
#'   \code{cond.growth} become \code{level_2} nodes. The first node of
#'   \code{level_1} is set as a \code{referenceNode} for \code{level_2}.
#'
#'   * The algorithm then evaluates if one or more cells adjacent to the first
#'   node in \code{level_2} evaluate positively to \code{cond.growth}. If any
#'   does, then the first node of \code{level_2} is set as \code{referenceNode}
#'   for \code{level_3}.
#'
#'   * The algorithm continues in this way until \code{level_n} where none of
#'   the adjacent cells evaluate positively to \code{cond.growth}.
#'
#'   * At this point the algorithm considers cells adjacent to the second node
#'   at \code{level_n-1}. If any cell has a positive evaluation, then second
#'   node in \code{level_n-1} is set as \code{referenceNode} for \code{level_n}
#'   and the algorithm continues as described above.
#'
#'   * On the contrary, if none of the cells evaluate positively to
#'   \code{cond.growth}, the algorithm will remain at at \code{level_n-1} and
#'   will consider the next node. If none of the nodes at \code{level_n-1} has
#'   any neighboring cell evaluating positively to \code{cond.growth}, then the
#'   algorithm will move down to \code{level_n-2}.
#'
#'   * The evaluation of \code{cond.growth} will stop once the algorithm returns
#'   to \code{level_0}.
#'
#'   **Relative conditions and evaluation lag**
#'   This implementation allows to track back to \code{level_0} all the
#'   connected cells having a positive evaluation. When conditions relative to
#'   the focal cell are considered (see \code{\link{conditions}}), the arguments
#'   \code{lag.growth} and \code{lag.isol} determine to which
#'   \code{referenceNode} newly evaluated cells have to be compared to. A lag of
#'   \code{1} indicates that the conditions at \code{level_n} have to be
#'   evaluated considering the \code{referenceNode} at \code{level_n-1}. A lag
#'   set as \code{Inf} indicates that at any level the conditions are evaluated
#'   considering the initial node cell at \code{level_0}. If the lag is greater
#'   than the level at which the evaluation occurs, then the conditions are
#'   evaluated considering the initial node cell at \code{level_0} (e.g.,
#'   \code{level_1}, \code{lag = 2}).
#'
#' @seealso [conditions()], [attTbl()], [ngbList()]
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
#' # SET FIGURE MARGINS
#' m <- c(4.2, 1, 1, 15)
#'
#' ############################################################################
#' # EXAMPLE PLOTS
#' ############################################################################
#' # 1.
#' as <- anchor.seed(attTbl = at, ngbList = nbs, rNumb = FALSE, class = NULL, silent = TRUE,
#'                   cond.filter = "dummy_var > 1", cond.seed = "dummy_var == max(dummy_var)")
#'
#' plot(cv.2.rast(r,classVector = as), type="classes", mar=m, col=terrain.colors(30)[-1],
#'      asp=NA, axes=FALSE, plg=list(x=1, y=1, cex=.75, title="Classes", y.intersp= 1.4))
#' text(r)
#' lines(r)
#' mtext(side=1, line=0, cex=0.9, font=2, adj=0, "cond.filter:")
#' mtext(side=1, line=0, cex=0.9, adj=1, "dummy_var > 1")
#' mtext(side=1, line=1, cex=0.9, font=2, adj=0, "cond.seed:")
#' mtext(side=1, line=1, cex=0.9, adj=1, "dummy_var == max(dummy_var)")
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth:")
#' mtext(side=1, line=2, cex=0.9, adj=1, "NULL")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "NULL")
#'
#' # 2.
#' as <- anchor.seed(attTbl = at, ngbList = nbs, rNumb = FALSE, class = NULL, silent = TRUE,
#'                   cond.filter = "dummy_var > 1", cond.seed = "dummy_var == max(dummy_var)",
#'                   cond.growth = "dummy_var<dummy_var[] & dummy_var>2")
#'
#' plot(cv.2.rast(r,classVector = as), type="classes", mar=m, asp=NA,
#'      axes=FALSE, plg=list(x=1, y=1, cex=.75, title="Classes", y.intersp= 1.4))
#' text(r)
#' lines(r)
#' mtext(side=1, line=0, cex=0.9, font=2, adj=0, "cond.filter:")
#' mtext(side=1, line=0, cex=0.9, adj=1, "dummy_var > 1")
#' mtext(side=1, line=1, cex=0.9, font=2, adj=0, "cond.seed:")
#' mtext(side=1, line=1, cex=0.9, adj=1, "dummy_var == max(dummy_var)")
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth:")
#' mtext(side=1, line=2, cex=0.9, adj=1, "dummy_var<dummy_var[] & dummy_var>2")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "NULL")
#'
#' # 3.
#' as <- anchor.seed(attTbl = at, ngbList = nbs, rNumb = FALSE, class = NULL, silent = TRUE,
#'                   cond.filter = "dummy_var > 1", cond.seed = "dummy_var == max(dummy_var)",
#'                   cond.growth = "dummy_var<dummy_var[] & dummy_var>2",
#'                   cond.isol = "dummy_var<dummy_var[]")
#'
#' plot(cv.2.rast(r,classVector = as), type="classes", mar=m, asp=NA,
#'      axes=FALSE, plg=list(x=1, y=1, cex=.75, title="Classes", y.intersp= 1.4))
#' text(r)
#' lines(r)
#' mtext(side=1, line=0, cex=0.9, font=2, adj=0, "cond.filter:")
#' mtext(side=1, line=0, cex=0.9, adj=1, "dummy_var > 1")
#' mtext(side=1, line=1, cex=0.9, font=2, adj=0, "cond.seed:")
#' mtext(side=1, line=1, cex=0.9, adj=1, "dummy_var == max(dummy_var)")
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth:")
#' mtext(side=1, line=2, cex=0.9, adj=1, "dummy_var<dummy_var[] & dummy_var>2")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "dummy_var<dummy_var[]")
#'
#' # 4.
#' as <- anchor.seed(attTbl = at, ngbList = nbs, rNumb = FALSE, class = NULL, silent = TRUE,
#'                   cond.seed = "dummy_var >= 5", cond.growth = "dummy_var >= 5")
#'
#' plot(cv.2.rast(r,classVector = as), type="classes", mar=m, asp=NA,
#'      axes=FALSE, plg=list(x=1, y=1, cex=.75, title="Classes", y.intersp= 1.4))
#' text(r)
#' lines(r)
#' mtext(side=1, line=0, cex=0.9, font=2, adj=0, "cond.filter:")
#' mtext(side=1, line=0, cex=0.9, adj=1, "NULL")
#' mtext(side=1, line=1, cex=0.9, font=2, adj=0, "cond.seed:")
#' mtext(side=1, line=1, cex=0.9, adj=1, "dummy_var >= 5")
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth:")
#' mtext(side=1, line=2, cex=0.9, adj=1, "dummy_var >= 5")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "NULL")
#'
#' # 5.
#' as <- anchor.seed(attTbl = at, ngbList = nbs, rNumb = FALSE, class = NULL, silent = TRUE,
#'                   cond.filter = "dummy_var > 1", cond.seed = "dummy_var == max(dummy_var)",
#'                   cond.growth = "dummy_var<dummy_var[]", lag.growth = Inf)
#'
#' plot(cv.2.rast(r,classVector = as), type="classes", mar=m, asp=NA,
#'      axes=FALSE, plg=list(x=1, y=1, cex=.75, title="Classes", y.intersp= 1.4))
#' text(r)
#' lines(r)
#' mtext(side=1, line=0, cex=0.9, font=2, adj=0, "cond.filter:")
#' mtext(side=1, line=0, cex=0.9, adj=1, "dummy_var > 1")
#' mtext(side=1, line=1, cex=0.9, font=2, adj=0, "cond.seed:")
#' mtext(side=1, line=1, cex=0.9, adj=1, "dummy_var == max(dummy_var)")
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth (LAG = INF):")
#' mtext(side=1, line=2, cex=0.9, adj=1, "dummy_var<dummy_var[]")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "NULL")
#'
#' # 6.
#' as <- anchor.seed(attTbl = at, ngbList = nbs, rNumb = FALSE, class = NULL, silent = TRUE,
#'                   cond.filter = "dummy_var > 1", cond.seed = "dummy_var == max(dummy_var)",
#'                   cond.growth = "dummy_var<dummy_var[]", lag.growth = 0)
#'
#' plot(cv.2.rast(r,classVector = as), type="classes", mar=m, asp=NA,
#'      axes=FALSE, plg=list(x=1, y=1, cex=.75, title="Classes", y.intersp= 1.4))
#' text(r)
#' lines(r)
#' mtext(side=1, line=0, cex=0.9, font=2, adj=0, "cond.filter:")
#' mtext(side=1, line=0, cex=0.9, adj=1, "dummy_var > 1")
#' mtext(side=1, line=1, cex=0.9, font=2, adj=0, "cond.seed:")
#' mtext(side=1, line=1, cex=0.9, adj=1, "dummy_var == max(dummy_var)")
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth (LAG = 0):")
#' mtext(side=1, line=2, cex=0.9, adj=1, "dummy_var<dummy_var[]")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "NULL")

anchor.seed <- function(attTbl,
                        ngbList,
                        rNumb = FALSE,
                        class = NULL,
                        cond.filter = NULL,
                        cond.seed,
                        sort.seed = NULL,
                        sort.col = NULL,
                        cond.growth  = NULL,
                        lag.growth = Inf,
                        cond.isol = NULL,
                        lag.isol = 1,
                        saveRDS = NULL,
                        overWrite = FALSE,
                        isolationBuff = FALSE,
                        silent = FALSE)
{

  cat("\n")
  timeStart <- Sys.time()

  # TEST IF FILENAMES ALREADY EXIST
  if(!overWrite){
    # RDS
    if(!is.null(saveRDS)){
      if(file.exists(saveRDS)) stop("RDS filename exists; use a different name")
    }
  }


  # TEST THAT CLASS IS NOT EQUAL TO -1 (ISOLATION CLASS)
  if(!is.null(class)){
    if(class == -1) stop("change class number; -1 is used for the isolation class")
  }


  # TEST FOR COLUMN CELL IN attTbl
  if(!("Cell" %in% names(attTbl))){
    stop(
      "the attribute table should have one column named 'Cell' with cell numbers
                                       that indicate the position of each row in the original Raster object"
    )
  }


  # TEST FOR CORRESPONDENCE attTbl, ngbList
  if(length(ngbList) != nrow(attTbl)){
    stop("ngbList and attTbl shoud have the same length/nrows")
  }

  # CONVERT ngbList FORM CELL IDS TO CELL INDECES
  if(!rNumb){
    fct     <- rep(seq_along(lengths(ngbList)), lengths(ngbList))
    ngbList <- match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- split(ngbList, fct)

    rm(fct, no_nas)
  }

  if(!is.null(cond.filter)){

    cond.filter <- paste0("(", cond.filter, ") & is.na(seedVector)" )

  } else {

    cond.filter <- "is.na(seedVector)"

  }

  c_all <- c(cond.seed, cond.growth, cond.isol)

  ind_nms <- c(TRUE, FALSE, FALSE)
  if(!is.null(cond.growth)){ind_nms[2] <- TRUE}
  if(!is.null(cond.isol))  {ind_nms[3] <- TRUE}

  c_nms <- c("cond.seed", "cond.growth", "cond.isol")[ind_nms]

  cond_list <- as.list(c_all)
  names(cond_list) <- c_nms

  v_list <- list()

  for(cond in 1:length(c_all)){

    nm <- c_nms[cond]
    v_list[[ nm ]] <- list()

    v_list[[ nm ]][["v_ab"]] <-
      names(attTbl)[stringr::str_detect(c_all[cond], paste0("\\b", names(attTbl), "(?!\\[|\\{)", "\\b"))]

    v_list[[ nm ]][["v_fc"]] <-
      names(attTbl)[stringr::str_detect(c_all[cond], paste0(names(attTbl), "\\[\\]"))]

  }

  for(cond in 1:length(c_all)){

    l <- v_list[[cond]]

    for(v in l$v_ab){cond_list[[cond]] <-
      stringr::str_replace_all(cond_list[[cond]], paste0("\\b", v, "(?!\\[|\\{)", "\\b"), paste0("l_ab$", v))}

    for(v in l$v_fc){cond_list[[cond]] <-
      stringr::str_replace_all(cond_list[[cond]], paste0(v, "\\[\\]" ), paste0("l_fc$", v))}

  }


  # CONDITION FILTER (TO APPLY GLOBALLY)
  cond_filter <- cond.filter

  vf <- names(attTbl)[stringr::str_detect(cond_filter, paste0("\\b", names(attTbl), "\\b"))]
  for(v in vf){cond_filter <-
    stringr::str_replace_all(cond_filter, v, paste0("attTbl$", v))
  }

  len_vf <- length(vf) != 0


  # CONDITION FILTER (TO APPLY LOCALLY)
  cf <- stringr::str_replace_all(cond_filter, "attTbl", "l_ab")
  cf <- stringr::str_replace_all(cf, "seedVector(?!\\[)", paste0("seedVector", "\\[nind\\]"))

  # CONDITION FILTER GROWTH, CAN OVERWERITE ISOLATION CLASS
  cfg<- stringr::str_replace_all(cf, "is.na\\(seedVector\\[nind\\]\\)",
                                 "\\(is.na\\(seedVector\\[nind\\]\\)\\|seedVector\\[nind\\] == -1\\)")


  # PARSE CONDITIONS
  cond_filter <- parse(text = cond_filter)
  cf          <- parse(text = cf)
  cfg         <- parse(text = cfg)
  cond_list   <- lapply(cond_list, function(x){parse(text = x)})


  # TEST CONDITIONS?
  gc_true <- !is.null(cond.growth)
  ic_true <- !is.null(cond.isol)

  ### INITIALIZE VARIABLES ###########################################################################
  seedVector <- rep(as.numeric(NA), NROW(attTbl))

  flt_ok <- which(eval(cond_filter))
  N      <- length(flt_ok)

  cnumb <- 0

  if(length(flt_ok) == 0){stop("\n No cell meeting filter conditions")}

  seeds <- TRUE
  while(length(flt_ok) > 0 & seeds){
    evaluate <- TRUE #test if any filter cond == T within growth and isolation
    cnumb <- cnumb + 1

    v    <- v_list$cond.seed
    l_ab <- lapply( as.list(v$v_ab), function(x) attTbl[[x]][flt_ok] )
    names(l_ab) <- v$v_ab

    # FOCAL CELL INDEX
    fc_ind <- which( eval(cond_list[["cond.seed"]]) )

    # STOP IF NO SEED CELL
    if(length(fc_ind) == 0){
      seeds <- FALSE

      if(cnumb == 1){stop("No cell meeting seed conditions")}

      next

    }

    fc_ind <- flt_ok[fc_ind]

    if(!is.null(sort.seed)){

      srt <- attTbl[[sort.col]][fc_ind]

      if(sort.seed == "max"){
        fc_ind <- fc_ind[match(max(srt), srt)][1]
      } else if(sort.seed == "min") {
        fc_ind <- fc_ind[match(min(srt), srt)][1]
      }

    } else {

      fc_ind <- fc_ind[1]

    }

    # CLASSIFY SEED
    seedVector[fc_ind] <- cnumb

    ## INITIALIZE ALGORITHM
    classification_t1 <- fc_ind

    if(gc_true){

      lag_fc0 <- fc_ind
      lag_fc  <- list()

      cgL     <- lag_fc0

      list_new_cell_ind <- list()

      v    <- v_list$cond.growth
      cond <- cond_list$cond.growth

    }

    k <- 1
    while(k > 0 & gc_true){
      ### DEFINE FOCAL CELL NEIGHBORHOOD AND TEST CONDITION ###############################################
      nind  <- ngbList[[ fc_ind ]]

      # EVALUATE FILTER
      if(len_vf){

        l_ab        <- lapply( as.list(vf), function(x) attTbl[[x]][nind] )
        names(l_ab) <- vf

      }

      i <- which(eval(cf))

      if(length(i) > 0){
        evaluate <- TRUE
        nind <- nind[i]
      }

      # EVALUATE GROWTH
      if(evaluate){

        fct <- 1:length(nind)

        l_fc <- lapply( as.list(v$v_fc), function(x) rep(attTbl[[x]][cgL], length(nind)) )
        names(l_fc) <- v$v_fc

        l_ab <- lapply( as.list(v$v_ab), function(x) attTbl[[x]][nind] )
        names(l_ab) <- v$v_ab

        i <- which(eval(cond))

      }

      ### AT LEAST ONE NEIGHBOUR CELL MEETING CONDITION ###################################################
      if(length(i) > 0){

        # CLASSIFY CELL AND SET NEXT NODE NEIGHBORHOOD
        nind <- nind[i]
        seedVector[ nind ] <- cnumb
        list_new_cell_ind[[ k ]] <- nind

        # CELLS TO BE USED TO TEST FOR ISOLATION CONDITIONS
        classification_t1 <- c(classification_t1, nind)

        # UPDATE FOCAL CELL & REMOVE IT FROM 'list_new_cell_ind'
        fc_ind                   <- list_new_cell_ind[[ k ]][1]
        list_new_cell_ind[[ k ]] <- list_new_cell_ind[[ k ]][-1]

        # DEFINE NODE CONDITION
        lag_fc[[k]] <- fc_ind

        # UPDATE LAG CONDITION TO TEST NEXT NODE CELLS
        if(k - lag.growth < 1){

          cgL <- lag_fc0

        } else {

          cgL <- lag_fc[[k - lag.growth]]

        }

        # SET NODE LEVEL
        k <- k + 1

        ### NO CELL MEETING ISOLATION CONDITION ###########################################################
      } else {

        # CHECK IF ANY NODE HAS UNEVALUATED CELLS
        node2eval <- TRUE
        while( node2eval ){
          k = k - 1

          if(k < 1) {
            node2eval <- FALSE
            next
          }

          if(length(list_new_cell_ind[[ k ]]) != 0) {

            # UPDATE FOCAL CELL FROM LOWER LEVEL NODES
            fc_ind                   <- list_new_cell_ind[[ k ]][1]
            list_new_cell_ind[[ k ]] <- list_new_cell_ind[[ k ]][-1]

            node2eval <- FALSE

            # UPDATE LOWER LOWER NODE CONDITION
            lag_fc[[k]] <- fc_ind

            # UPDATE LAG CONDITION TO TEST NEXT NODE CELLS
            if(k - lag.growth < 1){

              cgL <- lag_fc0

            } else {

              cgL <- lag_fc[[k - lag.growth]]

            }

            # SET NODE LEVEL
            k <- k + 1

          } # if(length(list_new_cell_ind[[ k ]]) != 0)

        } # while( node2eval )

      } # else of if(length(i) != 0)

      evaluate <- FALSE

    } # while(k > 0 & ic_true)


    ### APPLY ISOLATION CONDITION ###########################################################################
    if(ic_true){
      for(fc_ind in classification_t1){

        ### INITIALIZE ISOLATION ALGORITHM ####################################################################
        lag_fc0 <- fc_ind
        lag_fc  <- list()

        ciL     <- lag_fc0

        # new_lag_cond      <- T
        list_new_cell_ind <- list()

        v    <- v_list$cond.isol
        cond <- cond_list$cond.isol

        k = 1
        evaluate <- FALSE

        ### START ISOLATION CONDITION #########################################################################
        while(k > 0){

          ### DEFINE FOCAL CELL NEIGHBORHOOD AND TEST CONDITION ###############################################
          nind  <- ngbList[[ fc_ind ]]

          # EVALUATE FILTER
          if(len_vf){

            l_ab        <- lapply( as.list(vf), function(x) attTbl[[x]][nind] )
            names(l_ab) <- vf

          }

          i <- which(eval(cf))

          if(length(i) > 0){
            evaluate = TRUE
            nind <- nind[i]
          }

          # EVALUATE ISOLATION
          if(evaluate){

            fct <- 1:length(nind)

            l_fc <- lapply( as.list(v$v_fc), function(x) rep(attTbl[[x]][ciL], length(nind)) )
            names(l_fc) <- v$v_fc

            l_ab <- lapply( as.list(v$v_ab), function(x) attTbl[[x]][nind] )
            names(l_ab) <- v$v_ab

            i <- which(eval(cond))

          }

          ### AT LEAST ONE NEIGHBOUR CELL MEETING ISOLATION CONDITION #########################################
          if(length(i) != 0){

            # CLASSIFY CELL AND SET NEXT NODE NEIGHBORHOOD
            nind <- nind[i]
            seedVector[ nind ] <- -1
            list_new_cell_ind[[ k ]] <- nind

            # UPDATE FOCAL CELL
            fc_ind                   <- list_new_cell_ind[[ k ]][1]
            list_new_cell_ind[[ k ]] <- list_new_cell_ind[[ k ]][-1]

            # DEFINE NODE CONDITION
            lag_fc[[k]] <- fc_ind

            # UPDATE LAG CONDITION TO TEST NEXT NODE CELLS
            if(k - lag.isol < 1){

              ciL <- lag_fc0

            } else {

              ciL <- lag_fc[[k - lag.isol]]

            }

            # SET NODE LEVEL
            k <- k + 1

            ### NO CELL MEETING ISOLATION CONDITION ###########################################################
          } else {

            # CHECK IF ANY NODE HAS UNEVALUATED CELLS
            node2eval <- TRUE
            while( node2eval ){
              k = k - 1

              if(k < 1) {
                node2eval <- FALSE
                next
              }

              if(length(list_new_cell_ind[[ k ]]) != 0) {

                # UPDATE FOCAL CELL FROM LOWER LEVEL NODES
                fc_ind                   <- list_new_cell_ind[[ k ]][1]
                list_new_cell_ind[[ k ]] <- list_new_cell_ind[[ k ]][-1]

                node2eval <- FALSE

                # UPDATE LOWER LEVEL NODE CONDITION
                lag_fc[[k]] <- fc_ind

                # UPDATE LAG CONDITION TO TEST NEXT NODE CELLS
                if(k - lag.isol < 1){

                  ciL <- lag_fc0

                } else {

                  ciL <- lag_fc[[k - lag.isol]]

                }

                # SET NODE LEVEL
                k <- k + 1

              } # if(length(list_new_cell_ind[[ k ]]) != 0)

            } # while( node2eval )

          } # else of if(length(i) != 0)

          evaluate <- FALSE

        } # while(k > 0 & ic_true)

      } #for(i_ind in classification_t0)

    } #ic_true


    ### REINITIALIZE ALGORITHM ##############################################################################
    flt_ok <- which(eval(cond_filter))

    if(!silent){
      n      <- length(flt_ok)
      # p      <- round((1-n/N) * 100, 2)

      dtime <- round(difftime(Sys.time(), timeStart, units = "mins"), 2)

      # cat("\r", paste0(cnumb, ") ", p, "%"," complete, (", N-n, "/", N, ") cells classified, elapsed time ",
      #                  dtime, " mins"))

      cat("\r", paste0("Seeds Identified: ", cnumb, ", elapsed time ",
                       dtime, " mins"))
    }

  }#while

  ### FINALIZE FUNCTION #####################################################################################
  if(!is.null(class))seedVector[which( seedVector > 0 )] <- class

  if(!isolationBuff) seedVector[seedVector == -1] <- NA

  if(!is.null(saveRDS)) saveRDS(seedVector, saveRDS)

  if(!silent){cat("\n", paste0("Execution Time: ", round(difftime(Sys.time(), timeStart, units = 'mins') , 2), " minutes" ))}

  return(seedVector)

}
