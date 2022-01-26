#' Set Anchor Cells from Local Minima or Local Maxima
#'
#' Returns a vector of cell numbers at the locations of local minima or local
#' maxima. These cells can be used as anchor cells in other
#' \code{scapesClassification} functions.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{ngbList}}).
#' @param rNumb logic, \code{ngbList} contain the neighbors index position
#'   in the attribute table (see \code{\link{ngbList}}).
#' @param class numeric, the numeric class to attribute to local minima or
#'   maxima. If \code{NULL} then each local minima or maxima is classified with
#'   a different number starting from one. Buffers have the same classification
#'   number as the local minima or maxima to which they refer.
#' @param cond.filter character string, the conditions to define for what cells
#'   the arguments \code{cond.seed}, \code{cond.growth} and \code{cond.isol}
#'   have to be evaluated. It can be \code{NULL}. Absolute conditions can be
#'   used (see \code{\link{conditions}}).
#' @param cond.seed character string, the conditions to identify local maxima or
#'   minima. The condition should be similar to \code{"variable_x ==
#'   max(variable_x)"} or \code{"variable_x == min(variable_x)"}. It cannot be
#'   \code{NULL}.
#' @param sort.seed character, sort seeds based on column values. If
#'   \code{"max"} seeds are evaluated from the maximum to the minimum. If
#'   \code{"min"} seeds are evaluated from the minimum to the maximum.
#' @param sort.col character, the column name in the \code{attTbl} on which the
#'   sorting in based on.
#' @param cond.growth character string, the conditions to define a buffer around
#'   local maxima or minima. It can be \code{NULL}. Absolute and focal cell
#'   conditions can be used (see \code{\link{conditions}}).
#' @param lag.growth numeric, it defines the lag on which focal cell conditions
#'   in \code{cond.growth} are evaluated (see \code{\link{conditions}}).
#' @param cond.isol character string, the conditions to define how one local
#'   maxima or minima is isolated from another. It can be \code{NULL}. Absolute
#'   and focal cell conditions can be used (see \code{\link{conditions}}).
#' @param lag.isol numeric, it defines the lag on which focal cell conditions in
#'   \code{cond.isol} are evaluated (see \code{\link{conditions}}).
#' @param classVector numeric vector, if any of the condition arguments refers
#'   to a previous classification, then this classification have to be referred
#'   as \code{"classVector"}.
#' @param saveRDS filename, if a file name is provided save the anchor cell
#'   vector as an RDS file.
#' @param overWrite logic, if the RDS names already exist, existing files are
#'   overwritten.
#' @param isolationClass logic, return cells meeting isolation conditions
#'   (\code{isolator cells}) classified as \code{'-1'}.
#' @param silent logic, progress is not printed on the console.
#'
#' @return Classification vector.
#'
#' @details This function implements an algorithm to identify local maxima or
#'   minima based on a set of \code{\link{conditions}}. It is possible to
#'   identify local maxima or minima considering only a subset of cells using
#'   the argument \code{cond.filter}.
#'
#'   The algorithm starts identifying the cell(s) meeting the conditions of the
#'   argument \code{cond.seed}. This cell(s) can be referred to as \code{seed}
#'   cell(s). If a \code{cond.growth} argument is provided, all cells around the
#'   \code{seed} meeting the conditions of \code{cond.growth} are also
#'   classified as \code{seed}. The algorithm will continue looking for new
#'   \code{seed} cells as long as among the neighbors of the newly classified
#'   \code{seed} cells there are cells that meet the conditions of
#'   \code{cond.growth}. In the next step, cells adjacent to \code{seed} cells
#'   meeting the conditions of the argument \code{cond.isol} are classified as
#'   \code{isolator} cells. As for \code{seed} cells, the algorithm will
#'   continue looking for \code{isolator} cells as long as among the neighbors
#'   of the newly classified \code{isolator} cells there are cells that meet the
#'   conditions of \code{cond.isol}. When a cell is classified either as
#'   \code{seed} or as \code{isolator}, then it is ignored in all successive
#'   iterations of the algorithm. The \code{cond.seed}, \code{cond.growth},
#'   \code{cond.isol} process is repeated until all cells are classified either
#'   as \code{seed} or as \code{isolator}.
#'
#'   The evaluation of \code{cond.growth} and \code{cond.isol} follows a node
#'   structure. Let us consider \code{cond.growth} as an example and refer to
#'   'the conditions of the \code{cond.growth} argument' simply as \code{cond}.
#'   At \code{level_0} there is a single node, a cell previously classified as
#'   \code{seed}. Each cell adjacent to the node in \code{level_0} that evaluate
#'   positively to \code{cond} represents a node at \code{level_1}. Suppose that
#'   at \code{level_1} there are one or more nodes. The algorithm then considers
#'   the first of these nodes and looks for neighbor cells evaluating positively
#'   to \code{cond}. If one or more cells do, then the first node of
#'   \code{level_1} is set as \code{referenceNode} for \code{level_2}. The
#'   algorithm then evaluates if one or more cells adjacent to the first node in
#'   \code{level_2} evaluate positively to \code{cond}. If any does, then the
#'   first node of \code{level_2} is set as \code{referenceNode} for
#'   \code{level_3}. The algorithm continues in this way until \code{level_n}
#'   where none of the adjacent cells evaluate positively to \code{cond}. At
#'   this point the algorithm considers cells adjacent to the second node at
#'   \code{level_n-1}. If any cell has a positive evaluation, then second node
#'   in \code{level_n-1} is set as \code{referenceNode} for \code{level_n} and
#'   the algorithm continues as described above. On the contrary, if none of the
#'   cells evaluate positively to \code{cond}, the algorithm will remain at at
#'   \code{level_n-1} and will consider the next node. If none of the nodes at
#'   \code{level_n-1} has any neighboring cell evaluating positively to
#'   \code{cond}, then the algorithm will move down to \code{level_n-2}. The
#'   evaluation of \code{cond.growth} will stop once the algorithm returns to
#'   \code{level_0}.
#'
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
#'
#' @export

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
                        classVector = NULL,
                        saveRDS = NULL,
                        overWrite = FALSE,
                        isolationClass = FALSE,
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

  # Add classVector column if not null
  if(!is.null(classVector)){attTbl$classVector <- classVector}

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

  if(!isolationClass) seedVector[seedVector == -1] <- NA

  if(!is.null(saveRDS)) saveRDS(seedVector, saveRDS)

  if(!silent){cat("\n", paste0("Execution Time: ", round(difftime(Sys.time(), timeStart, units = 'mins') , 2), " minutes" ))}

  return(seedVector)

}
