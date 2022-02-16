#' Identify seed cells
#'
#' Returns a vector of cell numbers at the locations of seed cells and growth
#' buffers.
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
#'   meet the function conditions. If \code{NULL}, a new class number is
#'   assigned every time a new seed cell is identified. Growth buffers have the
#'   same classification number as the seed cell to which they refer.
#' @param cond.filter character string, defines for what cells the arguments
#'   \code{cond.seed}, \code{cond.growth} and \code{cond.isol} have to be
#'   evaluated. It can be \code{NULL}. Absolute conditions can be used (see
#'   \code{\link{conditions}}).
#' @param cond.seed character string, the conditions to identify seed cells.
#'   Absolute conditions can be used (see \code{\link{conditions}}). It cannot
#'   be \code{NULL}.
#' @param cond.growth character string, the conditions to define a growth buffer
#'   around seed cells. It can be \code{NULL}. Absolute and focal cell
#' conditions can be used (see \code{\link{conditions}}).
#' @param lag.growth numeric, it defines the lag on which focal cell conditions
#'   in \code{cond.growth} are evaluated.
#' @param cond.isol character string, the conditions to define an isolation
#'   buffer around seed cells and growth buffers. It can be \code{NULL}.
#'   Absolute and focal cell conditions can be used (see
#'   \code{\link{conditions}}).
#' @param lag.isol numeric, it defines the lag on which focal cell conditions in
#'   \code{cond.isol} are evaluated.
#' @param sort.col character, the column name in the \code{attTbl} on which the
#'   \code{sort.seed} is based on. It determines in what order seed buffers are
#'   computed.
#' @param sort.seed character, the order seed buffers are computed is based on
#'   the value seed cells have in the \code{sort.col}. If
#'   \code{sort.seed="max"}, buffers are computed from the seed cell having the
#'   maximum value to the seed cell having the minimum value. If
#'   \code{sort.seed="min"}, buffers are computed in the opposite order.
#' @param saveRDS filename, if a file name is provided save the class vector as
#'   an RDS file.
#' @param overWrite logic, if the RDS names already exist, existing files are
#'   overwritten.
#' @param isol.buff logic, return the isolation buffer (class = -999).
#' @param silent logic, progress is not printed on the console.
#'
#' @return Class vector. See \code{\link{conditions}} for more details about
#'   class vectors.
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
#'   The argument \code{cond.filter} defines the set of cells to be considered
#'   by the function.
#'
#'   1. A seed cell is identified based on \code{cond.seed} and receives a
#'   classification number as specified by the argument \code{class}. If
#'   \code{class=NULL}, then a new class is assigned to every new seed cell.
#'
#'   2. Cells contiguous and continuous to the seed cell meeting the conditions
#'   specified by \code{cond.growth} are assigned to the same class of the seed
#'   cell (growth buffer).
#'
#'   3. Cells contiguous and continuous to the seed cell (and its growth buffer)
#'   meeting the conditions specified by \code{cond.isol} are assigned to the
#'   isolation buffer (\code{class = -999}).
#'
#'   4. A new seed cell is identified and a new iteration starts. Seed, growth
#'   and isolation cells identified in previous iteration are ignored in
#'   successive iterations.
#'
#'   The function stops when it cannot identify any new seed cell.
#'
#'   \cr**Relative focal cell conditions and evaluation lag**
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
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth (lag.growth=Inf):")
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
#' mtext(side=1, line=2, cex=0.9, font=2, adj=0, "cond.growth (lag.growth=0):")
#' mtext(side=1, line=2, cex=0.9, adj=1, "dummy_var<dummy_var[]")
#' mtext(side=1, line=3, cex=0.9, font=2, adj=0, "cond.isol:")
#' mtext(side=1, line=3, cex=0.9, adj=1, "NULL")

anchor.seed <- function(attTbl,
                        ngbList,
                        rNumb = FALSE,
                        class = NULL,
                        cond.filter = NULL,
                        cond.seed,
                        cond.growth  = NULL,
                        lag.growth = Inf,
                        cond.isol = NULL,
                        lag.isol = 1,
                        sort.col = NULL,
                        sort.seed = "max",
                        saveRDS = NULL,
                        overWrite = FALSE,
                        isol.buff = FALSE,
                        silent = FALSE)
{

  if(!silent){cat("\n"); timeStart<-Sys.time()}

  # TEST ARGUMENTS ################################################################################

  # Filename already exists?
  if(!overWrite){
    if(!is.null(saveRDS)){
      if(file.exists(saveRDS))stop("RDS filename exists; use a different name")}}

  # Class -999 reserved to isolation buffer
  if(!is.null(class)){
    if(class == -999) stop("-999 isolation buffer number. Change `class` number")}

  # Column cell in `attTbl`
  if(!("Cell" %in% names(attTbl))){
    stop("attribute table mising 'Cell' column. Check ?attTbl")}

  # Length `ngbList` == nrows `attTbl`
  if(length(ngbList) != nrow(attTbl)){
    stop("ngbList and attTbl shoud have the same length/nrows")
  }

  # Sort.seed can be only 'max' or 'min'
  if(!sort.seed=="max"|sort.seed=="min"){
    stop("sort.seed can be either 'max' or 'min'")
  }

  # Lag on the evaluation of lag.growth and lag.isol
  if(!(lag.growth %in% c(0,Inf)|lag.isol %in% c(0,Inf))){
    stop("'lag.growth' and 'lag.isol' can be either '0' or 'Inf'")
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

  ### CONDITION STRINGS & CONDITION CONTROLS ######################################################
  v_list      <- list()
  cond_parsed <- list()

  # Filter
  if(!is.null(cond.filter)){
    cond_filter <- cond_parse(names(attTbl), cond.filter)

    ctype <- names(cond_filter[[2]])[lengths(cond_filter[[2]])>0]
    if(!all("v_ab" == ctype)){
      stop("Only absolute test cell conditions allowed for `cond.filter` (see ?conditions)")
    }

    cond_filter <- gsub("\\bl_ab\\b", "attTbl", cond_filter[[1]], perl = TRUE)
    cond_filter <- paste(cond_filter, "& is.na(seedVector)")
    cond_filter <- parse(text = cond_filter)
  }

  if(is.null(cond.filter)){
    cond_filter <- parse(text = "is.na(seedVector)")
  }

  # Seed
  cond <- cond_parse(names(attTbl), cond.seed)

  ctype <- names(cond[[2]])[lengths(cond[[2]])>0]
  if(!all("v_ab" == ctype)){
    stop("Only absolute test cell conditions allowed for `cond.seed`
         (see ?conditions)")
  }

  cond_parsed[["cond.seed"]] <- cond[[1]]
  v_list[["cond.seed"]]      <- cond[[2]]

  # Growth
  gc_true <- FALSE
  if(!is.null(cond.growth)){
    if(!is.null(cond.filter)){cond.growth <- paste0(cond.growth,"&",cond.filter)}
    cond <- cond_parse(names(attTbl), cond.growth)
    gc_true <- TRUE

    ctype <- names(cond[[2]])[lengths(cond[[2]])>0]
    if(any("v_nAB" == ctype | "v_n" == ctype)){
      stop("Only absolute test cell and focal cell conditions allowed for `cond.growth`
           (see ?conditions)")
    }

    cond_parsed[["cond.growth"]] <- cond[[1]]
    v_list[["cond.growth"]]      <- cond[[2]]

    if (length(cond[[2]][["v_ab"]]) > 0)  {Gfa = TRUE} else {Gfa = FALSE}
    if (length(cond[[2]][["v_fc"]]) > 0)  {Gfc = TRUE} else {Gfc = FALSE}
  }

  # Isol
  ic_true <- FALSE
  if(!is.null(cond.isol)){
    if(!is.null(cond.filter)){cond.isol <- paste0(cond.isol,"&",cond.filter)}
    cond <- cond_parse(names(attTbl), cond.isol)
    ic_true <- TRUE

    ctype <- names(cond[[2]])[lengths(cond[[2]])>0]
    if(any("v_nAB" == ctype | "v_n" == ctype)){
      stop("Only absolute test cell and focal cell conditions allowed for `cond.isol`
           (see ?conditions)")
    }

    cond_parsed[["cond.isol"]] <- cond[[1]]
    v_list[["cond.isol"]]      <- cond[[2]]

    if (length(cond[[2]][["v_ab"]]) > 0)  {Ifa = TRUE} else {Ifa = FALSE}
    if (length(cond[[2]][["v_fc"]]) > 0)  {Ifc = TRUE} else {Ifc = FALSE}
  }

  ### INITIALIZE VARIABLES ########################################################################
  seedVector <- rep(as.numeric(NA), NROW(attTbl))

  flt_ok <- which(eval(cond_filter))
  N      <- length(flt_ok)

  cnumb <- 0

  if(length(flt_ok) == 0){stop("\n No cell meeting filter conditions")}

  ### START ALGORITHM #############################################################################
  seeds <- TRUE
  while(length(flt_ok) > 0 & seeds){
    cnumb <- cnumb + 1

    v    <- v_list[["cond.seed"]]
    l_ab <- lapply( v$v_ab, function(x) attTbl[[x]][flt_ok] )
    names(l_ab) <- v$v_ab

    # FOCAL CELL INDEX
    fc_ind <- which( eval(cond_parsed[["cond.seed"]]) )

    # STOP IF NO SEED CELL
    if(length(fc_ind) == 0){
      if(cnumb == 1){stop("No cell meeting seed conditions")}

      seeds <- FALSE
      next
    }

    # REINDEX BASED ON FILTER
    fc_ind <- flt_ok[fc_ind]

    # SORT SEEDS BASED ON sort.seed
    if(!is.null(sort.col)){

      srt <- attTbl[[sort.col]][fc_ind]

      if(sort.seed == "max"){
        fc_ind <- fc_ind[match(max(srt), srt)][1]
      } else if(sort.seed == "min") {
        fc_ind <- fc_ind[match(min(srt), srt)][1]
      }
    }

    if(is.null(sort.col)){fc_ind <- fc_ind[1]}

    # CLASSIFY SEED
    seedVector[fc_ind] <- cnumb
    LS <- fc_ind

    # GROWTH ######################################################################################
    if(gc_true){

      list_new_cell <- integer()
      v    <- v_list[["cond.growth"]]
      cond <- cond_parsed[["cond.growth"]]

      # Isolation conditions
      classification_t1 <- fc_ind

    }

    evaluate <- FALSE
    continue <- TRUE
    while(continue & gc_true){

      continue <- FALSE
      while(length(fc_ind) > 0){

        # TEST CELLS
        nind  <- ngbList[[ fc_ind[1] ]]

        # EVALUATE FILTER
        i <- which(is.na(seedVector[nind]))

        if(length(i) > 0){
          evaluate <- TRUE
          nind <- nind[i]
        }

        # EVALUATE GROWTH
        if(evaluate){

          if(Gfc){ # Focal cell conditions

            if(is.infinite(lag.growth)){
              Lag <- LS
            } else {
              Lag <- fc_ind[1]
            }

            l_fc <- lapply( v$v_fc, function(x) rep(attTbl[[x]][Lag], length(nind)) )
            names(l_fc) <- v$v_fc
          }

          if(Gfa){ # Test cell conditions
            l_ab <- lapply( as.list(v$v_ab), function(x) attTbl[[x]][nind] )
            names(l_ab) <- v$v_ab
          }

          i <- which(eval(cond))

        }

        # AT LEAST ONE TEST CELL MEETS CONDITION
        if(length(i) > 0){

          # CLASSIFY CELL
          nind <- nind[i]
          seedVector[ nind ] <- cnumb
          list_new_cell <- c(list_new_cell, nind)

        }

        # REMOVE EVALUATED TEST CELL
        fc_ind <- fc_ind[-1]
        evaluate <- FALSE

      }

      # SET THE NEW GROUP OF TEST CELLS
      if(length(list_new_cell) > 0){
        fc_ind <- list_new_cell
        list_new_cell <- integer()
        continue <- TRUE

        if(ic_true){classification_t1<-c(classification_t1, fc_ind)}
      }
    }

    # ISOL ################################################################################################
    if(ic_true){

      list_new_cell <- integer()
      v    <- v_list[["cond.isol"]]
      cond <- cond_parsed[["cond.isol"]]

      if(gc_true){fc_ind <- classification_t1}

    }

    evaluate <- FALSE
    continue <- TRUE
    while(continue & ic_true){

      continue <- FALSE
      while(length(fc_ind) > 0){

        # TEST CELLS
        nind  <- ngbList[[ fc_ind[1] ]]

        # EVALUATE FILTER
        i <- which(is.na(seedVector[nind]))

        if(length(i) > 0){
          evaluate <- TRUE
          nind <- nind[i]
        }

        # EVALUATE GROWTH
        if(evaluate){

          if(Ifc){ # Focal cell conditions

            if(is.infinite(lag.isol)){
              Lag <- LS
            } else {
              Lag <- fc_ind[1]
            }

            l_fc <- lapply( v$v_fc, function(x) rep(attTbl[[x]][Lag], length(nind)) )
            names(l_fc) <- v$v_fc
          }

          if(Ifa){ # Test cell conditions
            l_ab <- lapply( as.list(v$v_ab), function(x) attTbl[[x]][nind] )
            names(l_ab) <- v$v_ab
          }

          i <- which(eval(cond)) # Evaluate conditions

        }

        # AT LEAST ONE TEST CELL MEETS CONDITION
        if(length(i) > 0){

          # CLASSIFY CELL
          nind <- nind[i]
          seedVector[ nind ] <- -999
          list_new_cell <- c(list_new_cell, nind)
        }

        # REMOVE EVALUATED TEST CELL
        fc_ind <- fc_ind[-1]
        evaluate <- FALSE

      }

      # SET THE NEW GROUP OF TEST CELLS
      if(length(list_new_cell) > 0){
        fc_ind <- list_new_cell
        list_new_cell <- integer()
        continue <- TRUE
      }
    }

    ### REINITIALIZE ALGORITHM ##############################################################################
    flt_ok <- which(eval(cond_filter))

    if(!silent){
      n      <- length(flt_ok)
      dtime  <- round(difftime(Sys.time(), timeStart, units = "mins"), 2)
      cat("\r", paste0("Evaluated cells: ", n, "/", N, "; Elapsed time: ", dtime, " mins"))
    }

  }#while (seed.cond)

  ### FINALIZE FUNCTION #####################################################################################
  if(!is.null(class))seedVector[which( seedVector > 0 )] <- class
  if(!isol.buff) seedVector[seedVector == -999] <- NA
  if(!is.null(saveRDS)) saveRDS(seedVector, saveRDS)
  if(!silent){cat("\n", paste0("Seed cells: ", cnumb,
                               "; Execution Time: ", round(difftime(Sys.time(), timeStart, units = 'mins') , 2), " mins" ))}

  return(seedVector)
}
