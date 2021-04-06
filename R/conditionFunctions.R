#' Test Conditions for Neighbors and Neighbors of Neighbors
#'
#' Evaluate conditions for cells neighboring specific classes and classify them
#' if conditions are \code{TRUE}.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{ngbList}}).
#' @param nbsIndex logic, \code{ngbList} contain the neighbors index position
#'   in the attribute table (see \code{\link{ngbList}}).
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified.
#' @param nbs_of numeric or numeric vector, indicates the neighbors of the
#'   classes that have to be evaluated. If it includes the class of the argument
#'   \code{class}, the function will evaluate the conditions also for the
#'   neighbors of the neighbors.
#' @param conditions character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{class}. Absolute, focal cell
#'   and focal neighborhood conditions can be used. Condition strings can
#'   include only one neighborhood condition (see \code{\link{conditions}}).
#' @param class numeric, the numeric class to attribute to the cells meeting
#'   conditions.
#' @param min.border numeric value between 0 and 1. It indicates the minimum
#'   percentage of cells adjacent to the cell in evaluation that have to belong
#'   to one of the classes in \code{nbs_of}. If one, then all neighboring cells
#'   have to belong to one of the classes in \code{nbs_of} for the cell in
#'   evaluation being classified as \code{class}.
#' @param overwrite_class logic, reclassify cells that have already been
#'   classified.
#' @param max.iter integer, the maximum number of iterations.
#' @param fn_perc numeric value between 0 and 1. If neighborhood conditions are
#'   considered, it determines the percentage of cells in the neighborhood for
#'   which the conditions have to be true in order to classify the cell being
#' evaluated as indicted by the argument \code{class} (see
#' \code{\link{conditions}}).
#' @param directional logic, only directional neighbors are considered to test
#'   neighborhood conditions. The argument \code{fn_perc} will also consider
#'   only directional neighbors (see \code{\link{conditions}}).
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function.
#'
#' @details The function evaluates the conditions of the argument
#'   \code{\link{conditions}} for all cells neighboring cells classified as one
#'   of the classes included in \code{nbs_of}. If the argument \code{nbs_of}
#'   includes the class of the argument \code{class}, then at each iteration the
#'   function will evaluate if among the neighbors of the newly classified cells
#'   there are cells meeting conditions and it will classify them accordingly.
#'
#' @seealso [conditions()], [attTbl()], [ngbList()]
#'
#' @export

cond.4.nofn <- function(attTbl,
                        ngbList,
                        nbsIndex = FALSE,
                        classVector,
                        nbs_of,
                        conditions,
                        class,
                        min.border = NULL,
                        overwrite_class = FALSE,
                        max.iter = +Inf,
                        fn_perc = 1,
                        directional = T) {

  # TEST FOR COLUMN CELL IN attTbl
  if (!("Cell" %in% names(attTbl))){
    stop("attribute table mising 'Cell' column. Check ?attTbl.")
  }


  # TEST FOR CORRESPONDENCE attTbl, ngbList
  if (length(ngbList) != nrow(attTbl)) {
    stop("ngbList and attTbl shoud have the same length/nrows")
  }


  ## HANDLE CONDITION STRING
  # REMOVE SPACES AND DOUBLES (&& ||)
  conditions <- stringr::str_replace_all(conditions, "\\&\\&", "\\&")
  conditions <- stringr::str_replace_all(conditions, "\\|\\|", "\\|")
  conditions <- stringr::str_replace_all(conditions, " ", "")

  # TEST FOR CONDITION INTEGRITY
  conditions(names(attTbl), conditions, silent = TRUE)

  # DECONSTRUCT SRING TO DETECT CONDITION TYPES
  cVect <- unlist(strsplit(conditions,"\\&|\\|"))
  cVect <- cVect[cVect != ""]

  vList <- list()

  for(x in seq_along(cVect)){

    v_ab <-
      names(attTbl)[stringr::str_detect(cVect[x], paste0(names(attTbl), "(?!\\[|\\{)"))]
    v_fc <-
      names(attTbl)[stringr::str_detect(cVect[x], paste0(names(attTbl), "\\[\\]"))]
    v_fn <-
      names(attTbl)[stringr::str_detect(cVect[x], paste0(names(attTbl), "\\{\\}"))]
    v_fnAB <- character()

    for (v in v_ab) {
      cVect[x] <-
        stringr::str_replace_all(cVect[x], paste0(v, "(?!\\[|\\{)"), paste0("l_ab$", v))
    }

    for (v in v_fc) {
      cVect[x] <-
        stringr::str_replace_all(cVect[x], paste0(v, "\\[\\]"), paste0("l_fc$", v))
    }

    if(length(v_ab) == 0 & length(v_fn) > 0){

      cnd    <- "l_fnAB$"
      v_fnAB <- v_fn
      v_fn   <- character()

      for (v in v_fnAB) {
        cVect[x] <-
          stringr::str_replace_all(cVect[x], paste0(v, "\\{\\}"), paste0("l_fnAB$", v))
      }

    } else {

      for (v in v_fn) {
        cVect[x] <-
          stringr::str_replace_all(cVect[x], paste0(v, "\\{\\}"), paste0("l_fn$", v))
      }

    }

    vList[["v_ab"]]   <- c(vList[["v_ab"]], v_ab)
    vList[["v_fc"]]   <- c(vList[["v_fc"]], v_fc)
    vList[["v_fn"]]   <- c(vList[["v_fn"]], v_fn)
    vList[["v_fnAB"]] <- c(vList[["v_fnAB"]], v_fnAB)

  }

  # RECONSTRUCT CONDITION STRINGS AND VARIABLES
  ands <- stringr::str_locate_all(conditions, "\\&")[[1]]
  ors  <- stringr::str_locate_all(conditions, c("\\|"))[[1]]

  lg_op      <- as.data.frame( rbind(ands, ors) )
  lg_op$type <- c(rep("&", nrow(ands)), rep("|", nrow(ors)))
  lg_op      <- lg_op[order(lg_op$start),]

  cond_parsed <- character()
  for(x in seq_along(cVect)){

    cond_parsed <- paste0(cond_parsed, cVect[x])

    if(x <= nrow(lg_op)){
      cond_parsed <- paste0(cond_parsed, lg_op$type[x])
    }

  }

  cond_parsed <- parse(text=cond_parsed)


  ## OVERWRITE CLASSES
  if (!overwrite_class | is.na(overwrite_class)) {
    flt <- parse(text = "is.na(classVector[n_ind])")
  } else {
    flt <-
      parse(text = "(classVector[n_ind] != class | is.na(classVector[n_ind]))")
  }


  ## CONVERT NBS FORM CELL IDS TO CELL INDECES
  if(!nbsIndex){
    fct     <- rep(seq_along(lengths(ngbList)), lengths(ngbList))
    ngbList <- match(unlist(ngbList), attTbl$Cell)
    no_nas  <- !is.na(ngbList)
    ngbList <- ngbList[no_nas]
    fct     <- fct[no_nas]

    ngbList <- split(ngbList, fct)

    rm(fct, no_nas)
  }

  ## CONDITIONS TYPE CONTROLS
  v_ab   <- vList$v_ab
  v_fc   <- vList$v_fc
  v_fn   <- vList$v_fn
  v_fnAB <- vList$v_fnAB

  fa   <- FALSE
  fc   <- FALSE
  fn   <- FALSE
  fnAB <- FALSE

  if (length(v_ab) > 0) {fa = TRUE}
  if (length(v_fc) > 0) {fc = TRUE}
  if (length(v_fn) > 0) {fn = TRUE}
  if (length(v_fnAB) > 0) {fnAB = TRUE}

  ###INITIALIZE ALGORITHM #########################################################################
  itr <- 0
  if (is.na(max.iter) | is.null(max.iter)) {
    max.iter <- +Inf
  }

  if (is.null(min.border)) {
    tb <- F
  } else {
    if (min.border > 1 |
        min.border < 0)
      stop("min.border have to be a numeric value between 0 and 1")
    tb <- T
  }

  if (is.na(fn_perc)) {
    fn_perc <- 1
  }

  continue          <- TRUE
  new_cell_id       <- which(classVector %in% nbs_of)
  classification_t0 <- new_cell_id
  conditions0       <- conditions


  ### RUN ALGORITHM #################################################################### while ####
  while (continue & itr <= max.iter) {
    continue <- FALSE
    itr      <- itr + 1

    k = 1
    list_new_cell_ind <- list()

    for (c in new_cell_id) {
      #

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
      fct <- 1:length(n_ind)
      if (fn) {
        # IF CONSIDERING FOCAL NEIGHBORHOOD
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

        l_fn <-
          lapply(as.list(v_fn), function(x)
            attTbl[[x]][unlist(fn_ind)])
        names(l_fn) <- v_fn
      }

      if (fnAB) {
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

        l_fnAB <-
          lapply(as.list(v_fnAB), function(x)
            attTbl[[x]][unlist(fn_ind)])
        names(l_fnAB) <- v_fnAB
      }

      # FOCAL CELL CONDITION
      if (fc) {
        l_fc <-
          lapply(as.list(v_fc), function(x)
            rep(attTbl[[x]][c], length(fct)))
        names(l_fc) <- v_fc
      }

      # ABSOLUTE CONDITION
      if (fa) {
        l_ab <-
          lapply(as.list(v_ab), function(x)
            attTbl[[x]][n_ind][fct])
        names(l_ab) <- v_ab
      }

      ##################################################### while//for//if//relative_condition ####

      ### TEST FOR CELLS MEETING CONDITIONS ########################### while//for//conditions ####
      ev_cond <- eval(cond_parsed)

      if (fn|fnAB) {

      # if (fn|fnAB) {
        ev_cond <-
          sapply(split(ev_cond, fct), function(x)
            sum(x) / length(x), USE.NAMES = F)
        i       <- which(ev_cond >= fn_perc)
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
        test_min_border <- rep(FALSE, length(i))
        for (mb in 1:length(i)) {
          nbg_index <- ngbList[[n_ind[mb]]]
          test_min_border[mb] <-
            sum(classVector[nbg_index] %in% nbs_of) / 8 >= min.border

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
    if (length(new_cell_id) != 0) {
      classification_t0 <- c(new_cell_id , classification_t0)
      continue          <- TRUE
    }

    ################################################################### while//if//new_cell_id ####

  } #WHILE ENDS

  return(classVector)

}


#' Test Conditions and Reclassify
#'
#' Evaluate conditions for cells of a class and reclassify them if conditions
#' are true.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param ngbList list, it has to contain the list of 8-neighbors of each cell
#'   in \code{attTbl$Cell} (see \code{\link{ngbList}}). If conditions do not
#'   include focal neighborohood conditions this argument can be \code{NULL}
#'   (see \code{\link{conditions}}).
#' @param nbsIndex logic, \code{ngbList} contain the neighbors index position
#'   in the attribute table (see \code{\link{ngbList}}).
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified.
#' @param conditions character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{reclass}. Absolute and focal
#'   neighborhood conditions can be used (see \code{\link{conditions}}).
#' @param class numeric or numeric vector, indicates the classes for which
#'   conditions have to be evaluated.
#' @param reclass numeric, the numeric class to attribute to the cells meeting
#'   conditions.
#' @param fn_perc numeric value between 0 and 1. If a focal neighborhood
#'   condition is considered, it determines the percentage of cells in the focal
#'   neighborhood for which the conditions have to be \code{TRUE} in order to
#'   classify the cell being evaluated as indicted by the argument
#'   \code{reclass} (see \code{\link{conditions}}).
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function.
#'
#' @details The function evaluates the conditions of the argument
#'   \code{\link{conditions}} for all cells in the classes of the argument
#'   \code{class}. The cells for which conditions are true are reclassified as
#'   indicted by the argument \code{reclass}.
#'
#' @seealso [conditions()], [attTbl()], [ngbList()]
#'
#' @export

cond.reclass <- function(attTbl,
                         ngbList = NULL,
                         nbsIndex = FALSE,
                         classVector,
                         conditions,
                         class,
                         reclass,
                         fn_perc = 1) {

  ### PARSE CONDITIONS
  v_ab <-
    names(attTbl)[stringr::str_detect(conditions, paste0(names(attTbl), "(?!\\[|\\{)"))]
  v_fn <-
    names(attTbl)[stringr::str_detect(conditions, paste0(names(attTbl), "\\{\\}"))]

  for (v in v_ab) {
    conditions <-
      stringr::str_replace_all(conditions, paste0(v, "(?!\\[|\\{)"), paste0("l_ab$", v))
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
  if(!nbsIndex){
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

    if (fn) {
      rc <- sum(ev_cond) / length(n_ind) >= fn_perc
    } else {
      rc <- ev_cond
    }

    if (rc) {
      classVector[c] <- reclass
    }

  } #FOR ENDS

  return(classVector)

}


#' Test Conditions for All Cells
#'
#' Evaluate conditions for cells that have not already been classified and
#' classify them if conditions are true.
#'
#' @param attTbl data.frame, the attribute table returned by the function
#'   \code{\link{attTbl}}.
#' @param conditions character string, the conditions a cell have to meet to be
#'   classified as indicated by the argument \code{class}. Absolute conditions
#'   can be used (see \code{\link{conditions}}).
#' @param classVector numeric vector, defines the cells in the attribute table
#'   that have already been classified.
#' @param class numeric, the numeric class to attribute to the cells meeting
#'   conditions.
#' @param overwrite_class logic, reclassify cells that have already been
#'   classified.
#'
#' @return Update \code{classVector} with the new cells that were classified by
#'   the function. If no \code{classVector} was provided, the function return a
#'   new classification vector.
#'
#' @details The function evaluates the conditions of the argument
#'   \code{\link{conditions}} for all cells that have not already been
#'   classified. The cells for which \code{\link{conditions}} are \code{TRUE}
#'   are classified as indicted by the argument \code{class}.
#'
#' @seealso [conditions()], [attTbl()]
#'
#' @export

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
    names(attTbl)[stringr::str_detect(conditions, names(attTbl))]
  for (v in v_ab) {
    conditions <-
      stringr::str_replace_all(conditions, v, paste0("attTbl$", v))
  }

  classVector[which(eval(parse(text = conditions)))] <- class

  return(classVector)

}



