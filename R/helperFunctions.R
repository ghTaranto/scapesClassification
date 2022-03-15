#' scapesClassification conditions
#'
#' Check for spelling and syntax errors in conditions (\code{cond} argument) and
#' detect the type of conditions being used.
#'
#' @param names_attTbl character vector, the column (i.e. variable) names of the
#'   attribute table returned by the function \code{\link{attTbl}}.
#' @param cond character string, the condition string used by the \code{cond.*}
#'   functions to classify raster cells (see \code{\link{conditions}}).
#' @param silent logic, when true, the function returns only error messages.
#'
#' @encoding UTF-8
#'
#' @return An error message if the function finds spelling or syntax errors or a
#'   string with the types of rules detected in the condition string.
#'
#' @seealso [cond.4.all()], [cond.4.nofn()], [anchor.seed()], [cond.reclass()],
#'   [conditions()]
#'
#' @details **Conditions (alias classification rules)**
#'
#'   * Classification rules evaluate either to true or false and determine what
#'   raster cells have to be classified.
#'
#'   * Conditions are passed to \code{scapesClassification} functions as a
#'   single character string. They can consist of combination of variables names
#'   (as named in the attribute table, see \code{\link{attTbl}}), arithmetic
#'   \code{(+|-|*|/|^|%%|%/%)}, relational \code{(>|<|>=|<=|==|!=|%/%)} and
#'   logic operators \code{(&||)} and base R functions (e.g.,
#'   \code{abs(variable_name)}).
#'
#'   * All variables included in an attribute table (see \code{\link{attTbl}})
#'   can be included in a condition string by name (e.g., var name =
#'   \code{"dummy_var"}; condition = \code{"dummy_var > 1"}).\cr
#'
#'   **Class vectors**
#'
#'   * Class vectors map raster cells to numeric classes.
#'
#'   * The n^th^ element of a class vector stores the class of the raster cell
#'   stored in the n^th^ row of the corresponding attribute table (see
#'   \code{\link{attTbl}}).
#'
#'   * Class vectors can serve also as a function input. As inputs, they provide
#'   information about the groups of cells that have already been classified.
#'
#'   * Every time a class vector is provided as a function input, it is
#'   _updated_ by assigning a numeric class to _unclassified cells_ that meet
#'   function conditions.
#'
#'   * Unclassified cells are represented as NA values.\cr
#'
#'   **Rule evaluation: Global evaluation**
#'
#'   * Classification rules are applied to all unclassified raster
#'   cells.
#'
#'   * Function using _global evaluation_: \code{\link{cond.4.all}}.
#'
#'   **Rule evaluation: Focal evaluation**
#'
#'   * Classification rules are applied only to raster cells at specific
#'   locations and are based on class (dis)contiguity and class continuity.
#'
#'   * __Class contiguity:__ \cr classification rules are applied only to raster
#'   cells contiguous to focal cells (identified in the \code{cond.*} functions
#'   by the argument \code{nbs_of}).
#'
#'   * __Class continuity:__ \cr join into the same class cells that respect the
#'   same rules and that are connected to the same focal cells. This means that,
#'   at each iteration, newly classified cells become focal cells and conditions
#'   are tested in their neighborhood.
#'
#'   * Function using _focal evaluation_: \code{\link{anchor.seed}},
#'   \code{\link{cond.4.nofn}}, \code{\link{cond.reclass}},
#'   \code{\link{reclass.nbs}} and \code{\link{classify.all}}.\cr
#'
#'   **Focal evaluation: Definitions**
#'
#'   * __Cell neighborhood:__ a cell with coordinates \code{(x, y)} has 8
#'   neighbors with coordinates: \code{(x±1, y)},  \code{(x, y±1)} and
#'   \code{(x±1, y±1)}. Cells on the edge of a raster have less than 8
#'   neighbors. See \code{\link{ngbList}}.
#'
#'   * __Focal cell:__ cells whose neighbors are evaluated against the
#'   classification rule(s). In the classification functions _focal cells_ are
#'   identified by the argument \code{nbs_of}.
#'
#'   * __Test cell:__ the cell in the neighborhood of the focal cell that is
#'   being tested. At turns all cells in the neighborhood of a focal cell are
#'   tested against the classification rule(s).
#'
#'   * __Directional neighborhood:__ it consists of the intersection between the
#'   focal and the test cell neighborhoods.\cr
#'
#'   **Condition type: Absolute conditions**
#'
#'   __1) Absolute test cell condition:__ compares cell values against a
#'   threshold value.
#'
#'   * This type of condition applies to all functions with a \code{conditions}
#'   argument.
#'
#'   * In global evaluations all cells meeting absolute conditions receive a
#'   classification number. In focal evaluations all \code{test cells} meeting
#'   absolute conditions receive a classification number.
#'
#'   * _Examples of valid conditions:_ \code{"variable_A > 1 & variable_B !=
#'   0"}; \code{"(variable_A^2 < 50 & variable_B == 0) | abs(variable_C) > 50"}.\cr
#'   _Functions:_ \code{\link{anchor.seed}}, \code{\link{cond.4.all}},
#'   \code{\link{cond.4.nofn}} and \code{\link{cond.reclass}}.\cr
#'
#'   __2) Absolute neighborhood condition:__ compares the values of the
#'   \code{test cell} and of its \code{neighborhood} against a threshold value.
#'
#'   * An absolute neighborhood condition is identified by a variable name
#'   followed by curly brackets (e.g., \code{"variable_name{}"}).
#'
#'   * A maximum of 9 evaluations are performed for each test cell (the test
#'   cell itself and the cells of its neighborhood are compared against a
#'   threshold value).
#'
#'   * Test cells receive a classification number if the rule is true for at
#'   least as many evaluations as the ones specified by the argument
#'   \code{peval}. The argument \code{peval} ranges from 0 to 1. When 9
#'   evaluations are performed, \code{peval = 1} means that all \code{9}
#'   evaluations have to be true; \code{peval = 0.5} means that at least
#'   \code{4.5} (rounded to \code{5}) evaluations have to be true.
#'
#'   * Only one neighborhood rule is allowed for each condition string (e.g., it
#'   is not possible to have a condition string like \code{"variable_A{} > 0 &
#'   variable_B{} > 1"}).
#'
#'   * The function \code{\link{cond.4.nofn}} can consider a \code{directional
#'   neighborhood} instead of the test cell neighborhood by setting the argument
#'   \code{directional = TRUE}.
#'
#'   * _Example of valid conditions:_ \code{"variable_A{} > 1 & abs(variable_B)
#'   != 0"}. \cr _Functions:_ \code{\link{cond.4.nofn}} and
#'   \code{\link{cond.reclass}}.\cr
#'
#'   **Condition type: Relative conditions**
#'
#'   __1) Relative focal cell condition:__ compares the \code{test cell} value
#'   against the \code{focal cell} value.
#'
#'   * A relative focal cell condition is identified by a variable name followed
#'   by square brackets (e.g., \code{"variable_name[]"}).
#'
#'   * Rules are defined repeating twice the same variable name, once with
#'   square brackets and once without. Square brackets indicate the focal cell
#'   value. As an example, the rule \code{"dummy_var < dummy_var[]"} compares
#'   the value of the the test cell (\code{"dummy_var"}) against the value of
#'   the focal cell (\code{"dummy_var[]"}).
#'
#'   * Test cells are classified if the rule is true.
#'
#'   * _Examples of valid conditions:_ \code{"variable_A > variable_A[]"};
#'   \code{"(variable_A > variable_A[] & variable_B{} < 10) | variable_C > 1"}.
#'   Note that the last example is a combination of absolute and focal cell
#'   conditions. \cr _Functions:_ \code{\link{anchor.seed}} and
#'   \code{\link{cond.4.nofn}}.
#'
#'   __2) Relative neighborhood condition:__ compares the values of the \code{test
#'   cell} against the values of the \code{test cell neighborhood}.
#'
#'   * A relative neighborhood condition is identified by a variable name
#'   followed by curly brackets (e.g., \code{"variable_name{}"}).
#'
#'   * Rules are defined repeating twice the same variable name, once with curly
#'   brackets and once without. Curly brackets indicate the test cell
#'   neighborhood. As an example, the rule \code{'dummy_var < dummy_var{}'}
#'   compares the value of the the test cell (\code{dummy_var}) against the
#'   values of cells included in the test cell neighborhood
#'   (\code{dummy_var{}}).
#'
#'   * A maximum of 8 evaluations are performed for each test cell (the test
#'   cell is compared against each cell included in its neighborhood).
#'
#'   * Test cells receive a classification number if the rule is true for at
#'   least as many evaluations as the ones specified by the argument
#'   \code{peval}. The argument \code{peval} ranges from 0 to 1. When 8
#'   evaluations are performed, \code{peval = 1} means that all \code{8}
#'   evaluations have to be true; \code{peval = 0.5} means that at least
#'   \code{4} evaluations have to be true.
#'
#'   * Only one neighborhood rule is allowed for each condition string (e.g., it
#'   is not possible to have a condition string like \code{"variable_A{} > 0 &
#'   variable_B{} > variable_B"}).
#'
#'   * The function \code{\link{cond.4.nofn}} can consider a \code{directional
#'   neighborhood} instead of the test cell neighborhood by setting the argument
#'   \code{directional = TRUE}.
#'
#'   * _Example of valid conditions:_ \code{"variable_A > variable_A{}"};
#'   \code{"(variable_A > variable_A{} & variable_B != variable_B[]) |
#'   variable_C > 1"}. Note that the last example is a combination of absolute
#'   and relative conditions. \cr _Functions:_ \code{\link{cond.4.nofn}} and
#'   \code{\link{cond.reclass}}.
#'
#' @seealso [anchor.seed()], [attTbl()], [cond.4.all()], [cond.4.nofn()],
#'   [cond.reclass()], [classify.all()]
#'
#' @export
#' @examples
#' # LOAD LIBRARIES
#' library(scapesClassification)
#' library(terra)
#'
#' ################################################################################
#' # TYPES OF CONDITIONS
#' ################################################################################
#'
#' # As an example consider an attribute with the following columns
#' names_attTbl <- c("bathymetry", "slope")
#'
#' # And the following conditions
#' cond <- "bathymetry>10"
#' conditions(names_attTbl, cond)
#'
#' cond <- "bathymetry[]>bathymetry | abs(slope{}) < 5"
#' conditions(names_attTbl, cond)
#'
#' cond <- "bathymetry[]>bathymetry | abs(slope{}) < slope"
#' conditions(names_attTbl, cond)
#'
#' ################################################################################
#' # FOCAL EVALUATION DEFINITIONS
#' ################################################################################
#'
#' # CELL NUMBERS OF A DUMMY RASTER (7X7)
#' r   <- terra::rast(matrix(1:49, nrow = 7, byrow = TRUE), extent=c(0,7,0,7))
#' nbs <- ngbList(r)
#'
#' # CLASS VECTOR WITH ONE TEST AND ONE FOCAL CELL
#' cv <- as.numeric(rep(NA, 49))
#' cv[c(32, 25)] <- c(1, 2) # tc (class 1), fc (class 2)
#' r_cv <- cv.2.rast(r, classVector = cv)
#'
#' # POLYGONS REPRESENTING NEIGHBORHOODS
#' fcn <- rbind(c(2,5), c(5,5), c(5,2), c(2,2), c(2,5))
#' fcn <- terra::vect(fcn, type="polygons")
#'
#' tcn <- rbind(c(2,4), c(5,4), c(5,1), c(2,1), c(2,4))
#' tcn <- terra::vect(tcn, type="polygons")
#'
#' # PLOT - FOCAL EVALUATION DEFINITIONS
#' m <- c(3.5, 8, 1.2, 8)
#' plot(r_cv, type = "class", asp = NA, legend = FALSE, axes = FALSE, mar = m,
#'      col = c("goldenrod3","#1088a0"), colNA= "grey90")
#' text(r)
#' mtext(side=3, line=0, adj=0, cex=1, font=2, "FOCAL EVALUATION")
#' mtext(side=1, line=0, adj=0, cex=0.9,
#'       "Focal cell neighborhood: 17, 18, 19, 24, 26, 31, 32, 33")
#' mtext(side=1, line=1, cex=0.9, adj=0,
#'       "Test cell neighborhood: 24, 25, 26, 31, 33, 38, 39, 40")
#' mtext(side=1, line=2, cex=0.9, adj=0,
#'      "Directional neighborhood: 24, 25, 26, 31, 32, 33")
#' lines(fcn, col="#1088a0", lwd=2)
#' lines(tcn, col="#cfad8999", lwd=2)
#' legend("bottomleft", ncol = 1, bg = "white", y.intersp= 1.3,
#'        legend = c("Focal cell", "Test cell"), fill = c("#1088a0", "goldenrod3"))

conditions <- function(names_attTbl,
                       cond,
                       silent = FALSE) {

  if(any(grepl("seedVector|classVector", names_attTbl))){
    nvn <- names_attTbl[grepl("seedVector|classVector", names_attTbl)]
    stop(paste0(paste0("'",unique(nvn),"'",collapse = " and "),
                ": not valid column name(s). Please change column name(s) in the 'attTbl'."))
  }

  cond <- gsub("\\&\\&", "\\&", cond)
  cond <- gsub("\\|\\|", "\\|", cond)
  cond <- gsub(" ", "", cond)

  # CHECK FOR DOUBLE TEST CELL NEIGHBORHOODS NOT ALLOWED
  doubleTNcond <- unlist(strsplit(cond,"\\&|\\|"))
  doubleTNcond <- doubleTNcond[doubleTNcond != ""]
  doubleTNcond <- sapply(as.list(doubleTNcond), function(x)grepl("\\{\\}", x, perl = TRUE))

  if(sum(doubleTNcond, na.rm = TRUE) > 1){stop("Only one 'neighborhood' condition allowed.")}

  # CHECK FOR SPELLING ERRORS
  match_pattern <-
    paste0(paste(names_attTbl, collapse = "|"), "|([0-9]+)")

  tc_spl <-
    unlist(strsplit(cond,
                    ",| |\\(|\\+|\\*|>|<|=|-|\\)|\\||\\&|\t|\\{|\\||\\!|\\^|\\%|\\{|\\}|\\%in\\%|\n"
    ))
  tc_spl <- tc_spl[tc_spl != ""]
  tc_spl <- tc_spl[which(!(tc_spl %in% c("between",ls("package:base"))))]
  tc_spl <- tc_spl[which(!(tc_spl %in% c("T", "TRUE", "F", "FALSE")))]
  tc_spl <- tc_spl[!grepl(match_pattern, tc_spl)]

  if (length(tc_spl > 0)) {
    stop(
      paste0(
        "Check: '",
        paste(unique(tc_spl), collapse = ", "),
        "' -- Possible spelling error."
      )
    )
  }

  # CHECK FOR SYNTAX ERRORS
  c_eval <- gsub(paste(c("classVector", names_attTbl), collapse = "\\b|\\b"), "1", cond)
  c_eval <- gsub("\\[|\\]|\\{|\\}", "", c_eval)
  c_eval <- try(eval(parse(text = c_eval)), silent = TRUE)

  if(!is.logical(c_eval)){stop(paste0("Check conditions, possible syntax error."))}

  # PRINT RESULTS IF NO ERRORS
  if(!silent){

    c_split <- unlist(strsplit(cond, "&|\\|"))
    counts  <-  sapply(c_split, function(x){length(gregexpr(paste0(names_attTbl,collapse = "|"), x)[[1]])})

    abs_cond <- grepl(paste0(names_attTbl, ("(?!\\[|\\{)"), collapse = "|"), c_split, perl=TRUE)
    abs_cond <- any(counts == 1 & abs_cond)

    fc_cond_r <- grepl("\\[\\]", c_split, perl = TRUE)
    fc_cond_a <- fc_cond_r
    fc_cond_r <- any(counts > 1 & fc_cond_r)
    fc_cond_a <- any(counts == 1 & fc_cond_a)

    tn_cond_r <- grepl("\\{\\}", c_split, perl = TRUE)
    tn_cond_a <- tn_cond_r
    tn_cond_r <- any(counts > 1 & tn_cond_r)
    tn_cond_a <- any(counts == 1 & tn_cond_a)

    verify_types <- c(abs_cond,
                      fc_cond_a,
                      fc_cond_r,
                      tn_cond_a,
                      tn_cond_r)

    cond_types <- c("'Absolute test cell'",
                    "'Absolute focal cell'",
                    "'Relative focal cell'",
                    "'Absolute neighborhood'",
                    "'Relative neighborhood'")

    detc_cond <- cond_types[verify_types]
    if(length(detc_cond) == 0){detc_cond <- "unknown"}

    print(paste0(paste0(detc_cond, collapse = " AND "), " condition type(s) detected."))
  }
}

#' Parse conditions
#'
#' Parse the condition string so that it can be evaluated by the `cond.*`
#' functions. Intended for internal use only.
#'
#' @param names_attTbl character vector, the column (i.e. variable) names of the
#'   attribute table returned by the function \code{\link{attTbl}}.
#' @param cond character string, the condition string used by the \code{cond.*}
#'   functions to classify raster cells (see \code{\link{conditions}}).
#'
#' @return The function returns a two-element list. The first element contains
#'   the parsed conditions to be evaluated by the \code{cond.*} functions. The
#'   second element defines the condition type each variable refers to.
#'
#' @seealso [cond.4.all()], [cond.4.nofn()], [anchor.seed()], [cond.reclass()],
#'   [conditions()]
#'
#' @export

cond.parse <- function(names_attTbl, cond){

  scapesClassification::conditions(names_attTbl, cond, silent = TRUE)

  cond <- gsub("\\&\\&", "\\&", cond)
  cond <- gsub("\\|\\|", "\\|", cond)
  cond <- gsub(" ", "", cond)

  cVect <- unlist(strsplit(cond,"\\&|\\|"))
  cVect <- cVect[cVect != ""]

  vList <- list()

  for(i in seq_along(cVect)){

    #ab
    ip <- sapply(paste0("\\b",names_attTbl,"(?!\\[|\\{)","\\b"),
                 grepl, x=cVect[i], perl = TRUE)
    v_ab <-names_attTbl[ip]

    #fc `[]`
    ip <- sapply(paste0("\\b",names_attTbl, "\\[\\]"),
                 grepl, x=cVect[i], perl = TRUE)
    v_fc <-names_attTbl[ip]

    #nbors `{}`
    ip <- sapply(paste0("\\b",names_attTbl, "\\{\\}"),
                 grepl, x=cVect[i], perl = TRUE)
    v_n <-names_attTbl[ip]

    v_nAB <- character()

    for (v in v_ab) {
      cVect[i] <- gsub(paste0("\\b", v, "(?!\\[|\\{)", "\\b"),
                       paste0("l_ab$", v), cVect[i], perl = TRUE)
    }

    for (v in v_fc) {
      cVect[i] <- gsub(paste0("\\b", v, "\\[\\]"),
                       paste0("l_fc$", v), cVect[i], perl = TRUE)
    }

    if(length(v_ab) == 0 & length(v_n) > 0){

      v_nAB <- v_n
      v_n   <- character()

      for (v in v_nAB) {
        cVect[i] <- gsub(paste0("\\b", v, "\\{\\}"),
                         paste0("l_nAB$", v), cVect[i], perl = TRUE)
      }

    } else {

      for (v in v_n) {
        cVect[i] <- gsub(paste0("\\b", v, "\\{\\}"),
                         paste0("l_n$", v), cVect[i], perl = TRUE)
      }

    }

    vList[["v_ab"]]  <- c(vList[["v_ab"]], v_ab)
    vList[["v_fc"]]  <- c(vList[["v_fc"]], v_fc)
    vList[["v_n"]]   <- c(vList[["v_n"]], v_n)
    vList[["v_nAB"]] <- c(vList[["v_nAB"]], v_nAB)

  }

  # RECONSTRUCT CONDITION STRINGS AND VARIABLES
  ands <- gregexpr("\\&", cond, perl = TRUE)[[1]]
  ors  <- gregexpr("\\|", cond, perl = TRUE)[[1]]

  lg_op <- data.frame(pos = c(ands, ors),
                      op = c(rep("&", length(ands)),
                             rep("|", length(ors))))

  lg_op <- lg_op[order(lg_op$pos),]
  lg_op <- lg_op[lg_op$pos > 0, ]

  cond_parsed <- character()
  for(i in seq_along(cVect)){
    cond_parsed <- paste0(cond_parsed, cVect[i])

    if(i <= nrow(lg_op)){
      cond_parsed <- paste0(cond_parsed, lg_op$op[i])
    }
  }

  l <- list()
  l[["cond_par"]] <- parse(text=cond_parsed)
  l[["vList"]]    <- vList

  return(l)
}


#' Class vector to raster
#'
#' Transform a class vector or a generic vector into a raster.
#'
#' @param r raster object.
#' @param classVector numeric vector, the values to be assigned to the cell
#'   numbers indicated by \code{index}.
#' @param index numeric vector, the cell numbers of the argument \code{r} to
#'   which assign the values of the argument \code{classVector}. If \code{NULL},
#'   the column \code{Cell} of the attribute table \code{attTbl(r)} is used (see
#'   \code{\link{attTbl}}).
#' @param plot logic, plot the raster.
#' @param type character, type of map/legend. One of "continuous", "classes", or
#'   "interval".
#' @param writeRaster filename, if a raster name is provided save the raster to
#'   a file.
#' @param overWrite logic, if the raster names already exist, the existing file
#'   is overwritten.
#'
#' @details The arguments \code{index} and \code{vector} need to have the same
#'   length. The function assign the values of \code{vector} at the positions of
#'   \code{index} to an empty raster having the same spatial properties of the
#'   raster \code{r}.
#'
#' @return A class vector or a generic vector transformed into a raster.
#'
#' @export
#' @examples
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
#' # Compute an example class vector
#' cv <- cond.4.all(attTbl = at, cond = "dummy_var > 1", class = 1)
#'
#' # Class vector to raster
#' cv.2.rast(r, cv, plot = TRUE)
#' text(r) # add raster values

cv.2.rast <- function(r, classVector, index = NULL, plot = FALSE, type = "classes",
                      writeRaster = NULL, overWrite = FALSE){

  if(is.null(index)){
    index <- attTbl(r)
    index <- index$Cell
  }

  r2 <- r[[1]]
  r2[] <- NA
  r2[index] <- classVector

  if (plot){
    terra::plot(r2, type = type)
  }

  if (!is.null(writeRaster)){
    terra::writeRaster(r2, writeRaster, overwrite = overWrite)
  }

  return(r2)

}
