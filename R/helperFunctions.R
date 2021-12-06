#' Check Conditions
#'
#' Check for spelling and syntax errors in \code{conditions} and evaluate the
#' type of conditions being used.
#'
#' @param names_attTbl character vector, the column (i.e. variable) names of the
#'   attribute table returned by the function \code{\link{attTbl}}.
#' @param cond character string, a set of conditions used to classify raster
#'   cells (see \code{\link{cond.4.nofn}}, \code{\link{cond.reclass}},
#'   \code{\link{cond.4.all}} and \code{\link{anchor.seed}}).
#' @param silent logic, returns only error messages.
#'
#' @return An error message as character string if spelling or syntax errors
#'   were found or a message detailing the types of conditions that were
#'   found.
#'
#' @details \cr **Conditions** \cr\cr Conditions are passed to
#'   \code{scapesClassification} functions as a single character string. They
#'   can consist of combination of arithmetic \code{(+|-|*|/|^|%%|%/%)},
#'   relational \code{(>|<|>=|<=|==|!=|%/%)} and logic operators \code{(&||)},
#'   base R functions (e.g., \code{abs(variable_name)}), variables names (i.e.,
#'   \code{names(attTbl)}) and class vectors (referred to as
#'   \code{"classVector"}). \cr\cr A combination of absolute and relative
#'   conditions can be used, but _**only one neighborhood condition per
#'   string**_ is allowed. In order to avoid errors, it should always be
#'   possible to evaluate the \code{conditions} as \code{TRUE} or \code{FALSE}
#'   when substituting variable names with numeric values. All variables of
#'   interest should be addressed in the \code{conditions} argument by name.
#'   \cr\cr\cr **Evaluation of conditions** \cr\cr The way conditions are
#'   evaluated depends on both the type of condition and the type of function
#'   used. When class contiguity is not being considered (i.e., functions
#'   without the argument \code{nbs_of}) absolute conditions can be used. Their
#'   evaluation is vectorized and exclude cells that have already been
#'   classified (i.e., \code{which(!is.na(classVector))}), unless the argument
#'   \code{overwrite_class} is \code{TRUE}. The function \code{cond.reclass} can
#'   evaluate absolute focal neighborhood conditions, but cannot evaluate
#'   directional neighborhood conditions (see functions \code{\link{cond.4.all}}
#'   and \code{\link{cond.reclass}}). \cr\cr When class contiguity is being
#'   considered (i.e., functions with the argument \code{nbs_of} and
#'   \code{anchor.seed}), then both absolute and relative conditions can be
#'   used. Their evaluation is iterative. At each iteration, conditions are
#'   evaluated for the 8-neighbors of one of the cells classified as
#'   \code{nbs_of}. Neighbors that have already been classified are excluded
#'   from the evaluation (i.e., \code{which(!is.na(classVector))}), unless the
#'   argument \code{overwrite_class} is \code{TRUE} (see functions
#'   \code{\link{anchor.seed}} and \code{\link{cond.4.nofn}}). Some
#'   classification functions do not have a \code{condition} argument. These
#'   functions consider only class contiguity to classify raster cells (see
#'   functions \code{\link{anchor.cell}}, \code{\link{anchor.svo}},
#'   \code{\link{reclass.nbs}} and \code{\link{classify.all}}). \cr\cr
#'   Hereinafter the following definitions are adopted (see also the example
#'   section): \itemize{ \item{**focal cell:** }{cell included in the classes of
#'   the argument \code{nbs_of} whose neighbors are in evaluation;} \item{**nbs
#'   cell**: }{one of the neighbors in evaluation;} \item{**focal
#'   neighborhood:** }{when _absolute focal neighborhood conditions_ are used,
#'   it includes \code{nbs} and its 8 neighbors; when _relative focal
#'   neighborhood conditions_ are used, it only includes the 8 neighbors of
#'   \code{nbs}.} \item{**directional neighborhood:** }{it consists of the
#'   intersection set of the \code{focal cell} neighbors with the \code{focal
#'   neighborhood}. When _absolute focal neighborhood conditions_ are used, it
#'   includes \code{nbs}, but it does not include the \code{focal cell}; when
#'   _relative focal neighborhood conditions_ are used, it includes the
#'   \code{focal cell}, but it does not include \code{nbs}.}} \cr **Absolute
#'   conditions** \cr This type of condition applies to all functions with a
#'   \code{conditions} argument. It compares between variables (including
#'   \code{"classVector"}) and numeric values. When class contiguity is not
#'   considered, all cells meeting absolute conditions receive a classification
#'   number. When class contiguity is considered, all \code{nbs} meeting
#'   absolute conditions receive a classification number. \cr\cr _Examples of
#'   valid conditions:_ \code{"variable_A > 1 & variable_B != 0"};
#'   \code{"(variable_A^2 < 50 & variable_B == 0) | abs(variable_C) > 50"}. \cr
#'   _Functions:_ \code{\link{anchor.seed}}, \code{\link{cond.4.all}},
#'   \code{\link{cond.4.nofn}} and \code{\link{cond.reclass}}. \cr\cr\cr
#'   **Absolute focal neighborhood conditions** \cr This type of condition
#'   applies to the functions \code{\link{cond.4.nofn}} and
#'   \code{\link{cond.reclass}}. It compares each cell of \code{focal
#'   neighborhood} with a numeric value. The argument \code{fn_perc} controls
#'   the percentage of evaluations that have to be \code{TRUE} in order to
#'   assign a classification number to \code{nbs}. The function
#'   \code{\link{cond.4.nofn}} can consider a \code{directional neighborhood}
#'   instead of a \code{focal neighborhood}. This type of condition is flagged
#'   by a variable name followed by curly brackets (i.e.,
#'   \code{"variable_name{}"}), where the curly brackets indicate the
#'   \code{focal} or \code{directional neighborhood} considered at each
#'   iteration. \cr\cr _Example of valid conditions:_ \code{"variable_A{} > 1 &
#'   abs(variable_B) != 0"}. \cr _Functions:_ \code{\link{cond.4.nofn}} and
#'   \code{\link{cond.reclass}}. \cr\cr\cr **Focal cell conditions.** \cr This
#'   type of condition applies only to functions considering class contiguity
#'   and with a \code{conditions} argument. It compares the value of \code{nbs}
#'   with the value of the \code{focal cell}. This type of condition is flagged
#'   by a variable name followed by square brackets (i.e.,
#'   \code{"variable_name[]"}), where the square brackets indicate the focal
#'   cell. \cr\cr _Examples of valid conditions:_ \code{"variable_A >
#'   variable_A[]"}; \code{"(variable_A > variable_A[] & variable_B{} < 10) |
#'   variable_C > 1"}.  Note that the last example is a combination of absolute
#'   and focal cell conditions. \cr _Functions:_ \code{\link{anchor.seed}} and
#'   \code{\link{cond.4.nofn}}. \cr\cr\cr **Relative focal neighborhood
#'   conditions.** \cr This type of condition applies only to the functions
#'   \code{\link{cond.4.nofn}} and \code{\link{cond.reclass}}. It compares each
#'   cell of \code{focal neighborhood} with \code{nbs}. The argument
#'   \code{fn_perc} controls the percentage of evaluations that have to be
#'   \code{TRUE} in order to assign a classification number to \code{nbs}. The
#'   function \code{\link{cond.4.nofn}} can consider a \code{directional
#'   neighborhood} instead of a \code{focal neighborhood}. This type of
#'   condition is flagged by a variable name followed by curly brackets (i.e.,
#'   \code{"variable_name{}"}), where the curly brackets indicate the
#'   \code{focal} or \code{directional neighborhood} considered at each
#'   iteration. \cr\cr _Example of valid conditions:_ \code{"variable_A >
#'   variable_A{}"}; \code{"(variable_A > variable_A{} & variable_B !=
#'   variable_B[]) | variable_C > 1"}. Note that the last example is a
#'   combination of absolute, focal cell and focal neighborhood conditions. \cr
#'   _Functions:_ \code{\link{cond.4.nofn}} and \code{\link{cond.reclass}}
#'
#' @seealso [anchor.seed()], [attTbl()], [cond.4.all()], [cond.4.nofn()],
#'   [cond.reclass()]
#'
#' @export
#' @examples
#' # EXAMPLES OF VALID AND INVALID CONDITIONS:
#'
#' names_attTbl <- c("bathymetry", "slope")
#
#' cond <- "bathymetry>10"
#' conditions(names_attTbl, cond)
#' cond <- "classVector != 1"
#' conditions(names_attTbl, cond)
#' cond <- "bathymetry[]>10 & abs(slope{}) < 5"
#' conditions(names_attTbl, cond)
#' \dontrun{cond <- "thymetry[]>10 & abs(slpe{}) < 5"
#' conditions(names_attTbl, cond)
#' cond <- "bathymetry[]>10 & | abs(slope{}) < 5"
#' conditions(names_attTbl, cond)
#' cond <- "bathymetry{}>10 & | abs(slope{}) < 5"
#' conditions(names_attTbl, cond)}
#'
#'
#' # EXAMPLES OF FOCAL CELL, NBS CELLS AND NEIGHBORHOODS:
#'
#' # Matrix m mocking a raster of 3 rows and 4 columns
#' m <- matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE)
#' m
#'
#' # FOCAL CELL
#' fc <- 6
#' fc
#'
#' # NBS CELLS
#' nbs <- nbg8(3, 4)[[6]]
#' nbs
#'
#' # CELL IN EVALUATION
#' nbs1 <- nbs[1]
#' nbs1
#'
#' # FOCAL NEIGHBORHOOD
#' # Absolute:
#' ab_fn_nbs1 <- c(nbs1, nbg8(3, 4)[[1]])
#' ab_fn_nbs1
#'
#' # Relative:
#' r_fn_nbs1  <- c(nbg8(3, 4)[[1]])
#' r_fn_nbs1
#'
#' # DIRECTIONAL FOCAL NEIGHBORHOOD
#' # Absolute:
#' ab_dfn_nbs1 <- c(nbs1, intersect(c(nbg8(3, 4)[[1]]) , nbg8(3, 4)[[6]]))
#' ab_dfn_nbs1
#'
#' # Relative:
#' r_dfn_nbs1  <- c(fc, intersect(c(nbg8(3, 4)[[1]]) , nbg8(3, 4)[[6]]))
#' r_dfn_nbs1


conditions <- function(names_attTbl,
                       cond,
                       silent = FALSE) {
  match_pattern <-
    paste0(paste(names_attTbl, collapse = "|"), "|classVector|([0-9]+)")

  # CHECK FOR DOUBLE FOCAL NEIGHBORHOOD CONDITION
  doubleFNcond <- stringr::str_replace_all(cond, "\\&\\&", "\\&")
  doubleFNcond <- stringr::str_replace_all(doubleFNcond, "\\|\\|", "\\|")
  doubleFNcond <- stringr::str_replace_all(doubleFNcond, " ", "")


  doubleFNcond <- unlist(strsplit(doubleFNcond,"\\&|\\|"))
  doubleFNcond <- doubleFNcond[doubleFNcond != ""]

  doubleFNcond <- sapply(as.list(doubleFNcond), function(x)stringr::str_detect(x, "\\{\\}"))

  if(sum(doubleFNcond, na.rm = TRUE) > 1){stop("Only one 'neighborhood' condition allowed.")}

  # CHECK FOR SPELLING ERRORS
  tc_spl <-
    unlist(strsplit(
      cond,
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
        paste(tc_spl, collapse = ", "),
        "' -- Possible spelling error."
      )
    )
  }

  # CHECK FOR SYNTAX ERRORS
  c_eval <- stringr::str_replace_all(cond, paste(c("classVector", names_attTbl), collapse = "\\b|\\b"), "1")
  c_eval <- stringr::str_replace_all(c_eval, "\\[|\\]|\\{|\\}", "")

  c_eval <- try(eval(parse(text = c_eval)), silent = TRUE)

  if(class(c_eval) == "try-error"){

    err <-as.character(stringr::str_match(c_eval, "\\'.+\\'"))
    stop(paste0("Check: ", err, " -- Possible syntax error."))

  }

  if(!silent){

    abs_cond <- stringr::str_detect(cond, paste0(names_attTbl, ("(?!\\[|\\{)"), collapse = "|"))
    fc_cond  <- stringr::str_detect(cond, "\\[\\]")
    fn_cond  <- stringr::str_detect(cond, "\\{\\}")
    cv_cond  <- stringr::str_detect(cond, "classVector")

    verify_types <- c(abs_cond, fc_cond, fn_cond, cv_cond)
    cond_types <- c("'Absolute'", "'Focal cell'", "'Focal neighborhood'", "'Class vector'")


    detc_cond <- cond_types[verify_types]
    if(length(detc_cond) == 0){detc_cond <- "unknown"}

    print(paste0(paste0(detc_cond, collapse = " AND "), " condition type(s) detected)."))
  }
}

#' Class Vector To Raster
#'
#' Transform a class vector into a raster.
#'
#' @param r A \code{Raster*} object on which the values of \code{classVector}
#'   are assigned at the raster cell numbers indicated by \code{index}. To all
#'   remaining raster cells are assigned \code{NA values}.
#' @param index The cell numbers to which each element of \code{classVector}
#'   correspond. If \code{classVector} was computed with
#'   \code{scapesClassificaton} functions and with an attribute table
#'   \code{\link{attTbl}}, this argument correspond to the column
#'   \code{attTbl$Cell} of the attribute table.
#' @param classVector numeric vector, the values to be assigned at the raster
#'   cell numbers indicated by \code{index}.
#' @param plot logic, plot anchor \code{classVector}.
#' @param writeRaster filename, if a raster name is provided save the
#'   \code{classVector} as a raster file.
#' @param overWrite logic, overwrite existing raster.
#'
#' @return The \code{classVector} transformed into a raster.
#'
#' @export
cv.2.rast <- function(r, index, classVector, plot = FALSE, writeRaster = NULL, overWrite = FALSE){

  r2 <- r[[1]]
  r2[] <- NA
  r2[index] <- classVector

  if (plot){
    raster::plot(r2)
  }

  if (!is.null(writeRaster)){
    raster::writeRaster(r2, writeRaster, overwrite = overWrite)
  }

  return(r2)

}

