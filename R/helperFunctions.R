#' scapesClassification conditions
#'
#' Check for spelling and syntax errors in conditions (\code{cond} argument) and
#' detect the type of conditions being used.
#'
#' @param names_attTbl character vector, the column (i.e. variable) names of the
#'   attribute table returned by the function \code{\link{attTbl}}.
#' @param cond character string, condition string used to classify raster cells
#'   (see \code{\link{cond.4.nofn}}, \code{\link{cond.reclass}},
#'   \code{\link{cond.4.all}} and \code{\link{anchor.seed}}).
#' @param silent logic, when true, the function returns only error messages.
#'
#' @encoding UTF-8
#'
#' @return An error message if the function finds spelling or syntax errors or a
#'   a string with the types of rules that were detected in the argument
#'   \code{cond}.
#'
#' @details **Conditions (or classification rules)**
#'
#'   * Classification rules evaluate either to true or false and determine what
#'   raster cells belong to a class.
#'
#'   * Conditions are passed to \code{scapesClassification} functions as a
#'   single character string. They can consist of combination of arithmetic
#'   \code{(+|-|*|/|^|%%|%/%)}, relational \code{(>|<|>=|<=|==|!=|%/%)} and
#'   logic operators \code{(&||)}, base R functions (e.g.,
#'   \code{abs(variable_name)}), variables names (as named in the attribute
#'   table, see \code{\link{attTbl}}) and previous classifications (either
#'   stored as \code{classVector} or as rasters).
#'
#'   * A combination of absolute and relative conditions can be used, but only
#'   _**one neighborhood condition per string**_ is allowed.
#'
#'   \cr **Rule evaluation**
#'
#'   * One of the arguments of the classification functions is the
#'   \code{classVector}, a numeric vector that identifies what raster cells have
#'   already been classified (non-NA values) and what have yet to be classified
#'   (NA values). Cells that have already been classified are excluded from the
#'   rule evaluation unless the argument \code{'overwrite_class = TRUE'}.
#'
#'   * __Global evaluation__ \cr Classification rules are applied to all raster
#'   cells (excluding the classified ones). This type of evaluation is common to
#'   classification functions that do not have the argument \code{nbs_of}. Only
#'   absolute conditions can have a global evaluation. See function
#'   \code{\link{cond.4.all}}.
#'
#'   * __Focal evaluation__ \cr Classification rules are applied only to raster
#'   cells contiguous to focal cells. This type of evaluation is common to
#'   classification functions that have the argument \code{nbs_of}. The argument
#'   \code{nbs_of} identifies the class(es) of the focal cells. See functions
#'   \code{\link{anchor.seed}}, \code{\link{cond.4.nofn}} and
#'   \code{\link{cond.reclass}}.
#'
#'       * Focal evaluation can take into account both absolute and relative rules;
#'
#'       * Some classification functions do not have a \code{condition} argument.
#'   Classifications performed by these functions are focal and take into
#'   account only the expected spatial relationships existing among classified
#'   and unclassified cells. See functions \code{\link{reclass.nbs}} and
#'   \code{\link{classify.all}}.
#'
#'   \cr **Focal evaluation, definitions**
#'
#'   * __Cell neighborhood:__ a cell with coordinates \code{(x, y)} has 8
#'   neighbors with coordinates: \code{(x±1, y)},  \code{(x, y±1)} and
#'   \code{(x±1, y±1)}. Cells on the edge of a raster have less than 8
#'   neighbors. See \code{\link{ngbList}}.
#'
#'   * __Focal cell:__ cell identified by one of the classes of the argument
#'   \code{nbs_of}.
#'
#'   * __Test cell:__ the cell in the neighborhood of the focal cell that is
#'   being tested. At turns all cells in the neighborhood of a focal cell are
#'   tested against the classification rule.
#'
#'   * __Directional neighborhood:__ it consists of the intersection between the
#'   focal and the test cell neighborhoods.
#'
#'   \cr**Absolute conditions**
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
#'   0"}; \code{"(variable_A^2 < 50 & variable_B == 0) | abs(variable_C) > 50"}.
#'   \cr _Functions:_ \code{\link{anchor.seed}}, \code{\link{cond.4.all}},
#'   \code{\link{cond.4.nofn}} and \code{\link{cond.reclass}}.
#'
#'   \cr __2) Absolute neighborhood condition:__ compares the values of the
#'   \code{test cell} and of its \code{neighborhood} against a threshold value.
#'
#'   * This type of condition applies to the functions \code{cond.4.nofn} and
#'   \code{cond.reclass}.
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
#'   \code{fn_perc}. The argument \code{fn_perc} ranges from 0 to 1. When 9
#'   evaluations are performed, \code{fn_perc = 1} means that all \code{9}
#'   evaluations have to be true; \code{fn_perc = 0.5} means that at least
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
#'   \code{\link{cond.reclass}}.
#'
#'   \cr **Relative conditions**
#'
#'   __1) Relative focal cell condition:__ compares the \code{test cell} value
#'   against the \code{focal cell} value.
#'
#'   * This type of condition applies only to functions performing focal
#'   evaluation (i.e. function with a \code{nbs_of} argument).
#'
#'   * It is identified by a variable name followed by square brackets (e.g.,
#'   \code{"variable_name[]"}).
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
#'   __2) Relative neighborhood rule:__ compares the values of the \code{test
#'   cell} against the values of the \code{test cell neighborhood}.
#'
#'   * This type of condition applies only to the functions
#'   \code{\link{cond.4.nofn}} and \code{\link{cond.reclass}}.
#'
#'   * It is identified by a variable name followed by curly brackets (e.g.,
#'   \code{"variable_name{}"}).
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
#'   \code{fn_perc}. The argument \code{fn_perc} ranges from 0 to 1. When 8
#'   evaluations are performed, \code{fn_perc = 1} means that all \code{8}
#'   evaluations have to be true; \code{fn_perc = 0.5} means that at least
#'   \code{5} evaluations have to be true.
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
#' # TYPES OF CONDITIONS
#'
#' # As an example consider an attribute with the following columns
#' names_attTbl <- c("bathymetry", "slope")
#'
#' # And the following conditions
#' cond <- "bathymetry>10"
#' conditions(names_attTbl, cond)
#'
#' cond <- "classVector != 1"
#' conditions(names_attTbl, cond)
#'
#' cond <- "bathymetry[]>bathymetry | abs(slope{}) < 5"
#' conditions(names_attTbl, cond)
#'
#' cond <- "bathymetry[]>bathymetry | abs(slope{}) < slope"
#' conditions(names_attTbl, cond)
#'
#' \dontrun{
#' # The function conditions detect syntax and spelling errors
#'
#' cond <- "bathymetry[]>10 & | abs(slope{}) < 5"
#' conditions(names_attTbl, cond)
#'
#' cond <- "baxxxthymetryxxx[]>10 &  abs(slope{}) < 5"
#' conditions(names_attTbl, cond)
#' }
#'
#' # PREPARE PLOT
#' library(scapesClassification)
#' library(raster)
#' library(ggplot2)
#' library(reshape2)
#' library(ggpubr)
#'
#' m     <- matrix(1:49, nrow = 7, ncol = 7, byrow = TRUE)
#'
#' r   <- raster(m)
#' nbs <- ngbList(r)
#'
#' m <- t(m)[,nrow(m):1]
#' m_long <- melt(m)
#'
#' m_long$tags[m_long$value == 32] <- "FC"
#' m_long$tags[m_long$value == 25] <- "TC"
#'
#' npool <- unique(unlist(nbs[c("32", "25")]))
#' m_long$whites[m_long$value %in% npool] <- m_long$value[m_long$value %in% npool]
#'
#' m_long$tags <- as.factor(m_long$tags)
#'
#' p1 <- ggplot(m_long, aes(x=Var1, y=Var2)) +
#'  geom_tile(aes(fill=tags), colour="gray90", lwd=1.5, show.legend = FALSE) +
#'  coord_fixed(ratio=1) +
#'  geom_text(aes(label=value), color = "grey50", family=c("serif"), size=8, na.rm=TRUE) +
#'  geom_text(aes(label=whites), color = "white", family=c("serif"), size=8, na.rm=TRUE) +
#'  scale_fill_manual(values = c("#1088a0", "goldenrod3"), na.value = "black") +
#'  geom_rect(aes(xmin = 2.5, xmax = 5.5, ymin = 1.5, ymax = 4.5),
#'  fill = "transparent", color="#1088a0", lwd = 1.5) +
#'  geom_rect(aes(xmin = 2.5, xmax = 5.5, ymin = 2.5, ymax = 5.5),
#'  fill = "transparent", color="goldenrod3", lwd = 1.5) +
#'  geom_rect(aes(xmin = 2.6, xmax = 5.4, ymin = 2.6, ymax = 4.4),
#'  fill = "transparent", color="red", lwd = 1.2) +
#'  theme_void()  +
#'  scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
#'  labs(title = "Cell numbers") +
#'  theme(plot.title=element_text(size = 16, face = "bold", color = "black", family = "serif"))
#'
#' # FOCAL EVALUATION DEFINITIONS
#'
#' p1
#' # FOCAL CELL
#' # Cell 32
#'
#' # TEST CELL
#' # Cell 25
#'
#' # FOCAL CELL NEIGHBORHOOD
#' # Cells 24, 25, 26, 31, 33, 38, 39, 40
#'
#' # TEST CELL NEIGHBORHOOD
#' # Cells 17, 18, 19, 24, 26, 31, 32, 33
#'
#' # DIRECTIONAL NEIGHBORHOOD
#' # Cells 24, 25, 26, 31, 32, 33

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

    c_split <- unlist(strsplit(cond, "&|\\|"))
    counts  <- stringr::str_count(c_split, paste0(names_attTbl,collapse = "|"))

    abs_cond <- stringr::str_detect(c_split, paste0(names_attTbl, ("(?!\\[|\\{)"), collapse = "|"))
    abs_cond <- any(counts == 1 & abs_cond)

    fc_cond_r <- stringr::str_detect(c_split, "\\[\\]")
    fc_cond_r <- any(counts > 1 & fc_cond_r)

    fc_cond_a <- stringr::str_detect(c_split, "\\[\\]")
    fc_cond_a <- any(counts == 1 & fc_cond_a)

    fn_cond_r <- stringr::str_detect(c_split, "\\{\\}")
    fn_cond_r <- any(counts > 1 & fn_cond_r)

    fn_cond_a <- stringr::str_detect(c_split, "\\{\\}")
    fn_cond_a <- any(counts == 1 & fn_cond_a)

    cv_cond  <- any(stringr::str_detect(c_split, "classVector"))

    verify_types <- c(abs_cond,
                      fc_cond_a,
                      fc_cond_r,
                      fn_cond_a,
                      fn_cond_r,
                      cv_cond)

    cond_types <- c("'Absolute test cell'",
                    "'Absolute focal cell'",
                    "'Relative focal cell'",
                    "'Absolute neighborhood'",
                    "'Relative neighborhood'",
                    "'Class vector'")


    detc_cond <- cond_types[verify_types]
    if(length(detc_cond) == 0){detc_cond <- "unknown"}

    print(paste0(paste0(detc_cond, collapse = " AND "), " condition type(s) detected."))
  }
}

#' Class vector to raster
#'
#' Transform a class vector or a generic vector into a raster.
#'
#' @param r The \code{Raster*} object used to compute the attribute table (see
#'   \code{\link{attTbl}}) or a generic \code{Raster*} object.
#' @param index numeric vector, the cell numbers of the argument \code{r} to
#'   which assign the values of the argument \code{classVector}. Within
#'   \code{scapesClassifications}, this argument normally correspond to the
#'   column \code{Cell} of the attribute table (see \code{\link{attTbl}}).
#' @param vector numeric vector, the values to be assigned to the cell numbers
#'   indicated by \code{index}.
#' @param plot logic, plot the raster.
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

