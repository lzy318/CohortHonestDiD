#' Paper color palette
#'
#' @param type Representation, procedure, or cohort colors.
#' @return A character vector of hexadecimal colors, named for representations
#'   and procedures and ordered for cohorts.
#' @export
paper_palette <- function(type = c("representation", "procedure", "cohort")) {
  type <- match.arg(type)
  switch(type,
         representation = c(cohort = "#0072B2", aggregate = "#D55E00"),
         procedure = c(FLCI = "#CC79A7", Conditional = "#E69F00", `C-LF` = "#009E73"),
         cohort = c("#56B4E9", "#009E73", "#CC79A7", "#E69F00", "#0072B2", "#D55E00", "#000000"))
}
