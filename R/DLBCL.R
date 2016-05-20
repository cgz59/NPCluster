#' Survival in diffuse large B-cell lymphoma
#'
#' A dataset containing survival times for 235 patients after chemotherapy and
#' their associated 3782 microarray expression measurements.   Five patients with 0
#' survival times have been removed and expression measurements were filtered via
#' [FILL-IN].
#'
#' @name dlbcl
#' @docType data
#'
#' @source Rosenwald, A. et al. (2002), “The Use of Molecular Profiling to Predict
#' Survival After Chemotherapy for Diffuse Large B-Cell Lymphoma,”
#' New England Journal of Medicine, 346, 1937–1947.
#'
#' @format A data frame with 235 rows and 3784 variables:
#' \describe{
#'   \item{years}{observation time}
#'   \item{death}{0/1 indicator of death at end of observation}
#'   ...
#' }
NULL
