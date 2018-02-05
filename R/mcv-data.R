#' Collaborative Cross QTL data for mean red blood cell volume (MCV)
#'
#' Data for mean red blood cell volume (MCV) QTL idenfied in Kelada et al. (2012).
#'
#' @docType data
#'
#' @usage data(mcv.data)
#'
#' @format A list of inputs formatted for use with TIMBR.
#'
#' @keywords datasets
#'
#' @references Kelada et al. (2012) G3 2:157-165
#' (\href{http://www.ncbi.nlm.nih.gov/pubmed/22384394}{PubMed})
#'
#' @source MCV data from Kelada et al. (2012); recomputed diplotype state probabilities using similar methods
#'
#' @examples
#' data(mcv.data)
#' str(mcv.data)
#' TIMBR.results <- TIMBR(mcv.data$y, mcv.data$prior.D, mcv.data$prior.M$crp)
"mcv.data"