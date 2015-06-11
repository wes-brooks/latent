#' Analysis results for 55 water samples from Milwaukee and Madison storm sewers
#'
#' A dataset containing the analysis results and metadata for 55 samples of water
#' from storm sewers in Milwaukee and Madison, Wisconsin. The analytic results
#' include counts of several fecal indicator bacteria as well as summary optical
#' properties of the water sample. Optical summary values are calculated from 
#' the full absorption and excitation/emission spectra
#'
#' @format A data frame with 55 rows and 399 variables:
#' \describe{
#'   \item{mei}{count of \emph{Enterococci} in CFU per 100mL}
#'   \item{modmtec}{count of \emph{E. coli} in CFU per 100mL}
#'   \item{FC}{count of fecal coliforms in CFU per 100mL}
#'   \item{Bac.human}{count of human \emph{Bacteroides} in CFU per 100mL}
#'   \item{Lachno.2}{count of human \emph{Lachnobacteria} in CFU per 100mL}
#'   ...
#' }
#' @source Samples collected by the U.S. Geological Survey Wisconsin Water Science Center under Steve Corsi (srcorsi@usgs.gov)
"dfOptAnalysisDataSSJan2015"