#' @title df.absorbance
#' @description Raw Absorbance data for 12 RNA helices. Note the UAUAUAUA helix was not included in the final analysis because it did not have well defined lower absorbance baselines.
#' @format A data frame with 1428 rows and 10 variables:
#' \describe{
#'   \item{\code{Temperature}}{double Temperature in degCelcius}
#'   \item{\code{Absorbance}}{double The absorbance measured using a spectrometer}
#'   \item{\code{Experiment}}{character The experiment name}
#'   \item{\code{RNA}}{character The RNA sequence}
#'   \item{\code{Buffer}}{logical The buffer}
#'   \item{\code{File}}{character The file path}
#'   \item{\code{Cell}}{character The cell in a run}
#'   \item{\code{Blank}}{character The blank the spectrometer uses}
#'   \item{\code{Pathlength}}{double The pathlength of the cuvette in cm}
#'   \item{\code{Wavelength}}{double The wavelength the absorbance was measured at in nm} 
#'}
#' @source \url{http://somewhere.important.com/}
"df.absorbance"

#' @title df.Meltwin
#' @description Meltwin results for 12 RNA helix data sets. Note the UAUAUAUA helix was not included in the final analysis because it did not have well defined lower absorbance baselines.
#' @format A data frame with 24 rows and 9 variables:
#' \describe{
#'   \item{\code{Sequence}}{character The RNA the data was collected on}
#'   \item{\code{Method}}{character The Meltwin method}
#'   \item{\code{dH}}{num The enthalpy of helix formation in kcal/mol}
#'   \item{\code{SE.dH}}{num The error in the enthalpy of helix formation in kcal/mol}
#'   \item{\code{dS}}{num The entropy of helix formation in cal/mol/K}
#'   \item{\code{SE.dS}}{num The error in the entropy of helix formation in cal/mol/K}
#'   \item{\code{dG}}{num The Gibbs free energy of helix formation in kcal/mol}
#'   \item{\code{SE.dG}}{num The error in the Gibbs free energy of helix formation in kcal/mol}
#'   \item{\code{Tm}}{num The Tm assuming the Ct is 0.1 mM in degC}
#'}
#' @source \url{http://somewhere.important.com/}
"df.absorbance"