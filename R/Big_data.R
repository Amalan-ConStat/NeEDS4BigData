#' Electric consumption data
#'
#' Hebrail and Berard (2012) described data which contains 2,049,280 completed measurements for a house
#' located at Sceaux, France between December 2006 and November 2010.
#' The log scale minute-averaged current intensity is selected as the response and
#' the covariates are active electrical energy (watt-hour) in the kitchen, the laundry room,
#' and electric water-heater and an air-conditioner.
#'
#' @format A data frame with 4 columns and 2,049,280 rows.
#' \describe{
#' \item{\code{Intensity}}{Minute-averaged current intensity}
#' \item{\code{EE_Kitchen}}{Active electrical energy (watt-hour) in the kitchen}
#' \item{\code{EE_Laundry}}{Active electrical energy (watt-hour) in the laundry room}
#' \item{\code{EE_WH_AC}}{Active electrical energy (watt-hour) of electric water-heater and an air-conditioner}
#' }
#'
#' @examples
#' nrow(Electric_consumption)
#'
#' @source
#' Extracted from
#'
#' Hebrail G, Berard A (2012) Individual Household Electric Power Consumption.
#' UCI Machine Learning Repository.
#'
#' Available at: \doi{10.24432/C58K54}
#'
"Electric_consumption"

#' Skin segmentation data
#'
#' Rajen and Abhinav (2012) addressed the challenge of detecting skin-like
#' regions in images as a component of the intricate process of facial recognition.
#' To achieve this goal, they curated the “Skin segmentation” data set, comprising
#' RGB (R-red, G-green, B-blue) values of randomly selected pixels from
#' N = 245,057 facial images, including 50,859 skin samples and 194,198 nonskin
#' samples, spanning diverse age groups, racial backgrounds, and genders.
#'
#' @format A data frame with 4 columns and 245,057 rows.
#' \describe{
#' \item{\code{Skin_presence}}{Skin presence in the randomly selected pixels}
#' \item{\code{Red}}{Red values in the randomly selected pixels}
#' \item{\code{Green}}{Green values in the randomly selected pixels}
#' \item{\code{Blue}}{Blue values in the randomly selected pixels}
#' }
#'
#' @examples
#' nrow(Skin_segmentation)
#'
#' @source
#' Extracted from
#'
#' Rajen B, Abhinav D (2012) Skin segmentation. UCI Machine Learning Repository.
#'
#' Available at: \doi{10.24432/C5T30C}
#'
"Skin_segmentation"

#' Bike sharing data
#'
#' Fanaee-T (2013) collected data to understand the bike sharing demands under
#' the rental and return process. The data in total contains 17,379 observations
#' where we consider the covariates temperature, humidity and windspeed
#' to model the response, the number of bikes rented hourly.
#'
#' @format A data frame with 4 columns and 17,379 rows.
#' \describe{
#' \item{\code{Rented_Bikes}}{Number of bikes rented hourly}
#' \item{\code{Temperature}}{Hourly temperature}
#' \item{\code{Humidity}}{Hourly humidity}
#' \item{\code{Windspeed}}{Hourly windspeed}
#' }
#'
#' @examples
#' nrow(Bike_sharing)
#'
#' @source
#' Extracted from
#'
#' Fanaee-T H (2013) Bike Sharing. UCI Machine Learning Repository.
#'
#' Available at: \doi{10.24432/C5W894}
#'
"Bike_sharing"
