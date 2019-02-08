#' @title Raw Accelerometry Data Sample
#'
#' @description Sample raw accelerometry data collected during 25 minutes of an
#' outdoor run. Data were collected at frequency 100 Hz with two
#' ActiGraph GT9X Link sensors located at left hip and left ankle.
#'
#' @format A \code{data.frame} with 300000 observations of 5 variables:
#' \itemize{
#'   \item \code{x} - acceleration measurement time series collected from a "x" axis of the sensor accelerometer,
#'   \item \code{y} - acceleration measurement time series collected from a "y" axis of the sensor accelerometer,
#'   \item \code{z} - acceleration measurement time series collected from a "z" axis of the sensor accelerometer,
#'   \item \code{date_time} - date and time of acceleration measurement collection stored as \code{POSIXct},
#'   \item \code{sensor_location} - sensor location label, one of: \code{"left_hip"}, \code{"left_ankle"}.
#' }
#'
#' @details
#' Sample raw accelerometry data were collected during 25 minutes of an
#' outdoor run performed by an adult healthy female, 180 cm tall and of 67 kg weight.
#'
#' Based on a mobile tracking application output, the ground elevation difference between
#' start and end point of the data collection is approximately 36 m (17 m at
#' the start point, 53 m at the finish point). The distance covered is approximately 3.35 km.
#'
#' The data were collected at frequency 100 Hz with two
#' ActiGraph GT9X Link sensors. One of the sensors was attached to the shoe with a clip
#' on the outside side of a left foot, just below the left ankle. The other sensor
#' was attached to the elastic belt located at hip, on the left side of a hip.
#'
#'
#' The person from which the data were collected is Marta Karas, a co-author of the package.
#' The IRB Office Determination Request Form for Primary (New) Data Collection
#' request form was submitted in regard to the collection and further publishing of
#' these data. Based on preliminary review of the
#' request form submitted, it was determined that the data collection and further data publishing
#' activity described in the determination request does not qualify as human subjects research as
#' defined by DHHS regulations 45 CFR 46.102, and does not require IRB oversight.
#'
"acc_running"
