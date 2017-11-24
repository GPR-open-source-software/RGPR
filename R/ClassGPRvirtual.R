
#------------------------------------------#
#----------- CLASS DEFINITION -------------#

#' An S4 class to represent a ground-penetrating radar (GPR) data.
#'
#' @slot version      A length-one character vector indicating the version of 
#'                    RGPR
#' @slot data         A \eqn{m \times n} numeric matrix consiting of a 
#'                    cross-section of signal amplitudes as a function of the 
#'                    GPR position. The columns of \code{data} correspond to 
#'                    the GPR traces and the row of \code{data} to the
#'                    time/depth samples.
#' @slot traces       A length-m numeric vector corresponding to the trace 
#'                    number.
#' @slot z            A length-n numeric vector indicating the sampling time or
#'                    the vertical position of the trace samples.
#' @slot x            A length-m numeric vector indicating the relative 
#'                    position of the trace along the survey profile.
#' @slot time0        A length-m numeric vector containing the 'time-zero' of
#'                    every trace.
#' @slot time         A length-m numeric vector containing the recording time 
#'                    of  every trace.
#' @slot fid          A length-m character vector containing fiducial markers 
#'                    associated with the traces.
#' @slot ann          A length-m character vector containing annotations 
#'                    associated  with the traces.
#' @slot coord        A \eqn{m \times 3} matrix containing the (x, y, z) 
#'                    positions of every trace.
#' @slot rec          A \eqn{m \times 3} matrix containing the (x, y, z) 
#'                    positions of the receiver for every trace.
#' @slot trans        A \eqn{m \times 3} matrix containing the (x, y, z) 
#'                    positions of the transmitter for every trace.
#' @slot coordref     A length-3 numeric vector containing the coordinates of a
#'                    local reference.
#' @slot freq         A length-one numeric vector corresponding to the GPR 
#'                    antennae frequency (in MHz).
#' @slot dz           A length-one numeric vector corresponding to the time or
#'                    depth sampling step.
#' @slot dx           A length-one numeric vector corresponding to the trace
#'                    step.
#' @slot antsep       A length-one numeric vector corresponding to the antenna
#'                    separation.
#' @slot name         A length-one character vector containing the name of the
#'                    GPR data.
#' @slot description  A length-one character vector containing the description
#'                    of the GPR data.
#' @slot filepath     A length-one character vector containing the file path
#'                    of the original GPR data.
#' @slot zunit        A length-one character vector corresponding to the 
#'                    time/depth unit (e.g., "ns", "m").
#' @slot xunit        A length-one character vector corresponding to the 
#'                    (x, y)-unit (e.g., "m").
#' @slot surveymode   A length-one character vector containing the survey mode
#'                    (e.g., "Reflection", "CMP")
#' @slot date         A length-one character vector containing the date of the
#'                    survey in the format "yyyy-mm-dd".
#' @slot crs          A length-one character vector containing the coordinate
#'                    reference system following the R notation of proj4string
#'                    from the PROJ.4 library. 
#' @slot proc         A length-varying character vector whose each element
#'                    correspond to a processing step applied to the data.
#' @slot vel          A list containing the velocity model.
#' @slot delineations A list containing delineated structures.
#' @slot hd           A list of less relevant additional informations.
setClass(
  Class = "GPRvirtual",
  contains = "VIRTUAL",
  slots = c(
    version      = "character",  # class version
    data         = "ANY",        # 
    z            = "numeric",    # depth/time position
    x            = "numeric",    # trace position
    time0        = "numeric",    # time-zero (first air-wave arrival)
    time         = "numeric",    # time of the trace recording
    fid          = "character",  # fiducial marks/markers
    ann          = "character",  # annotation (e.g. intersections)
    coord        = "matrix",     # trace coordinates (x,y,z)
    rec          = "matrix",     # receiver coordinates (x,y,z)
    trans        = "matrix",     # transmitter coordinates (x,y,z)
    coordref     = "numeric",    # coordinates references
    freq         = "numeric",    # antenna frequency
    dz           = "numeric",    # time/depth sampling
    dx           = "numeric",    # spatial trace sampling
    antsep       = "numeric",    # antenna separation
    name         = "character",  # data name
    description  = "character",  # data description
    filepath     = "character",  # data filepath
    zunit    = "character",  # time/depth unit
    xunit      = "character",  # spatial unit
    surveymode   = "character",  # survey mode (reflection/CMP)
    date         = "character",  # survey date (format %Y-%m-%d)
    crs          = "character",  # coordinate reference system of coord
    proc         = "character",  # processing steps
    vel          = "list",       # velocity model
    delineations = "list",       # delineations
    hd           = "list"        # data from header file/meta-data
  )
)

