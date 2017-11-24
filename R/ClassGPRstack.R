#------------------------------------------#
#----------- CLASS DEFINITION -------------#

#' An S4 class to represent a ground-penetrating radar (GPR) data.
#'
#' @slot version      A length-one character vector indicating the version of 
#'                    RGPR
#' @slot data         A \eqn{m \times n \time p} numeric arry consiting of
#'                    several two-dimensional ground-penetrating (GPR) data 
#'                    having the same size (often results of GPR data 
#'                    decomposition through transform functions.
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
  Class = "GPRstack",
  contains = "GPRvirtual",
  slots = c( data = "array")
)


#' @export
setMethod(f = "stack", signature = "GPR", definition = function(x, ...){
  rlist <- list(...)
  if ( length(rlist) == 1 && inherits(rlist[[1]], "GPR")) {
    y <- rlist[[1]]
    if( all.equal(dim(x), dim(y)) ){
      newData <- array(dim = c(dim(x), 2))
      newData[,,1] <- x@data
      newData[,,2] <- y@data
    }
  }else if( length(rlist) == 0){
    newData <- array(dim = c(dim(x), 1))
    newData[,,1] <- x@data
  }
  y <- new("GPRstack",
            version      = x@version,
            data         = newData,
            z            = x@z,
            x            = x@x,
            time0        = x@time0,
            time         = x@time,
            fid          = x@fid,
            ann          = x@ann,
            coord        = x@coord,
            rec          = x@rec,
            trans        = x@trans,
            coordref     = x@coordref,
            freq         = x@freq,
            dz           = x@dz,
            dx           = x@dx,
            antsep       = x@antsep,
            name         = x@name,
            description  = x@description,
            filepath     = x@filepath,
            zunit        = x@zunit,
            xunit        = x@xunit,
            surveymode   = x@surveymode,
            date         = x@date,
            crs          = x@crs,
            proc         = x@proc,
            vel          = x@vel,
            delineations = x@delineations,
            hd           = x@hd
          )
  return(y)
})    
    
    
#' Length of an object of the class \code{GPR}
#'
#' The length of an object of the class \code{GPR} is equal to the number of
#' traces.
#' @examples
#' data(frenkeLine00)
#' length(frenkeLine00)
#' @export
setMethod(f = "length", signature(x = "GPRstack"), definition = function(x){
  return( prod(dim(x@data)[2:3]) )
})
    
    
    
    
    
    
    
    
    
    