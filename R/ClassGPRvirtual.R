
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


#' @export
setAs(from = "GPRvirtual", to = "vector", def = function(from){
  return( as.vector(from@data) )
})

#' @rdname coerceFromGPR
#' @export
setMethod(f = "as.vector", signature(x = "GPRvirtual"),
          definition = function(x, mode = "any"){
  return( as.vector(x@data) )
})

#' @rdname coerceFromGPR
#' @export
setMethod(f = "as.numeric", signature(x = "GPRvirtual"),
          definition = function(x, ...){
  return( as.numeric(x@data) )
})




#' Number of rows of an object of the class \code{GPR}
#'
#' = number of time samples
#' @examples
#' data(frenkeLine00)
#' nrow(frenkeLine00)
#' @export
setMethod(f = "nrow", signature(x = "GPRvirtual"), definition = function(x){
    nrow(x@data)
})

#' Number of columns of an object of the class \code{GPR}
#'
#' = number of traces.
#' @examples
#' data(frenkeLine00)
#' ncol(frenkeLine00)
#' @export
setMethod(f = "ncol", signature(x = "GPRvirtual"), definition = function(x){
    ncol(x@data)
})

#' Dimenstionf an object of the class \code{GPR}
#'
#' = number of time samples, number of traces.
#' @examples
#' data(frenkeLine00)
#' dim(frenkeLine00)
#' @export
setMethod(f = "dim", signature = "GPRvirtual", definition = function(x){
    dim(x@data)
})

#' Summary of GPR traces
#'
#' Summary of GPR traces
#' @examples
#' data(frenkeLine00)
#' summary(frenkeLine00)
#' @export
setMethod(f = "summary", signature(object = "GPRvirtual"),
          definition = function(object){
  return( summary(object@data) )
})

#' Mean of GPR data
#'
#' Mean of GPR data
#' @examples
#' data(frenkeLine00)
#' mean(frenkeLine00)
#' @export
setMethod(f = "mean", signature(x = "GPRvirtual"), definition = function(x){
  return( mean(x@data) )
})

#' Median of GPR data
#'
#' Median of GPR data
#' @examples
#' data(frenkeLine00)
#' median(frenkeLine00)
#' @export
setMethod(f = "median", signature(x = "GPRvirtual"), definition = function(x){
  return( median(x@data) )
})

#' Range of GPR data
#'
#' Range of GPR data
#' @examples
#' data(frenkeLine00)
#' range(frenkeLine00)
#' @export
setMethod(f= "range", signature(x = "GPRvirtual"), 
          definition = function(x, ..., na.rm = FALSE){
  return( range(as.matrix(x), ..., na.rm = na.rm) )
})

#' Apply Functions Over Array Margins
#'
#' Returns a vector or array or list of values obtained by applying a function 
#' to margins of an array or matrix.
#' @examples
#' data(frenkeLine00)
#' apply(frenkeLine00, 2, mean)
#' @export
setMethod(f = "apply", signature( X = "GPRvirtual"), 
          definition = function(X, MARGIN, FUN, ...){
  return(apply(X@data, MARGIN, FUN, ...))
})


#' Mininum value of GPR data
#'
#' Mininum value of GPR data
#' @examples
#' data(frenkeLine00)
#' min(frenkeLine00)
#' @export
setMethod(f = "min", signature = "GPRvirtual", 
          definition = function(x, ..., na.rm = FALSE){
  min(x@data,na.rm=na.rm)
})

#' Maximum value of GPR data
#'
#' Maximum value of GPR data
#' @examples
#' data(frenkeLine00)
#' min(frenkeLine00)
#' @export
setMethod(f = "max", signature = "GPRvirtual", 
          definition = function(x, ..., na.rm = FALSE){
  max(x@data,na.rm=na.rm)
})



#' Basic mathematical functions
#'
#' Methods for the base Math methods \link[methods]{S4groupGeneric}
#' @param x An object of the class RGPR.
#' @details Currently implemented methods include:
#' \itemize{
#'  \item{"abs", "sign", "sqrt", "ceiling", "floor", "trunc",
#'        "exp", "expm1", "log", "log10", "log2", "log1p", "cos",
#'        "cosh", "sin", "sinh", "tan", "tanh"}
#'  }
#' @examples
#' data(frenkeLine00)
#' A <- exp(frenkeLine00)
#' @rdname Math-methods
#' @aliases Math-GPR-method
#' @export
# getGroupMembers("Math")
setMethod(
  f="Math",
  signature="GPRvirtual",
  definition=function(x){
    switch(.Generic,
      abs     = abs(x@data),
      sign    = sign(x@data),
      sqrt    = sign(x@data)*sqrt(abs(x@data)),
      ceiling = ceiling(x@data),
      floor   = floor(x@data),
      trunc   = trunc(x@data),
      cummax  = paste("not allowed"),  
      cumprod = paste("not allowed"),  
      cumsum  = paste("not allowed"),  
      exp     = exp(x@data),
      expm1   = expm1( x@data),
      log     = sign(x@data) *log(abs(x@data)),
      log10   = sign(x@data) *log10(abs(x@data)),
      log2    = sign(x@data) *log2(abs(x@data)), 
      log1p   = sign(x@data) *log1p(abs(x@data)), 
      cos     = cos(x@data),
      cosh    = cosh(x@data),
      sin     = sin(x@data),
      sinh    = sinh(x@data),
      tan     = tan(x@data),
      tanh    = tanh(x@data),
        # "acos" "acosh" "asin" "asinh" "atan" "atanh" "gamma" "lgamma"
        # "digamma" "trigamma"
      stop(paste(.Generic, "not allowed on GPR objects"))
    )
    proc(x) <- getArgs()
    return(x)
  }
)

# > getGroupMembers("Arith")
# [1] "+"   "-"   "*"   "^"   "%%"  "%/%" "/" 
.GPR.add <- function(a, b){
  #FIXME #TODO: case where a (or b) is a vector trace
  if(is(b,"GPRvirtual")){
    x <- b
    b <- b@data
  }
  if(is(a,"GPRvirtual")){
    x <- a
    a <- a@data
  }
  x@data <- a + b
  return(x)
}
.GPR.sub <- function(a, b){
  if(is(b,"GPRvirtual")){
    x <- b
    b <- b@data
  }
  if(is(a,"GPRvirtual")){
    x <- a
    a <- a@data
  }
  if(!is.null(dim(a)) && !is.null(dim(b))){
    if(dim(a)[2] == 1 && dim(b)[2] > 1){
      a <- as.vector(a)
    }else if(dim(b)[2] == 1 && dim(a)[2] > 1){
      b <- as.vector(b)
    }
  }
  x@data <- a - b
  return(x)
}
.GPR.mul <- function(a, b){
  if(is(b,"GPRvirtual")){
    x <- b
    b <- b@data
  }
  if(is(a,"GPRvirtual")){
    x <- a
    a <- a@data
  }
  x@data <- a * b
  return(x)
}
.GPR.div <- function(a, b){
  if(is(b,"GPRvirtual")){
    x <- b
    b <- b@data
  }
  if(is(a,"GPRvirtual")){
    x <- a
    a <- a@data
  }
  x@data <- a / b
  return(x)
}
.GPR.pow <- function(a, b){
  if(is(b, "GPRvirtual")){
    x <- b
    b <- b@data
  }
  if(is(a,"GPRvirtual")){
    x <- a
    a <- a@data
  }
  x@data <- a ^ b
  return(x)
}

.GPR.arith <- function(e1,e2){
  switch(.Generic,
    "+" = .GPR.add(e1, e2),
    "-" = .GPR.sub(e1, e2),
    "*" = .GPR.mul(e1, e2),
    "/" = .GPR.div(e1, e2),
    "^" = .GPR.pow(e1, e2),
    stop(paste("binary operator \"", .Generic, "\" not defined for GPR"))
  )
}



#' Basic arithmetical functions
#'
#' 
#' @param e1 An object of the class GPR, GPRstack or GPRcube.
#' @param e2 An object of the class GPR, GPRstack or GPRcube
#' @examples
#' data(frenkeLine00)
#' A <- exp(frenkeLine00)
#' B <- A + frenkeLine00
#' @name Arith
#' @rdname Arith-methods
#' @export
setMethod(f = "Arith", signature(e1 = "GPRvirtual", e2 = "ANY"),
          definition = .GPR.arith)

#' @name Arith
#' @rdname Arith-methods
#' @export
setMethod(f = "Arith", signature(e1 = "GPRvirtual", e2 = "GPRvirtual"),
          definition = .GPR.arith)

#' @name Arith
#' @rdname Arith-methods
#' @export
setMethod(f = "Arith", signature(e1 = "ANY", e2 = "GPRvirtual"),
          definition = .GPR.arith)

          
##================================ GETTER/SETTER =============================##
##============================== SETTER/GETTER ===============================##
#' @export
setMethod("gethd", "GPRvirtual", function(x,hd=NULL){
    if(is.null(hd)){
      return(x@hd)
    }else{
      if(!is.null(x@hd[[hd]])){
        if(is.character(x@hd[[hd]])){
          as.character(x@hd[[hd]])
        }else{
          as.numeric(x@hd[[hd]])
        }
      }
    }
  } 
)

#' Name of the GPR data
#' 
#' @name name
#' @rdname name
#' @export
setMethod("name", "GPRvirtual", function(x){
  return(x@name)
} 
)

#' @name name<-
#' @rdname name
#' @export
setReplaceMethod(
  f="name",
  signature="GPRvirtual",
  definition=function(x,value){
    value <- as.character(value)[1]
    x@name <- value
    x@proc <- c(x@proc, "name<-")
    return(x)
  }
)

#' Depth unit of the GPR data
#' 
#' @name trZUnit
#' @rdname trZUnit
#' @export
setMethod("trZUnit", "GPRvirtual", function(x){
  return(x@zunit)
} 
)

#' @name trZUnit<-
#' @rdname trZUnit
#' @export
setReplaceMethod(
  f="trZUnit",
  signature="GPRvirtual",
  definition=function(x,value){
    value <- as.character(value)[1]
    x@zunit <- value
    x@proc <- c(x@proc, "trZUnit<-")
    return(x)
  }
)

#' Position unit of the GPR data
#' 
#' @name tpUnit
#' @rdname tpUnit
#' @export
setMethod("tpUnit", "GPRvirtual", function(x){
  return(x@xunit)
} 
)

#' @name tpUnit<-
#' @rdname tpUnit
#' @export
setReplaceMethod(
  f="tpUnit",
  signature="GPRvirtual",
  definition=function(x,value){
    value <- as.character(value)[1]
    x@xunit <- value
    x@proc <- c(x@proc, "tpUnit<-")
    return(x)
  }
)

#' Description of the GPR data
#' 
#' @name description
#' @rdname description
#' @export
setMethod("description", "GPRvirtual", function(x){
  return(x@description)
} 
)

#' @name description<-
#' @rdname description
#' @export
setReplaceMethod(
  f="description",
  signature="GPRvirtual",
  definition=function(x, value){
    x@description <- as.character(value)[1]
    x@proc <- c(x@proc, "description<-")
    return(x)
  }
)

#' Filepath of the GPR data
#' 
#' @name filepath
#' @rdname filepath
#' @export
setMethod("filepath", "GPRvirtual", function(x){
    return(x@filepath)
  } 
)

#' @name filepath<-
#' @rdname filepath
#' @export
setReplaceMethod(
  f="filepath",
  signature="GPRvirtual",
  definition=function(x,value){
    x@filepath <- as.character(value)[1]
    x@proc <- c(x@proc, "filepath<-")
    return(x)
  }
)

#' Survey date
#' 
#' Return NULL if no date exists, else an object of the class 'Date'
#' @name svDate
#' @rdname svDate
#' @export
setMethod("svDate", "GPRvirtual", function(x){
    if(length(x@date) > 0){
    return(as.Date(x@date))
    }else{
    return(NULL)
  }
  } 
)

#' @param x An object of the class 'GPR'
#' @param value An object of the class 'Date'
#' @name svDate<-
#' @rdname svDate
#' @export
setReplaceMethod(
  f="svDate",
  signature="GPRvirtual",
  definition=function(x,value){
    value <- value[1]
    if(class(value) == "Date"){
      x@date <- as.character(value)
      x@proc <- c(x@proc, "svDate<-")
      return(x)
    }else{
      stop("'value' must be from the class 'Date'!")
  }
  }
)

#' Coordinate reference system (CRS) of the GPR data
#' 
#' @name crs
#' @rdname crs
#' @export
setMethod("crs", "GPRvirtual", function(x){
    return(x@crs)
  } 
)

#' @name crs<-
#' @rdname crs
#' @export
setReplaceMethod(
  f="crs",
  signature="GPRvirtual",
  definition=function(x,value){
    value <- as.character(value)[1]
    x@crs <- value
    x@proc <- c(x@proc, "crs<-")
    return(x)
  }
)

#' Velocity model of the GPR data
#' 
#' @name vel
#' @rdname vel
#' @export
setMethod("vel", "GPRvirtual", function(x){
  if(length(x@vel) == 1){
    return(x@vel[[1]])
  }else{
    return(x@vel)
  }
} 
)

#' @name vel<-
#' @rdname vel
#' @export
setReplaceMethod(
  f="vel",
  signature="GPRvirtual",
  definition=function(x, value){
    if(typeof(value) != "list"){
      value <- list(as.numeric(value))
    }
    x@vel <- value
    x@proc <- c(x@proc, "vel<-")
    return(x)
  }
)


#' Depth/time of the GPR data
#' 
#' @name trZ
#' @rdname trZ
#' @export
setMethod("trZ", "GPRvirtual", function(x){
  return(x@z)
} 
)
#' @name trZ<-
#' @rdname trZ
#' @export
setReplaceMethod(
  f = "trZ",
  signature = "GPRvirtual",
  definition = function(x, value){
    if(length(value) == length(x@z)){
      x@z <- value
    }else{
      stop("length(value) != length(x@z)")
    }
    x@proc <- c(x@proc, "trZ<-")
    return(x)
  }
)

#' Values of the GPR data
#' 
#' @name values
#' @rdname values
#' @export
setMethod("values", "GPRvirtual", function(x){
    return(x@data)
  } 
)

#' @name values<-
#' @rdname values
#' @export
setReplaceMethod(
  f="values",
  signature="GPRvirtual",
  definition=function(x,value){
    if(all(dim(value)==dim(x))){
      x@data <- value
      x@proc <- c(x@proc, "values<-")
      return(x)
    }else{
      stop("x [",nrow(x),"x",ncol(x),"] and A [",nrow(value),"x",ncol(value),
            "] should have the same size\n")
    }
  } 
)




#' Processing steps applied to the data
#' 
#' \code{processing} returns all the processing steps applied to the data.
#' 
#' @param x An object of the class GPR.
#' @return A character vector whose elements contain the name of the 
#' processing functions with their arguments applied previously on the
#' GPR data.
#' @examples
#' data(frenkeLine00)
#' A <- dewow(frenkeLine00, type = "Gaussian")
#' proc(A)
#' @name proc
#' @rdname proc
#' @export
setMethod("proc", "GPRvirtual", function(x){
    return(x@proc)
  } 
)

#' Add a processing step
#' 
#' @name proc
#' @rdname proc
#' @export
setReplaceMethod(
  f="proc",
  signature="GPRvirtual",
  definition=function(x,value){
    value <- as.character(value)
    x@proc <- c(x@proc, value)
    return(x)
  }
)