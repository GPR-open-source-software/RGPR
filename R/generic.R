
#--------------------------------------------#
#---------------- SETGENERIC ----------------#
# setGenericVerif <- function(x,y){
#   if(!isGeneric(x)){
#     setGeneric(x,y)
#   }else{
#     cat("setGeneric", x,"already exists!\n")
#   }
# }
setGenericVerif <- function(x,y){setGeneric(x,y)}


# # #' @name readGPR
# # #' @rdname readGPR
# # #' @export

#' @title Read a GPR data file
#'
#' @description Function to read GPR data
#'
#' @details
#' Supported format: 
#' \itemize{
#'   \item Sensors and Software (.DT1, .HD)
#'   \item MALA (.rd3, .rad)
#'   \item SEG-Y for RadSys Zond GPR device (.sgy)
#'   \item R (.rds)
#' }
#' If the antenna separation distance is missing, it will be estimated from
#' the GPR frequency.
#'
#' @param fPath Filepath (character).
#' @param desc  Short description of the file (character).
#' @return An object of the class GPR.
#' @examples
#' NULL
#' @name readGPR
#' @rdname readGPR
#  @aliases readGPR-methods
#' @export
setGenericVerif("readGPR", function(fPath, desc = "") standardGeneric("readGPR")
)


#------------------------------

setGenericVerif("as.SpatialPoints", function(x) 
standardGeneric("as.SpatialPoints"))

setGenericVerif("as.SpatialLines", function(x) 
standardGeneric("as.SpatialLines"))


#------------------------------
#' @name coordref
#' @rdname coordref
#' @export
setGenericVerif("coordref", function(x) standardGeneric("coordref"))

#' @name coordref<-
#' @rdname coordref
#' @export
setGenericVerif("coordref<-", function(x, value) standardGeneric("coordref<-"))

setGenericVerif("intersections", function(x) standardGeneric("intersections"))

#' @name filepath
#' @rdname filepath
#' @export
setGenericVerif("filepath", function(x) standardGeneric("filepath"))

#' @name filepath<-
#' @rdname filepath
#' @export
setGenericVerif("filepath<-", function(x, value) standardGeneric("filepath<-"))

setGenericVerif("coords", function(x,i) standardGeneric("coords"))
setGenericVerif("coords<-",function(x,value){standardGeneric("coords<-")})

#' @name coord
#' @rdname coord
#' @export
setGenericVerif("coord", function(x, i, ...) standardGeneric("coord"))

#' @name coord<-
#' @rdname coord
#' @export
setGenericVerif("coord<-",function(x,value){standardGeneric("coord<-")})


#' @name svDate
#' @rdname svDate
#' @export
setGenericVerif("svDate", function(x, i, ...) standardGeneric("svDate"))

#' @name svDate<-
#' @rdname svDate
#' @export
setGenericVerif("svDate<-",function(x,value){standardGeneric("svDate<-")})


#' @name vel
#' @rdname vel
#' @export
setGenericVerif("vel", function(x) standardGeneric("vel"))

#' @name vel<-
#' @rdname vel
#' @export
setGenericVerif("vel<-",function(x,value){standardGeneric("vel<-")})

#' @name ann
#' @rdname ann
#' @export
setGenericVerif("ann", function(x) standardGeneric("ann"))

#' @name ann<-
#' @rdname ann
#' @export
setGenericVerif("ann<-",function(x,value){standardGeneric("ann<-")})

#' @name name
#' @rdname name
#' @export
setGenericVerif("name", function(x) standardGeneric("name"))

#' @name name<-
#' @rdname name
#' @export
setGenericVerif("name<-",function(x,value){standardGeneric("name<-")})

#' @name trZUnit
#' @rdname trZUnit
#' @export
setGenericVerif("trZUnit", function(x) standardGeneric("trZUnit"))

#' @name trZUnit<-
#' @rdname trZUnit
#' @export
setGenericVerif("trZUnit<-",function(x,value){standardGeneric("trZUnit<-")})

#' @name tpUnit
#' @rdname tpUnit
#' @export
setGenericVerif("tpUnit", function(x) standardGeneric("tpUnit"))

#' @name tpUnit<-
#' @rdname tpUnit
#' @export
setGenericVerif("tpUnit<-",function(x,value){standardGeneric("tpUnit<-")})

#' @name crs
#' @rdname crs
#' @export
setGenericVerif("crs", function(x) standardGeneric("crs"))

#' @name crs<-
#' @rdname crs
#' @export
setGenericVerif("crs<-",function(x,value){standardGeneric("crs<-")})

#' @name trZ
#' @rdname trZ
#' @export
setGenericVerif("trZ", function(x) standardGeneric("trZ"))

#' @name trZ<-
#' @rdname trZ
#' @export
setGenericVerif("trZ<-", function(x,value) standardGeneric("trZ<-"))

#' @name tpPos
#' @rdname tpPos
#' @export
setGenericVerif("tpPos", function(x) standardGeneric("tpPos"))

#' @name tpPos<-
#' @rdname tpPos
#' @export
setGenericVerif("tpPos<-", function(x,value) standardGeneric("tpPos<-"))

#' @name time0
#' @rdname time0
#' @export
setGenericVerif("time0", function(x) standardGeneric("time0"))

#' @name time0<-
#' @rdname time0
#' @export
setGenericVerif("time0<-",function(x,value){standardGeneric("time0<-")})

#' Time of data collection for each trace
#'
#' @name trTime
#' @rdname trTime
#' @export
setGenericVerif("trTime", function(x) standardGeneric("trTime"))

#' @name fid
#' @rdname fid
#' @export
setGenericVerif("fid", function(x) standardGeneric("fid"))

#' @name fid<-
#' @rdname fid
#' @export
setGenericVerif("fid<-",function(x,value){standardGeneric("fid<-")})

#' @name values
#' @rdname values
#' @export
setGenericVerif("values", function(x) standardGeneric("values"))

#' @name values<-
#' @rdname values
#' @export
setGenericVerif("values<-", function(x,value) standardGeneric("values<-"))


#' @name proc
#' @rdname proc
#' @export
setGenericVerif("proc", function(x) standardGeneric("proc"))


#' @name proc<-
#' @rdname proc
#' @export
setGenericVerif("proc<-",function(x,value){standardGeneric("proc<-")})

#' @name description
#' @rdname description
#' @export
setGenericVerif("description", function(x) standardGeneric("description"))

#' @name description<-
#' @rdname description
#' @export
setGenericVerif("description<-", function(x, value) 
standardGeneric("description<-"))
#####
setGenericVerif("papply", function(x, prc = NULL) standardGeneric("papply"))
##########
                  
#------------------------------GPR
setGenericVerif("gethd", function(x,hd=NULL) standardGeneric("gethd"))

#' @name plotAmpl
#' @rdname plotAmpl
#' @export
setGenericVerif("plotAmpl", function(x, FUN = mean, add = FALSE, 
                all = FALSE,...) standardGeneric("plotAmpl"))
setGenericVerif("ampl", function(x, FUN=mean, ...) standardGeneric("ampl"))

#' @name interpPos
#' @rdname interpPos
#' @export
setGenericVerif("interpPos", function(x, topo, plot = FALSE,
                                      r = NULL, ...) 
    standardGeneric("interpPos"))

#' @name regInterpPos
#' @rdname regInterpPos
#' @export
setGenericVerif("regInterpPos", function(x, type = c("linear", "cosine"), 
          dx = NULL)  standardGeneric("regInterpPos"))

#' @name relPos
#' @rdname relPos
#' @export
setGenericVerif("relPos", function(x) 
    standardGeneric("relPos"))
    


#' @name writeGPR
#' @rdname writeGPR
#' @export
setGeneric("writeGPR", function(x, fPath = NULL, 
                type = c("DT1", "rds", "ASCII", "xyzv"),
                overwrite = FALSE, ...){ standardGeneric("writeGPR")})

#' @name exportCoord
#' @rdname exportCoord
#' @export
setGenericVerif("exportCoord",  
          function(x, type = c("SpatialPoints", "SpatialLines", "ASCII"),
  fPath = NULL, driver = "ESRI Shapefile", ...) standardGeneric("exportCoord"))
#' @name exportFid
#' @rdname exportFid
#' @export
setGenericVerif("exportFid", function(x, fPath = NULL) 
                  standardGeneric("exportFid"))

#' @name exportProc
#' @rdname exportProc
#' @export
setGenericVerif("exportProc",  function(x,fPath=NULL,sep="\t", row.names=FALSE,
              col.names=FALSE, ...) standardGeneric("exportProc"))

#' @name reverse
#' @rdname reverse
#' @export
setGenericVerif("reverse", function(x, id = NULL,  tol = 0.3) 
                standardGeneric("reverse"))
  
#' @name shiftEst
#' @rdname shiftEst
#' @export
setGenericVerif("shiftEst", function(x, y = NULL, 
                method=c("phase", "WSSD"), dxy = NULL, ...) 
                standardGeneric("shiftEst"))

setGenericVerif("NMOCor", function(x, v = NULL, asep = NULL) 
  standardGeneric("NMOCor"))
setGenericVerif("CMPAnalysis", function(x, method = c("semblance", 
               "winsemblance", "wincoherence", "wincoherence2"), v = NULL, 
               asep = NULL, w = NULL) standardGeneric("CMPAnalysis"))

setGenericVerif("migration", function(x,type=c("static","kirchhoff"), ...) 
standardGeneric("migration"))
setGenericVerif("upsample", function(x,n) standardGeneric("upsample"))
setGenericVerif("timeCorOffset", function(x, t0 = NULL, c0 = 0.299) 
  standardGeneric("timeCorOffset"))

#' @name filter1D
#' @rdname filter1D
#' @export
setGenericVerif("filter1D", function(x, type = c("median", "hampel", 
                "Gaussian"), ...) standardGeneric("filter1D"))

#' @name filter2D
#' @rdname filter2D
#' @export
setGenericVerif("filter2D", function(x, type=c("median3x3", "adimpro"), ...) 
                standardGeneric("filter2D"))
                
setGenericVerif("dewow", function(x, type=c("MAD", "Gaussian"), w ) 
                standardGeneric("dewow"))
                
setGenericVerif("gain", function(x, type=c("power", "exp", "agc"), ...) 
                standardGeneric("gain")) 

setGenericVerif("trAmplCor", 
                function(x, type=c("power", "exp", "agc"),  ...) 
                standardGeneric("trAmplCor"))
                
setGenericVerif("dcshift", function(x, u=1:10, FUN=mean) 
                standardGeneric("dcshift"))
                
setGenericVerif("firstBreak", function(x, method = c("coppens", "coppens2",
                "threshold",  "MER"), thr = 0.12, w = 11, ns = NULL, 
                bet = NULL)
                standardGeneric("firstBreak"))

setGenericVerif("clip", function(x, Amax=NULL,Amin=NULL) 
                standardGeneric("clip"))
                
setGenericVerif("gammaCorrection", function(x, a = 1, b = 1) 
                standardGeneric("gammaCorrection"))
                
setGenericVerif("traceScaling", function(x, 
                  type = c("stat", "min-max", "95", "eq", "sum", "rms", 
                           "mad", "invNormal")) 
                  standardGeneric("traceScaling"))

setGenericVerif("spec", function(x, type = c("f-x", "f-k"), plotSpec = TRUE, 
                unwrapPhase = TRUE, ...) standardGeneric("spec"))
                
setGenericVerif("fFilter", function(x, f = 100, 
                type = c('low', 'high', 'bandpass'),
                L = 257, plotSpec = FALSE) standardGeneric("fFilter"))
                
setGenericVerif("fkFilter", function(x, fk = NULL, L = c(5, 5), npad = 1) 
                standardGeneric("fkFilter"))

setGenericVerif("traceShift", function(x,  ts, method = c("spline", "linear", 
                "nearest", "pchip", "cubic", "none"), crop = TRUE) 
                standardGeneric("traceShift"))
                
setGenericVerif("traceAverage", function(x, w = NULL, FUN = mean, ...) 
                standardGeneric("traceAverage"))

setGenericVerif("time0Cor",  function(x, t0 = NULL, 
                method = c("spline", "linear", "nearest", "pchip", "cubic", 
                "none"), crop = TRUE, keep = 0) 
                standardGeneric("time0Cor"))

setGenericVerif("deconv", function(x, method=c("spiking", "wavelet",
                "min-phase", "mixed-phase"), ...) standardGeneric("deconv"))
setGenericVerif("conv1D", function(x, w) standardGeneric("conv1D"))
setGenericVerif("conv2D", function(x, w) standardGeneric("conv2D"))
setGenericVerif("rotatePhase", function(x, phi) standardGeneric("rotatePhase"))


#------------------------------GRPsurvey
setGenericVerif("getGPR", function(x,id) standardGeneric("getGPR"))
setGenericVerif("surveyIntersect", function(x) 
                standardGeneric("surveyIntersect"))
setGenericVerif("writeSurvey", function(x, fPath, overwrite=FALSE){ 
                standardGeneric("writeSurvey")})

#' Georeferencing
#'
#' Perform on a set of x,y coordinates
#' (1) a translation by \code{-cloc}, then
#' (2) a rotation by \code{alpha} (radian), and (3)
#' a translation by \code{creg}. If \code{creg}
#' is \code{NULL}, then \code{creg} is set equal
#' to \code{cloc}.
#' @param x A matrix with the first two columns corresponding
#'          to coordinates.
#' @param alpha A length-one numeric vector corresponding to 
#'              the rotation angle in radians. If \code{alpha = NULL},
#'              \code{alpha} is estimated from the pairs of points in
#'              the local reference system (\code{ploc}) and in the
#'              regional reference system (\code{preg}).
#' @param cloc A length-two numeric vector corresponding to the coordinate
#'             center of the local reference system
#' @param creg A length-two numeric vector corresponding to the coordinate
#'             center of the regional reference system. Setting 
#'             \code{creg = NULL} (default) is equivalent to apply a rotation
#'             of angle \code{alpha} and center \code{cloc}.
#' @param ploc A matrix with the first two columns corresponding
#'             to coordinates in the local reference system.
#' @param preg A matrix with the first two columns corresponding
#'             to coordinates in the regional reference system.
#' @param FUN If \code{alpha = NULL}, a function to estimate the rotation angle
#'            from the angles computed for each pairs of coordinates of
#'            \code{ploc}-\code{preg}.
#' @export
#' @name georef
#' @rdname georef
setGenericVerif("georef", function(x, alpha = NULL, cloc = c(0,0), creg = NULL,
                   ploc = NULL, preg = NULL, FUN = mean){ 
                standardGeneric("georef")})
                
#------------------------------BOTH
setGenericVerif("plot3DRGL", 
          function(x, addTopo = FALSE, clip = NULL, normalize = NULL, 
                  nupspl = NULL, add = TRUE, xlim = NULL, ylim = NULL, 
                  zlim = NULL,...) 
standardGeneric("plot3DRGL"))

setGenericVerif("exportPDF", function(x, fPath = NULL, addTopo = FALSE, 
                clip = NULL, normalize = NULL, nupspl = NULL, ...) 
standardGeneric("exportPDF"))

#setGenericVerif("adimproSmooth", function(x,hmax=2,...) standardGeneric("
# adimproSmooth"))

#---------------------- DELINEATIONS ---------------------#
#' @name delineate
#' @rdname delineation
#' @export
setGenericVerif("delineate", function(x, name = NULL,
                type = c("raster", "wiggles"),
                addTopo = FALSE, nupspl = NULL, n = 10000, ...) 
                  standardGeneric("delineate"))
#' @name rmDelineations<-
#' @rdname delineation
#' @export
setGenericVerif("rmDelineations<-", function(x,value=NULL) 
                  standardGeneric("rmDelineations<-"))
#' @name delineations
#' @rdname delineation
#' @export
setGenericVerif("delineations", function(x,sel=NULL,...) 
                  standardGeneric("delineations"))
#' @name addDelineation
#' @rdname delineation
#' @export
setGenericVerif("addDelineation", function(x,...) 
                  standardGeneric("addDelineation"))
setGenericVerif("showDelineations", function(x,sel=NULL,...) 
                  standardGeneric("showDelineations"))
#' @name exportDelineations
#' @rdname delineation
#' @export
setGenericVerif("exportDelineations", function(x, dirpath="") 
                  standardGeneric("exportDelineations"))
#' @name plotDelineations3D
#' @rdname delineation
#' @export
setGenericVerif("plotDelineations3D", 
                function(x,sel=NULL,col=NULL,add=TRUE,...)
                standardGeneric("plotDelineations3D"))
#' @name plotDelineations
#' @rdname delineation
#' @export
setGenericVerif("plotDelineations", function(x,sel=NULL,col=NULL,...) 
                  standardGeneric("plotDelineations"))
#' @name identifyDelineation
#' @rdname delineation
#' @export
setGenericVerif("identifyDelineation", function(x,sel=NULL,...) 
                  standardGeneric("identifyDelineation"))

#' Structure tensor field of GPR data 
#'
#' @name strTensor
#' @rdname strTensor
#' @export
setGenericVerif("strTensor", function(x,  blksze = c(2, 4),
                        kBlur   = list(n = 1, m = 1, sd = 1), 
                        kEdge   = list(n = 5, m = 5, sd = 1), 
                        kTensor = list(n = 5, m = 5, sd = 1),
                        thresh = 0.02, what = c("tensor", "mask"), ...)
                        standardGeneric("strTensor"))
                  
                  
                  
if (!isGeneric("stack")) {
  setGeneric("stack", function(x, ...)
    standardGeneric("stack"))
} 