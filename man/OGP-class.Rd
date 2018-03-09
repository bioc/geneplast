\name{OGP-class}
\docType{class}
\alias{OGP-class}
\alias{gplast,OGP-method}
\alias{gplast.get,OGP-method}

\title{Class \code{"OGP"}: an S4 class for genetic plasticity analysis.}

\description{
This S4 class includes methods to access the genetic plasticity of orthologous groups.
}

\section{Objects from the Class}{
Objects can be created by calls to the \code{"gplast.preprocess"} constructor.
}

\section{Slots}{
  \describe{
    \item{\code{sspids}:}{Object of class \code{"data.frame"},
    a data frame with species annotation.
    }
    \item{\code{cogids}:}{Object of class \code{"data.frame"},
    a data frame with COG's data.
    }
    \item{\code{orthodist}:}{Object of class \code{"matrix"},
    a matrix with COG's information (see return values in the OGP methods).
    }
    \item{\code{abundance}:}{Object of class \code{"numeric"},
    a numeric vector with results from the \code{gplast} function (see return values in the OGP methods).
    }
    \item{\code{diversity}:}{Object of class \code{"numeric"},
    a numeric vector with results from the \code{gplast} function (see return values in the OGP methods).
    }
    \item{\code{plasticity}:}{Object of class \code{"numeric"},
    a numeric vector with results from the \code{gplast} function (see return values in the OGP methods). 
    }
    \item{\code{status}:}{Object of class \code{"character"},
    a character value specifying the status of the OGP object 
    based on the available methods.
    }
  }
}

\section{Methods}{
  \describe{
    \item{gplast}{\code{signature(object = "OGP")}: see  \code{\link[geneplast:gplast]{gplast}}}
    \item{gplast.get}{\code{signature(object = "OGP")}: see  \code{\link[geneplast:gplast.get]{gplast.get}}}   
	 }
}

\author{
Rodrigo Dalmolin, Mauro Castro
}

\seealso{
  \code{\link[geneplast:gplast.preprocess]{gplast.preprocess}}
}

\keyword{classes}
