\name{OGR-class}
\Rdversion{1.1}
\docType{class}
\alias{OGR-class}
\alias{groot,OGR-method}
\alias{groot.get,OGR-method}

\title{Class \code{"OGR"}: an S4 class for rooting analysis.}

\description{
This S4 class includes methods to do inferential analysis of evolutionary roots in a given species tree.
}

\section{Objects from the Class}{
Objects can be created by calls to the \code{"groot.preprocess"} constructor.
}

\section{Slots}{
  \describe{
    \item{\code{cogids}:}{Object of class \code{"data.frame"},
    a data frame with COG's data.
    }
    \item{\code{tree}:}{Object of class \code{"phylo"},
    a given species tree.
    }    
    \item{\code{spbranches}:}{Object of class \code{"data.frame"},
    a data frame listing branches of a given species tree.
    }
    \item{\code{orthoroot}:}{Object of class \code{"data.frame"},
    a data.frame with results from the 'groot' function (see return values in the OGR methods).
    }
    \item{\code{orthoct}:}{Object of class \code{"data.frame"},
    a data.frame with results from the 'groot.preprocess' function 
    (see return values in the OGR methods).
    }    
    \item{\code{status}:}{Object of class \code{"character"},
    a character value specifying the status of the OGR object 
    based on the available methods.
    }
  }
}

\section{Methods}{
  \describe{
    \item{groot}{\code{signature(object = "OGR")}: see  \code{\link[geneplast:groot]{groot}}}
    \item{groot.get}{\code{signature(object = "OGR")}: see  \code{\link[geneplast:groot.get]{groot.get}}} 
	 }
}

\author{
Rodrigo Dalmolin, Mauro Castro
}

\seealso{
  \code{\link[geneplast:groot.preprocess]{groot.preprocess}}
}

\keyword{classes}
