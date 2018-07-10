\name{FPCA3D-package}
\alias{FPCA3D-package}
\alias{FPCA3D}
\docType{package}
\title{
\packageTitle{FPCA3D}
}
\description{
\packageDescription{FPCA3D}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{FPCA3D}
\packageIndices{FPCA3D}
data_in = array(runif(4000,0,1),dim=c(10,10,10,4))
test = FPCA_3D_score(data_in,0.8)
}
\author{
\packageAuthor{FPCA3D}

Maintainer: \packageMaintainer{FPCA3D}
}
\references{
Lin N, Jiang J, Guo S, Xiong M. Functional Principal Component Analysis and Randomized Sparse Clustering Algorithm for Medical Image Analysis. PLOS ONE. 2015;10(7):e0132945.
}
\keyword{ package }
\seealso{
\code{\link{FFT2FS_3D}}, \code{\link{FPCA_3D_score}}
}
\examples{
data_in = array(runif(4000,0,1),dim=c(10,10,10,4))
test = FPCA_3D_score(data_in,0.8)
}
