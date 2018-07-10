\name{FFT2FS_3D}
\alias{FFT2FS_3D}
\title{
  Three dimensional Fourier Series
}
\description{
  Calculate the three dimensional Fourier series coeffiients of the input three dimensional array.
}
\usage{
FFT2FS_3D(A)
}

\arguments{
  \item{A}{
A three dimensional numerical data array.  For example, \emph{A} can be the data array of an three dimensional image.
}
}
\details{
Calcualte the three dimensional numerical data array.  The input \emph{A} array can be any three dimensional data array.  For image input data, the input should be data array only without any header information.
}
\value{
A three dimensional Fourier series coefficients array of the input \emph{A} data array.
}
\references{
Lin N, Jiang J, Guo S, Xiong M. Functional Principal Component Analysis and Randomized Sparse Clustering Algorithm for Medical Image Analysis. PLOS ONE. 2015;10(7):e0132945.
}

\examples{
test_data = array(runif(1000,0,1),dim = c(10,10,10))
rlt = FFT2FS_3D(test_data)
}

