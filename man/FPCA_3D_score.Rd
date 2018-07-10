\name{FPCA_3D_score}
\alias{FPCA_3D_score}
\title{
Three Dimensional Functional Component Analysis
}
\description{
Calculation of three dimensional functional principal component scores for a series of three dimensional array data.
}
\usage{
FPCA_3D_score(X, prop)
}
\arguments{
  \item{X}{
The input data array.  \emph{X} is a four dimensional data array.  The first three dimensional data represents the three dimensional data array for each observation.  The fourth dimention represents the observations.
}
  \item{prop}{
The prespecified proportion of variance the calcuatled functional principal component scores can explain in the functional domain.
}
}
\details{
Calculate the three dimensional functional principal component scores for a series of three dimensional data.
}
\value{
A two dimensional score matrix.  The row of the score matrix represents each individual and the column of the score matrix represent each component score.
}
\references{
Lin N, Jiang J, Guo S, Xiong M. Functional Principal Component Analysis and Randomized Sparse Clustering Algorithm for Medical Image Analysis. PLOS ONE. 2015;10(7):e0132945.
}

\examples{
data_in = array(runif(4000,0,1),dim=c(10,10,10,4))
test = FPCA_3D_score(data_in,0.8)
}

