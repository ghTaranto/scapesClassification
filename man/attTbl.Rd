% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/formatInputs.R
\encoding{UTF-8}
\name{attTbl}
\alias{attTbl}
\title{Attribute table}
\usage{
attTbl(r, var_names = NULL)
}
\arguments{
\item{r}{single or multi-layer raster of the class \code{SpatRaster} (see
\code{help("rast", terra)}).}

\item{var_names}{character vector, raster layers' names in the attribute
table. If \code{NULL}, then the original layers' names are used.}
}
\value{
data.frame
}
\description{
Converts a single or a multi-layer raster into an attribute table
(\code{data.frame}).
}
\details{
Attribute tables come with a column named \strong{\code{"Cell"}} which
stores raster cell numbers and associate each row of the attribute table
with a cell of the raster object. The remaining columns of the attribute
table store the data contained in the raster layers. Note that only raster
cells having no missing value in no layer (\strong{complete cases}) are included
in the attribute table.
}
\note{
\strong{Attribute table contains only complete cases}, i.e., raster cells
having a value for every layer of the stack.
}
\examples{
library(scapesClassification)
library(terra)

## CREATE A DUMMY RASTER ##
r <- terra::rast(matrix(c(NA,100,100,NA,100,100,0,0,0),
                        nrow = 3,
                        ncol = 3,
                        byrow = TRUE))

## RASTER CELL NUMBERS ##
rcn <- r; rcn[] <- 1:9

## PLOT DATA AND CELL NUMBERS ##
oldpar <- par(mfrow = c(1,2))
m <- c(4, 1, 4, 1)

plot(r, col="grey90", colNA="red3", mar=m, asp = NA, axes=FALSE, legend=FALSE)
text(r)
lines(r)
mtext(side=3, line=0.2, adj=0, cex=1.5, font=2, "Dummy_var")
legend("bottomright", ncol=1, bg="white", fill=c("red3"),
       legend = c("NA cells (1 and 4)"))

plot(rcn, col="grey90", mar=m, asp=NA, axes=FALSE, legend=FALSE)
text(rcn)
lines(rcn)
mtext(side=3, line=0.2, adj=0, cex=1.5, font=2, "Cell numbers")
par(oldpar)

## VISUALIZE ATTRIBUTE TABLE ##

at <- attTbl(r, var_names = c("dummy_var"))
at

# Note that cells 1 and 4 have missing values and therefore are not included in the table
any(at$Cell \%in\% c(1,4))
}
