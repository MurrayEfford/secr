\name{writeGPS}
\alias{writeGPS}

\title{
Upload to GPS
}
\description{

Upload a set of point locations as waypoints to a GPS unit connected
by USB or via a serial port. Intended primarily for detector locations
in a traps object. Uses the GPSBabel package which must have been
installed. Coordinates are first inverse-projected to latitude and
longitude using function \code{\link[sf]{st_transform}} from \pkg{sf}.

}
\usage{

writeGPS(xy, o = "garmin", F = "usb:", proj = "+proj=nzmg")

}

\arguments{
  \item{xy}{2-column matrix or dataframe of x-y coordinates}
  \item{o}{character output format (see GPSBabel documentation)}
  \item{F}{character for destination (see Details)}
  \item{proj}{character string describing projection}
}
\details{
  
This function is derived in part from \code{readGPS} in \pkg{maptools}.

For users of Garmin GPS units, useful values of \code{o} are "garmin"
for direct upload via USB or serial ports, and "gdb" for a file in
Mapsource database format.

\code{F} may be "usb:" or "com4:" etc. for upload via USB or serial
ports, or the name of a file to create.

% 2020-04-05 https://www.nceas.ucsb.edu/scicomp/recipes/projections not working

The \code{proj} argument may be complex. For further information see the
Examples and the vignette 
\href{https://www.otago.ac.nz/density/pdfs/secr-spatialdata.pdf}{secr-spatialdata.pdf}.
If \code{proj} is an empty string then coordinates are assumed already to
be latitudes (column 1) and longitudes (column 2).

Waypoint names are derived from the rownames of \code{xy}.

}

\value{
  
No value is returned. The effect is to upload waypoints to an attached
  GPS or file. 

}

\seealso{

  \code{\link{make.systematic}}
 
}

\note{
  
  GPSBabel is available free from
\url{https://www.gpsbabel.org/}. Remember to add it to the Path. On
Windows this means following something like Settings > Control panel >
System > Advanced settings > Environment variables > (select Path) Edit
and adding ";C:/Program Files (x86)/gpsbabel" to the end (without the
quotes). Or ";C:/Program Files/gpsbabel" on 32-bit systems.
  
}

\examples{

## Example using shapefile "possumarea.shp" in
## "extdata" folder. As 'cluster' is not specified,
## the grid comprises single multi-catch detectors.

\dontrun{

## test for availability of GPSBabel

if (nzchar(Sys.which("gpsbabel"))) {

library(sf)
shpfilename <- system.file("extdata/possumarea.shp", package = "secr")
possumarea <- st_read(shpfilename)

possumgrid <- make.systematic(spacing = 100, region = possumarea, 
    plt = TRUE)

## May upload directly to GPS...
# writeGPS(possumgrid, proj = "+proj=nzmg")

## ...or save as Mapsource file
writeGPS(possumgrid, o = "gdb", F = "tempgrid.gdb",
    proj = "+proj=nzmg")

## If `region' had been specified in another projection we
## would need to specify this as in Proj.4. Here is a
## hypothetical example for New Zealand Transverse Mercator
## with datum NZGD2000 (EPSG:2193)

NZTM <- paste("+proj=tmerc +lat_0=0 +lon_0=173 +k=0.9996",
    "+x_0=1600000 +y_0=10000000 +ellps=GRS80",
    " +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")

# writeGPS(possumgridNZTM, o = "gdb", F = "tempNZTM.txt", 
#    proj = NZTM)

## Or to upload coordinates from UTM Zone 18 in eastern
## Maryland, USA...

# writeGPS(MarylandUTMgrid, proj = 
#    "+proj=utm +zone=18 +ellps=WGS84")

}

}

}

\keyword{ IO }


