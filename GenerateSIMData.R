

## ------------------------------------------------------------------------------------ I.Set simulation parameters ----

habitat.resolution <- 1
buffer <- 0
grid.size <- 50
detector.resolutions <- 1.2
# detector.divisions <- 1
detector.resolutions.cam <- 1
N <- 50
a0 <- 1
a0.cam <- 0.75
sigma <- 1
aug.factor <- 1
a <- 0.5

## ------------------------------------------------------------------------------------ II.Load required libraries ----- 

library(raster)
library(sp)
library(rgdal)
library(rgeos)

## ------------------------------------------------------------------------------------ III.Generate data set ----------

### ==== 1.CREATE A SQUARE SPATIAL DOMAIN WHERE DETECTORS WILL BE PLACED ====

coords <- matrix(c(0, 0,
                   grid.size, 0,
                   grid.size, grid.size,
                   0, grid.size,
                   0, 0), ncol = 2, byrow = TRUE)
P1 <- Polygon(coords)
myStudyArea <- SpatialPolygons(list(Polygons(list(P1), ID = "a")),
                               proj4string=CRS("+proj=utm +zone=33 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

### ==== 2. HABITAT OBJECTS ====

r <- raster(nrow=50, ncol=50, xmn=0, xmx=grid.size,
            ymn=0, ymx=grid.size)

r[] <- 1  ##--we assume a homogeneous habitat quality

IDCells.r <- r
IDCells.r[] <- 1:length(IDCells.r)							         
IDCells.mx <- as.matrix(IDCells.r)

habitat.xy <- xyFromCell(r, 1:ncell(r))
dimnames(habitat.xy) <- list(1:length(habitat.xy[,"x"]), c("x","y"))
habitat.sp <- SpatialPointsDataFrame(data.frame(habitat.xy[,c("x","y")]), data=data.frame(habitat.xy), 
                                     proj4string=CRS(projection(myStudyArea)))
x.max <- round(max(habitat.sp$x))
y.max <- round(max(habitat.sp$y))

### ==== 3.GENERATE Primary DETECTORS: NGS ====

co <- seq(0,grid.size, by = detector.resolutions)
x <- rep(co, length(co))

y <- sort(rep(co, length(co)), decreasing = T)
detectors.xy <- cbind(x, y)
detectors.sp <- SpatialPoints(detectors.xy,
                              proj4string = CRS(proj4string(myStudyArea)))

### ==== 4.GENERATE Secondary DETECTORS: CAM ====
co <- seq(1,grid.size, by = detector.resolutions.cam)
x <- rep(co, length(co))

y <- sort(rep(co, length(co)), decreasing = T)
detectors.xy.cam <- cbind(x, y)
detectors.sp.cam <- SpatialPoints(detectors.xy.cam,
                                  proj4string = CRS(proj4string(myStudyArea)))

### ==== 5.SIMULATE INDIVIDUAL AC LOCATIONS ====

mu <- rep(1,length(habitat.sp))
temp <- as.numeric(rmultinom(1, N, mu))
index <- rep(1:length(temp), temp)
AC.sp <- habitat.sp[index, ]

### ==== 6.SIMULATE DETECTIONS BY NGS ====

D <- gDistance(detectors.sp, AC.sp, byid=TRUE)

pZero <- 1-exp(-(a0/(2*pi*sigma^2)))              ##--equation (8) in the manuscript
p <- pZero * exp(-D*D/(2*sigma*sigma))

y <- apply(p, c(1,2), function(x) rbinom(1, 1, x))
table(y)

### ==== 7.SIMULATE DETECTIONS BY CAM ====

D.cam <- gDistance(detectors.sp.cam, AC.sp, byid = TRUE)

pZero.cam <- 1-exp(-(a0.cam/(2*pi*sigma^2)))     ##--equation (8) in the manuscript 
p.cam <- pZero.cam * exp(-D.cam*D.cam/(2*sigma*sigma))

y.cam <- apply(p.cam, c(1,2), function(x) rbinom(1, 1, x))
table(y.cam)

### ==== 8.GENERATE THE DIFFERENT TYPES OF DETECTIONS ====

y.ORIGINAL.ngs <- y
identifiedDets <- apply(y.ORIGINAL.ngs, c(1,2), function(x)rbinom(1, 1, a))

y.IDENTIFIED.ngs <- ifelse(y.ORIGINAL.ngs == 1 & identifiedDets == 1, 1, 0)
y.POOLED.ngs <- as.numeric(apply(y.ORIGINAL.ngs, 2, function(x)any(x >= 1)))
y.UNIDENTIFIED.ngs <- ifelse(y.ORIGINAL.ngs == 1 & identifiedDets == 0, 1, 0)
y.UNIDENTIFIED.ngs <- as.numeric(apply(y.UNIDENTIFIED.ngs, 2, function(x)any(x >= 1)))  

y.POOLED.cam <- as.numeric(apply(y.cam, 2, function(x)any(x >= 1)))


### ==== 9.AUGMENT DATA ====

dimnames(y.IDENTIFIED.ngs) <- list(1:dim(y.IDENTIFIED.ngs)[1], 1:dim(y.IDENTIFIED.ngs)[2])
n.tot <- round(dim(y.IDENTIFIED.ngs)[1]*(1 + aug.factor))
y.aug <- matrix(NA, n.tot, dim(y.IDENTIFIED.ngs)[2])
y.aug[1:dim(y.IDENTIFIED.ngs)[1], ] <- y.IDENTIFIED.ngs
dimnames(y.aug) <- list(c( dimnames(y.IDENTIFIED.ngs)[[1]], rep("Augmented", n.tot - dim(y.IDENTIFIED.ngs)[1])),
                        dimnames(y.IDENTIFIED.ngs)[[2]])

y <- y.aug

### ==== 10.CONSTRUCT z ====

z <- z.init <- rep(NA, dim(y)[1])
detected <- apply(y, 1, function(x)any(x >= 1))
z[detected] <- 1
z.init[!detected] <- 1
