# Impediments to Understanding Seagrasses Response to Global Change
## Brianna M. Rock & Barnabas H. Daru
#### Department of Life Sciences, Daru Lab, Texas A&M University-Corpus Christi, Corpus Christi, Texas, USA 

![alt text](https://barnabasdaru.files.wordpress.com/2018/05/title-picture.png?w=300&h=300g)

## 1. Spatial sampling trends

First, download the needed packages. 
```
rm(list = ls())
library(data.table)
library(rgdal)
library(devtools)
library(bioregion)
library(raster)
library(sp)
library(ggplot2)
library(maptools)
library(phytools)
data("wrld_simpl")
```
Then read in the IUCN spatial data and a file with the desired grid size. These should also be set to an approprite coordinate system.
```
ss <- readOGR(dsn = "/Users/darulab/Desktop/BriannaR/Review/Data/ShapeFiles/Seagrasses_SHP_raw/", layer = "SEAGRASSES")
proj4string(ss) <- CRS("+proj=longlat +datum=WGS84")
S <- as.character(unique(ss$binomial))
s1 <- readOGR(dsn = "/Users/darulab/Desktop/BriannaR/Review/Data/ShapeFiles/grids", layer = "grids_100km")
proj4string(s1) <- CRS("+proj=longlat +datum=WGS84")

```
To plot seagrass species richness from the IUCN polygons, we utilized the following for-loop which determines the number of species of seagrasses occuring in each grid cell.
```
out <- NULL
for (i in 1:length(S)) {
  x1 <- subset(ss, ss$binomial %in% S[i])
  proj4string(x1) <- proj4string(s1)
  x2 <- over(s1, x1)
  x3 <- cbind(data.frame(s1), x2)
  M <- as.data.frame.matrix(table(x3$grids, x3$binomial))
  M[M>0] <- 1
  M1 <- picante::matrix2sample(M)
  w <- M1[,c(1,3)]
  names(w) <- c("grids", "species")
  out <- rbind(out, w)
  print(i)
}
r <- data.frame(out)
```
Next, we write the output to file, and read it in for subsequent analysis.
```
write.csv(r, "/Users/darulab/Desktop/Brianna R (SPD)/Data/CSVs/PRESAB_100km.csv", row.names = FALSE)
r <-fread("/Users/darulab/Desktop/BriannaR/Review/Data/CSVs/PRESAB_100km.csv")
```
To plot in geographic space, we convert the data to a dataframe, layer the grid cells, and select for species richness (SR) while replacing NAs with 0.
```
mm <- data.frame(table(r$grids))
names(mm) <- c("grids", "SR")
index1 <- match(s1$grids, mm$grids)
zm <- cbind(s1, mm$SR[index1])
names(zm)[3] <- "SR"
zd <- as.data.frame(zm)
zd[is.na(zd)] <- 0
zm1 <- zd[zd$SR>0, ]
k=10
COLOUR <- hcl.colors(k, palette = "Zissou 1")
y = choropleth(zm1, values=zm1$SR, k)
plot(wrld_simpl, col="white", border="grey")
plot(y, layer="SR", col=COLOUR[y$values], border = NA, add=T)
add.color.bar(cols = COLOUR, lims=c(1,22), digits=1, prompt=TRUE,title = NA,
               lwd=4, outline=TRUE)

# See below for resulting plot:
```
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure4a_seagrass_sampling_density.png)

To produce a plot representing the density of seagrass sampling by species, we first import the dataset we created during the previous analysis.
```
dat <- fread("/Users/darulab/Desktop/BriannaR/Review/Data/CSVs/PRESAB_100km.csv")
```
We then proceed to edit the dataset for the desired traits, in this case the number of point records for each seagrass species, and then combine these data with a data layer consistng of 100km grid cells.
```
dd <- data.frame(table(dat$species))
names(dd) <- c("species", "range")
trait_range <- dd$range
names(trait_range) <- dd$species
gg <- mapTraits(dat, trait = trait_range, FUN = sd)
colnames(gg) <- c("grids","traits")
s1 <- readOGR(dsn = "/Users/darulab/Desktop/BriannaR/Review/Data/ShapeFiles/grids", layer = "grids_100km")
proj4string(s1) <- CRS("+proj=longlat +datum=WGS84")
index <- match(s1$grids, gg$grids)
z <- cbind(s1, gg$traits[index])
names(z)[3] <- "traits"
z <- z[complete.cases(z@data$traits),]
z1 <- z[z@data$traits>0, ]
```
Finally, we plot the data in geographic space.
```
k=10
COLOUR <- hcl.colors(k, palette = "Zissou 1")
y = choropleth(z1, values=z1$traits, k)
data("wrld_simpl")
plot(wrld_simpl, fill=NA, border="grey")
plot(z1, col=COLOUR[y$values], border = NA, add=TRUE)
add.color.bar(cols = COLOUR, lims=c(1,532), digits=1, prompt=TRUE,title = NA,
               lwd=4, outline=TRUE)

# See below for resulting plot:
```
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure4a_seagrass_sampling_density.png)

## 2. Temporal increase in sampling
The goal of this analysis is to visualize how the collection of seagrass records has changed across time. First, we begin by reading in the needed packages.
```
library(lubridate)
library(circular)
library(tidyverse)
library(phyloregion)
library(raster)
library(data.table)
library(colorRamps)
```
Then, we read in the datasets needed for the analysis, and perform data manipulation to isolate temporal information provided in the point records and convert these to Julian day of year format.
```
w <- shapefile("/Users/darulab/Desktop/BriannaR/Review/Data/MEOW/meow_dissolved.shp")
d <- fread("/Users/darulab/Desktop/BriannaR/Review/Data/SeagrassGBIF/GBIF Seagrasses/seagrass occurrances jan_17 2020.csv")
d <- d[, c("family", "species", "decimalLongitude", "decimalLatitude", "year", "month", "day")]
names(d)[c(3,4)] <- c("lon", "lat")
d <- d[complete.cases(d),]
d <- d %>% mutate(julian = as.Date(paste(year, month, day, sep="-")))
d <- d %>% mutate(julian = yday(as.Date(paste(year, month, day, sep = "-"))))
d <- d[!is.na(d$julian),]
coordinates(d) <- ~ lon+lat
proj4string(w) <- proj4string(d)
x <- over(d, w)
y <- cbind(as.data.table(d), x)
y <- y[complete.cases(y),]
```
We next identify month names as its own variable, then create a data frame that assigns a number to each day of the month in a given year, and finally create variables for the Julian day data and years within the data.
```
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
            "Jul", "Aug", "Sep", "Oct","Nov", "Dec")
mos <- data.frame(month = fct_inorder(months), day=(cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[1:12])+15))
Julian_Day <- y$julian
Year <- y$year
```
Finally, we plot the point records across time based off of the year and Julian day of year assigned to that point record, and save the resulting plot as a PNG.
```
p <- ggplot() +
  geom_point(data=y, aes(Julian_Day, Year, color=year), size=1) +
  geom_text(data=mos, aes(day, 2019, label=month), size=7, angle = 45) +
  ylim(1700, NA) +
  theme_minimal()+
  scale_color_gradientn(colours=rev(blue2red(10)))+
  theme(axis.text = element_text(colour = "black", size = 20))

# save as .png
png("/Users/darulab/Desktop/temp2.png", width=8, height=12, units = "in", res=300)
p
dev.off()

# See below for resulting plot: 
```
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%204.png)

## 3. Temporal sampling trends 
## 4. Taxonomic sampling trends
## 5. Family rank correlations
## 6. Extinction risk
## 7. Quantifying increase of sampling
## References
## Acknowledgments
