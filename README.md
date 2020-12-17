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
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%204%20(SR).png)

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
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%204%20(sampling).png)

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
We next identify months as a variable, then create a data frame that accumulates the total number of days per month and assigns a number to each day within in a given year, then finally create variables for the Julian day and year.
```
months <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
            "Jul", "Aug", "Sep", "Oct","Nov", "Dec")
mos <- data.frame(month = fct_inorder(months), day=(cumsum(c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[1:12])+15))
Julian_Day <- y$julian
Year <- y$year
```
Lastly, we plot the point records across time based off of the year and Julian day of year assigned to that point record, and save the resulting plot as a PNG.
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
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%205.png)
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/MEOW_plot.png)
## 4. Taxonomic sampling trends
This analysis sets out to determine the presence of a phylogenetic signal in the sampling of seagrasses. First, we read in the needed packages.
```
library(ape) 
library(data.table) 
library(phyloregion) 
library(picante) 
library(phytools)
library(phylobase) 
library(adephylo) 
library(ggtree)
```
Then, we read in the downloaded seagrass occurrence data from GBIF as well as a dataset of seagrass phylogeny. These are then subjected to data manipulation, and phylogenetic data is matched with the point records in an attempt to track phylogenetic signal in species sampling.
```
d <- fread("/Users/darulab/Desktop/BriannaR/Review/Data/SeagrassGBIF/GBIF Seagrasses/seagrass occurrances jan_17 2020.csv", stringsAsFactors = FALSE)
viewd <- d[, c("genus", "species")]
d$species <- gsub(" ", "_", d$species)
tree <- read.nexus("/Users/darulab/Desktop/BriannaR/Review/Data/phylogeny/FINAL_tree.tre")
tree <- phylobuilder(d$species, tree = tree)
s <- data.frame(table(d$species))
z <- as.data.frame(s)
rownames(z) <- z[,1]
z[,1] <- NULL
index <- intersect(tree$tip.label, rownames(z))
tr1 <- keep.tip(tree, index)
M <- subset(z, rownames(z) %in% index)
subphy <- match.phylo.data(tree, M)$phy
subdat <- match.phylo.data(tree, M)$data
names(subdat)[1] <- "x"
subdat$x <- log(subdat$x+1)
```
Next, we visualize the matched data by plotting.
```
p <- ggtree(subphy, lwd=0.13, layout="fan", open.angle=180) + 
  geom_tiplab(aes(label=label, angle=angle), size=2, align=TRUE, linesize=.025, offset=100) + 
  theme_tree2()
pp <- (p + scale_y_continuous(expand=c(0, 0.93))) %>%
  gheatmap(subdat, offset=3, width=0.5, colnames=TRUE, colnames_angle=-45, hjust=-0.5, low="#F6FBF4", high="#004616", color = "grey", font.size=2) %>%
  scale_x_ggtree()
pp + theme(legend.position="right")
vp2 <- open_tree(pp, 180

# See below for resulting plot: 
```
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%206.png)
We also performed statistcal analyses (i.e., Moran's I, Abouheif's mean, Pagel's λ, and Blomberg's K) to test for a signifcant phylogenetic signal in point records.
```
# 1. Moran's I
# First, you need one file that combines both the phylogenetic and trait data
phylotraits <- phylo4d(subphy, subdat)
# Then, we do moran test using some Monte Carlo simulations (default is 999 randomizations)
moran.test <- abouheif.moran(phylotraits, method="Abouheif")
plot(moran.test)


# 2. Abouheif's mean
# First, you need one file that combines both the phylogenetic and trait data
phylotraits <- phylo4d(subphy, subdat)
abouheif.test <- abouheif.moran(phylotraits,method="oriAbouheif")
abouheif.test
plot(abouheif.test)


# 3. Pagel's λ
# First, you need to define which trait you want to test and give names to each value according to species
trait<-subdat[,1]
names(trait)<-rownames(subdat)
phylosig(subphy, trait, method="lambda", test=TRUE, nsim=999)


# 4. Blomberg's K
# First, you need to define which trait you want to test and give names to each value according to species
trait<-subdat[,1]
names(trait)<-rownames(subdat)
#Then, we do the test with 999 randomizations:
phylosig(subphy, trait, method="K", test=TRUE, nsim=999)
```


## 5. Family rank correlations
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%20MEOWs_seagrass_families3.png)
## 6. Extinction risk
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%203.png)
## 7. Quantifying increase of sampling
![alt text](https://github.com/brirock35/Impediments-to-Understanding-Seagrasses-Response-to-Global-Change/blob/main/Figure%209.png)

## References
## Acknowledgments
