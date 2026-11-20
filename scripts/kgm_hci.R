rm(list = ls())
library(lubridate)
library(fields)
library(ncdf4)
library(viridis)
library(cmocean)

load('c:/Users/brendan.turley/Documents/R_projects/KGM-environ-sensitivity/data/kmk_smdat.RData')
quantile(kmk_pos$TEMPSURF,c(.025,.975),na.rm=T)


setwd("c:/Users/brendan.turley/Documents/R_projects/esr-indicator-scratch/data/intermediate_files")

sst <- readRDS('ct_sst_2025')
i <- 300
sst$time[i]
imagePlot(sst$sst[,, i])
contour(sst$sst[,, i], levels = c(20, 30), add = T)

f_ls <- list.files()
fi <- which(grepl('ct_sst', f_ls))[-c(42:44)]
hci_ts <- matrix(NA, length(fi), 366)
hci_mts <- array(NA, dim = c(41, dim(sst$sst)[1:2]))
for (f in fi) {
  print(which(fi == f) / length(fi))
  sst <- readRDS(f_ls[f])
  ### matrix of HCI time series for each file
  hr1 <- apply(sst$sst, 3, function(x) sum(which(x > 20 & x < 30)))
  ht1 <- apply(sst$sst, 3, function(x) sum(which(!is.na(x))))
  hci_ts[which(fi == f), 1:length(hr1)] <- hr1 / ht1

  hr2 <- apply(sst$sst, c(1, 2), function(x) sum(which(x > 20 & x < 30)))
  ht2 <- apply(sst$sst, c(1, 2), function(x) sum(which(!is.na(x))))
  hci_mts[which(fi == f), , ] <- hr2 / ht2
}

# cols <- viridis(nrow(hci_ts))
cols <- cmocean::cmocean('thermal')(nrow(hci_ts))
matplot(
  t(hci_ts),
  type = 'l',
  lty = 1,
  col = cols,
  xlab = 'Time',
  ylab = 'HCI',
  main = 'HCI Time Series'
)

plot(1985:2025, apply(hci_ts, 1, mean, na.rm = TRUE), type = 'b')
plot(1985:2025, apply(hci_ts, 1, sd, na.rm = TRUE), type = 'b')

xs <- c(1:366, 366:1)
yiqr <- c(
  apply(hci_ts, 2, quantile, probs = 0.25, na.rm = TRUE),
  rev(apply(hci_ts, 2, quantile, probs = 0.75, na.rm = TRUE))
)
ymx <- c(
  apply(hci_ts, 2, min, na.rm = TRUE),
  rev(apply(hci_ts, 2, max, na.rm = TRUE))
)
plot(1:366, apply(hci_ts, 2, median, na.rm = TRUE), type = 'l', ylim = c(0, 1))
polygon(xs, ymx, col = rgb(1, 0, 0, 0.3), border = NA)
polygon(xs, yiqr, col = rgb(0, 0, 1, 0.3), border = NA)

plot(1985:2025, apply(hci_ts, 1, function(x) length(which(x == 1))), type = 'b')
plot(1985:2025, apply(hci_ts, 1, function(x) which(x == 1)[1]), type = 'b')

# extract linear slope from matrix of HCI time series
slopes <- apply(hci_ts, 2, function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  lm(x ~ c(1:length(x))) |> summary() |> coef()
})
slopes[2, which(slopes[8, ] < .05)] |> hist()
plot(
  1:366,
  slopes[2, ],
  type = 'l',
  xlab = 'Year',
  ylab = 'Slope',
  main = 'Linear Slope of HCI Time Series'
)
abline(h = 0, col = 'red', lty = 2)
points(
  c(1:366)[which(slopes[8, ] < .05)],
  slopes[2, which(slopes[8, ] < .05)],
  col = 'blue',
  pch = 16
)


slopes_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  coef(lm(x ~ c(1:length(x))))[2]
}) ### this should return a matrix of 360x260

imagePlot(slopes_m)


### terra

library(terra)

burl <- 'https://www.ngdc.noaa.gov/thredds/dodsC/global/ETOPO2022/60s/60s_bed_elev_netcdf/ETOPO_2022_v1_60s_N90W180_bed.nc'
bdat <- nc_open(burl)
# crs <- 'EPSG:4326'

ln <- ncvar_get(bdat, 'lon')
ln_i <- which(ln>=min_lon & ln<=max_lon)
lt <- ncvar_get(bdat, 'lat')
lt_i <- which(lt>=min_lat & lt<=max_lat)

bathy <- ncvar_get(bdat, 'z', 
                   start = c(ln_i[1],lt_i[1]),
                   count = c(length(ln_i), length(lt_i)))
# bathy <- ncvar_get(bdat, 'Band1', 
#                    start = c(ln_i[1],lt_i[1]),
#                    count = c(length(ln_i), length(lt_i)))
nc_close(bdat)

bathy[which(bathy > -10 | bathy < -100)] <- NA
bathy[which(!is.na(bathy))] <- 1

bathy_r <- rast(
  t(bathy[,ncol(bathy):1]),
  crs = "EPSG:4326"
) #|> as.polygons()
ext(bathy_r) <- c(min(ln[ln_i]), max(ln[ln_i]), min(lt[lt_i]), max(lt[lt_i]))
bathy <- as.polygons(bathy_r)

setwd("~/data/shapefiles/gulf_iho")
iho <- vect('iho.shp') |> makeValid()

gulf_shelf <- terra::intersect(bathy, iho)
plot(gulf_shelf)

# define years  --------------------------------
styear <- 1985
enyear <- 2025

# define spatial domain  --------------------------------
min_lon <- -98
max_lon <- -80
min_lat <- 18
max_lat <- 31

url <- 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/crw/5km/v3.1/nc/v1.0/daily/sst/1985/coraltemp_v3.1_19850101.nc'
dat <- nc_open(url)

lon <- ncvar_get(dat, 'lon')
lat <- ncvar_get(dat, 'lat')

i_lon <- which(lon >= (min_lon) & lon <= (max_lon))
i_lat <- which(lat <= (max_lat) & lat >= (min_lat))
lons <- lon[i_lon]
lats <- lat[i_lat]

setwd("c:/Users/brendan.turley/Documents/R_projects/esr-indicator-scratch/data/intermediate_files")

# raster_list <- list()
hci_ts2 <- matrix(NA, length(styear:enyear), 366)
for (i in styear:enyear) {
  cat(i, '\n')
  tmp <- paste0('ct_sst_', i) |> readRDS()
  # sst_r <- rast(tmp$sst, crs="EPSG:4326")
  sst_r <- rast(aperm(tmp$sst[, 260:1, ], c(2, 1, 3)), crs = "EPSG:4326")
  ext(sst_r) <- c(range(lons), range(lats))
  # sst_r <- project(sst_r, "EPSG:4326")
  time(sst_r) <- as.Date(tmp$time)

  sst_r <- crop(sst_r, gulf_shelf) |> mask(gulf_shelf)

  ### find number of cells per time step x > 20 & x < 30 in raster sst_r
  hot_cells <- sst_r > 20 & sst_r < 30
  valid_cells <- !is.na(sst_r)
  ### sum per time step
  hot_count <- global(hot_cells, sum, na.rm = TRUE)
  valid_count <- global(valid_cells, sum, na.rm = TRUE)
  hci_ts2[which(i == (styear:enyear)), 1:nrow(hot_count)] <- unlist(hot_count / valid_count)

  # Store in list
  # raster_list[[which(i == (styear:enyear))]] <- sst_r
}
# sst_r <- rast(raster_list)

### decadal medians
yrs <- styear:enyear
decadal_medians <- apply(hci_ts2, 2, function(x) {
  tapply(x, (styear:enyear), median, na.rm = TRUE)
})
yr_80 <- subset(hci_ts2, yrs<1990) |>
  apply(2, median, na.rm = TRUE)
yr_90 <- subset(hci_ts2, yrs>=1990 & yrs<2000) |>
  apply(2, median, na.rm = TRUE)
yr_00 <- subset(hci_ts2, yrs>=2000 & yrs<2010) |>
  apply(2, median, na.rm = TRUE)
yr_10 <- subset(hci_ts2, yrs>=2010 & yrs<2020) |>
  apply(2, median, na.rm = TRUE)
yr_20 <- subset(hci_ts2, yrs>=2020) |>
  apply(2, median, na.rm = TRUE)
dcols <- cmocean::cmocean('amp')(5)

jday <- data.frame(date = seq(as.Date("2020-01-01"), as.Date("2020-12-31"), by = "day"))
# cols <- viridis(nrow(hci_ts))
cols <- cmocean::cmocean('thermal')(nrow(hci_ts2))
matplot(jday$date,
  t(hci_ts2),
  type = 'l',
  lty = 1,
  col = cols,
  xlab = 'Time',
  ylab = 'HCI',
  main = 'HCI Time Series',
  xaxt = 'n'
)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")

plot(1985:2025, apply(hci_ts2, 1, mean, na.rm = TRUE), type = 'b')
plot(1985:2025, apply(hci_ts2, 1, sd, na.rm = TRUE), type = 'b')

xs <- c(jday$date, rev(jday$date))
yiqr <- c(
  apply(hci_ts2, 2, quantile, probs = 0.25, na.rm = TRUE),
  rev(apply(hci_ts2, 2, quantile, probs = 0.75, na.rm = TRUE))
)
ymx <- c(
  apply(hci_ts2, 2, min, na.rm = TRUE),
  rev(apply(hci_ts2, 2, max, na.rm = TRUE))
)

plot(jday$date, apply(hci_ts2, 2, median, na.rm = TRUE), type = 'l', ylim = c(0, 1), xaxt = 'n', lwd = 2)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
polygon(xs, ymx, col = rgb(1, 0, 0, 0.3), border = NA)
polygon(xs, yiqr, col = rgb(0, 0, 1, 0.3), border = NA)
points(jday$date, yr_80, type = 'l', col = dcols[1], lwd = 2)
points(jday$date, yr_90, type = 'l', col = dcols[2], lwd = 2)
points(jday$date, yr_00, type = 'l', col = dcols[3], lwd = 2)
points(jday$date, yr_10, type = 'l', col = dcols[4], lwd = 2)
points(jday$date, yr_20, type = 'l', col = dcols[5], lwd = 2)

plot(1985:2025, apply(hci_ts2, 1, function(x) length(which(x == 1))), type = 'b')
plot(1985:2025, apply(hci_ts2, 1, function(x) which(x == 1)[1]), type = 'b')

# extract linear slope from matrix of HCI time series
slopes <- apply(hci_ts2, 2, function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  lm(x ~ c(1:length(x))) |> summary() |> coef()
})
# slopes[2, which(slopes[8, ] < .05)] |> hist()
plot(
  jday$date, 
  slopes[2, ] * 41,
  xlab = 'Year',
  ylab = 'Slope',
  main = 'Linear Slope of HCI Time Series',
  type = 'h', 
  xaxt = 'n'
)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
abline(h = 0, col = 'red', lty = 2, lwd = 2)
points(
  jday$date[which(slopes[8, ] < .05)],
  slopes[2, which(slopes[8, ] < .05)] * 41,
  col = 'blue',
  pch = 16
)


plot(jday$date, apply(hci_ts2, 2, median, na.rm = TRUE), type = 'l', ylim = c(0, 1.2), xaxt = 'n', lwd = 2)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
# polygon(xs, ymx, col = rgb(1, 0, 0, 0.3), border = NA)
# polygon(xs, yiqr, col = rgb(0, 0, 1, 0.3), border = NA)
arrows(jday$date, 
  apply(hci_ts2, 2, median, na.rm = TRUE), 
  jday$date, 
  apply(hci_ts2, 2, median, na.rm = TRUE) + slopes[2, ] * 41, 
  length = 0.02, col = 1)

