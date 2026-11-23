rm(list = ls())
library(lubridate)
library(fields)
library(ncdf4)
library(viridis)
library(cmocean)
library(terra)

load('c:/Users/brendan.turley/Documents/R_projects/KGM-environ-sensitivity/data/kmk_smdat.RData')
quantile(kmk_pos$TEMPSURF,c(.025,.975),na.rm=T)
# lower <- round(quantile(kmk_pos$TEMPSURF,c(.025,.975),na.rm=T)[1])
### but fishing forums and articles state 68-70 is the lower limit 
# https://www.floridasportsman.com/editorial/tracking-the-spring-king-mackerel-run/545891
# https://oldsaltfishing.org/kingfishing-in-tampa-bay/
lower <- 20
upper <- round(quantile(kmk_pos$TEMPSURF,c(.025,.975),na.rm=T)[2])

setwd("c:/Users/brendan.turley/Documents/R_projects/esr-indicator-scratch/data/intermediate_files")
sst <- readRDS('ct_sst_2025')
i <- 300
sst$time[i]
imagePlot(sst$sst[,, i])
contour(sst$sst[,, i], levels = c(lower, upper), add = T)

f_ls <- list.files()
fi <- which(grepl('ct_sst', f_ls))[-c(42:44)]
hci_ts <- matrix(NA, length(fi), 366)
hci_mts <- array(NA, dim = c(41, dim(sst$sst)[1:2]))
for (f in fi) {
  print(which(fi == f) / length(fi))
  sst <- readRDS(f_ls[f])
  ### matrix of HCI time series for each file
  hr1 <- apply(sst$sst, 3, function(x) sum(which(x >= lower & x <= upper)))
  ht1 <- apply(sst$sst, 3, function(x) sum(which(!is.na(x))))
  hci_ts[which(fi == f), 1:length(hr1)] <- hr1 / ht1

  hr2 <- apply(sst$sst, c(1, 2), function(x) sum(which(x >= lower & x <= upper)))
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
# slopes[2, which(slopes[8, ] <= .05)] |> hist()

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
  c(1:366)[which(slopes[8, ] <= .05)],
  slopes[2, which(slopes[8, ] <= .05)],
  col = 'blue',
  pch = 16
)


slopes_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  coef(lm(x ~ c(1:length(x))))[2]
}) ### this should return a matrix of 360x260
hist(slopes_m, main = 'Histogram of Slopes', xlab = 'Slope', breaks = 30)
range(slopes_m, na.rm = TRUE)
b_brks <- seq(-.009, .009, by = 0.0006) #length.out = 30)
b_col <- cmocean::cmocean('balance')(length(b_brks)-1)


min_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  min(x, na.rm = TRUE)
}) ### this should return a matrix of 360x260
hist(min_m, main = 'Histogram of Min', xlab = 'HCI', breaks = 30)
range(min_m, na.rm = TRUE)
m_brks <- seq(.4, 1, by = .02) #length.out = 30)
m_col <- cmocean::cmocean('matter')(length(m_brks)-1)


max_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  max(x, na.rm = TRUE)
}) ### this should return a matrix of 360x260
hist(max_m, main = 'Histogram of Max', xlab = 'HCI', breaks = 30)
range(max_m, na.rm = TRUE)
mx_brks <- seq(.75, 1, by = .01) #length.out = 30)
mx_col <- cmocean::cmocean('speed')(length(mx_brks)-1)


mean_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  mean(x, na.rm = TRUE)
}) ### this should return a matrix of 360x260
hist(mean_m, main = 'Histogram of Mean', xlab = 'HCI', breaks = 30)
range(mean_m, na.rm = TRUE)
mn_brks <- seq(.55, 1, by = .05) #length.out = 30)
mn_col <- cmocean::cmocean('rain')(length(mn_brks)-1)


sd_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  sd(x, na.rm = TRUE)
}) ### this should return a matrix of 360x260
hist(sd_m, main = 'Histogram of SD', xlab = 'HCI', breaks = 30)
range(sd_m, na.rm = TRUE)
sd_brks <- seq(0, 0.14, by = 0.01) #length.out = 30)
sd_col <- cmocean::cmocean('dense')(length(sd_brks)-1)


setwd("c:/Users/brendan.turley/Documents/R_projects/kgm-environ-sensitivity/figures")
png('hci_plots.png', width = 7, height = 9, units = 'in', res = 300)
par(mfrow = c(3, 2), mar = c(4, 4, 2, 2))
imagePlot(slopes_m, breaks = b_brks, col = b_col, asp = 1)
mtext('Slope of HCI')
contour(slopes_m, add = TRUE, col = 'gray20', levels = 0)

imagePlot(min_m, breaks = m_brks, col = m_col, asp = 1)
mtext('Min of HCI')
contour(min_m, add = TRUE, col = 'gray20', levels = seq(.4, 1, by = .1))

imagePlot(max_m, breaks = mx_brks, col = mx_col, asp = 1)
mtext('Max of HCI')
contour(max_m, add = TRUE, col = 'gray50', levels = seq(.75, 1, by = .05))

imagePlot(mean_m, breaks = mn_brks, col = mn_col, asp = 1)
mtext('Mean of HCI')
contour(mean_m, add = TRUE, col = 'gray50', levels = seq(.55, 1, by = .05))

imagePlot(sd_m, breaks = sd_brks, col = sd_col, asp = 1)
mtext('SD of HCI')
contour(sd_m, add = TRUE, col = 'gray20', levels = seq(0, 0.14, by = 0.01))
dev.off()



median_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  median(x, na.rm = TRUE)
}) ### this should return a matrix of 360x260
hist(median_m, main = 'Histogram of Median', xlab = 'HCI', breaks = 30)
range(median_m, na.rm = TRUE)
med_brks <- seq(.55, 1, by = .05) #length.out = 30)
med_col <- cmocean::cmocean('rain')(length(med_brks)-1)


iqr_m <- apply(hci_mts, c(2, 3), function(x) {
  if (all(is.na(x))) {
    return(NA)
  }
  IQR(x, na.rm = TRUE)
}) ### this should return a matrix of 360x260
hist(iqr_m, main = 'Histogram of IQR', xlab = 'HCI', breaks = 30)
range(iqr_m, na.rm = TRUE)
iqr_brks <- seq(0, 0.25, by = 0.01) #length.out = 30)
iqr_col <- cmocean::cmocean('dense')(length(iqr_brks)-1)

par(mfrow = c(2, 2), mar = c(4, 4, 2, 2))
imagePlot(mean_m, breaks = mn_brks, col = mn_col)
mtext('Mean of HCI')
contour(mean_m, add = TRUE, col = 'gray50', levels = seq(.55, 1, by = .05))

imagePlot(sd_m, breaks = sd_brks, col = sd_col)
mtext('SD of HCI')
contour(sd_m, add = TRUE, col = 'gray20', levels = seq(0, 0.14, by = 0.01))

imagePlot(median_m, breaks = med_brks, col = med_col)
mtext('Median of HCI')
contour(median_m, add = TRUE, col = 'gray50', levels = seq(.55, 1, by = .05))

imagePlot(iqr_m, breaks = iqr_brks, col = iqr_col)
mtext('IQR of HCI')
contour(iqr_m, add = TRUE, col = 'gray20', levels = seq(0, 0.25, by = 0.05))

### terra

# define years  --------------------------------
styear <- 1985
enyear <- 2025

# define spatial domain  --------------------------------
min_lon <- -98
max_lon <- -80
min_lat <- 18
max_lat <- 31

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
# plot(gulf_shelf)

### load coral watch temperature data lon/lat
url <- 'https://www.ncei.noaa.gov/thredds-ocean/dodsC/crw/5km/v3.1/nc/v1.0/daily/sst/1985/coraltemp_v3.1_19850101.nc'
dat <- nc_open(url)

lon <- ncvar_get(dat, 'lon')
lat <- ncvar_get(dat, 'lat')

i_lon <- which(lon >= (min_lon) & lon <= (max_lon))
i_lat <- which(lat <= (max_lat) & lat >= (min_lat))
lons <- lon[i_lon]
lats <- lat[i_lat]
### can use these if server is down
lons <- min_lon: max_lon
lats <- min_lat: max_lat

setwd("c:/Users/brendan.turley/Documents/R_projects/esr-indicator-scratch/data/intermediate_files")
# raster_list <- list()
hci_ts2 <- matrix(NA, length(styear:enyear), 366)
for (i in styear:enyear) {
  cat(i, '\n')
  tmp <- paste0('ct_sst_', i) |> readRDS()
  sst_r <- rast(aperm(tmp$sst[, 260:1, ], c(2, 1, 3)), crs = "EPSG:4326")
  ext(sst_r) <- c(range(lons), range(lats))
  time(sst_r) <- as.Date(tmp$time)

  sst_r <- crop(sst_r, gulf_shelf) |> mask(gulf_shelf)

  ### find number of cells per time step x > 20 & x < 30 in raster sst_r
  hot_cells <- sst_r >= lower & sst_r <= upper
  valid_cells <- !is.na(sst_r)
  ### sum per time step
  hot_count <- global(hot_cells, sum, na.rm = TRUE)
  valid_count <- global(valid_cells, sum, na.rm = TRUE)
  hci_ts2[which(i == (styear:enyear)), 1:nrow(hot_count)] <- unlist(hot_count / valid_count)
}
setwd("c:/Users/brendan.turley/Documents/R_projects/kgm-environ-sensitivity/data")
save(hci_ts, hci_mts, hci_ts2, file = "kgm_habitat_compression.RData") ### need to rerun when coraltemp server is available
load("kgm_habitat_compression.RData")

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
dcols <- cmocean::cmocean('algae')(5)

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

xs <- c(jday$date, rev(jday$date))
yiqr <- c(
  apply(hci_ts2, 2, quantile, probs = 0.25, na.rm = TRUE),
  rev(apply(hci_ts2, 2, quantile, probs = 0.75, na.rm = TRUE))
)
ymx <- c(
  apply(hci_ts2, 2, min, na.rm = TRUE),
  rev(apply(hci_ts2, 2, max, na.rm = TRUE))
)


setwd("c:/Users/brendan.turley/Documents/R_projects/kgm-environ-sensitivity/figures")
png('kgm_hci_stats.png', width = 7, height = 5, units = 'in', res = 300)
par(mfrow = c(2,2), mar = c(3, 4, 2, 1))
plot(1985:2025, apply(hci_ts2, 1, mean, na.rm = TRUE), 
  type = 'o', pch =16, col = 4, las = 1,
  ylab = expression(paste('HCI (', bar(x), ')')), xlab = '')
plot(1985:2025, apply(hci_ts2, 1, sd, na.rm = TRUE), 
  type = 'o', pch =16, col = 3, las = 1,
  ylab = 'HCI (std dev)', xlab = '')
plot(1985:2025, apply(hci_ts2, 1, function(x) length(which(x == 1))), 
  type = 'o', pch =16, col = 2, las = 1,
  main = 'Number of Days HCI = 1', ylab = 'Days')
plot(1985:2025, apply(hci_ts2, 1, function(x) which(x == 1)[1]), 
  type = 'o', pch =16, col = 1, las = 1,
  main = 'First Jday HCI = 1', ylab = 'Julian Day', xlab = '')
dev.off()

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


setwd("c:/Users/brendan.turley/Documents/R_projects/kgm-environ-sensitivity/figures")
png('kgm_hci_climatology.png', width = 12, height = 6, units = 'in', res = 300)
par(mfrow = c(2,2), mar = c(3, 4, 2, 1))
plot(jday$date, apply(hci_ts2, 2, median, na.rm = TRUE), 
  type = 'n', ylim = c(0, 1), xaxt = 'n', xlab = '', ylab = 'HCI', las = 1)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
polygon(xs, ymx, col = rgb(t(col2rgb('gold'))/255, alpha = 0.3), border = NA)
polygon(xs, yiqr, col = rgb(t(col2rgb('purple3'))/255, alpha = 0.4), border = NA)
points(jday$date, apply(hci_ts2, 2, median, na.rm = TRUE), type = 'l', lwd = 2)
legend('bottomright', legend = c('Median', 'IQR', 'Min/Max'), 
  fill = c('black', 'purple3', 'gold'), bty = 'n', cex = .8)

plot(jday$date, apply(hci_ts2, 2, median, na.rm = TRUE), 
  type = 'n', ylim = c(0, 1), xaxt = 'n', xlab = '', ylab = 'HCI', las = 1)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
polygon(xs, ymx, col = rgb(t(col2rgb('purple3'))/255, alpha = 0.3), border = NA)
# polygon(xs, yiqr, col = rgb(t(col2rgb('purple3'))/255, alpha = 0.3), border = NA)
points(jday$date, yr_80, type = 'l', col = dcols[1], lwd = 2)
points(jday$date, yr_90, type = 'l', col = dcols[2], lwd = 2)
points(jday$date, yr_00, type = 'l', col = dcols[3], lwd = 2)
points(jday$date, yr_10, type = 'l', col = dcols[4], lwd = 2)
points(jday$date, yr_20, type = 'l', col = dcols[5], lwd = 2)
legend('bottomright', legend = c('1980s', '1990s', '2000s', '2010s', '2020s'), 
  col = dcols, lwd = 2, bty = 'n', cex = .8)

# plot(jday$date, slopes[2, ] * 41,
#   main = 'Linear Slope of HCI Time Series',
#   type = 'h', xaxt = 'n', xlab = '', ylab = 'Slope', las = 1)
# axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
# abline(h = 0, col = 'red', lty = 2, lwd = 2)
# points(jday$date[which(slopes[8, ] < .05)],
#   slopes[2, which(slopes[8, ] < .05)] * 41,
#   col = 'blue', pch = 16)
# legend('bottomright', legend = c('Significant'), col = 'blue', pch = 16, bty = 'n', cex = .8)

plot(jday$date, apply(hci_ts2, 2, median, na.rm = TRUE),
  type = 'l', ylim = c(0, 1.2), xaxt = 'n', lwd = 2,
  main = 'HCI slope', xlab = '', ylab = 'HCI', las = 1)
axis.Date(1, at = jday$date[which(format(jday$date, "%d") == "01")], format = "%b")
# polygon(xs, ymx, col = rgb(1, 0, 0, 0.3), border = NA)
# polygon(xs, yiqr, col = rgb(0, 0, 1, 0.3), border = NA)
arrows(jday$date, 
  apply(hci_ts2, 2, median, na.rm = TRUE), 
  jday$date, 
  apply(hci_ts2, 2, median, na.rm = TRUE) + slopes[2, ] * 41, 
  length = 0.05, lwd = .7, col = ifelse(slopes[2, ] > 0, 'green3', 'orange'))
points(jday$date[which(slopes[8, ] < .05)],
  apply(hci_ts2, 2, median, na.rm = TRUE)[which(slopes[8, ] < .05)],
  bg = 'purple', pch = 21, cex = .7)
legend('bottomright', legend = c('Significant'), 
  pt.bg = 'purple', pch = 21, bty = 'n', cex = .8)

dev.off()

