
library(dplyr)
library(readxl)
library(viridisLite)
library(fishgrowth)
library(FSA)
library(nlstools)
library(vctrs)
library(minpack.lm)
library(MASS)
library(fields)

setwd('c:/Users/brendan.turley/Documents/CMP/data/kmk_lengths')
kgm_laa_s99 <- read_xlsx("KMK_age_ALL_20260323C.xlsx")

kgm_dat <- subset(kgm_laa_s99, State_Landed %in% c('TX','LA','MS','AL','FL')) |>
  subset(Stock == 'Gulf of America') |>
  subset(Gear_Group_Code == 'HL')

colnames(kgm_dat) <- tolower(colnames(kgm_dat))

with(aggregate(final_age ~ year, data = kgm_dat,
               function(x) length(x)),barplot(final_age,names.arg=year,las=2))
with(aggregate(final_age ~ state_landed, data = kgm_dat,
               function(x) length(x)),barplot(final_age,names.arg=state_landed,las=2))

all_yrs <- expand.grid(year = 1986:2024, final_age = 0:24)
gulf_yrage <- aggregate(final_length_mm ~ year + final_age, data = kgm_dat,
                        function(x) length(x)) |>
  merge(all_yrs, by = c('year','final_age'), all = T) |>
  arrange(year, final_age) |>
  reshape(idvar = 'final_age', timevar = 'year', direction = 'wide')
gulf_pro <- t(as.matrix(gulf_yrage[,-1]))/apply(t(as.matrix(gulf_yrage[,-1])),1,sum,na.rm=T)

imagePlot(1986:2024,0:24,gulf_pro,asp=1,col=rev(mako(60)))
boxplot(gulf_pro, names = gulf_yrage[,1])

gulf_yrsd <- aggregate(final_length_mm ~ year + final_age, data = kgm_dat,
                        sd, na.rm=T) |>
  merge(all_yrs, by = c('year','final_age'), all = T) |>
  arrange(year, final_age) |>
  reshape(idvar = 'final_age', timevar = 'year', direction = 'wide')

image(1986:2024,0:24,t(as.matrix(gulf_yrsd[,-1])),asp=1,col=rev(mako(60)))
boxplot(t(as.matrix(gulf_yrsd[,-1])), names = gulf_yrsd[,1])

# Growth per year ---------------------------------------------------------
kgm_test <- subset(kgm_dat, final_age>0 & final_age<14)
aggregate(final_length_mm ~ final_age + year, 
  data = kgm_test, 
  quantile, probs =c(.001,.5,.999), na.rm = TRUE)

kgm_filtered  <- kgm_test |>
  group_by(final_age, year, final_sex) |>
  filter(
    final_length_mm >= quantile(final_length_mm, 0.001, na.rm = TRUE),
    final_length_mm <= quantile(final_length_mm, 0.999, na.rm = TRUE)
  ) |>
  ungroup()

kgm_n_yrm <- aggregate(final_length_mm ~ year + final_age, 
  data = kgm_filtered,
  length) |>
  merge(expand.grid(year = 1986:2024, final_age = 1:13), by = c('year','final_age'), all = T) |>
  arrange(year, final_age) |>
  reshape(idvar = 'final_age', timevar = 'year', direction = 'wide')

imagePlot(1986:2024,1:13,t(as.matrix(kgm_n_yrm[,-1])),asp=1,col=rev(mako(60)))
boxplot(t(as.matrix(kgm_n_yrm[,-1])), names = kgm_n_yrm[,1])

vb2 <- makeGrowthFun(type="von Bertalanffy",pname="Original")
vb2_form <- formula(final_length_mm ~ Linf - (Linf - L0) * exp(-K * final_age))

yrs <- sort(unique(kgm_dat$year))
# output <- list()
output <- matrix(NA, length(yrs), 3) |> as.data.frame()
r <- rep(NA, length(yrs))
for(i in yrs){
  tmp <- subset(kgm_filtered, year==i)

  tmp_filtered <- tmp |>
  dplyr::group_by(final_age) |>
  dplyr::filter(dplyr::n() >= 5) |>
  dplyr::ungroup()

  resv1 <- nlsLM(final_length_mm ~ vb2(final_age, Linf, K, L0), data=tmp_filtered, 
                 start=c(Linf=1300, K=.2, L0=10),
                 control = nls.lm.control(maxiter = 1024, maxfev = 10000))
  output[i-1985,] <- rbind(coef(resv1))
  
}
output <- cbind(yrs, output) |>
  setNames(c('year','Linf','K','L0'))
# k_ts <- unlist(output)[seq(2,9*32,9)]
# l_ts <- unlist(output)[seq(1,9*32,9)]

# plot(yrs, log(k_ts), typ = 'o')
# plot(yrs, l_ts, typ = 'o')

plot(output$year, output$K, typ = 'o')
plot(output$year, output$Linf, typ = 'o')
plot(output$year, output$L0, typ = 'o')
# plot(yrs, r, typ = 'o')
plot(kgm_filtered$final_age, kgm_filtered$final_length_mm)

kgm_laa_yrm <- aggregate(final_length_mm ~ year + final_age + final_sex, 
  data = kgm_filtered,
  median, na.rm = TRUE)

with(subset(kgm_laa_yrm, final_sex == 'F'),
  plot(year, final_length_mm, typ = 'n'))
for(i in unique(subset(kgm_laa_yrm, final_sex == 'F')$final_age)){
  tmp <- subset(subset(kgm_laa_yrm, final_sex == 'F'), final_age==i)
  lines(tmp$year, tmp$final_length_mm, typ = 'o', col = i)
}

with(subset(kgm_laa_yrm, final_sex == 'M'),
  plot(year, final_length_mm, typ = 'n'))
for(i in unique(subset(kgm_laa_yrm, final_sex == 'M')$final_age)){
  tmp <- subset(subset(kgm_laa_yrm, final_sex == 'M'), final_age==i)
  lines(tmp$year, tmp$final_length_mm, typ = 'o', col = i)
}


