# Mathematical function for CHLA profiles
# Following Mignot et al. 2011 and Carranza et al. 2018
fgauss_expo <- function(z, Fsurf, Zdemi, Fmax, Zmax, dz){
  Fsurf*exp((-log(2)/Zdemi)*z) + Fmax*exp(-(z-Zmax)^2/dz^2)
}

fsigmoid <- function(z, Fsurf, Zdemi, s){
  Fsurf*(1/(1+exp((Zdemi-z)*s)))
}

fexpo <- function(z, Fsurf, Zdemi){
  Fsurf*exp((-log(2)/Zdemi)*z)
}

fgauss <- function(z, Fmax, Zmax, dz){
  Fmax*exp(-(z-Zmax)^2/dz^2)
}

fgauss_sigmoid <- function(z, Fsurf, Zdemi, s, Fmax, Zmax, dz){
  Fmax*exp(-(z-Zmax)^2/dz^2) + Fsurf*(1/(1+exp((Zdemi-z)*s)))
}

#########################
##### BEGIN THE FIT #####
#########################

fit <- ldply(as.list(unique(profiles$juld)), function(i){
  
  tmp <- profiles[profiles$juld==i,]
  
  #data needed for R^2 and R^2 adjusted
  mean_fluo <- mean(tmp$fluo, na.rm=T) # mean fluo over all profile
  ss_tot <- sum((tmp$fluo - mean_fluo)^2, na.rm=T) # for R^2 computation
  n <- length(tmp$depth) # number of data points
  # Arthur Capet's suggestion => Fsurf should be a mean of the fluo in the MLD
  index_MLD <- which.min(tmp$depth <= MLDdf$MLD[which(MLDdf$juld == i)]) 
  
  ##############
  # PARAMETERS #
  ##############
  
  #Parameters estimation
  Fsurf <- mean(tmp$fluo[1:index_MLD], na.rm = T)
  Zdemi <- tmp$depth[which.max(tmp$fluo <=  Fsurf/2)]
  if(MLDdf$MLD[which(MLDdf$juld == i)] > Zdemi){
    Zdemi <- MLDdf$MLD[which(MLDdf$juld == i)]
  }
  maxindex <- which.max(tmp$fluo)
  Zmax <- tmp$depth[maxindex]
  Fmax <- tmp$fluo[maxindex]
  dz <- 5 # First guess
  s <- -0.01 # First guess without explication like dz, it is fixed and adapted if initial values do not work
 
  #############################
  # GAUSSIAN-EXPONENTIAL (GE) #
  #############################
  
  res_GE <- tryCatch(nlsLM(fluo ~ fgauss_expo(depth, Fsurf, Zdemi, Fmax, Zmax, dz),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  
  if(inherits(res_GE, "error")){# handle error if initial dz is not good
    dz <- 10
    res_GE <- tryCatch(nlsLM(fluo ~ fgauss_expo(depth, Fsurf, Zdemi, Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, Fmax = Fmax, Zmax = Zmax, dz = dz),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  #R^2 and R^2 adjusted
  ge <- fgauss_expo(tmp$depth,coef(res_GE)["Fsurf"], coef(res_GE)["Zdemi"], coef(res_GE)["Fmax"], coef(res_GE)["Zmax"], coef(res_GE)["dz"])
  ss_res_ge <- sum((ge - tmp$fluo)^2, na.rm=T)  
  rcoef_ge <- 1-(ss_res_ge/ss_tot)
  rcoef_ge <- 1 - (1-rcoef_ge)*(n-1)/(n-5-1) # final R^2 ADJUSTED
  
  ###############
  # SIGMOID (S) #
  ###############
  
  res_S <- tryCatch(nlsLM(fluo ~ fsigmoid(depth, Fsurf, Zdemi, s),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  if(inherits(res_S, "error")){# added due to a bug
    Zdemi <- 30
    #s <- -0.01
    res_S <- tryCatch(nlsLM(fluo ~ fsigmoid(depth, Fsurf, Zdemi, s),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  
  sigmoid <- fsigmoid(tmp$depth, coef(res_S)["Fsurf"], coef(res_S)["Zdemi"], coef(res_S)["s"])
  ss_res_sigmoid <- sum((sigmoid - tmp$fluo)^2, na.rm = T)
  rcoef_sigmoid <- 1-(ss_res_sigmoid/ss_tot)               
  rcoef_sigmoid <- 1 - (1-rcoef_sigmoid)*(n-1)/(n-3-1)
  
  ################
  # GAUSSIAN (G) #
  ################
  
  res_G <- tryCatch(nlsLM(fluo ~ fgauss(depth, Fmax, Zmax, dz),
                        start = c(Fmax = Fmax, Zmax = Zmax, dz = dz),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  gauss <- fgauss(tmp$depth, coef(res_G)["Fmax"], coef(res_G)["Zmax"], coef(res_G)["dz"])
  ss_res_gauss <- sum((gauss - tmp$fluo)^2, na.rm=T)                
  rcoef_gauss <- 1-(ss_res_gauss/ss_tot)  
  rcoef_gauss <- 1 - (1-rcoef_gauss)*(n-1)/(n-3-1)
  
  ###################
  # EXPONENTIAL (E) #
  ###################
  
  res_E <- tryCatch(nlsLM(fluo ~ fexpo(depth, Fsurf, Zdemi),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  expo <- fexpo(tmp$depth, coef(res_E)["Fsurf"], coef(res_E)["Zdemi"])
  ss_res_expo <- sum((expo - tmp$fluo)^2, na.rm=T)
  rcoef_expo <- 1-(ss_res_expo/ss_tot)
  rcoef_expo <- 1 - (1-rcoef_expo)*(n-1)/(n-2-1)
  
  #########################
  # GAUSSIAN-SIGMOID (GS) #
  #########################
  
  res_GS <- tryCatch(nlsLM(fluo ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                           Fmax, Zmax, dz),
                        start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                        data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                        minFactor = 1/2048, warnOnly=T)),
                  error=function(e) e)
  
  if(inherits(res_GS, "error")){
    Zdemi <- 30
    s <- -0.1 
    res_GS <- tryCatch(nlsLM(fluo ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                             Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  if(inherits(res_GS, "error")){
    dz <- 10
    s <- -0.1
    res_GS <- tryCatch(nlsLM(fluo ~ fgauss_sigmoid(depth, Fsurf, Zdemi, s,
                                             Fmax, Zmax, dz),
                          start = c(Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, Zmax=  Zmax, dz = dz),
                          data=tmp[,c(1,3)], control = nls.control(maxiter=1000,
                                                          minFactor = 1/2048, warnOnly=T)),
                    error=function(e) e)
  }
  
  gs <- fgauss_sigmoid(tmp$depth, coef(res_GS)["Fsurf"], coef(res_GS)["Zdemi"], coef(res_GS)["s"], coef(res_GS)["Fmax"], coef(res_GS)["Zmax"], coef(res_GS)["dz"])
  ss_res_gs <- sum((gs - tmp$fluo)^2, na.rm=T)
  rcoef_gs <- 1-(ss_res_gs/ss_tot)
  rcoef_gs <- 1 - (1-rcoef_gs)*(n-1)/(n-6-1)
  
  #######################
  # SHAPE DETERMINATION #
  #######################
  
  # Which is the best fit (the order was chosen so that if multiple values are equal, the first is taken)
  bestfit <- which.max(c(rcoef_sigmoid, rcoef_expo, rcoef_gauss, rcoef_ge, rcoef_gs))

  # Rejection criteria 
  Rcrit <- 0.9
  
  if(c(rcoef_sigmoid, rcoef_expo, rcoef_gauss, rcoef_ge, rcoef_gs)[bestfit] < Rcrit){
    Zmax <- NA
    r2adj <- NA
    shape <- "O" # "O" for "Other"
    chla_at_dcm <- NA
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- NA
    Zmax <- NA
    dz <- NA
  }else if(bestfit == 1){
    Zmax <- NA # no Zmax
    r2adj <- rcoef_sigmoid
    shape <- "S"
    chla_at_dcm <- NA
    Fsurf <- coef(res_S)["Fsurf"]
    Zdemi <- coef(res_S)["Zdemi"]
    s <- coef(res_S)["s"]
    Fmax <- NA
    Zmax <- NA
    dz <- NA
  }else if(bestfit == 2){
    Zmax <- NA
    r2adj <- rcoef_expo
    shape <- "E"
    chla_at_dcm <- NA
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- NA
    Zmax <- NA
    dz <- NA
  }else if(bestfit == 3){
    Zmax <- coef(res_G)["Zmax"]
    r2adj <- rcoef_gauss
    shape <- "G"
    chla_at_dcm <- tmp$fluo[which.min(abs(tmp$depth - Zmax))] 
    Fsurf <- NA
    Zdemi <- NA
    s <- NA
    Fmax <- coef(res_G)["Fmax"]
    Zmax <- coef(res_G)["Zmax"]
    dz <- coef(res_G)["dz"]
  }else if(bestfit == 4){
    Zmax <- coef(res_GE)["Zmax"]
    r2adj <- rcoef_ge
    shape <- "GE"
    chla_at_dcm <- tmp$fluo[which.min(abs(tmp$depth - Zmax))]
    Fsurf <- coef(res_GE)["Fsurf"]
    Zdemi <- coef(res_GE)["Zdemi"]
    s <- NA
    Fmax <- coef(res_GE)["Fmax"]
    Zmax <- coef(res_GE)["Zmax"]
    dz <- coef(res_GE)["dz"]
  }else{
    Zmax <- coef(res_GS)["Zmax"]
    r2adj <- rcoef_gs
    shape <- "GS"
    chla_at_dcm <- tmp$fluo[which.min(abs(tmp$depth - Zmax))]
    Fsurf <- coef(res_GS)["Fsurf"]
    Zdemi <- coef(res_GS)["Zdemi"]
    s <- coef(res_GS)["s"]
    Fmax <- coef(res_GS)["Fmax"]
    Zmax <- coef(res_GS)["Zmax"]
    dz <- coef(res_GS)["dz"]
  }
  
  show(i)
  
  data.frame(shape = shape, r2adj = r2adj, Zmax = Zmax, juld = tmp$juld[1],
             lon = tmp$lon[1], lat = tmp$lat[1], day = tmp$day[1],
             month = tmp$month[1], year = tmp$year[1], DOY = tmp$DOY[1],
             platform = tmp$platform[1], chla_surf = tmp$fluo[1], #chla_surf is the chla value at the surface of the profile
             day_night = tmp$day_night[1], MLD = MLDdf$MLD[which(MLDdf$juld == i)], chla_dcm = chla_at_dcm, Fsurf = Fsurf, Zdemi = Zdemi, s = s, Fmax = Fmax, dz = dz,
             date = tmp$date_long[1], chla_dcm_chla_surf = chla_at_dcm/tmp$fluo[1]) 
})

# save fit
init_fit <- fit

# Let's keep only profiles fitted as "G", "GE" or "GS"
profiles <- filter(profiles, juld %in% unique(fit[fit$shape %in% c("G", "GE", "GS"),]$juld))

# clean data according to JULD in PROFILES 
fit <- filter(fit, juld %in% unique(profiles$juld))
density_profiles <- filter(density_profiles, juld %in% unique(profiles$juld))
MLDdf <- filter(MLDdf, juld %in% unique(profiles$juld))
doxy <- filter(doxy, juld %in% unique(profiles$juld))
PAR <- filter(PAR, juld %in% unique(profiles$juld))

# Compute the SIGMA@MAX MLD for each year and each platform (based on INIT MLD, not filtered!!)
MLDdf_info <- ddply(init_MLD, ~year~platform, summarize,
                    max = max(MLD, na.rm = T),
                    sigmaMaxMLD = sigmaMaxMLD[which.max(MLD)])

################
# DATA SORTING #
################

# Additional criteria for elimination of false DCMs
profiles <- filter(profiles, juld %in% fit[which(fit$chla_dcm > 1/0.75 * fit$chla_surf),]$juld) # 33% criteria

# clean data according to JULD in PROFILES 
fit <- filter(fit, juld %in% unique(profiles$juld))
density_profiles <- filter(density_profiles, juld %in% unique(profiles$juld))
MLDdf <- filter(MLDdf, juld %in% unique(profiles$juld))
doxy <- filter(doxy, juld %in% unique(profiles$juld))
PAR <- filter(PAR, juld %in% unique(profiles$juld))

# for a "boxplot", add ratio chla_dcm/chla_surf
ratio_chla_dcm_chla_surf <- ldply(as.list(fit$juld), function(i){
  tmp <- fit[fit$juld == i,]
  ratio <- tmp$chla_dcm/tmp$chla_surf
  data.frame(ratio = ratio)
})

fit$ratio_chla_dcm_chla_surf <- ratio_chla_dcm_chla_surf$ratio
rm(ratio_chla_dcm_chla_surf)

sigma_ratio <- ldply(as.list(MLDdf$juld), function(i){
  tmp <- MLDdf[MLDdf$juld == i,]
  tmp2 <- fit[fit$juld == i,]
  tmp3 <- profiles[profiles$juld == i,]
  sigmaMaxMLD <- MLDdf_info$sigmaMaxMLD[which(MLDdf_info$year == tmp$year & MLDdf_info$platform == tmp$platform)]
  sigmaDCM <- approx(tmp3$depth, tmp3$sigma, tmp2$Zmax)$y #interpolation
  ratio = sigmaDCM/sigmaMaxMLD
  data.frame(ratio = ratio, juld = tmp$juld[1], month = tmp$month[1], sigmaMaxMLD = sigmaMaxMLD, sigmaDCM = sigmaDCM)
})

rm(MLDdf_info)

# add some final touch to create the Rdata file

# ADD SOME DATA FOR THE .Rmd file
# get unique positions of CHLA profiles
lat_lon <- subset(profiles, select = c("lat","lon","juld"))
lat_lon <- lat_lon[!duplicated(lat_lon),]
lat_lon <- transform(lat_lon, id=as.numeric(factor(juld)))
rownames(lat_lon) <- NULL

#same for profiles (based on the .Rmd script)
rownames(profiles) <- NULL
profiles <- transform(profiles, id = as.numeric(factor(juld)))

save.image(file='DataPaper.RData') # save workspace directly
