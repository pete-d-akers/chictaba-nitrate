# chictaba-nitrate
This repository contains the R script for analyzing and plotting the nitrate (NO3-) data from the CHICTABA traverse, Antarctica.  
The script uses data from four CSV files, available both here and also on PANGAEA database at https://doi.org/10.1594/PANGAEA.948355.  

Units and column information for:  
chictaba_alldata.csv  
Nitrate samples were originally collected as 0.2–1.5 kg snow samples in the field along the CHICTABA traverse, Antarctica. Snow was collected in three methods: snow pits, skin layer, and 1 m depth layer. For snow pits, snow was collected in increments of 3 cm thick layers. Skin layer samples were collected as the loose snow grains 2–6 mm deep on the surface. The 1 m depth layer samples were taken as a 5–10 cm thick layer surrounding 1 m depth, which was then thoroughly mixed. Samples were melted and nitrate concentrated with a ion exchange resin at Concordia Station, Antarctica. Isotopic analysis performed with MAT253 mass spectrometer at IGE, Grenoble, France.  

id: The unique sample identifier for each NO3- sample.  
depth.mean: The average depth of the snow sample's location in the snowpack, in cm.  
depth.top: The top depth of the snow sample's location in the snowpack, in cm. For skin layer and 1 m depth layer samples, this was simply set as 0 and 100, respectively.  
depth.btm: The bottom depth of the snow sample's location in the snowpack, in cm. For skin layer and 1 m depth layer samples, this was simply set as 0 and 100, respectively.  
lat: Latitude of sampling site.  
long: Longitude of sampling site.  
NO3: Mass fraction of nitrate in ng g-1, measured by colorimetry method.  
d17O: δ17O isotopic ratio of nitrate, measured with MAT253 mass spectrometer.  
d18O: δ18O isotopic ratio of nitrate, measured with MAT253 mass spectrometer.  
D17O: Δ17O oxygen isotope anomaly (17-O excess) of nitrate, calculated from δ17O and δ18O.  
d15N: δ15N isotopic ratio of nitrate, measured with MAT253 mass spectrometer.  
NO3err: Standard error of NO3 measurements.  
d17Oerr: Standard error of δ17O measurements.  
d18Oerr: Standard error of δ18O measurements.  
D17Oerr: Standrad error of Δ17O calculations.  
d15Nerr: Standard error of δ15N measurements.  
group: Sampling method/origin group, as pits (p1–p5), skin layer (skl) or 1 m depth layer (dl).  
site: Name of stop along CHICTABA route where samples were taken.  
pseudo: 1 = the row of data is a "pseudo" sample duplicated or extrapolated from other samples for analytical and/or plotting purposes.  
smb.1979_2021: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, 1979–2021 mean.  
smb.1979_2021.mean.ci: 95% confidence interval of the mean SMB calculated from each annual SMB from 1979–2021.  
smb.1979_2021.MARavgerr: Mean error of annual SMB estimated by comparing MARv3.12.1 SMB output to in situ SMB observations for 1979–2021.  
smb.2011_2013: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, 2011–2013 mean.  
smb.2011_2013.mean.ci: 95% confidence interval of the mean SMB calculated from each annual SMB from 2011–2013.  
smb.2011_2013.MARavgerr: Mean error of annual SMB estimated by comparing MARv3.12.1 SMB output to in situ SMB observations for 2011–2013.  
smb.2013: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, 2013 alone.  
smb.2013.MARavgerr: Mean error of annual SMB estimated by comparing MARv3.12.1 SMB output to in situ SMB observations for 2013 alone.  

chictaba_rema_profile.csv  
distance.km: Distance along CHICTABA transect, in km.  
plot.dist.mean: Mean of distance.km of same row + previous row, used for accurate stepped line plotting.  
elevation.rema: Elevation, in m above sea level, at listed distance along CHICTABA. Data extracted from the Reference Elevation Model of Antarctica (https://www.pgc.umn.edu/data/rema/) with QGIS.  

chictaba_SMBMAR_profile.csv  
distance.km: Distance along CHICTABA transect, in km.  
plot.dist.mean: Mean of distance.km of same row + previous row, used for accurate stepped line plotting.  
smb.mar.1979_2021: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, 1979–2021 mean. Data extracted using QGIS.  
smb.mar.error.1979_2021: Mean error of annual SMB estimated by comparing MARv3.12.1 SMB output to in situ SMB observations for 1979–2021. Data extracted using QGIS.  
smb.mar.2011_2013: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, 2011–2013 mean. Data extracted using QGIS.  
smb.mar.error.2011_2013: Mean error of annual SMB estimated by comparing MARv3.12.1 SMB output to in situ SMB observations for 2011–2013. Data extracted using QGIS.  
smb.mar.2013: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, 2013 alone. Data extracted using QGIS.  
smb.mar.error.2013: Mean error of annual SMB estimated by comparing MARv3.12.1 SMB output to in situ SMB observations for 2013 alone. Data extracted using QGIS.  

chictaba_SMBMAR_profile_allyears.csv  
distance.km: Distance along CHICTABA transect, in km.  
plot.dist.mean: Mean of distance.km of same row + previous row, used for accurate stepped line plotting.  
Remaining columns: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.12.1 using ERA5 data, each column representing one year of accumulation from 1979 through 2021. Data extracted using QGIS.  

