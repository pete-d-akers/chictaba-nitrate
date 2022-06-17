# chictaba-nitrate
This repository contains the R scripts for analyzing and plotting the nitrate (NO3-) data from the CHICTABA traverse, Antarctica.
The script uses data from three CSV files, available both here and also on PANGAEA database (permanent links to come).

Units and column information for:
chictaba_NO3data.csv

Nitrate samples were originally collected as 0.2–1.5 kg snow samples in the field along the CHICTABA traverse, Antarctica. Snow was collected in three methods: snow pits, skin layer, and 1 m depth layer. For snow pits, snow was collected in increments of 3 cm thick layers. Skin layer samples were collected as the loose snow grains 2–6 cm deep on the surface. The 1 m depth layer samples were taken as a 5–10 cm thick layer surrounding 1 m depth, which was then thoroughly mixed. Samples were melted and nitrate concentrated with a ion exchange resin at Concordia Station, Antarctica. Isotopic analysis performed with MAT253 mass spectrometer at IGE, Grenoble, France.

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
smb: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.6.4 using ERA-Interim data 1979-2017.
pseudo: 1 = the row of data is a "pseudo" sample duplicated or extrapolated from other samples for analytical and/or plotting purposes.

chictaba_rema_profile.csv
distance.km: Distance along CHICTABA transect, in km.
plot.dist.mean: Mean of distance.km of same row + following row, used for accurate stepped line plotting.
elevation.rema: Elevation, in m above sea level, at listed distance along CHICTABA. Data extracted from the Reference Elevation Model of Antarctica (https://www.pgc.umn.edu/data/rema/) with QGIS.

chictaba_SMBMAR_profile.csv
distance.km: Distance along CHICTABA transect, in km.
plot.dist.mean: Mean of distance.km of same row + following row, used for accurate stepped line plotting.
smb.mar: Surface mass balance in kg m-2 a-1 output from from MAR regional climate model version 3.6.4 using ERA-Interim data 1979-2017. Data extracted using QGIS.
