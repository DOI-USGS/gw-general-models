#'''
#Ancillary data for use with NAWQA General Model (Genmod1.0) notebooks. These notebooks also can
#be used prior to using a separate series of notebooks that determine groundwater residence time distributions (GRTD) for any model. 
#
#This file contains pathnames that are common to all general model notebooks.
#
#Links for downloading the input datasets (in compressed format and opened using the Windows "extract here" command) follow:
#---------------------------------------------------------------------------------------------------------------
#Top-level NHDPlus website
#    http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php
#select a vector processing unit (VPU), for example “Northeast 01” for the example model Assabet
#    http://www.horizon-systems.com/NHDPlus/NHDPlusV2_01.php
#Create a folder named "NHDPlusV2Data" so the following path will exist C:\General_Models_User_Guide\input_data\NHDPlusV2Data
#From this page download the following five files: “NEDSnapshot”, “NHDPlusAttributes”, “NHDPlusCatchment”, 
#“NHDSnapshot”, and “WBDSnapshot”, which are compressed files,
#into the above directory. Use Windows "extract here" on each file and the result will be automaically placed into a folder 
#named "NHDPlusNE". A similar directory structure will be automatically created for other VPU locations. See the Assabet example 
#model dictionary in gen_mod_dict.py for specifying directory path locations.
#
#Path and actual file name in the "/input-data" folder for the NHDPlus datasets used are below:
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NEDSnapshot/Ned01a/elev_cm
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NHDPlusCatchment/cat
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NHDSnapshot/Hydrography/NHDFlowline.shp
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NHDSnapshot/Hydrography/NHDWaterbody.shp
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NHDPlusAttributes/PlusFlowlineVAA.dbf
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NHDSnapshot/NHDFCode.dbf
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/NHDPlusAttributes/ElevSlope.dbf
#       /NHDPlusV2Data/NHDPlusNE/NHDPlus01/WBDSnapshot/WBD/WBD_Subwatershed.shp
#
#---------------------------------------------------------------------------------------------------------------
#National Land Cover Database (NLCD)
#Create a folder named "NLCD" so the following path will exist C:\General_Models_User_Guide\input_data\NLCD
#    Available at: https://www.mrlc.gov/nlcd11_data.php
#click on "NLCD 2011 Land Cover" to download the compressed file to the above directory.
#Use Windows "extract here" command and the result will be automaically placed into a folder 
#named C:\General_Models_User_Guide\input_data\NLCD\nlcd_2011_landcover_2011_edition_2014_10_10
#
#Path and actual file name in the "/input-data" folder for the NLCD dataset used is below:
#        /NLCD/nlcd_2011_landcover_2011_edition_2014_10_10/nlcd_2011_landcover_2011_edition_2014_10_10.img
#
#---------------------------------------------------------------------------------------------------------------
#Datasets used for bedrock surface elevation and surficial thickness are part of a single report found at:
#    https://doi.org/10.3133/sim3392
#Create a folder named "Soller" so the following path will exist C:\General_Models_User_Guide\input_data\Soller
#Download the "Spatial Data" compressed file from the above link to the "Soller" folder.
#Use Windows "extract here" command and the result will be automaically placed into a folder 
#named C:\General_Models_User_Guide\input_data\Soller\sim3392_spatialdata\
#The dataset used for surficial thickness is then found at:
#C:\General_Models_User_Guide\input_data\Soller\sim3392_spatialdata\sim3392_sheet1
#The dataset for bedrock surface elevation is then found at:
#C:\General_Models_User_Guide\input_data\Soller\sim3392_spatialdata\sim3392_sheet2
#
#Path and actual file name in the "/input-data" folder for the bedrock surface elevation and surficial thickness datasets
#used are below:
#        /Soller/sim3392_spatialdata/sim3392_sheet1/sim3392_sheet1_driftthickness.img
#        /Soller/sim3392_spatialdata/sim3392_sheet2/sim3392_sheet2_bedrocktopo.img
#
#---------------------------------------------------------------------------------------------------------------
#Reitz recharge
#    Available at: https://www.sciencebase.gov/catalog/item/55d383a9e4b0518e35468e58 or
#    https://doi.org/10.5066/F7PN93P0
#Create a folder named "recharge_Reitz" so the following path will exist C:\General_Models_User_Guide\input_data\recharge_Reitz
#Download the file for effective recharge for 2013 (EffRecharge_2013.zip) to the above folder
#Use Windows "extract here" command and the result will be automaically placed into a folder named
#C:\General_Models_User_Guide\input_data\recharge_Reitz\2013
#
#Path and actual file name in the "/input-data" folder for the Reitz recharge dataset used is below:
#        /recharge_Reitz/2013/RC_eff_2013.tif
#
#---------------------------------------------------------------------------------------------------------------
#Wolock recharge
#    Available at: https://water.usgs.gov/lookup/getspatial?rech48grd
#Create a folder named "recharge_Wolock" so the following path will exist C:\General_Models_User_Guide\input_data\recharge_Wolock
#Click on ( Get this dataset here! ), then click on  "https://water.usgs.gov/GIS/dsdl/rech48grd.zip"
#and download to C:\General_Models_User_Guide\input_data\recharge_Wolock
#Use Windows "extract here" command and the result will be automaically placed into the above folder.
#
#Path and actual file name in the "/input-data" folder for the Wolock recharge dataset used is below:
#        /recharge_Wolock/rech48grd
#
#---------------------------------------------------------------------------------------------------------------
#
#Names of files used from above downloads are specified in JN1.
#Paths to the files are specified below.
#
#ft2m : float
#    Converts U.S. survey feet to meters by multiplication.
#
#
#'''

ft2m = 0.3048006096012192


# set pathnames for data sources

base_dir = 'C:/General_Models_User_Guide'
proj_dir = base_dir + '/subprojects/siteGeneral'

#this is the current means of identifying surficial geology       not yet published
qa_dir = base_dir + '/input_data/Geology'

#Used for many things. Published and available at:  http://www.horizon-systems.com/NHDPlus/NHDPlusV2_data.php 
nhd_dir = base_dir + '/input_data/NHDPlusV2Data'

#Used for landcover in JN5 in a separate series of notebooks that determine groundwater residence time distributions (GRTD)       
#            Homer, C.G., Dewitz, J.A., Yang, L., Jin, S., Danielson, P., Xian, G., Coulston, J., Herold, N.D., 
#            Wickham, J.D., and Megown, K., 2015, Completion of the 2011 National Land Cover Database for the 
#            conterminous United States-Representing a decade of land cover change information. Photogrammetric 
#            Engineering and Remote Sensing, v. 81, no. 5, p. 345-354. Data accessed online 
#            at: https://www.mrlc.gov/nlcd11_data.php
nlcd_dir = base_dir + '/input_data/NLCD/nlcd_2011_landcover_2011_edition_2014_10_10'

#Used for bedrock surface elevation
#           Soller, D.R., and Garrity, C.P., 2018, Quaternary sediment thickness and bedrock topography of the 
#           glaciated United States east of the Rocky Mountains: U.S. Geological Survey Scientific 
#           Investigations Map SIM-3392, doi 10.3133/sim3392.
#           https://doi.org/10.3133/sim3392
soller_bedrock_topo_dir = base_dir + '/input_data/Soller/sim3392_spatialdata/sim3392_sheet2'

#Used for surficial thickness
#           Soller, D.R., and Garrity, C.P., 2018, Quaternary sediment thickness and bedrock topography of the 
#           glaciated United States east of the Rocky Mountains: U.S. Geological Survey Scientific 
#           Investigations Map SIM-3392, doi 10.3133/sim3392.
#           https://doi.org/10.3133/sim3392
soller_surf_thick_dir = base_dir + '/input_data/Soller/sim3392_spatialdata/sim3392_sheet1'   

#Two sources for recharge values:

#Used for recharge, Reitz published values
#            Reitz, Meredith, Sanford, W.E., Senay, G.B., and Cazenas, Jeffrey, 2017, Annual estimates of recharge, 
#            quick-flow runoff, and ET for the contiguous US using empirical regression equations, 2000-2013: 
#            U.S. Geological Survey data release, https://doi.org/10.5066/F7PN93P0
#            or
#            https://www.sciencebase.gov/catalog/item/55d383a9e4b0518e35468e58
recharge_reitz = base_dir + '/input_data/recharge_Reitz/2013'

#The Wolock recharge below is included as an alternative to Reitz for comparison purposes; user will have to specify in JN3
#            Wolock, D.M., 2003, Estimated mean annual natural ground-water recharge in the 
#            conterminous United States: U.S. Geological Survey Open-File Report 03-311, 
#            available at https://water.usgs.gov/lookup/getspatial?rech48grd
alt_recharge_wolock = base_dir + '/input_data/recharge_Wolock'



# paths to executables
mfpth = base_dir + '/executables/MODFLOW-NWT_1.0.9/bin/MODFLOW-NWT_64.exe'
mp_exe_name = base_dir + '/executables/modpath.6_0/bin/mp6.exe'


