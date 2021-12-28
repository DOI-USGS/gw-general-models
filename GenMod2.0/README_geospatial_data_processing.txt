Five data sets are needed prior to construction General Simulation Models (GSMs). Other data sets also are required, but these will be downloaded automatically in the GSM process. For the pre-specified pathnames to work, the following data should be placed in a directory called national_spatial_data directly under the directory that contains the GSM notebooks. 

1. National-Scale Grid to Support Regional Groundwater Availability Studies and a National Hydrogeologic Database

    Clark, B.R., Barlow, P.M., Peterson, S.M., Hughes, J.D., Reeves, H.W., and Viger, R.J., 2018, National-scale grid to support regional groundwater availability studies and a national hydrogeologic database: U.S. Geological Survey data release, https://doi.org/10.5066/F7P84B24.
    
2. Average annual rates of evapotranspiration, quick-flow runoff, and recharge for the CONUS, 2000-2013

    Reitz, M., Sanford, W. E., Senay, G., & Cazenas, J. G. (2017). Annual estimates of recharge, quick-flow runoff, and ET for the contiguous US using empirical regression equations, 2000-2013 [Data set]. U.S. Geological Survey. https://doi.org/10.5066/F7PN93P0
    
3. K_b_dry_mod_md.tif
4. K_s_dry_mod_md.tif

    where 
    
    K is hydraulic conductivity,
    b is bedrock,
    s is surficial,
    mod is modified based on the following notebook, and
    and md is meters / day.
    
    These data layers were constructed using the notebook "create K datasets.ipynb" and the following data
    sources:
    
    Output_CONUS_trans_dtw, simulated depth to water and estimated transmissivity, https://water.usgs.gov/GIS/metadata/usgswrd/XML/zell2020_wrr.xml (also at https://doi.org/10.5066/P91LFFN1)
    
    which is documented in 

    Zell, W. O., Sanford, W. E. (2020). Calibrated simulation of the long-term average surficial groundwater system and derived spatial distributions of its characteristics for the contiguous United States. Water Resources Research, 55, e2019WR026724. https://doi.org/10.1029/2019WR026724
    
    and 
    
    Tashie, A., Pavelsky, T., Band, L., & Topp, S. (2021). Watershed-scale effective hydraulic properties of 
    the continental United States. Journal of Advances in Modeling Earth Systems, 13, e2020MS002440. https://doi.org/10.1029/2020MS002440
    
    The latter data set is in tabular form by HUC12. It was joined to the NHD Watershed Boundary Dataset (WBD)
    to produce the data set used here. The notebook for this processing is called "feather_reader.ipynb"
    
5. b_r_Shan.tif 

    Shangguan, W., T. Hengl, J. Mendes de Jesus, H. Yuan, and Y. Dai (2017), Mapping the global depth to 
    bedrock for land surface modeling, J. Adv. Model. Earth Syst., 9, 65–88, doi:10.1002/2016MS000686.

