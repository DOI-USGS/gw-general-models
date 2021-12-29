This file gives an overview of the General Simulation Model (GSM) process. More detail on the process 
is provided in text inside the notebooks. GSMs are discussed in 

    Starn, J. J., & Belitz, K. (2018).Regionalization of groundwater residence time using metamodeling. Water Resources Research, 54, 6357–6373. https://doi.org/10.1029/2017WR021531

    Starn, J. J., Kauffman, L. J., Carlson, C. S., Reddy, J. E., & Fienen, M. N. (2021). Three-dimensional distribution of groundwater residence time metrics in the glaciated United States using metamodels trained on general numerical simulation models. Water Resources Research, 57, e2020WR027335. https://doi.org/10.1029/2020WR027335

Three preparatory steps are required before the GSM process is started.
    1. Create a Python environment as described in ReadMe_create_python_environment.txt
    2. Copy the national_spatial_data directory as described in README_geospatial_data_processing.txt
    3. Have available executable files for MODFLOW 6 (6.1.1 was used in the example) and MODPATH 7 (7.2.001 was used in the example). The location of these files is specified in the notebooks.
    
Once the environment has been activated, GSMs are created using a sequence of Jupyter notebooks. Static 
views of the notebooks, including selected output, are provided in pdf renderings of the notebooks.

General_Simulation_Model_0_Download_NHD_NGWM.ipynb
    This notebook will help you select the area to be modeled using a Leaflet map. It will download 
    selected data from the National Hydrography Dataset High Resolution (NHDPlus HR) and the National Groundwater Model. It will create two directories, called "/gis" and "/downloads", that will be used 
    in subsequent notebooks.
    
General_Simulation_Model_1_Geospatial_processing.ipynb
    This notebook will resample NHD and national spatial data onto a finite difference model grid. The user
    only needs to specify the grid cell size. Spatial data will be created in the "/gis" directory.
    
General_Simulation_Model_2_Create_GSM_NGWM.ipynb
    This notebook creates a MODFLOW6 finite difference model. Some model properties can be specified by the 
    user near the top of the notebook, but defaults are provided. The GSM is created in a directory called "/model_ws" (model workspace).
    
General_Simulation_Model_3_Calibrate_Model_Bayes.ipynb
    This notebook calibrates the GSM. There is user interaction available throughout the notebook, but defaults are provided that should work in most cases. The calibrated model is placed in a directory 
    called "optimal_model".
    
General_Simulation_Model_4_Multiple_runs_and_plots.ipynb
    This notebook is for evaluating uncertainty. It runs the model for multiple parameter sets that lie on the Pareto front, but it is just a prototype.  More work needs to be done to evaluate the results of multiple parameter sets on travel time distributions.
    
General_Simulation_Model_5_ParticleTrackingfromStreams.ipynb
    This notebook creates a MODFLOW7 simulation. Particles are placed in a line at the location of streams in 
    proportion to their groundwater inflow. Particles are tracked backward to their point of entry into the aquifer and their travel times are recorded in a modified endpoint file. MODPATH7 output is placed in the "/optimal_model" directory.
    
General_Simulation_Model_6_Baseflow_Age_Distributions.ipynb
    This notebook produces a spatial dataset of travel time characteristics by stream reach (nhd_age shapefile) and several images (png) of their spatial distribution.   
    
Notebooks using to process particle travel time data for water wells rather streams can be found at

     Starn, J.J., Kauffman, L.J., Carlson, C.S., Reddy, J.E., and Fienen, M.N., 2020, Data for three-dimensional distribution of groundwater residence time metrics in the glaciated United States using metamodels trained on general numerical simulation models: U.S. Geological Survey data release, https://doi.org/10.5066/P9BNWWCU.
    
