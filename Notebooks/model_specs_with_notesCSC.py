'''
Ancillary data for use with NAWQA General Model notebooks.  These notebooks also can
be used to generate age distributions for any model. 

This file contains variable values and pathnames that are common to all general model notebooks.
Final values are in meters; conversion factors can be added.  
 
The variables are roughly in the order of how often the user might want to change them.
NB# indicates which notebooks the variable is used

m2ft : constants
    Converts U.S. survey feet to meters by multiplication (should have been named ft2m though).
scenario_dir : string 
    Model files will be created under the model work space directory (see below) with this name.  (NB3)
add_bedrock : boolean
    If true adds a bottom layer below the top of bedrock and assigns properties as defined herein. (NB3)
num_surf_layers : integer
    Number of layers by which to divide the surficial (glacial) aquifer. (NB3)
scen : integer
    Selects one of the predefined scenarios.  Others can be added. (NB3)
K_dict : dictionary; values are floats
    The keys shouldn't be changed, but the values can be changed. (NB3)
    Hydraulic conductivity in m/d.
L : float
    Cell size in m (NB1, NB3)
min_thk : float
    Minimum thickness of each layer into which the surficial aquifer is divided. (NB3)
stream_width : float
    Width assigned to river or drain in m. (NB3)
stream_bed_thk : float
    Streambed thickness assigned to river or drain in m. (NB3)
bedrock_thk : float
    Thickness assigned to bottom bedrock layer if present, in m. (NB3)
nhd_dir : string
    Pathname to NHD directory. The default directory structure defined in the NHDPlus is assumed. (NB1)
geol_dir : string
    Pathname to USGS_DS_656/Shapefiles. (NB1)
thickness_dir : string
    Pathname to  BedrockTopo_DRS_11-13. (NB1)
recharge_dir : string
    Pathname to rchg_mm_geotiff. (NB1)
mohp_dir : string
    Pathname to  HydroPosition (optional) (NB1)
mohp2_dir : string
    Pathname to  Hydroposition/MOHP_TopoCatchments_CleanedVersion. (optional) (NB1)
thies2_dir : string
    Pathname to  Hydroposition/MOHP_ThiessCatchments_CleanedVerion. (optional) (NB1)
pour_dir : string
    Pathname to  Merged. (optional) (NB1)
swb_dir : string
    Pathname to  SWB_RechargeGrid. (NB1)
surg_dir : string
    Pathname to  SSURGO. (optional) (NB1)
qa_dir : string
    Pathname to  Quat Atlas. (optional) (NB1)
soller_thick_dir : string
    Pathname to  SollerProvisionalDriftThickness. (NB1)
model_dict : nested dictionary
    Outer dictionary has one key per model; inner dictionary, has the following:
        model_ws : string
            The model workspace directory. (all NB)
        vpu : string
            Path to vector processing unit. vpu is part of the default NHD directory structure. 
            Appended to nhd_dir. (NB1)
        rpu : string
            Path to raster processing unit. rpu is part of the default NHD directory structure.
            Appended to nhd_dir and vpu. (NB1)
        df : string
            Domain file.  The name of the shapefile containg a simple outline of the model domain. (NB1)
            Must have an attribute called "ibound"
        ib_filter : integer
            The number of filter to use to eiminate isolated cells. see NB1 for documentation. (NB1)
        bedrock : float.
            Hydraulic conductivity of the bedrock (m/d).
        NROW, NCOL : integers (optional)
            Values to use for number of rows and columns if they are available from an existing model. (NB1)
    mfpth : string
        Path to MODFLOW executable. (NB3, NB4)
    mp_exe_name : string
        Path to MODPATH (version 6) executable. (NB5)
    NPER : integer
        Number of stress periods. currently set up for one steady state, but not hard to add transient capability. (NB1)
    hnoflo : integer
        Code to use for inactive cells in MODFLOW. (all notebooks)
    hdry : integer
        Code to use for cells that convert to dry in MODFLOW. (all notebooks)

'''
ft2m = 0.3048006096012192


# set pathnames for data sources
base_dir = 'C:/General_Models_WRR'
proj_dir = base_dir + '/subprojects/siteGeneral'

#used for many things                 published
nhd_dir = base_dir + '/input_data/NHDPlusV2Data'

#used for landcover used in nb5       published
nlcd_dir = base_dir + '/input_data/NLCD'


#this is the current means of identifying surficial geology     factor_added_QA_from_Diick_Yager.shp     not yet published
qa_dir = base_dir + '/input_data/Geology'

#used for surficial thickness, I assume this one-Soller and Garrity is not yet published?  and is different than USGS DS 656?
soller_thick_dir = base_dir + '/input_data/Soller_Provisional_DriftThickness'   


#used for bedrock surface elevation       can not find it online, not yet published?
thickness_dir = base_dir + '/input_data/BedrockTopo_DRS_11-13'

#used for multi order hydrologic position          can not find it online, not yet published?
#finds tif files in all subdirectories in each below
mohp_dir = base_dir + '/input_data/HydroPosition'
mohp2_dir = base_dir + '/input_data/Hydroposition2/MOHP_TopoCatchments_CleanedVersion'
thies2_dir = base_dir + '/input_data/Hydroposition2/MOHP_ThiessCatchments_CleanedVerion'

# is pour_dir part of mohp download?
pour_dir = base_dir + '/input_data/Merged'

#used for recharge.  Is there a link or reference for this? Published yet?
swb_dir = base_dir + '/input_data/SWB_RechargeGrid'

#The following seem to not be used. Could they be removed from model_specs and the notebooks?

#not used anymore, superceded by qa_dir    but geol_dir includes     factor_added_stack_map.shp
geol_dir = base_dir + '/input_data/USGS_DS_656/shapefiles'

#not used, contains Reitz and rchg_mm, but there are two paths used for Reitz
recharge_dir = base_dir + '/input_data/recharge'

#not used for recharge, Reitz is not used
alt_recharge_file = '/input_data/recharge/rch_mm_Reitz.tif'


#perhaps used only if available?
surg_dir = base_dir + '/input_data/SSURGO'

#Are these used?
#well_pth = base_dir + '/input_data/all_glac_wells.shp'
well_pth = base_dir + '/input_data/glac_wells_PAS2.shp'
pumping_rate_pth = base_dir + '/input_data/PumpingRates.txt'
lpm_dir = base_dir + '/input_data/LPM_AgeDists'

#not used, but included for now beacuse it has to be defined for nb to run
gage_file = base_dir + '/GenModData/gages/GagesRunoff_w_stats.shp'



# paths to executables
mfpth = base_dir + '/executables/MODFLOW-NWT_1.0.9/bin/MODFLOW-NWT_64.exe'
mp_exe_name = base_dir + '/executables/modpath.6_0/bin/mp6.exe'


