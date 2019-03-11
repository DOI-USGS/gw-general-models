'''
Ancillary data for use with NAWQA General Model notebooks (Genmod1.0).  These notebooks also can
be used prior to using a separate series of notebooks that determine groundwater residence time distributions (GRTD) for any model. 

This file contains variable values that are common to all general model notebooks.
Final values are in meters; conversion factors can be added.  
 
The variables are roughly in the order of how often the user might want to change them.
JN# indicates which jupyter notebooks (JN) the variable is used

ft2m : float
    Converts U.S. survey feet to meters by multiplication.
scenario_dir : string 
    Model files will be created under the model work space directory (see below) with this name.  (JN3)
add_bedrock : boolean
    If true adds a bottom layer below the top of bedrock and assigns properties as defined herein. (JN3)
num_surf_layers : integer
    Number of layers by which to divide the surficial (glacial) aquifer. (JN3)
scen : integer
    Selects one of the predefined scenarios.  Others can be added. (JN3)
K_dict : dictionary; values are floats
    The keys shouldn't be changed, but the values can be changed. (JN3)
    Hydraulic conductivity in m/d.
rech_fact : float
    ratio of recharge coarse to fine materials. default = 2.7
L : float
    Cell size in m (JN1, JN3)
min_thk : float
    Minimum total thickness into which the surficial aquifer is divided. K is assigned cell-by-cell-by-cbased
    based on whether each layer is in glacial sediments or bedrock. (JN3)
stream_width : float
    Width assigned to river or drain in m. (JN3)
stream_bed_thk : float
    Streambed thickness assigned to river or drain in m. (JN3)
stream_bed_kadjust : float
    Factor to adjust K used for computing the drain/river conductance, dimensionless. (JN3)
bedrock_thk : float
    Thickness assigned to bottom bedrock layer if present, in m. (JN3)
model_dict : nested dictionary
    Any of the above default model characteristics can be overwritten by including that characteristic in the model_dict using the same keyword. Outer dictionary has one key per model; inner dictionary has the following:
        model_ws : string
            The model workspace directory. (all JN)
        vpu : string
            Path to vector processing unit. vpu is part of the default NHD directory structure. 
            Appended to nhd_dir. (JN1)
        rpu : string
            Path to raster processing unit. rpu is part of the default NHD directory structure.
            Appended to nhd_dir and vpu. (JN1)
        df : string
            Domain file.  The name of the shapefile containg a simple outline of the model domain. (JN1)
            Must have an attribute called "ibound"
        ib_filter : integer
            The number of filter to use to eiminate isolated cells. see JN1 for documentation. (JN1)
        K_bedrock : float.
            Hydraulic conductivity of the bedrock (m/d).
        NROW, NCOL : integers (optional)
            Values to use for number of rows and columns if they are available from an existing model. (JN1)
        ibound_src : string
            If the model is from an inset or other detailed model, the ibound array can be read from this model_dict model instead.
        mfnam : string
            If the model is from an inset or other detailed model, the name file can be read from this model_dict model instead.
    NPER : integer
        Number of stress periods. currently set up for one steady state, but not hard to add transient capability. (JN1)
    hnoflo : integer
        Code to use for inactive cells in MODFLOW. (all notebooks)
    hdry : integer
        Code to use for cells that convert to dry in MODFLOW. (all notebooks)
    ratio_2_mean : float 
        Used to determine the number of particles per cell in relation to the mean volume of all cells. (JN5)
    err_tol : float
        +/- this value is the range within which the difference between head and model top is not considered an error.
        In other words, heads that are within this distance of the model top are "close enough"  (model distance units)
        (JN3 and JN4)

'''
ft2m = 0.3048006096012192
err_tol = 2.
ratio_2_mean = 0.25

# set up model scenario
# choose one of the following or make a new one

scen = 5

if scen == 1:
    scenario_dir = 'base'
    add_bedrock = False
    num_surf_layers = 1
    GHB = False
    GHB_sea = False

if scen == 2:
    scenario_dir = 'bedrock'
    add_bedrock = True
    num_surf_layers = 1
    GHB = False
    GHB_sea = False

if scen == 3:
    scenario_dir = 'layers'
    add_bedrock = True
    num_surf_layers = 3
    GHB = False
    GHB_sea = False

if scen == 4:
    scenario_dir = 'layers_GHB'
    add_bedrock = True
    num_surf_layers = 3
    GHB = True
    GHB_sea = False
	
if scen == 5:
    scenario_dir = 'layers_GHB_sea'
    add_bedrock = True
    num_surf_layers = 3
    GHB = False
    GHB_sea = True
    
# make initial guesses for K
# new values will happen because of the JN4 in which parameter estimation occurs
K_dict = {}
K_dict['K_fine'] = 10.0 * ft2m
K_dict['K_coarse'] = 100.0 * ft2m
K_dict['K_streambed'] = 100.0 * ft2m
K_dict['K_lakes'] = 10000.0 * ft2m
K_dict['K_bedrock'] = 10. * ft2m

# recharge factor--ratio of recharge volume from coarse material to fine material
rech_fact = 2.7

# set up grid cell size
L = 1000. * ft2m

# bedrock and stream geometry
rock_riv_dict = {}
rock_riv_dict['min_thk'] = 15. * ft2m
rock_riv_dict['stream_width'] = 10. * ft2m
rock_riv_dict['stream_bed_thk'] = 1. * ft2m
rock_riv_dict['river_depth'] = 10. * ft2m
rock_riv_dict['bedrock_thk'] = 100.0 * ft2m
rock_riv_dict['stream_bed_kadjust'] = 1.0
rock_riv_dict['coastal_sed_thk'] = 5.0
rock_riv_dict['coastal_sed_kadjust'] = 50.0

#variables for sea-level
sea_level = 0 #elevation of sea level (m)
den_salt = 1022 #density of saltwater in kg/m3, this is based on 29 ppt salinity at 10C
den_fresh = 1000 #density of freshwater in kg/m3

# Create dictionary for model location information. Example dictionaries, commented-out, for other model areas follow the
# example model, "Assabet".
model_dict = {}

model_dict['Assabet'] = {'ws':'Assabet','vpu':'NHDPlusNE/NHDPlus01','rpu':'ned01a','df':'Concord_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
model_dict['CoastalCT'] = {'ws':'CoastalCT','vpu':'NHDPlusNE_MA/NHDPlus01_02','rpu':'ned01a_02b','df':'CoastalCT_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['Board2'] = {'ws':'Board2','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Board_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
#model_dict['CONN'] = {'ws':'CONN','vpu':'NHDPlusNE/NHDPlus01','rpu':'ned01a','df':'Pomperaug_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07030001_domain'] = {'ws':'huc_07030001_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07c','df':'huc_07030001_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07050002_domain'] = {'ws':'huc_07050002_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07c','df':'huc_07050002_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07070001_domain'] = {'ws':'huc_07070001_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'huc_07070001_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
#model_dict['huc_07070003_domain'] = {'ws':'huc_07070003_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'huc_07070003_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07080205_domain'] = {'ws':'huc_07080205_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'huc_07080205_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07090006_domain'] = {'ws':'huc_07090006_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'huc_07090006_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07120003_domain'] = {'ws':'huc_07120003_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'huc_07120003_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07120004_domain'] = {'ws':'huc_07120004_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'huc_07120004_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['huc_07130009_domain'] = {'ws':'huc_07130009_domain','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'huc_07130009_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['IA_Willow_02'] = {'df': 'IA_Willow_02.shp', 'rpu': 'ned10d', 'vpu': 'NHDPlusMS/NHDPlus10L', 'ws': 'IA_Willow_02', 'ib_filter':0} 
#model_dict['IL_West_Fork_Mazon_03'] = {'df': 'IL_West_Fork_Mazon_03.shp', 'rpu': 'ned07a', 'vpu': 'NHDPlusMS/NHDPlus07', 'ws': 'IL_West_Fork_Mazon_03', 'ib_filter':0}
#model_dict['Kala2'] = {'ws':'Kala2','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Kala_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
#model_dict['MN_Talcot_Lake-Des_Moines_04'] = {'df': 'MN_Talcot_Lake-Des_Moines_04.shp', 'rpu': 'ned07b', 'vpu': 'NHDPlusMS/NHDPlus07', 'ws': 'MN_Talcot_Lake-Des_Moines_04', 'ib_filter':0} 
#model_dict['MO_U16_Mississippi'] = {'ws':'MO_U16_Mississippi','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'MO_U16_Mississippi.shp', 'ib_filter':1, 'K_bedrock':.8 * ft2m, 'K_fine':.1 * ft2m, 'K_coarse':100 * ft2m}
#model_dict['MO_Wildcat_02'] = {'df': 'MO_Wildcat_02.shp', 'rpu': 'ned10a', 'vpu': 'NHDPlusMS/NHDPlus10L', 'ws': 'MO_Wildcat_02', 'ib_filter':0} 
#model_dict['NE_Upper_Logan_Creek_02'] = {'df': 'NE_Upper_Logan_Creek_02.shp', 'rpu': 'ned10d', 'vpu': 'NHDPlusMS/NHDPlus10L', 'ws': 'NE_Upper_Logan_Creek_02', 'ib_filter':0}
#model_dict['NorthSkunk'] = {'ws':'NorthSkunk','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'North Skunk_domain.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
#model_dict['NY_Ramapo_02'] = {'df': 'NY_Ramapo_02.shp', 'rpu': 'ned02b', 'vpu': 'NHDPlusMA/NHDPlus02', 'ws': 'NY_Ramapo_02', 'ib_filter':0} 
#model_dict['Oconto'] = {'ws':'Oconto','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Oconto_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['OH_Three_Brothers_Creek-Grand_05'] = {'df': 'OH_Three_Brothers_Creek-Grand_05.shp', 'rpu': 'ned04c', 'vpu': 'NHDPlusGL/NHDPlus04', 'ws': 'OH_Three_Brothers_Creek-Grand_05', 'ib_filter':0} 
#model_dict['Racoon'] = {'ws':'Racoon','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'IA_U01_Racoon_alb.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
#model_dict['SD_Willow_07'] = {'df': 'SD_Willow_07.shp', 'rpu': 'ned10e', 'vpu': 'NHDPlusMS/NHDPlus10U', 'ws': 'SD_Willow_07', 'ib_filter':0}
#model_dict['SugarCreek'] = {'ws':'SugarCreek','vpu':'NHDPlusMS/NHDPlus05','rpu':'ned05b','df':'SugarCreekModel.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
#model_dict['Tomorrow'] = {'ws':'Tomorrow','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Tomorrow_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
#model_dict['Upper_fox'] = {'ws':'Upper_Fox','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'Upper_Fox_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
#model_dict['Whitedam3'] = {'ws':'Whitedam3','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'White_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
#model_dict['WI_Waumaundee_04'] = {'df': 'WI_Waumaundee_04.shp', 'rpu': 'ned07c', 'vpu': 'NHDPlusMS/NHDPlus07', 'ws': 'WI_Waumaundee_04', 'ib_filter':0}


#
# specify steady state
NPER = 1

# constants
hnoflo = -9999.
hdry = -8.89E+03

