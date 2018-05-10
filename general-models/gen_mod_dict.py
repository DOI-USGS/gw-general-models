'''
Ancillary data for use with NAWQA General Model notebooks.  These notebooks also can
be used to generate age distributions for any model. 

This file contains variable values and pathnames that are common to all general model notebooks.
Final values are in meters; conversion factors can be added.  
 
The variables are roughly in the order of how often the user might want to change them.
NB# indicates which notebooks the variable is used

ft2m : float
    Converts U.S. survey feet to meters by multiplication.
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
rech_fact : float
    ratio of recharge coarse to fine materials. default = 2.7
L : float
    Cell size in m (NB1, NB3)
min_thk : float
    Minimum total thickness into which the surficial aquifer is divided. K is assigned cell-by-cell-by-cbased
    based on whether each layer is in glacial sediments or bedrock. (NB3)
stream_width : float
    Width assigned to river or drain in m. (NB3)
stream_bed_thk : float
    Streambed thickness assigned to river or drain in m. (NB3)
stream_bed_kadjust : float
    Factor to adjust K used for computing the drain/river conductance, dimensionless. (NB3)

bedrock_thk : float
    Thickness assigned to bottom bedrock layer if present, in m. (NB3)
model_dict : nested dictionary
    Any of the above default model characteristics can be overwritten by including that characteristic in the model_dict using the same keyword. Outer dictionary has one key per model; inner dictionary has the following:
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
        K_bedrock : float.
            Hydraulic conductivity of the bedrock (m/d).
        NROW, NCOL : integers (optional)
            Values to use for number of rows and columns if they are available from an existing model. (NB1)
        ibound_src : string
            If the model is from an inset or other detailed model, the ibound array can be read from this model_dict model instead.
        mfnam : string
            If the model is from an inset or other detailed model, the name file can be read from this model_dict model instead.
    NPER : integer
        Number of stress periods. currently set up for one steady state, but not hard to add transient capability. (NB1)
    hnoflo : integer
        Code to use for inactive cells in MODFLOW. (all notebooks)
    hdry : integer
        Code to use for cells that convert to dry in MODFLOW. (all notebooks)
    ratio_2_mean : float 
        Used to determine the number of particles per cell in relation to the mean volume of all cells. (NB5)
    err_tol : float
        +/- this value is the range within which the difference between head and model top is not considered an error.
        In other words, heads that are within this distance of the model top are "close enough"  (model distance units)
        (NB3 and NB4)

'''
ft2m = 0.3048006096012192
err_tol = 2.
ratio_2_mean = 0.25

# set up model scenario
# choose one of the following or make a new one

scen = 3

if scen == 1:
    scenario_dir = 'base'
    add_bedrock = False
    num_surf_layers = 1
    GHB = False

if scen == 2:
    scenario_dir = 'bedrock'
    add_bedrock = True
    num_surf_layers = 1
    GHB = False

if scen == 3:
    scenario_dir = 'layers'
    add_bedrock = True
    num_surf_layers = 3
    GHB = False
   
if scen == 4:
    scenario_dir = 'layers_GHB'
    add_bedrock = True
    num_surf_layers = 3
    GHB = True
    
# make initial guesses for K
# new values will happen because of the NB4 in which parameter estimation occurs
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

# create dictionary for model location information.
model_dict = {}
model_dict['CONN'] = {'ws':'CONN','vpu':'NHDPlusNE/NHDPlus01','rpu':'ned01a','df':'Pomperaug_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
# model_dict['WHMI'] = {'ws':'WHMI','vpu':'NHDPlusMS/NHDPlus05','rpu':'ned05c','df':'Miami_domain.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
model_dict['Upper_fox'] = {'ws':'Upper_Fox','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'Upper_Fox_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
# model_dict['Wolf'] = {'ws':'Wolf','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Wolf_domain.shp', 'ib_filter':2, 'K_bedrock':1 * ft2m}
#model_dict['fwp2'] = {'ws':'fwp2','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'grid_aea.shp', 'ib_filter':2, 'K_bedrock':1 * ft2m, 
#'NROW':930, 'NCOL':650, 'ibound_src':'FWP', 'mf_src':'../Subprojects/GLAC/FoxWolf/fwp2/2_input/A_modflow', 'mfnam':'FWP_BR.nam'}
# model_dict['FWP'] = {'ws':'FWP','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Fox_domain.shp', 'ib_filter':2, 'K_bedrock':1 * ft2m}
model_dict['Oconto'] = {'ws':'Oconto','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Oconto_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
model_dict['Tomorrow'] = {'ws':'Tomorrow','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Tomorrow_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
model_dict['SugarCreek'] = {'ws':'SugarCreek','vpu':'NHDPlusMS/NHDPlus05','rpu':'ned05b','df':'SugarCreekModel.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
model_dict['NorthSkunk'] = {'ws':'NorthSkunk','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'North Skunk_domain.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
model_dict['Racoon'] = {'ws':'Racoon','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'IA_U01_Racoon_alb.shp', 'ib_filter':1, 'K_bedrock':100 * ft2m}
#model_dict['Assabet'] = {'ws':'Assabet','vpu':'NHDPlusNE/NHDPlus01','rpu':'ned01a','df':'Assabet_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
model_dict['Assabet'] = {'ws':'Assabet','vpu':'NHDPlusNE/NHDPlus01','rpu':'ned01a','df':'Concord_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}
# model_dict['Mass_MTBE_Front'] = {'ws':'Mass_MTBE_Front','vpu':'NHDPlusNE/NHDPlus01','rpu':'ned01a','df':'North River_domain.shp', 'ib_filter':0, 'K_bedrock':1 * ft2m}

model_dict['Kala2'] = {'ws':'Kala2','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Kala_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
# model_dict['Kala2LMB'] = {'ws':'Kala2LMB','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'grid_aea.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m, 
# 'NROW':730, 'NCOL':1230, 'ibound_src':'Kala2', 'mf_src':'../Subprojects/GLAC/LMB inset/Kala2-INPUT&LST', 'mfnam':'kala2-postGK-BLUR-SFR_DRIV-DRN-SS-MNW2SKINhi-NWT-MP.nam'}
model_dict['Board2'] = {'ws':'Board2','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Board_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
# model_dict['Board2LMB'] = {'ws':'Board2LMB','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'grid_aea.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m, 
# 'NROW':830, 'NCOL':750, 'ibound_src':'Board2', 'mf_src':'../Subprojects/GLAC/LMB inset/Board2-INPUT&LST', 'mfnam':'board2-postGK-SFR_DRIV-DRN-SS-MNW2SKINHI-NWT.nam'}
model_dict['Mani3'] = {'ws':'Mani3','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'Mani_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
# model_dict['Mani3LMB'] = {'ws':'Mani3LMB','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'grid_aea.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m, 
# 'NROW':910, 'NCOL':530, 'ibound_src':'Mani3', 'mf_src':'../Subprojects/GLAC/LMB inset/Mani3-INPUT&LST', 'mfnam':'Mani3-postGK-SFR_DRIV-DRN-SS-MNW2SKINHI-NWT.nam'}
model_dict['Whitedam3'] = {'ws':'Whitedam3','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'White_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
# model_dict['Whitedam3LMB'] = {'ws':'Whitedam3LMB','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'grid_aea.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m, 
# 'NROW':570 , 'NCOL':670, 'ibound_src':'Whitedam3', 'mf_src':'../Subprojects/GLAC/LMB inset/Whitedam3-INPUT&LST', 'mfnam':'whitedam3-postGK-SFR_DRIV-DRN-SS-MNW2SKINHI-NWT.nam'}
# model_dict['Uofox3'] = {'ws':'Uofox3','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'Upper_Fox_domain.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m}
# model_dict['Uofox3LMB'] = {'ws':'Uofox3LMB','vpu':'NHDPlusGL/NHDPlus04','rpu':'ned04d','df':'grid_aea.shp', 'ib_filter':0, 'K_bedrock':100 * ft2m, 
# 'NROW':1050 , 'NCOL':430, 'ibound_src':'Uofox3', 'mf_src':'../Subprojects/GLAC/LMB inset/Uofox3-INPUT&LST', 'mfnam':'upfox3-postGK-SFR_DRIV-DRN-SS-MNW2SKINHI-NWT.nam'}

# model_dict['IA_U77_NorthSkunk'] = {'ws':'IA_U77_NorthSkunk','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'IA_U77_NorthSkunk.shp', 'ib_filter':1, 'K_bedrock':.1 * ft2m}
# model_dict['IA_U78_Nodaway'] = {'ws':'IA_U78_Nodaway','vpu':'NHDPlusMS/NHDPlus10L','rpu':'ned10b','df':'IA_U78_Nodaway.shp', 'ib_filter':1, 'K_bedrock':.1 * ft2m}
# model_dict['MO_U14_Grand'] = {'ws':'MO_U14_Grand','vpu':'NHDPlusMS/NHDPlus10L','rpu':'ned10a','df':'MO_U14_Grand.shp', 'ib_filter':1, 'K_bedrock':.8 * ft2m, 'K_fine':5 * ft2m}
# model_dict['MO_U15_Weldon'] = {'ws':'MO_U15_Weldon','vpu':'NHDPlusMS/NHDPlus10L','rpu':'ned10a','df':'MO_U15_Weldon.shp', 'ib_filter':1, 'K_bedrock':.05 * ft2m}
# model_dict['MO_U17_Missouri'] = {'ws':'MO_U17_Missouri','vpu':'NHDPlusMS/NHDPlus10L','rpu':'ned10b','df':'MO_U17_Missouri.shp', 'ib_filter':1, 'K_bedrock':.8 * ft2m, 'K_fine':.1 * ft2m, 'K_coarse':100 * ft2m}
model_dict['MO_U16_Mississippi'] = {'ws':'MO_U16_Mississippi','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07a','df':'MO_U16_Mississippi.shp', 'ib_filter':1, 'K_bedrock':.8 * ft2m, 'K_fine':.1 * ft2m, 'K_coarse':100 * ft2m}
# model_dict['IA_U01_Racoon'] = {'ws':'IA_U01_Racoon','vpu':'NHDPlusMS/NHDPlus07','rpu':'ned07b','df':'IA_U01_Racoon.shp', 'ib_filter':1, 'K_bedrock':.1 * ft2m}

# specify steady state
NPER = 1

# constants
hnoflo = -9999.
hdry = -8888.

