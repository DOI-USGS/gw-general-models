import os
import pandas as pd
import flopy as fp
import datetime as dt
import numpy as np
import scipy.stats as ss
import scipy.optimize as so
import geopandas as gpd
from shapely.geometry import LineString, MultiLineString, Point

# import warnings
# warnings.filterwarnings("ignore", message="converting a masked element to nan")

class RTD_util(object):
    '''Class to perform various functions from setting up an age-based backtracking MODPATH simulation
    and analyzing the results'''
    
    def __init__(self, sim, weight_label, group):
        # import various variables from the MODFLOW model
        self.ml = sim.get_model()
        self.sim_name = 'mfsim.nam'
        self.model_ws = self.ml.model_ws
        self.oc = self.ml.get_package('OC')
        self.dis = self.ml.get_package('DIS')
        self.npf = self.ml.get_package('NPF')
        self.tdis = sim.get_package('TDIS')
        self.namefile = self.ml.namefile
        self.prng = np.random.RandomState(9591029)
        self.delr = self.dis.delr.array
        self.delc = self.dis.delc.array
        self.nlay = self.dis.nlay.array
        self.nrow = self.dis.nrow.array
        self.ncol = self.dis.ncol.array
        self.l, self.r, self.c = np.indices((self.nlay, self.nrow, self.ncol))
        self.bot = self.dis.botm.array
        self.top = self.dis.top.array
        self.seqnums = np.arange(self.nlay * self.nrow * self.ncol)
        
        # self.hnoflo = self.bas.hnoflo
        # self.hdry = self.upw.hdry
        self.ibound = np.asarray(self.dis.idomain.array)
        self.hk = np.asarray(self.npf.k.array)
        self.vka = np.asarray(self.npf.k33.array)

        self.weight_label = weight_label
        self.group = group
        self.mpname = '{}_{}_{}'.format(self.ml.name, self.weight_label, self.group)
        
        self._len_mult()
        
        # Create dictionary of multipliers for converting model time units to years
        time_dict = dict()
        time_dict['unknown'] = 1.0 # undefined assumes days
        time_dict['seconds'] = 24 * 60 * 60 * 365.25
        time_dict['minutes'] = 24 * 60 * 365.25
        time_dict['hours'] = 24 * 365.25
        time_dict['days'] = 365.25
        time_dict['years'] = 1.0
        self.time_dict = time_dict
    
    def get_node(self, lrc_list):
        """
        Get node number from a list of MODFLOW layer, row, column tuples.

        Returns
        -------
        v : list of MODFLOW nodes for each layer (k), row (i),
            and column (j) tuple in the input list
        """
        if not isinstance(lrc_list, list):
            lrc_list = [lrc_list]
        nrc = self.nrow * self.ncol
        v = []
        for [k, i, j] in lrc_list:
            node = int(((k) * nrc) + ((i) * self.ncol) + j)
            v.append(node)
        return v

    def _get_output_dfs(self):
        # Make dataframes of budget information
        src = os.path.join(self.model_ws, '{}.cbb'.format(self.ml.name))
        self.bud_obj = fp.utils.CellBudgetFile(src, precision='double')
        all_bud_df = pd.DataFrame(self.bud_obj.recordarray)
        # convert to zero base
        all_bud_df['kper'] -= 1
        all_bud_df['kstp'] -= 1
        self.all_bud_df = all_bud_df
        
        headfile = '{}.hds'.format(self.ml.name)
        src = os.path.join(self.model_ws, headfile)
        self.hds = fp.utils.binaryfile.HeadFile(src, precision='double')       

    def _get_kstpkper(self, mf_start_date_str = '01/01/1900', mp_release_date_str = '01/01/2018' ):   
        # Use calendar release date and MODFLOW start date to pick out head and budget
        # items from transient model output
        self._get_output_dfs()
        
        # convert string representation of dates into Python datetime objects
        self.mf_start_date = dt.datetime.strptime(mf_start_date_str , '%m/%d/%Y')
        self.mp_release_date = dt.datetime.strptime(mp_release_date_str , '%m/%d/%Y')
    
        # check to make sure they are valid
        assert self.mf_start_date < self.mp_release_date, 'The particle release date has\
to be after the start of the MODFLOW simulation'
    
        # group by period and step
        kdf = self.all_bud_df.groupby(['kper', 'kstp']).median()
        kdf = kdf[['pertim', 'totim']]
    
        # make a datetime series for timesteps starting with 0
        # totim is elapsed time in simulation time
        units = self.tdis.time_units.array.lower()
        if units == 'days':
            units = 'D'
        end_date = self.mf_start_date + pd.to_timedelta(np.append(0, kdf.totim), unit=units)
        end_date = end_date.map(lambda t: t.strftime('%Y-%m-%d %H:%M'))
        kdf.loc[:, 'start_date'] = end_date[0:-1]
        kdf.loc[:, 'end_date'] = end_date[1:]
    
        # make a datetime series for timesteps starting with 0
        # totim is elapsed time in simulation time
        # reformat the dates to get rid of seconds
        end_date = self.mf_start_date + pd.to_timedelta(np.append(0, kdf.totim), unit=units)
        kdf.loc[:, 'start_date'] = end_date[0:-1].map(lambda t: t.strftime('%Y-%m-%d %H:%M'))
        kdf.loc[:, 'end_date'] = end_date[1:].map(lambda t: t.strftime('%Y-%m-%d %H:%M'))
    
        # reference time and date are set to the end of the last stress period
        self.ref_time = kdf.totim.max()
        self.ref_date = end_date.max()
    
        # release time is calculated in tracking time (for particle release) and 
        # in simulation time (for identifying head and budget components)
        self.release_time_trk = np.abs((self.ref_date - self.mp_release_date).days)
        self.release_time_sim = (self.mp_release_date - self.mf_start_date).days
    
        # find the latest group index that includes the release date
        idx = (kdf.totim >= self.release_time_sim).idxmax()
        kdf.loc[idx, 'particle_release'] = True
    
        # switch period and step 
        self.kstpkper = (idx[1], idx[0])
        
        assert self.ref_date > self.mp_release_date, 'The reference date has \
to be after the particle release'
    
    def get_heads(self):
        # Get the highest non-dry head in the 2D representation of the MODFLOW model
        # in each vertical stack of cells
        self._get_kstpkper()
        heads = self.hds.get_data(kstpkper=self.kstpkper)
        hd = heads.copy()
        hd[self.dis.idomain.array != 1] = np.nan
        # hd[np.isclose(self.upw.hdry, hd, atol=10)] = np.nan
        self.hd = hd
    
    def get_watertable(self):
        # Get the highest non-dry head in the 2D representation of the MODFLOW model
        # in each vertical stack of cells
        self.get_heads()
        hin = np.argmax(np.isfinite(self.hd), axis=0)    
        self.water_table =  np.squeeze(self.hd[hin, self.r[0,:,:], self.c[0,:,:]])
    
    def make_particle_array(self, parts_per_cell):  
        # Given the number of desired particles per cell, return an array in the 
        # format of MODPATH starting location information
        if not hasattr(self, 'release_time_trk'):
            self._get_kstpkper()
            
        self.parts_per_cell = parts_per_cell
        lg = self.l.ravel()
        rg = self.r.ravel()
        cg = self.c.ravel()
        label = parts_per_cell
    
        lrep = np.repeat( lg, parts_per_cell.ravel() )
        rrep = np.repeat( rg, parts_per_cell.ravel() )
        crep = np.repeat( cg, parts_per_cell.ravel() )
        label = np.repeat( label, parts_per_cell.ravel() )
    
        self.num_parts = lrep.shape[0]
    
        # generate random relative coordinates within a cell in 3D
        cell_coords = self.prng.rand( self.num_parts, 3 )
    
        grp = 1
    
        particles = np.zeros( ( self.num_parts, 11 ) )
        particles[:, 0] = np.arange( 1, self.num_parts + 1 )
        particles[:, 1] = grp
        particles[:, 2] = 1
        particles[:, 3] = lrep + 1
        particles[:, 4] = rrep + 1
        particles[:, 5] = crep + 1
        particles[:, 6:9] = cell_coords
        particles[:, 9] = self.release_time_trk
        particles[:, 10] = label
        
        return particles
        
    def make_arbitrary_particle_array(self, seqnum, label, parts_per_cell=1000, top_face=False):
        # Given the number of desired particles per cell, return an array in the 
        # format of MODPATH starting location information
        if not hasattr(self, 'release_time_trk'):
            self._get_kstpkper()
            
        self.parts_per_cell = parts_per_cell
        lg = self.l.ravel()[seqnum]
        rg = self.r.ravel()[seqnum]
        cg = self.c.ravel()[seqnum]
    
        lrep = np.repeat( lg, parts_per_cell)
        rrep = np.repeat( rg, parts_per_cell)
        crep = np.repeat( cg, parts_per_cell)
        label = np.repeat( label, parts_per_cell)
    
        self.num_parts = lrep.shape[0]
    
        # generate random relative coordinates within a cell in 3D
        cell_coords = self.prng.rand( self.num_parts, 3 )
        if top_face:
            cell_coords[:, 2] = 6
    
        grp = 1
    
        particles = np.zeros( ( self.num_parts, 11 ) )
        particles[:, 0] = np.arange( 1, self.num_parts + 1 )
        particles[:, 1] = grp
        particles[:, 2] = 1
        particles[:, 3] = lrep + 1
        particles[:, 4] = rrep + 1
        particles[:, 5] = crep + 1
        particles[:, 6:9] = cell_coords
        particles[:, 9] = self.release_time_trk
        particles[:, 10] = label
        
        return particles
    
    def write_starting_locations_file(self, particles):
        # Given a particle starting array, write a MODPATH starting location file with 
        # header information
        line = '{:5d}\n{:5d}\n'.format(1, 1)
        line = line + 'group_{}\n'.format(1)
        npart = particles.shape[0]
        line = line + '{:6d}'.format(npart)
        self.ep_file_name = os.path.join(self.model_ws, '{}_{}_{}'.format(self.ml.name, self.weight_label, self.group))
        form = '%6d %6d %3d %3d %3d %3d %12.9f %12.9f %12.9f %12.9e %15d'
        np.savetxt(self.ep_file_name+'.loc', particles, delimiter=' ', fmt=form, header=line, comments='')
        
    def run_MODPATH(self, por, mp_exe_name):
        # Run backtracking MODPATH simulation using a starting locations file
        # prepare Modpath files   
        SimulationType = 1              # 1 endpoint; 2 pathline; 3 timeseries
        TrackingDirection = 2           # 1 forward; 2 backward
        WeakSinkOption = 1              # 1 pass; 2 stop
        WeakSourceOption = 1            # 1 pass; 2 stop
        ReferemceTimeOption = 1         # 1 time value; 2 stress period, time step, relative offset
        StopOption = 2                  # 1 stop with simulation 2; extend if steady state 3; specify time
        ParticleGenerationOption = 2    # 1 automatic; 2 external file
        TimePointOption = 1             # 1 none; 2 number at fixed intervals; 3 array
        BudgetOutputOption = 3          # 1 none; 2 summary; 3 list of cells; 4 trace mode
        ZoneArrayOption = 1             # 1 none; 2 read zone array(s) 
        RetardationOption = 1           # 1 none; 2 read array(s) 
        AdvectiveObservationsOption = 1 # 1 none; 2 saved for all time pts 3; saved for final time pt

        options = [SimulationType, TrackingDirection, WeakSinkOption, WeakSourceOption, ReferemceTimeOption, 
                   StopOption, ParticleGenerationOption, TimePointOption, BudgetOutputOption, ZoneArrayOption, 
                   RetardationOption, AdvectiveObservationsOption]

        mpnf = '{}_{}_{}.mpnam'.format(self.ml.name, self.weight_label, self.group)
        mplf = '{}_{}_{}.mplst'.format(self.ml.name, self.weight_label, self.group)

        mp = fp.modpath.Modpath(modelname=self.mpname, modflowmodel=self.ml, dis_file=self.dis.file_name[0], exe_name=mp_exe_name,
                                model_ws=self.model_ws, simfile_ext='mpsim', dis_unit=self.dis.unit_number[0])

        mpsim = fp.modpath.ModpathSim(mp, mp_name_file=mpnf, 
                                      mp_list_file=mplf, 
                                      option_flags=options,
                                      ref_time=self.ref_time,
                                      cell_bd_ct=0, 
        #                               bud_loc=bud_chk_dict[group].loc[:, ('Grid', 'Layer', 'Row', 'Column')].values.tolist(),
                                      extension='mpsim')

        mpbas = fp.modpath.ModpathBas(mp, hnoflo=self.bas.hnoflo, hdry=self.upw.hdry, 
                                      def_face_ct=1, bud_label=['RECHARGE'], def_iface=[6], 
                                      laytyp=self.upw.laytyp.get_value(), ibound=self.bas.ibound.array, 
                                      prsity=por, prsityCB=0.20)    

        mp.write_input()
        success, msg = mp.run_model(silent=True, report=False)

        #     delete starting locations to save space--this information is now in the endpoint file
        if success:
            dst_pth = os.path.join(self.model_ws, '{}_{}_{}.loc'.format(self.ml.name, self.weight_label, self.group))
            os.remove(dst_pth)
            
    def modify_endpoint_file(self, ep_data_, write=False):  
        if not hasattr(self, 'water_table'):
            self.get_watertable()
        ep_data_ = ep_data_.copy()
        # Clean up and enhance an MODPATH endpoint file
        # set the Z coordinate for particles that end in dry cells to the 
        # head of the nearest non-dry cell below the dry cell.
        # ind = np.isclose(ep_data_.loc[:, 'Final Global Z'], self.upw.hdry, atol=100)
        # ep_data_.loc[:, 'Final Global Z'] = np.where(ind, self.water_table[ep_data_.loc[:, 'Final Row'] - 1, 
                                            # ep_data_.loc[:, 'Final Column']-1], ep_data_.loc[:, 'Final Global Z'])

        # eliminate particles that start in dry cells
        # ind = np.isclose(ep_data_.loc[:, 'Initial Global Z'], self.npf.hdry, rtol=0.99999)
        # self.ep_data = ep_data_.loc[~ind, :]

        # calculate approximate linear path distances
        x_dist = ep_data_.loc[:, 'Final global x'] - ep_data_.loc[:, 'Initial global x']
        y_dist = ep_data_.loc[:, 'Final global y'] - ep_data_.loc[:, 'Initial global y']
        z_dist = ep_data_.loc[:, 'Final global z'] - ep_data_.loc[:, 'Initial global z']
        ep_data_.loc[:, 'xy_path_len'] = np.sqrt(x_dist**2 + y_dist**2)
        ep_data_.loc[:, 'xyz_path_len'] = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)

        mendpoint_file = '{}_mod.{}'.format(self.mpname, 'mpend')
        mendpoint_file = os.path.join(self.model_ws, mendpoint_file)
        if write:
            ep_data_.to_csv(mendpoint_file)
            endpoint_file = '{}.{}'.format(self.mpname, 'mpend')
            endpoint_file = os.path.join(self.model_ws, endpoint_file)
            if os.path.exists(endpoint_file):
                os.remove(endpoint_file)
        self.ep_data = ep_data_

    def get_budget(self, text):
        # Get the MODFLOW budget file for the time period specified by the MODPATH release date
        # and the MODFLOW start date. 
        self._get_kstpkper()
        # budget = self.bud_obj.get_data(kstpkper=self.kstpkper, text=text, full3D=True)[0]
        budget = self.bud_obj.get_data(kstpkper=self.kstpkper, text=text)
        self.budget = budget
        
    def _len_mult(self):
        # the database values are in feet; if the model is in meters, 
        # provide a multiplier to convert database values to match the model
        lenuni_dict = {0: 'undefined units', 1: 'feet', 2: 'meters', 3: 'centimeters'}
        if self.dis.length_units.array == 'meters':
            self.len_mult = 0.3048006096012192
        elif self.dis.length_units.array == 'feet':
            self.len_mult = 1.0
        else:
            print('unknown length units')
            self.len_mult = 1.0

    def read_endpoints(self, endpoint_file):
        # read MODPATH 6 endpoint file
        # count the number of header lines
        i = 0
        with open(endpoint_file) as f:
            while True:
                line = f.readline()
                i += 1
                if 'END HEADER' in line:
                    break
                elif not line:
                    break

        # columns names from MP6 docs 
        cols = ['Sequence number', 'Particle Group', 'Particle ID', 'Status', 'Initial tracking time', 'Final tracking time', 'Initial cell number', 'Initial layer', 'Initial local x', 'Initial local y', 'Initial local z', 'Initial global x', 'Initial global y', 'Initial global z', 'Initial zone', 'Initial face', 'Final cell number','Final layer', 'Final local x', 'Final local y', 'Final local z', 'Final global x', 'Final global y', 'Final global z', 'Final zone', 'Final face']            

        # read the endpoint data
        ep_data = pd.read_csv(endpoint_file, names=cols, header=None, skiprows=i, delim_whitespace=True)

        # select only 'Normally Terminated' particles; status code = 2
        ep_data = ep_data.loc[ep_data.Status == 2, :]

        # tmp = ep_data[['Initial layer', 'Initial Row', 'Initial Column']] - 1
        ep_data['initial_node_num'] = ep_data['Initial cell number'] - 1

        # tmp = ep_data[['Final Layer', 'Final Row', 'Final Column']] - 1
        ep_data['final_node_num'] = ep_data['Final cell number'] - 1
        
        # calculate particle travel time in years
        ep_data['rt'] = (ep_data['Final tracking time'] - ep_data['Initial tracking time']) / self.time_dict[self.tdis.time_units.array]
        ep_data.set_index('initial_node_num', drop=True, inplace=True)
        return ep_data   
        
    def _distfit(self, t, sh_e, sc_e):
        first = t.min()
        _cdf = self.dist.cdf(t, sh_e, first, sc_e)
        return  _cdf

    def _explicit(self, t, sh_e, sc_e, sh_l, sc_l, fy):
        first = t.min()
        _cdf_e = self.dist.cdf(t, sh_e, first, sc_e)
        _cdf_l = self.dist.cdf(t, sh_l, first, sc_l)
        return fy * _cdf_e + (1 - fy) * _cdf_l

    def _fit(self, lprt, ly, func, components, bnds):
        for dist in self.dist_list:
            self.dist = dist
            label = '{}_{}'.format(components, dist.name)
            if components == 'one':
                numpar = 2
            if components == 'two':
                numpar = 5

            try:
                # with np.errstate(divide='ignore',invalid='ignore'):
                par1, cov1 = so.curve_fit(func, lprt, ly, bounds = bnds, method='trf')
                cdf1 = func(lprt, *par1)
                resid1 = ly - cdf1
                ssr1 = resid1.T.dot(resid1)
                error_message1 = 'Solved using TRF -- no warnings'
            except Exception as e: 
                ssr1 = np.inf
                error_message1 = 'TRF with warning: {}'.format(e.args)
                cov1 = np.nan

            try:
                # with np.errstate(divide='ignore',invalid='ignore'):
                par2, cov2 = so.curve_fit(func, lprt, ly, bounds = bnds, method='dogbox')
                cdf2 = func(lprt, *par2)
                resid2 = ly - cdf2
                ssr2 = resid2.T.dot(resid2)
                error_message2 = 'Solved using Dogbox -- no warnings'
            except Exception as e: 
                ssr2 = np.inf
                error_message2 = 'Dogbox with warning: {}'.format(e.args)
                cov2 = np.nan

            if ssr1 < ssr2:
                par = par1
                cov = cov1
                cdf = cdf1
                rmse = np.sqrt(ssr1 / self.s)
                method = 'TRF'
                error_message = error_message1
                
            elif ssr1 > ssr2:
                par = par2
                cov = cov2
                cdf = cdf2
                rmse = np.sqrt(ssr2 / self.s)
                method = 'Dogbox'
                error_message = error_message2
                
            else:
                par = np.zeros((numpar))
                cov = np.zeros((numpar, numpar))
                cdf = np.zeros((lprt.shape[0]))
                rmse = np.nan
                method = 'TRF and Dogbox did not solve'   
                error_message = 'No solution with TRF or Dogbox'
                
            return par, cdf, cov, rmse, method, label, error_message

    def fit_dists(self, ly, lprt, dist_list, fit_one=True, fit_two=False):
        # fit a list of distributions in either one or two components
        # return a dictionary of results
        self.dist_list = dist_list
        param_dict = {}
        cdf_dict = {}
        error_dict = {}
        tt_dict = {}

        # dist_list = [ss.invgauss, ss.gamma, ss.weibull_min] 
        first = lprt.min()
        s = ly.shape[0]
        tt_dict['rt'] = lprt
        tt_dict['rt_cdf'] = ly
        self.s = ly.shape[0]

        if fit_one:
            bnds = (0, [+np.inf, +np.inf])
            func = self._distfit
            
            par_one, cdf_one, cov_one, rmse_one, method_one, label, error_message_one = self._fit(lprt, ly, func, 'one', bnds)
            
            par_one = par_one[0:2]
            par_one = np.insert(par_one, 1, first)
            param_dict[label] = par_one
            error_dict[label] = (rmse_one, cov_one, error_message_one, method_one)
            cdf_dict[label] = cdf_one

        if fit_two:
            bnds = (0, [+np.inf, +np.inf, +np.inf, +np.inf, 1.0]) 
            func = self._explicit
            
            par_two, cdf_two, cov_two, rmse_two, method_two, label, error_message_two = self._fit(lprt, ly, func, 'two', bnds)
            
            par_two = np.insert(par_two, 1, first)
            par_two = np.insert(par_two, 4, first)
            param_dict[label] = par_two
            error_dict[label] = (rmse_two, cov_two, error_message_two, method_two)
            cdf_dict[label] = cdf_two

        return {'cdf' : cdf_dict, 'par' : param_dict, 'err' : error_dict, 'tt' : tt_dict}

    def make_stream_particle_array(self, shp_drn_df, data, particles_per_flow, seg_ref='comid'):
        '''
        Creates local cell stream coordinates given a stream shapefile that has been intersected with the model grid.

        Parameters:
            shp_drn_df: GeoDataFrame (from reading a shapefile with geopandas)
                    Each segment of a reach within a cell has to be attributed with the flow (q),
                    length (lengths), and the node number (node). 
                    Multiple segments withing a cell are allowed, including multipart streams.
            data: DataFrame
                    Created in General Simulation Model notebook #1, or can be constructed manually, as for the test 
                    problem. Has one row for each model grid cell in 2D (plan view), so the number of rows is 
                    nrow * ncol. Has to be attributed with node number (node) and cell boundary coordinates: 
                    x coord of left edge (xvert_l) , x coord of right edge (xvert_r),
                    y coord of bottom edge (yvert_b), and y coord of top edge (yvert_t).
            total_number_of_particles: integer
                    The total number of particles desired for the simulation
            total_flow: float
                    The total flow rate for the drain package
            seg_ref: int32
                    A field in shp_drn_df that contains a unique identifier for stream reaches.
                    In the NHDPlus 2.0, this is "comid". In the NHDPlus High Resolution, it
                    could be reachcode or NHDPlusID or something else. It has to be no larger
                    than about 9-10 digits for FloPy to accept it, and both these choices have
                    more digits. Some of tdigits aren't needed for the identifier to be unique
                    and can be truncated.

        Returns:
            x_partloc, y_partloc: lists
                    Relative cell coordinates (range 0-1) in x and y corresponding to equally spaced
                    particles along streams. Each cell receives a number of particles in proportion to 
                    the surface-water flow out of the cell (in this case drains).
        '''
        x_partloc = list()
        y_partloc = list()
        node_list = list()
        label_list = list()
        
        shp_grp = shp_drn_df.groupby('node')

        def make_particle_locations(k, data, cell_q, node, total_length, comid):
            x, y = k.xy
            rel_x = (x - data.loc[node, 'xvert_l']) / \
                (data.loc[node, 'xvert_r'] - data.loc[node, 'xvert_l'])
            rel_y = (y - data.loc[node, 'yvert_b']) / \
                (data.loc[node, 'yvert_t'] - data.loc[node, 'yvert_b'])
            xp = list()
            yp = list()
            num_lines = rel_x.shape[0] - 1
            for n in range(num_lines):
                seg_len = np.hypot(x[n+1] - x[n], y[n+1] - y[n]) / total_length
                nprt = np.rint(particles_per_flow * cell_q * seg_len).astype(np.int32)
                xp.extend(np.linspace(
                    rel_x[n], rel_x[n+1], nprt, endpoint=True).tolist())
                yp.extend(np.linspace(
                    rel_y[n], rel_y[n+1], nprt, endpoint=True).tolist())
                node_list.extend(np.repeat(node, nprt))
                label_list.extend(np.repeat(comid, nprt))
            return xp, yp

        for node, segments_in_node in shp_grp:
            total_length = segments_in_node.lengths.sum()
            cell_q = segments_in_node.q.mean()

            for p, segment_geometry in segments_in_node.iterrows():
                comid = segment_geometry[seg_ref]
                line_geometry = segment_geometry.geometry

                if line_geometry.geom_type == 'LineString':
                    line_coords = line_geometry.coords
                    xp, yp = make_particle_locations(
                        line_coords, data, cell_q, node, total_length, comid)
                    x_partloc.extend(xp)
                    y_partloc.extend(yp)

                elif line_geometry.geom_type == 'MultiLineString':
                    for line_part in line_geometry:
                        line_coords = line_part.coords
                        xp, yp = make_particle_locations(
                            line_coords, data, cell_q, node, total_length, comid)
                        x_partloc.extend(xp)
                        y_partloc.extend(yp)

                else:
                    'non-line geometry'

        return x_partloc, y_partloc, node_list, label_list

    def run_test(self, num_particles):
        L1 = LineString(((1000, 550, 1), (1050, 550, 1), (1050, 600, 1)))
        L2 = LineString(((1050, 525, 1), (1050, 550, 1)))
        L3a = LineString(((1075, 500, 1), (1100, 525, 1)))
        L3b = LineString(((1100, 575, 1), (1075, 575, 1), (1075, 600, 1)))
        L3 = MultiLineString((L3a, L3b))
        L4 = LineString(((1100, 525, 1), (1150, 575, 1), (1100, 575, 1)))
        gs = gpd.GeoSeries((L1, L2, L3, L4))

        test_shp = gpd.GeoDataFrame({'geometry': gs,
                                     'comid': [0, 1, 2, 2],
                                     'q': [100, 100, 100, 200],
                                     'node': [1000, 1000, 1000, 2000]})
        test_shp['lengths'] = test_shp.length

        test_data = pd.DataFrame({'xvert_l': [1000, 1100], 'xvert_r': [
                                 1100, 1200], 'yvert_b': [500, 500], 'yvert_t': [600, 600], 'node': [1000, 2000]})
        test_data.set_index('node', inplace=True)
            
        x_partloc, y_partloc, node_list, comid = self.make_stream_particle_array(test_shp, test_data, num_particles)

        xx = list()
        yy = list()


        for i in range(len(node_list)):
            n = node_list[i]
            xx.append(x_partloc[i] * 100 + test_data.loc[n, 'xvert_l'])
            yy.append(y_partloc[i] * 100 + test_data.loc[n, 'yvert_b'])

        test_points = list()
        for i in (zip(xx, yy)):
                  test_points.append(Point(i))

        points_gdf = gpd.GeoDataFrame({'geometry': test_points})

        ax = points_gdf.plot(markersize=40, marker='+', c='k', zorder=2, linewidth=1)
        test_shp.plot(ax=ax, column='comid', zorder=1, alpha=0.8) 

        for l in [1000, 1100, 1200]:
            ax.axvline(l, c='k', ls='dashed', lw=0.8)
            
        for l in [500, 600]:
            ax.axhline(l, c='k', ls='dashed', lw=0.8)
