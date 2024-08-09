#!/usr/bin/env python

import argparse
import logging
import f90nml
import os
from netCDF4 import Dataset
import numpy as np
from datetime import datetime, timedelta


###############################################################################
# Global settings                                                             #
###############################################################################

# Path to the directory containing processed case input files
CASE_NML_DIR = '../case_config'

# Path to the directory containing processed case input files
PROCESSED_CASE_DIR = '../../data/processed_case_input'

DEFAULT_MISSING_VALUE = -9999.0
DEFAULT_NUDGING_TIMESCALE = 7200.0 #s

###############################################################################
# Command line arguments                                                      #
###############################################################################

parser = argparse.ArgumentParser()
parser.add_argument('-n',   '--case_name',       help='name of case',         required=True)
parser.add_argument('-a',   '--use_area',        help='use column_area namelist attribute as forcing_scale',  action='store_true')
parser.add_argument('-d',   '--debug',           help='enable debugging output', action='store_true')


########################################################################################
#
########################################################################################
def parse_arguments():
    """Parse command line arguments"""
    args           = parser.parse_args()
        
    return (args.case_name, args.use_area, args.debug)

########################################################################################
#
########################################################################################
def setup_logging(debug):
    """Sets up the logging module."""
    if debug:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)
    else:
        logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.DEBUG)


def get_case_nml(case_name):
    """Returns case configuration Fortran namelist"""
    
    filename = os.path.join(CASE_NML_DIR, case_name + '.nml')
    
    print(filename)
    
    error = False
    nml = ''
    if (os.path.exists(filename)):
        nml = f90nml.read(filename)
    else:
        error = True
    
    return (nml, error)

def get_case_data(case_name):
    """Returns proprietery CCPP SCM case data in NetCDF Dataset format"""
    
    #TODO: need to handle LSM ICs
    
    filename = os.path.join(PROCESSED_CASE_DIR, case_name + '.nc')
    
    error = False
    try:
        nc_fid = Dataset(filename , 'r')
    except:
        error = True
        
    if (not error):
        # Store case data in a dict for easy loop assignment
        case_data_dict = dict()

        #read global variable(s)
        try:
            case_data_dict['missing_value'] = nc_fid.getncattr('missing_value')
        except:
            case_data_dict['missing_value'] = DEFAULT_MISSING_VALUE
        
        case_data_dict['time']   = nc_fid.variables['time'][:]
        case_data_dict['levels'] = nc_fid.variables['levels'][:]
        try:
            case_data_dict['soil_depth'] = nc_fid.variables['soil_depth'][:]
        except KeyError:
            case_data_dict['soil_depth'] = case_data_dict['missing_value']
        
        #read variables from scalar group; lat and lon can not be missing
        scalars_grp = nc_fid.groups['scalars']
        case_data_dict['lat'] = scalars_grp.variables['lat'][:]
        case_data_dict['lon'] = scalars_grp.variables['lon'][:]


        scalars = [ 'slmsk',
                    'vegsrc',
                    'vegtyp',
                    'soiltyp',
                    'scolor',
                    'slopetyp',
                    'tsfco',
                    'vegfrac',
                    'shdmin',
                    'shdmax',
                    'canopy',
                    'hice',
                    'fice',
                    'tisfc',
                    'snowd',
                    'snoalb',
                    'tg3',
                    'uustar',
                    'alvsf',
                    'alnsf',
                    'alvwf',
                    'alnwf',
                    'facsf',
                    'facwf',
                    'weasd',
                    'f10m',
                    't2m',
                    'q2m',
                    'ffmm',
                    'ffhh',
                    'tprcp',
                    'srflag',
                    'sncovr',
                    'tsfcl',
                    'zorl',
                    'zorll',
                    'zorli',
                    'zorlw',
                    'zorlwav',
                    'tvxy',
                    'tgxy',
                    'tahxy',
                    'canicexy',
                    'canliqxy',
                    'eahxy',
                    'cmxy',
                    'chxy',
                    'fwetxy',
                    'sneqvoxy',
                    'alboldxy',
                    'qsnowxy',
                    'wslakexy',
                    'taussxy',
                    'waxy',
                    'wtxy',
                    'zwtxy',
                    'xlaixy',
                    'xsaixy',
                    'lfmassxy',
                    'stmassxy',
                    'rtmassxy',
                    'woodxy',
                    'stblcpxy',
                    'fastcpxy',
                    'smcwtdxy',
                    'deeprechxy',
                    'rechxy',
                    'snowxy',
                    'wetness',
                    'clw_surf_land',
                    'clw_surf_ice',
                    'qwv_surf_land',
                    'qwv_surf_ice',
                    'tsnow_land',
                    'tsnow_ice',
                    'snowfallac_land',
                    'snowfallac_ice',
                    'sncovr_ice',
                    'sfalb_lnd',
                    'sfalb_lnd_bck',
                    'emis_ice',
                    'lai',
                    'area',
                    'stddev',
                    'convexity',
                    'oa1',
                    'oa2',
                    'oa3',
                    'oa4',
                    'ol1',
                    'ol2',
                    'ol3',
                    'ol4',
                    'theta_oro',
                    'gamma',
                    'sigma',
                    'elvmax',
                    'oro',
                    'oro_uf',
                    'landfrac',
                    'lakefrac',
                    'lakedepth',
                    'tref',
                    'z_c',
                    'c_0',
                    'c_d',
                    'w_0',
                    'w_d',
                    'xt',
                    'xs',
                    'xu',
                    'xv',
                    'xz',
                    'zm',
                    'xtts',
                    'xzts',
                    'd_conv',
                    'ifd',
                    'dt_cool',
                    'qrains']


        for var in scalars:
            try:
                case_data_dict[var] = scalars_grp.variables[var][:]
            except KeyError:
                case_data_dict[var] = case_data_dict['missing_value']

        # snowd is a special case

        try:
            case_data_dict['snowd'] = scalars_grp.variables['snowd'][:]
        except KeyError:
            try:
                case_data_dict['snowd'] = scalars_grp.variables['snwdph'][:]
            except KeyError:
                case_data_dict['snowd'] = case_data_dict['missing_value']


        #read variables from initial group
        initial_grp = nc_fid.groups['initial']

        # Some variables can't be missing

        initial_nomissing = ['qt','ql','qi','u','v','tke','ozone']

        for var in initial_nomissing:
            case_data_dict[var] = initial_grp.variables[var][:]

        initial = ['height',
                   'thetail',
                   'temp',
                   'stc',
                   'smc',
                   'slc',
                   'snicexy',
                   'snliqxy',
                   'tsnoxy',
                   'smoiseq',
                   'zsnsoxy',
                   'tiice',
                   'tslb',
                   'smois',
                   'sh2o',
                   'smfr',
                   'flfr']


        for var in initial_nomissing:
            try:
                case_data_dict[var] = initial_grp.variables[var][:]
            except KeyError:
                case_data_dict[var] = case_data_dict['missing_value']

        
        #read variables from forcing group
        forcing_grp = nc_fid.groups['forcing']

        case_data_dict['p_surf'] = forcing_grp.variables['p_surf'][:]
        case_data_dict['T_surf'] = forcing_grp.variables['T_surf'][:]
        case_data_dict['w_ls'] = forcing_grp.variables['w_ls'][:]
        case_data_dict['omega'] = forcing_grp.variables['omega'][:]
        case_data_dict['u_g'] = forcing_grp.variables['u_g'][:]
        case_data_dict['v_g'] = forcing_grp.variables['v_g'][:]
        case_data_dict['u_nudge'] = forcing_grp.variables['u_nudge'][:]
        case_data_dict['v_nudge'] = forcing_grp.variables['v_nudge'][:]
        case_data_dict['T_nudge'] = forcing_grp.variables['T_nudge'][:]
        case_data_dict['thil_nudge'] = forcing_grp.variables['thil_nudge'][:]
        case_data_dict['qt_nudge'] = forcing_grp.variables['qt_nudge'][:]
        case_data_dict['dT_dt_rad'] = forcing_grp.variables['dT_dt_rad'][:]
        case_data_dict['h_advec_thil'] = forcing_grp.variables['h_advec_thetail'][:]
        case_data_dict['v_advec_thil'] = forcing_grp.variables['v_advec_thetail'][:]
        case_data_dict['h_advec_qt'] = forcing_grp.variables['h_advec_qt'][:]
        case_data_dict['v_advec_qt'] = forcing_grp.variables['v_advec_qt'][:]


        # Only two forcing fields can be missing
        try:
            case_data_dict['sh_flux_sfc'] = forcing_grp.variables['sh_flux_sfc'][:]
        except KeyError:
            case_data_dict['sh_flux_sfc'] = ''    
        try:
            case_data_dict['lh_flux_sfc'] = forcing_grp.variables['lh_flux_sfc'][:]
        except KeyError:
            case_data_dict['lh_flux_sfc'] = ''
        
        nc_fid.close()
    
    return(case_data_dict, error)

def write_SCM_case_file(case_nml, case_data, use_area):
    """Write all data to a netCDF file in the DEPHY-SCM format"""
    
    #TODO: need to handle LSM ICs
    
    # Working types
    wp = np.float64
    wi = np.int32
    
    # Local switches
    forcing_on  = 1
    forcing_off = 0
    
    nml_keys = case_nml['case_config'].todict().keys()
    nml_filename = os.path.join(CASE_NML_DIR, case_nml['case_config']['case_name'] + '.nml')
    
    # Output file
    com = 'mkdir -p ' + PROCESSED_CASE_DIR
    logging.info(com)
    os.system(com)
    fileOUT = os.path.join(PROCESSED_CASE_DIR, case_nml['case_config']['case_name'] + '_dephy' + '_SCM_driver.nc')
    
    nc_file = Dataset(fileOUT, 'w', format='NETCDF3_CLASSIC')
    nc_file.description = "Case data for {} from CCPP SCM".format(case_nml['case_config']['case_name'])

    nc_file.missing_value   = case_data['missing_value']
    
    #not all namelists will have minutes, set to 0 if nml doesn't have
    try:
        minute = case_nml['case_config']['minute']
    except KeyError:
        minute = 0    
    
    start_date = datetime(case_nml['case_config']['year'],case_nml['case_config']['month'],case_nml['case_config']['day'],case_nml['case_config']['hour'],minute,0)
    start_date_string = start_date.strftime("%Y-%m-%d %H:%M:%S")
    runtime = case_nml['case_config']['runtime']
    delta = timedelta(seconds=runtime)
    end_date = start_date + delta
    end_date_string   = end_date.strftime("%Y-%m-%d %H:%M:%S")
    loc_string  = str(case_data['lon']) + "E" + str(case_data['lat']) + "N"
    case_string = case_nml['case_config']['case_name'] + '_' + start_date_string + '_' + loc_string
    
    logging.debug('Case string: {}'.format(case_string))
    logging.debug('Case start date: {}'.format(start_date))
    logging.debug('Case duration: {}'.format(delta))
    logging.debug('Case end date: {}'.format(end_date))
    
    if (case_nml['case_config']['sfc_type'] > 1.5):
        surface_string = 'ice'
    elif (case_nml['case_config']['sfc_type'] > 0.5):
        surface_string = 'land'
    else:
        surface_string = 'ocean'
    
    #override case nml with LSM/model data
    if ('lsm_ics' in nml_keys):
        if (case_nml['case_config']['lsm_ics']):
            if (case_data['slmsk'] > 1.5):
                surface_string = 'ice'
            elif (case_data['slmsk'] > 0.5):
                surface_string = 'land'
            else:
                surface_string = 'ocean'
    
    #DEPHY v1 format specifies the global attributes in this order. Some attributes are rewritten below after the order is established in the file.
    nc_file.case              = case_string
    nc_file.title             = 'Forcing and Initial Conditions for ' + case_string
    nc_file.reference         = 'https://dtcenter.org/sites/default/files/paragraph/scm-ccpp-guide-v6-0-0.pdf'
    nc_file.author            = 'Grant J. Firl and Dustin Swales'
    nc_file.version           = 'Created on ' + datetime.today().strftime('%Y-%m-%d %H:%M:%S')
    nc_file.format_version    = 'DEPHY SCM format version 1'
    nc_file.modifications     = ''
    nc_file.script            = os.path.basename(__file__)
    nc_file.comment           = 'converted from ' + case_nml['case_config']['case_name'] + '.nc'
    nc_file.start_date        = start_date_string
    nc_file.end_date          = end_date_string
    
    if (use_area and case_nml['case_config']['column_area']):
        nc_file.forcing_scale     = case_nml['case_config']['column_area']
    else:    
        nc_file.forcing_scale     = -1
        
    nc_file.adv_ta            = forcing_off
    nc_file.adv_qv            = forcing_off
    nc_file.adv_ua            = forcing_off  #no mom_forcing_type implemented this (providing pre-calculated advective terms for u)
    nc_file.adv_va            = forcing_off  #no mom_forcing_type implemented this (providing pre-calculated advective terms for v)
    nc_file.adv_theta         = forcing_off
    nc_file.adv_thetal        = forcing_off
    nc_file.adv_qt            = forcing_off
    nc_file.adv_rv            = forcing_off 
    nc_file.adv_rt            = forcing_off
    nc_file.radiation         = "on"         #not implemented in CCPP SCM - controlled by CCPP SDF and/or namelist
    nc_file.forc_wap          = forcing_off
    nc_file.forc_wa           = forcing_off
    nc_file.forc_geo          = forcing_off
    nc_file.nudging_ua        = forcing_off
    nc_file.nudging_va        = forcing_off
    nc_file.nudging_ta        = forcing_off
    nc_file.nudging_theta     = forcing_off
    nc_file.nudging_thetal    = forcing_off
    nc_file.nudging_qv        = forcing_off
    nc_file.nudging_qt        = forcing_off
    nc_file.nudging_rv        = forcing_off
    nc_file.nudging_rt        = forcing_off
    nc_file.zh_nudging_ta     = forcing_off
    nc_file.zh_nudging_theta  = forcing_off
    nc_file.zh_nudging_thetal = forcing_off
    nc_file.zh_nudging_qv     = forcing_off
    nc_file.zh_nudging_qt     = forcing_off
    nc_file.zh_nudging_rv     = forcing_off
    nc_file.zh_nudging_rt     = forcing_off
    nc_file.zh_nudging_ua     = forcing_off
    nc_file.zh_nudging_va     = forcing_off
    nc_file.pa_nudging_ta     = forcing_off
    nc_file.pa_nudging_theta  = forcing_off
    nc_file.pa_nudging_thetal = forcing_off
    nc_file.pa_nudging_qv     = forcing_off
    nc_file.pa_nudging_qt     = forcing_off
    nc_file.pa_nudging_rv     = forcing_off
    nc_file.pa_nudging_rt     = forcing_off
    nc_file.pa_nudging_ua     = forcing_off
    nc_file.pa_nudging_va     = forcing_off
    #
    nc_file.surface_type      = surface_string
    nc_file.surface_forcing_temp     = 'none'
    nc_file.surface_forcing_moisture = 'none'
    nc_file.surface_forcing_wind     = 'none'
    nc_file.surface_forcing_lsm      = 'none'
    
    #rewrite forc_wa, forc_wap, forc_geo, nudging_ua, nudging_va depending on mom_forcing_type provided in case_config nml    
    if (case_nml['case_config']['mom_forcing_type'] == 2):
        #CCPP SCM proprietery forcing interprets mom_forcing_type = 2 as calculating vertical advective terms from provided vertical velocity AND applying geostrophic winds
        
        #CCPP SCM proprietery cases could have either w or omega available (or both); use omega by default?
        w_ls_avail = True if np.any(case_data['w_ls'][:,:]) else False
        omega_avail = True if np.any(case_data['omega'][:,:]) else False
        if omega_avail:
            nc_file.forc_wap = forcing_on
            nc_file.forc_wa  = forcing_off
        elif w_ls_avail:
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the momentum variables should be forced using vertical velocity (through mom_forcing = 2), but neither w_ls or omega have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_off
        
        geostrophic_avail = True if (np.any(case_data['u_g'][:,:]) or np.any(case_data['v_g'][:,:])) else False
        if geostrophic_avail:
            nc_file.forc_geo = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the momentum variables should be forced using geostrophic winds (through mom_forcing = 2), but neither u_g or v_g have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_geo = forcing_off
        nc_file.nudging_ua = forcing_off
        nc_file.nudging_va = forcing_off
                        
    elif (case_nml['case_config']['mom_forcing_type'] == 3):
        #CCPP SCM proprietery forcing interprets mom_forcing_type = 3 as calculating momentum forcing as nudging toward u and v profiles (only)
        
        nc_file.forc_wa  = forcing_off
        nc_file.forc_wap = forcing_off
        nc_file.forc_geo = forcing_off
        
        u_nudge_avail = True if np.any(case_data['u_nudge'][:,:]) else False
        v_nudge_avail = True if np.any(case_data['v_nudge'][:,:]) else False
        relax_time_avail = True if case_nml['case_config']['relax_time'] else False
        if relax_time_avail:
            relax_time = case_nml['case_config']['relax_time']
        else:
            relax_time = DEFAULT_NUDGING_TIMESCALE
            message = 'The case namelist ({0}) specifies the momentum variables should be forced using nudging (through mom_forcing = 3), but relax_time was not provided -- using default value of {1}s '.format(nml_filename, DEFAULT_NUDGING_TIMESCALE)
            logging.info(message)
        
        nc_file.nudging_ua = forcing_on*relax_time
        nc_file.nudging_va = forcing_on*relax_time
                
    if (case_nml['case_config']['thermo_forcing_type'] == 1):
        #total advective forcing + radiative heating
        nc_file.adv_thetal        = forcing_on
        nc_file.adv_qt            = forcing_on
        #nc_file.radiation         = 'tend'  #radiation isn't turned off in the CCPP SCM through the forcing
    elif (case_nml['case_config']['thermo_forcing_type'] == 2):
        #horizontal advective forcing + (radiative heating) + vertical velocity
        nc_file.adv_thetal        = forcing_on
        nc_file.adv_qt            = forcing_on
        #nc_file.radiation         = 'tend'  #radiation isn't turned off in the CCPP SCM through the forcing
        
        w_ls_avail = True if np.any(case_data['w_ls'][:,:]) else False
        omega_avail = True if np.any(case_data['omega'][:,:]) else False
        if omega_avail:
            nc_file.forc_wap = forcing_on
            nc_file.forc_wa  = forcing_off
        elif w_ls_avail:
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the thermo variables should be forced using vertical velocity (through thermo_forcing = 2), but neither w_ls or omega have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_off
            
    elif (case_nml['case_config']['thermo_forcing_type'] == 3):
        #nudging + vertical velocity
        
        T_nudge_avail = True if np.any(case_data['T_nudge'][:,:]) else False
        qt_nudge_avail = True if np.any(case_data['qt_nudge'][:,:]) else False
        relax_time_avail = True if case_nml['case_config']['relax_time'] else False
        if relax_time_avail:
            relax_time = case_nml['case_config']['relax_time']
        else:
            relax_time = DEFAULT_NUDGING_TIMESCALE
            message = 'The case namelist ({0}) specifies the thermo variables should be forced using nudging (through thermo_forcing = 3), but relax_time was not provided -- using default value of {1}s '.format(nml_filename, DEFAULT_NUDGING_TIMESCALE)
            logging.info(message)
        
        nc_file.nudging_ta = forcing_on*relax_time
        nc_file.nudging_qt = forcing_on*relax_time
        
        w_ls_avail = True if np.any(case_data['w_ls'][:,:]) else False
        omega_avail = True if np.any(case_data['omega'][:,:]) else False
        if omega_avail:
            nc_file.forc_wap = forcing_on
            nc_file.forc_wa  = forcing_off
        elif w_ls_avail:
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_on
        else:
            message = 'The case namelist ({0}) specifies the thermo variables should be forced using vertical velocity (through thermo_forcing = 2), but neither w_ls or omega have nonzero values'.format(nml_filename)
            logging.info(message)
            nc_file.forc_wap = forcing_off
            nc_file.forc_wa  = forcing_off
    
    time_dim   = nc_file.createDimension('time', case_data['time'].shape[0])
    timei_dim  = nc_file.createDimension('t0',    1)
    lev_dim    = nc_file.createDimension('lev',  case_data['levels'].shape[0])
    
    if ('lsm_ics' in nml_keys or 'model_ics' in nml_keys):
        nc_file.surface_forcing_lsm      = 'lsm'
        if (('lsm_ics' in nml_keys and case_nml['case_config']['lsm_ics']) or ('model_ics' in nml_keys and case_nml['case_config']['model_ics'])):
            if (case_data['soil_depth'][0] != case_data['missing_value']):
                soil_dim   = nc_file.createDimension('nsoil', case_data['soil_depth'].shape[0])
            else:
                message = 'LSM ICs are expected from the case_nml file, but no soil depth is provided.'
                logging.critical(message)
                raise Exception(message)
            
            if (case_data['snicexy'][0] != case_data['missing_value']):
                snow_dim   = nc_file.createDimension('nsnow', case_data['snicexy'].shape[0])
                nslsnw_dim = nc_file.createDimension('nsoil_plus_nsnow',case_data['snicexy'].shape[0] + case_data['soil_depth'].shape[0])
            
            if (case_data['tiice'][0] != case_data['missing_value']):
                ice_dim    = nc_file.createDimension('nice',  case_data['tiice'].shape[0])
    
    #
    timei_var                    = nc_file.createVariable('t0', wp, ('t0'))
    timei_var.units              = 'seconds since ' + start_date_string
    timei_var.standard_name      = 'Initial time'
    timei_var.calendar           = 'gregorian'
    timei_var[:]                 = 0.0
    #
    timef_var                    = nc_file.createVariable('time', wp, ('time'))
    timef_var.units              = 'seconds since ' + start_date_string
    timef_var.standard_name      = 'Forcing time'
    timef_var.calendar           = 'gregorian'
    timef_var[:]                 = case_data['time'][:]
    #
    lev_var                      = nc_file.createVariable('lev', wp, ('lev'))
    lev_var.units                = 'Pa'
    lev_var.standard_name        = 'pressure'
    lev_var[:]                   = case_data['levels'][:]
    
    if (case_data['soil_depth'][0] != case_data['missing_value']):
        soil_depth_var               = nc_file.createVariable('soil_depth', wp, ('nsoil'))
        soil_depth_var.units         = 'm'
        soil_depth_var.standard_name = 'depth of bottom of soil layers'
        soil_depth_var[:]            = case_data['soil_depth'][:]
    
    #
    lon_var                      = nc_file.createVariable('lon', wp, ('time'))
    lon_var.units                = 'degrees_east'
    lon_var.standard_name        = 'longitude'
    lon_var[:]                   = case_data['lon']

    #
    lat_var                      = nc_file.createVariable('lat', wp, ('time'))
    lat_var.units                = 'degrees_north'
    lat_var.standard_name        = 'latitude'
    lat_var[:]                   = case_data['lat']
    
    if (case_data['slmsk'] != case_data['missing_value']):
        slmsk_var               = nc_file.createVariable('slmsk', wp)
        slmsk_var.units         = 'none'
        slmsk_var.standard_name = 'land_sea_ice_mask'
        slmsk_var[:]            = case_data['slmsk']
    
    if (case_data['vegsrc'] != case_data['missing_value']):
        vegsrc_var               = nc_file.createVariable('vegsrc', wp)
        vegsrc_var.units         = 'none'
        vegsrc_var.standard_name = 'vegetation source (1-2)'
        vegsrc_var[:]            = case_data['vegsrc']
    
    if (case_data['vegtyp'] != case_data['missing_value']):
        vegtyp_var               = nc_file.createVariable('vegtyp', wp)
        vegtyp_var.units         = 'none'
        vegtyp_var.standard_name = 'vegetation type (1-12)'
        vegtyp_var[:]            = case_data['vegtyp']
    
    if (case_data['soiltyp'] != case_data['missing_value']):
        soiltyp_var               = nc_file.createVariable('soiltyp', wp)
        soiltyp_var.units         = 'none'
        soiltyp_var.standard_name = 'soil type (1-12)'
        soiltyp_var[:]            = case_data['soiltyp']
    
    if (case_data['scolor'] != case_data['missing_value']):
        scolor_var               = nc_file.createVariable('scolor', wp)
        scolor_var.units         = 'none'
        scolor_var.standard_name = 'soil color'
        scolor_var[:]            = case_data['scolor']
    
    if (case_data['slopetyp'] != case_data['missing_value']):
        slopetyp_var               = nc_file.createVariable('slopetyp', wp)
        slopetyp_var.units         = 'none'
        slopetyp_var.standard_name = 'slope type (1-9)'
        slopetyp_var[:]            = case_data['slopetyp']
    
    if (case_data['tsfco'] != case_data['missing_value']):
        tsfco_var               = nc_file.createVariable('tsfco', wp)
        tsfco_var.units         = 'none'
        tsfco_var.standard_name = 'slope type (1-9)'
        tsfco_var[:]            = case_data['tsfco']
    
    if (case_data['vegfrac'] != case_data['missing_value']):
        vegfrac_var               = nc_file.createVariable('vegfrac', wp)
        vegfrac_var.units         = 'none'
        vegfrac_var.standard_name = 'slope type (1-9)'
        vegfrac_var[:]            = case_data['vegfrac']
    
    if (case_data['shdmin'] != case_data['missing_value']):
        shdmin_var               = nc_file.createVariable('shdmin', wp)
        shdmin_var.units         = 'none'
        shdmin_var.standard_name = 'slope type (1-9)'
        shdmin_var[:]            = case_data['shdmin']
    
    if (case_data['shdmax'] != case_data['missing_value']):
        shdmax_var               = nc_file.createVariable('shdmax', wp)
        shdmax_var.units         = 'none'
        shdmax_var.standard_name = 'slope type (1-9)'
        shdmax_var[:]            = case_data['shdmax']
    
    if (case_data['canopy'] != case_data['missing_value']):
        canopy_var               = nc_file.createVariable('canopy', wp)
        canopy_var.units         = 'kg m-2'
        canopy_var.standard_name = 'amount of water stored in canopy'
        canopy_var[:]            = case_data['canopy']
    
    if (case_data['hice'] != case_data['missing_value']):
        hice_var               = nc_file.createVariable('hice', wp)
        hice_var.units         = 'm'
        hice_var.standard_name = 'sea ice thickness'
        hice_var[:]            = case_data['hice']
    
    if (case_data['fice'] != case_data['missing_value']):
        fice_var               = nc_file.createVariable('fice', wp)
        fice_var.units         = 'none'
        fice_var.standard_name = 'ice fraction'
        fice_var[:]            = case_data['fice']
    
    if (case_data['tisfc'] != case_data['missing_value']):
        tisfc_var               = nc_file.createVariable('tisfc', wp)
        tisfc_var.units         = 'K'
        tisfc_var.standard_name = 'ice surface temperature'
        tisfc_var[:]            = case_data['tisfc']
    
    if (case_data['snowd'] != case_data['missing_value']):
        snowd_var               = nc_file.createVariable('snowd', wp)
        snowd_var.units         = 'mm'
        snowd_var.standard_name = 'water equivalent snow depth'
        snowd_var[:]            = case_data['snowd']
    
    if (case_data['snoalb'] != case_data['missing_value']):
        snoalb_var               = nc_file.createVariable('snoalb', wp)
        snoalb_var.units         = 'none'
        snoalb_var.standard_name = 'maximum snow albedo'
        snoalb_var[:]            = case_data['snoalb']
    
    if (case_data['tg3'] != case_data['missing_value']):
        tg3_var               = nc_file.createVariable('tg3', wp)
        tg3_var.units         = 'K'
        tg3_var.standard_name = 'deep soil temperature'
        tg3_var[:]            = case_data['tg3']
    
    if (case_data['uustar'] != case_data['missing_value']):
        uustar_var               = nc_file.createVariable('uustar', wp)
        uustar_var.units         = 'm s-1'
        uustar_var.standard_name = 'surface_friction_velocity'
        uustar_var[:]            = case_data['uustar']
    
    if (case_data['alvsf'] != case_data['missing_value']):
        alvsf_var               = nc_file.createVariable('alvsf', wp)
        alvsf_var.units         = 'none'
        alvsf_var.standard_name = '60 degree vis albedo with strong cosz dependency'
        alvsf_var[:]            = case_data['alvsf']
    
    if (case_data['alnsf'] != case_data['missing_value']):
        alnsf_var               = nc_file.createVariable('alnsf', wp)
        alnsf_var.units         = 'none'
        alnsf_var.standard_name = '60 degree nir albedo with strong cosz dependency'
        alnsf_var[:]            = case_data['alnsf']
    
    if (case_data['alvwf'] != case_data['missing_value']):
        alvwf_var               = nc_file.createVariable('alvwf', wp)
        alvwf_var.units         = 'none'
        alvwf_var.standard_name = '60 degree vis albedo with weak cosz dependency'
        alvwf_var[:]            = case_data['alvwf']
    
    if (case_data['alnwf'] != case_data['missing_value']):
        alnwf_var               = nc_file.createVariable('alnwf', wp)
        alnwf_var.units         = 'none'
        alnwf_var.standard_name = '60 degree nir albedo with weak cosz dependency'
        alnwf_var[:]            = case_data['alnwf']
    
    if (case_data['facsf'] != case_data['missing_value']):
        facsf_var               = nc_file.createVariable('facsf', wp)
        facsf_var.units         = 'none'
        facsf_var.standard_name = 'fractional coverage with strong cosz dependency'
        facsf_var[:]            = case_data['facsf']
    
    if (case_data['facwf'] != case_data['missing_value']):
        facwf_var               = nc_file.createVariable('facwf', wp)
        facwf_var.units         = 'none'
        facwf_var.standard_name = 'fractional coverage with weak cosz dependency'
        facwf_var[:]            = case_data['facwf']
    
    if (case_data['weasd'] != case_data['missing_value']):
        weasd_var               = nc_file.createVariable('weasd', wp)
        weasd_var.units         = 'mm'
        weasd_var.standard_name = 'water equivalent accumulated snow depth'
        weasd_var[:]            = case_data['weasd']
    
    if (case_data['f10m'] != case_data['missing_value']):
        f10m_var               = nc_file.createVariable('f10m', wp)
        f10m_var.units         = 'none'
        f10m_var.standard_name = 'ratio of sigma level 1 wind and 10m wind'
        f10m_var[:]            = case_data['f10m']
    
    if (case_data['t2m'] != case_data['missing_value']):
        t2m_var               = nc_file.createVariable('t2m', wp)
        t2m_var.units         = 'K'
        t2m_var.standard_name = '2-meter absolute temperature'
        t2m_var[:]            = case_data['t2m']
    
    if (case_data['q2m'] != case_data['missing_value']):
        q2m_var               = nc_file.createVariable('q2m', wp)
        q2m_var.units         = 'kg kg-1'
        q2m_var.standard_name = '2-meter specific humidity'
        q2m_var[:]            = case_data['q2m']
    
    if (case_data['ffmm'] != case_data['missing_value']):
        ffmm_var               = nc_file.createVariable('ffmm', wp)
        ffmm_var.units         = 'none'
        ffmm_var.standard_name = 'Monin-Obukhov similarity function for momentum'
        ffmm_var[:]            = case_data['ffmm']
    
    if (case_data['ffhh'] != case_data['missing_value']):
        ffhh_var               = nc_file.createVariable('ffhh', wp)
        ffhh_var.units         = 'none'
        ffhh_var.standard_name = 'Monin-Obukhov similarity function for heat'
        ffhh_var[:]            = case_data['ffhh']
    
    if (case_data['tprcp'] != case_data['missing_value']):
        tprcp_var               = nc_file.createVariable('tprcp', wp)
        tprcp_var.units         = 'm'
        tprcp_var.standard_name = 'instantaneous total precipitation amount'
        tprcp_var[:]            = case_data['tprcp']
    
    if (case_data['srflag'] != case_data['missing_value']):
        srflag_var               = nc_file.createVariable('srflag', wp)
        srflag_var.units         = 'none'
        srflag_var.standard_name = 'snow/rain flag for precipitation'
        srflag_var[:]            = case_data['srflag']
    
    if (case_data['sncovr'] != case_data['missing_value']):
        sncovr_var               = nc_file.createVariable('sncovr', wp)
        sncovr_var.units         = 'none'
        sncovr_var.standard_name = 'surface snow area fraction'
        sncovr_var[:]            = case_data['sncovr']
    
    if (case_data['tsfcl'] != case_data['missing_value']):
        tsfcl_var               = nc_file.createVariable('tsfcl', wp)
        tsfcl_var.units         = 'K'
        tsfcl_var.standard_name = 'surface skin temperature over land'
        tsfcl_var[:]            = case_data['tsfcl']
    
    if (case_data['zorl'] != case_data['missing_value']):
        zorl_var               = nc_file.createVariable('zorl', wp)
        zorl_var.units         = 'cm'
        zorl_var.standard_name = 'surface roughness length'
        zorl_var[:]            = case_data['zorl']
    
    if (case_data['zorll'] != case_data['missing_value']):
        zorll_var               = nc_file.createVariable('zorll', wp)
        zorll_var.units         = 'cm'
        zorll_var.standard_name = 'surface roughness length over land'
        zorll_var[:]            = case_data['zorll']
    
    if (case_data['zorli'] != case_data['missing_value']):
        zorli_var               = nc_file.createVariable('zorli', wp)
        zorli_var.units         = 'cm'
        zorli_var.standard_name = 'surface roughness length over ice'
        zorli_var[:]            = case_data['zorli']
    
    if (case_data['zorlw'] != case_data['missing_value']):
        zorlw_var               = nc_file.createVariable('zorlw', wp)
        zorlw_var.units         = 'cm'
        zorlw_var.standard_name = 'surface roughness length over ocean'
        zorlw_var[:]            = case_data['zorlw']
    
    if (case_data['zorlwav'] != case_data['missing_value']):
        zorlwav_var               = nc_file.createVariable('zorlwav', wp)
        zorlwav_var.units         = 'cm'
        zorlwav_var.standard_name = 'surface_roughness_length_from_wave_model'
        zorlwav_var[:]            = case_data['zorlwav']
    
    if (case_data['tvxy'] != case_data['missing_value']):
        tvxy_var               = nc_file.createVariable('tvxy', wp)
        tvxy_var.units         = 'K'
        tvxy_var.standard_name = 'vegetation temperature for NoahMP'
        tvxy_var[:]            = case_data['tvxy']
    
    if (case_data['tgxy'] != case_data['missing_value']):
        tgxy_var               = nc_file.createVariable('tgxy', wp)
        tgxy_var.units         = 'K'
        tgxy_var.standard_name = 'ground temperature for NoahMP'
        tgxy_var[:]            = case_data['tgxy']
    
    if (case_data['tahxy'] != case_data['missing_value']):
        tahxy_var               = nc_file.createVariable('tahxy', wp)
        tahxy_var.units         = 'K'
        tahxy_var.standard_name = 'canopy air temperature for NoahMP'
        tahxy_var[:]            = case_data['tahxy']
    
    if (case_data['canicexy'] != case_data['missing_value']):
        canicexy_var               = nc_file.createVariable('canicexy', wp)
        canicexy_var.units         = 'mm'
        canicexy_var.standard_name = 'canopy intercepted ice mass for NoahMP'
        canicexy_var[:]            = case_data['canicexy']
    
    if (case_data['canliqxy'] != case_data['missing_value']):
        canliqxy_var               = nc_file.createVariable('canliqxy', wp)
        canliqxy_var.units         = 'mm'
        canliqxy_var.standard_name = 'canopy intercepted liquid water for NoahMP'
        canliqxy_var[:]            = case_data['canliqxy']
    
    if (case_data['eahxy'] != case_data['missing_value']):
        eahxy_var               = nc_file.createVariable('eahxy', wp)
        eahxy_var.units         = 'Pa'
        eahxy_var.standard_name = 'canopy air vapor pressure for NoahMP'
        eahxy_var[:]            = case_data['eahxy']
    
    if (case_data['cmxy'] != case_data['missing_value']):
        cmxy_var               = nc_file.createVariable('cmxy', wp)
        cmxy_var.units         = 'none'
        cmxy_var.standard_name = 'surface drag coefficient for momentum for NoahMP'
        cmxy_var[:]            = case_data['cmxy']
    
    if (case_data['chxy'] != case_data['missing_value']):
        chxy_var               = nc_file.createVariable('chxy', wp)
        chxy_var.units         = 'none'
        chxy_var.standard_name = 'surface exchange coeff heat & moisture for NoahMP'
        chxy_var[:]            = case_data['chxy']
    
    if (case_data['fwetxy'] != case_data['missing_value']):
        fwetxy_var               = nc_file.createVariable('fwetxy', wp)
        fwetxy_var.units         = 'none'
        fwetxy_var.standard_name = 'area fraction of canopy that is wetted/snowed for NoahMP'
        fwetxy_var[:]            = case_data['fwetxy']
    
    if (case_data['sneqvoxy'] != case_data['missing_value']):
        sneqvoxy_var               = nc_file.createVariable('sneqvoxy', wp)
        sneqvoxy_var.units         = 'mm'
        sneqvoxy_var.standard_name = 'snow mass at previous time step for NoahMP'
        sneqvoxy_var[:]            = case_data['sneqvoxy']
    
    if (case_data['alboldxy'] != case_data['missing_value']):
        alboldxy_var               = nc_file.createVariable('alboldxy', wp)
        alboldxy_var.units         = 'none'
        alboldxy_var.standard_name = 'snow albedo at previous time step for NoahMP'
        alboldxy_var[:]            = case_data['alboldxy']
    
    if (case_data['qsnowxy'] != case_data['missing_value']):
        qsnowxy_var               = nc_file.createVariable('qsnowxy', wp)
        qsnowxy_var.units         = 'mm s-1'
        qsnowxy_var.standard_name = 'snow precipitation rate at surface for NoahMP'
        qsnowxy_var[:]            = case_data['qsnowxy']
    
    if (case_data['wslakexy'] != case_data['missing_value']):
        wslakexy_var               = nc_file.createVariable('wslakexy', wp)
        wslakexy_var.units         = 'mm'
        wslakexy_var.standard_name = 'lake water storage for NoahMP'
        wslakexy_var[:]            = case_data['wslakexy']
    
    if (case_data['taussxy'] != case_data['missing_value']):
        taussxy_var               = nc_file.createVariable('taussxy', wp)
        taussxy_var.units         = 'none'
        taussxy_var.standard_name = 'non-dimensional snow age for NoahMP'
        taussxy_var[:]            = case_data['taussxy']
    
    
    if (case_data['waxy'] != case_data['missing_value']):
        waxy_var               = nc_file.createVariable('waxy', wp)
        waxy_var.units         = 'mm'
        waxy_var.standard_name = 'water storage in aquifer for NoahMP'
        waxy_var[:]            = case_data['waxy']
    
    if (case_data['wtxy'] != case_data['missing_value']):
        wtxy_var               = nc_file.createVariable('wtxy', wp)
        wtxy_var.units         = 'mm'
        wtxy_var.standard_name = 'ater storage in aquifer and saturated soil for NoahMP'
        wtxy_var[:]            = case_data['wtxy']
    
    if (case_data['zwtxy'] != case_data['missing_value']):
        zwtxy_var               = nc_file.createVariable('zwtxy', wp)
        zwtxy_var.units         = 'm'
        zwtxy_var.standard_name = 'water table depth for NoahMP'
        zwtxy_var[:]            = case_data['zwtxy']
    
    if (case_data['xlaixy'] != case_data['missing_value']):
        xlaixy_var               = nc_file.createVariable('xlaixy', wp)
        xlaixy_var.units         = 'none'
        xlaixy_var.standard_name = 'leaf area index for NoahMP'
        xlaixy_var[:]            = case_data['xlaixy']
    
    if (case_data['xsaixy'] != case_data['missing_value']):
        xsaixy_var               = nc_file.createVariable('xsaixy', wp)
        xsaixy_var.units         = 'none'
        xsaixy_var.standard_name = 'stem area index for NoahMP'
        xsaixy_var[:]            = case_data['xsaixy']
    
    if (case_data['lfmassxy'] != case_data['missing_value']):
        lfmassxy_var               = nc_file.createVariable('lfmassxy', wp)
        lfmassxy_var.units         = 'g m-2'
        lfmassxy_var.standard_name = 'leaf mass for NoahMP'
        lfmassxy_var[:]            = case_data['lfmassxy']
    
    if (case_data['stmassxy'] != case_data['missing_value']):
        stmassxy_var               = nc_file.createVariable('stmassxy', wp)
        stmassxy_var.units         = 'g m-2'
        stmassxy_var.standard_name = 'stem mass for NoahMP'
        stmassxy_var[:]            = case_data['stmassxy']
    
    if (case_data['rtmassxy'] != case_data['missing_value']):
        rtmassxy_var               = nc_file.createVariable('rtmassxy', wp)
        rtmassxy_var.units         = 'g m-2'
        rtmassxy_var.standard_name = 'fine root mass for NoahMP'
        rtmassxy_var[:]            = case_data['rtmassxy']
    
    if (case_data['woodxy'] != case_data['missing_value']):
        woodxy_var               = nc_file.createVariable('woodxy', wp)
        woodxy_var.units         = 'g m-2'
        woodxy_var.standard_name = 'wood mass including woody roots for NoahMP'
        woodxy_var[:]            = case_data['woodxy']
    
    if (case_data['stblcpxy'] != case_data['missing_value']):
        stblcpxy_var               = nc_file.createVariable('stblcpxy', wp)
        stblcpxy_var.units         = 'g m-2'
        stblcpxy_var.standard_name = 'stable carbon in deep soil for NoahMP'
        stblcpxy_var[:]            = case_data['stblcpxy']
    
    if (case_data['fastcpxy'] != case_data['missing_value']):
        fastcpxy_var               = nc_file.createVariable('fastcpxy', wp)
        fastcpxy_var.units         = 'g m-2'
        fastcpxy_var.standard_name = 'short-lived carbon in shallow soil for NoahMP'
        fastcpxy_var[:]            = case_data['fastcpxy']
    
    if (case_data['smcwtdxy'] != case_data['missing_value']):
        smcwtdxy_var               = nc_file.createVariable('smcwtdxy', wp)
        smcwtdxy_var.units         = 'm3 m-3'
        smcwtdxy_var.standard_name = 'oil water content between the bottom of the soil and the water table for NoahMP'
        smcwtdxy_var[:]            = case_data['smcwtdxy']
    
    if (case_data['deeprechxy'] != case_data['missing_value']):
        deeprechxy_var               = nc_file.createVariable('deeprechxy', wp)
        deeprechxy_var.units         = 'm'
        deeprechxy_var.standard_name = 'echarge to or from the water table when deep for NoahMP'
        deeprechxy_var[:]            = case_data['deeprechxy']
    
    if (case_data['rechxy'] != case_data['missing_value']):
        rechxy_var               = nc_file.createVariable('rechxy', wp)
        rechxy_var.units         = 'm'
        rechxy_var.standard_name = 'recharge to or from the water table when shallow for NoahMP'
        rechxy_var[:]            = case_data['rechxy']
    
    if (case_data['snowxy'] != case_data['missing_value']):
        snowxy_var               = nc_file.createVariable('snowxy', wp)
        snowxy_var.units         = 'none'
        snowxy_var.standard_name = 'number of snow layers for NoahMP'
        snowxy_var[:]            = case_data['snowxy']
    
    if (case_data['wetness'] != case_data['missing_value']):
        wetness_var               = nc_file.createVariable('wetness', wp)
        wetness_var.units         = 'none'
        wetness_var.standard_name = 'normalized soil wetness for RUC LSM'
        wetness_var[:]            = case_data['wetness']
        
    if (case_data['clw_surf_land'] != case_data['missing_value']):
        clw_surf_land_var               = nc_file.createVariable('clw_surf_land', wp)
        clw_surf_land_var.units         = 'kg kg-1'
        clw_surf_land_var.standard_name = 'cloud condensed water mixing ratio at surface over land for RUC LSM'
        clw_surf_land_var[:]            = case_data['clw_surf_land']
        
    if (case_data['clw_surf_ice'] != case_data['missing_value']):
        clw_surf_ice_var               = nc_file.createVariable('clw_surf_ice', wp)
        clw_surf_ice_var.units         = 'kg kg-1'
        clw_surf_ice_var.standard_name = 'cloud condensed water mixing ratio at surface over ice for RUC LSM'
        clw_surf_ice_var[:]            = case_data['clw_surf_ice']
    
    if (case_data['qwv_surf_land'] != case_data['missing_value']):
        qwv_surf_land_var               = nc_file.createVariable('qwv_surf_land', wp)
        qwv_surf_land_var.units         = 'kg kg-1'
        qwv_surf_land_var.standard_name = 'water vapor mixing ratio at surface over land for RUC LSM'
        qwv_surf_land_var[:]            = case_data['qwv_surf_land']
    
    if (case_data['qwv_surf_ice'] != case_data['missing_value']):
        qwv_surf_ice_var               = nc_file.createVariable('qwv_surf_ice', wp)
        qwv_surf_ice_var.units         = 'kg kg-1'
        qwv_surf_ice_var.standard_name = 'water vapor mixing ratio at surface over ice for RUC LSM'
        qwv_surf_ice_var[:]            = case_data['qwv_surf_ice']
    
    if (case_data['tsnow_land'] != case_data['missing_value']):
        tsnow_land_var               = nc_file.createVariable('tsnow_land', wp)
        tsnow_land_var.units         = 'K'
        tsnow_land_var.standard_name = 'snow temperature at the bottom of the first snow layer over land for RUC LSM'
        tsnow_land_var[:]            = case_data['tsnow_land']
    
    if (case_data['tsnow_ice'] != case_data['missing_value']):
        tsnow_ice_var               = nc_file.createVariable('tsnow_ice', wp)
        tsnow_ice_var.units         = 'K'
        tsnow_ice_var.standard_name = 'snow temperature at the bottom of the first snow layer over land for RUC LSM'
        tsnow_ice_var[:]            = case_data['tsnow_ice']
    
    if (case_data['snowfallac_land'] != case_data['missing_value']):
        snowfallac_land_var               = nc_file.createVariable('snowfallac_land', wp)
        snowfallac_land_var.units         = 'kg m-2'
        snowfallac_land_var.standard_name = 'run-total snow accumulation on the ground over land for RUC LSM'
        snowfallac_land_var[:]            = case_data['snowfallac_land']
    
    if (case_data['snowfallac_ice'] != case_data['missing_value']):
        snowfallac_ice_var               = nc_file.createVariable('snowfallac_ice', wp)
        snowfallac_ice_var.units         = 'kg m-2'
        snowfallac_ice_var.standard_name = 'run-total snow accumulation on the ground over land for RUC LSM'
        snowfallac_ice_var[:]            = case_data['snowfallac_ice']
        
    if (case_data['sncovr_ice'] != case_data['missing_value']):
        sncovr_ice_var               = nc_file.createVariable('sncovr_ice', wp)
        sncovr_ice_var.units         = 'none'
        sncovr_ice_var.standard_name = 'surface snow area fraction over ice'
        sncovr_ice_var[:]            = case_data['sncovr_ice']
    
    if (case_data['sfalb_lnd'] != case_data['missing_value']):
        sfalb_lnd_var               = nc_file.createVariable('sfalb_lnd', wp)
        sfalb_lnd_var.units         = 'none'
        sfalb_lnd_var.standard_name = 'surface albedo over land for RUC LSM'
        sfalb_lnd_var[:]            = case_data['sfalb_lnd']
    
    if (case_data['sfalb_lnd_bck'] != case_data['missing_value']):
        sfalb_lnd_bck_var               = nc_file.createVariable('sfalb_lnd_bckß', wp)
        sfalb_lnd_bck_var.units         = 'none'
        sfalb_lnd_bck_var.standard_name = 'surface snow-free albedo over land for RUC LSM'
        sfalb_lnd_bck_var[:]            = case_data['sfalb_lnd_bck']
    
    if (case_data['emis_ice'] != case_data['missing_value']):
        emis_ice_var               = nc_file.createVariable('emis_ice', wp)
        emis_ice_var.units         = 'none'
        emis_ice_var.standard_name = 'surface emissivity over ice for RUC LSM'
        emis_ice_var[:]            = case_data['emis_ice']
    
    if (case_data['lai'] != case_data['missing_value']):
        lai_var               = nc_file.createVariable('lai', wp)
        lai_var.units         = 'none'
        lai_var.standard_name = 'leaf area index for RUC LSM'
        lai_var[:]            = case_data['lai']
    
    area_var               = nc_file.createVariable('area', wp)
    area_var.units         = 'm2'
    area_var.standard_name = 'grid cell area'
    if ('column_area' in nml_keys and case_nml['case_config']['column_area']):
        area_var[:]            = case_nml['case_config']['column_area']
        message = 'Since column_area was supplied in the case namelist, it will be used instead of the data from the case data file (if it exists).'
        logging.info(message)
    elif (case_data['area'] != case_data['missing_value']):
        area_var[:]            = case_data['area']
    
    if (case_data['stddev'] != case_data['missing_value']):
        stddev_var               = nc_file.createVariable('stddev', wp)
        stddev_var.units         = 'm'
        stddev_var.standard_name = 'standard deviation of subgrid orography'
        stddev_var[:]            = case_data['stddev']
        
    if (case_data['convexity'] != case_data['missing_value']):
        convexity_var               = nc_file.createVariable('convexity', wp)
        convexity_var.units         = 'none'
        convexity_var.standard_name = 'convexity of subgrid orography'
        convexity_var[:]            = case_data['convexity']
    
    if (case_data['oa1'] != case_data['missing_value']):
        oa1_var               = nc_file.createVariable('oa1', wp)
        oa1_var.units         = 'none'
        oa1_var.standard_name = 'assymetry of subgrid orography 1'
        oa1_var[:]            = case_data['oa1']
    
    if (case_data['oa2'] != case_data['missing_value']):
        oa2_var               = nc_file.createVariable('oa2', wp)
        oa2_var.units         = 'none'
        oa2_var.standard_name = 'assymetry of subgrid orography 2'
        oa2_var[:]            = case_data['oa2']
    
    if (case_data['oa3'] != case_data['missing_value']):
        oa3_var               = nc_file.createVariable('oa3', wp)
        oa3_var.units         = 'none'
        oa3_var.standard_name = 'assymetry of subgrid orography 3'
        oa3_var[:]            = case_data['oa3']
    
    if (case_data['oa4'] != case_data['missing_value']):
        oa4_var               = nc_file.createVariable('oa4', wp)
        oa4_var.units         = 'none'
        oa4_var.standard_name = 'assymetry of subgrid orography 4'
        oa4_var[:]            = case_data['oa4']
        
    if (case_data['ol1'] != case_data['missing_value']):
        ol1_var               = nc_file.createVariable('ol1', wp)
        ol1_var.units         = 'none'
        ol1_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 1'
        ol1_var[:]            = case_data['ol1']
        
    if (case_data['ol2'] != case_data['missing_value']):
        ol2_var               = nc_file.createVariable('ol2', wp)
        ol2_var.units         = 'none'
        ol2_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 2'
        ol2_var[:]            = case_data['ol2']
        
    if (case_data['ol3'] != case_data['missing_value']):
        ol3_var               = nc_file.createVariable('ol3', wp)
        ol3_var.units         = 'none'
        ol3_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 3'
        ol3_var[:]            = case_data['ol3']
        
    if (case_data['ol4'] != case_data['missing_value']):
        ol4_var               = nc_file.createVariable('ol4', wp)
        ol4_var.units         = 'none'
        ol4_var.standard_name = 'fraction of grid box with subgrid orography higher than critical height 4'
        ol4_var[:]            = case_data['ol4']
        
    if (case_data['theta_oro'] != case_data['missing_value']):
        theta_oro_var               = nc_file.createVariable('theta_oro', wp)
        theta_oro_var.units         = 'deg'
        theta_oro_var.standard_name = 'angle with respect to east of maximum subgrid orographic variations'
        theta_oro_var[:]            = case_data['theta_oro']
        
    if (case_data['gamma'] != case_data['missing_value']):
        gamma_var               = nc_file.createVariable('gamma', wp)
        gamma_var.units         = 'none'
        gamma_var.standard_name = 'anisotropy of subgrid orography'
        gamma_var[:]            = case_data['gamma']
        
    if (case_data['sigma'] != case_data['missing_value']):
        sigma_var               = nc_file.createVariable('sigma', wp)
        sigma_var.units         = 'none'
        sigma_var.standard_name = 'slope of subgrid orography'
        sigma_var[:]            = case_data['sigma']
    
    if (case_data['elvmax'] != case_data['missing_value']):
        elvmax_var               = nc_file.createVariable('elvmax', wp)
        elvmax_var.units         = 'm'
        elvmax_var.standard_name = 'maximum of subgrid orography'
        elvmax_var[:]            = case_data['elvmax']
        
    if (case_data['oro'] != case_data['missing_value']):
        oro_var               = nc_file.createVariable('oro', wp)
        oro_var.units         = 'm'
        oro_var.standard_name = 'orography'
        oro_var[:]            = case_data['oro']
        
    if (case_data['oro_uf'] != case_data['missing_value']):
        oro_uf_var               = nc_file.createVariable('oro_uf', wp)
        oro_uf_var.units         = 'm'
        oro_uf_var.standard_name = 'unfiltered orography'
        oro_uf_var[:]            = case_data['oro_uf']
        
    if (case_data['landfrac'] != case_data['missing_value']):
        landfrac_var               = nc_file.createVariable('landfrac', wp)
        landfrac_var.units         = 'none'
        landfrac_var.standard_name = 'fraction of horizontal grid area occupied by land'
        landfrac_var[:]            = case_data['landfrac']
        
    if (case_data['lakefrac'] != case_data['missing_value']):
        lakefrac_var               = nc_file.createVariable('lakefrac', wp)
        lakefrac_var.units         = 'none'
        lakefrac_var.standard_name = 'fraction of horizontal grid area occupied by lake'
        lakefrac_var[:]            = case_data['lakefrac']
        
    if (case_data['lakedepth'] != case_data['missing_value']):
        lakedepth_var               = nc_file.createVariable('lakedepth', wp)
        lakedepth_var.units         = 'none'
        lakedepth_var.standard_name = 'lake depth'
        lakedepth_var[:]            = case_data['lakedepth']
    
    if (case_data['tref'] != case_data['missing_value']):
        tref_var               = nc_file.createVariable('tref', wp)
        tref_var.units         = 'K'
        tref_var.standard_name = 'sea surface reference temperature for NSST'
        tref_var[:]            = case_data['tref']
    
    if (case_data['z_c'] != case_data['missing_value']):
        z_c_var               = nc_file.createVariable('z_c', wp)
        z_c_var.units         = 'm'
        z_c_var.standard_name = 'sub-layer cooling thickness for NSST'
        z_c_var[:]            = case_data['z_c']
        
    if (case_data['c_0'] != case_data['missing_value']):
        c_0_var               = nc_file.createVariable('c_0', wp)
        c_0_var.units         = 'none'
        c_0_var.standard_name = 'coefficient 1 to calculate d(Tz)/d(Ts) for NSST'
        c_0_var[:]            = case_data['c_0']
        
    if (case_data['c_d'] != case_data['missing_value']):
        c_d_var               = nc_file.createVariable('c_d', wp)
        c_d_var.units         = 'none'
        c_d_var.standard_name = 'coefficient 2 to calculate d(Tz)/d(Ts) for NSST'
        c_d_var[:]            = case_data['c_d']
        
    if (case_data['w_0'] != case_data['missing_value']):
        w_0_var               = nc_file.createVariable('w_0', wp)
        w_0_var.units         = 'none'
        w_0_var.standard_name = 'coefficient 3 to calculate d(Tz)/d(Ts) for NSST'
        w_0_var[:]            = case_data['w_0']
        
    if (case_data['w_d'] != case_data['missing_value']):
        w_d_var               = nc_file.createVariable('w_d', wp)
        w_d_var.units         = 'none'
        w_d_var.standard_name = 'coefficient 4 to calculate d(Tz)/d(Ts) for NSST'
        w_d_var[:]            = case_data['w_d']
        
    if (case_data['xt'] != case_data['missing_value']):
        xt_var               = nc_file.createVariable('xt', wp)
        xt_var.units         = 'K m'
        xt_var.standard_name = 'heat content in diurnal thermocline layer for NSST'
        xt_var[:]            = case_data['xt']
        
    if (case_data['xs'] != case_data['missing_value']):
        xs_var               = nc_file.createVariable('xs', wp)
        xs_var.units         = 'ppt m'
        xs_var.standard_name = 'salinity content in diurnal thermocline layer for NSST'
        xs_var[:]            = case_data['xs']
        
    if (case_data['xu'] != case_data['missing_value']):
        xu_var               = nc_file.createVariable('xu', wp)
        xu_var.units         = 'm2 s-1'
        xu_var.standard_name = 'u-current in diurnal thermocline layer for NSST'
        xu_var[:]            = case_data['xu']
        
    if (case_data['xv'] != case_data['missing_value']):
        xv_var               = nc_file.createVariable('xv', wp)
        xv_var.units         = 'm2 s-1'
        xv_var.standard_name = 'v-current in diurnal thermocline layer for NSST'
        xv_var[:]            = case_data['xv']
        
    if (case_data['xz'] != case_data['missing_value']):
        xz_var               = nc_file.createVariable('xz', wp)
        xz_var.units         = 'm'
        xz_var.standard_name = 'thickness of diurnal thermocline layer for NSST'
        xz_var[:]            = case_data['xz']
        
    if (case_data['zm'] != case_data['missing_value']):
        zm_var               = nc_file.createVariable('zm', wp)
        zm_var.units         = 'm'
        zm_var.standard_name = 'thickness of ocean mixed layer for NSST'
        zm_var[:]            = case_data['zm']
        
    if (case_data['xtts'] != case_data['missing_value']):
        xtts_var               = nc_file.createVariable('xtts', wp)
        xtts_var.units         = 'm'
        xtts_var.standard_name = 'sensitivity of diurnal thermocline layer heat content to surface temperature [d(xt)/d(ts)] for NSST'
        xtts_var[:]            = case_data['xtts']
    
    if (case_data['xzts'] != case_data['missing_value']):
        xzts_var               = nc_file.createVariable('xzts', wp)
        xzts_var.units         = 'm K-1'
        xzts_var.standard_name = 'sensitivity of diurnal thermocline layer thickness to surface temperature [d(xz)/d(ts)] for NSST'
        xzts_var[:]            = case_data['xzts']
        
    if (case_data['d_conv'] != case_data['missing_value']):
        d_conv_var               = nc_file.createVariable('d_conv', wp)
        d_conv_var.units         = 'm'
        d_conv_var.standard_name = 'thickness of free convection layer for NSST'
        d_conv_var[:]            = case_data['d_conv']
        
    if (case_data['ifd'] != case_data['missing_value']):
        ifd_var               = nc_file.createVariable('ifd', wp)
        ifd_var.units         = 'none'
        ifd_var.standard_name = 'index to start DTM run for NSST'
        ifd_var[:]            = case_data['ifd']
        
    if (case_data['dt_cool'] != case_data['missing_value']):
        dt_cool_var               = nc_file.createVariable('dt_cool', wp)
        dt_cool_var.units         = 'K'
        dt_cool_var.standard_name = 'sub-layer cooling amount for NSST'
        dt_cool_var[:]            = case_data['dt_cool']
        
    if (case_data['qrain'] != case_data['missing_value']):
        qrain_var               = nc_file.createVariable('qrain', wp)
        qrain_var.units         = 'W m-2'
        qrain_var.standard_name = 'sensible heat due to rainfall for NSST'
        qrain_var[:]            = case_data['qrain']
              
    if (case_data['theta_il'][0] != case_data['missing_value']):
        thetal_var                   = nc_file.createVariable('thetal', wp, ('t0','lev'))
        thetal_var.units             = 'K'
        thetal_var.standard_name     = 'air_liquid_potential_temperature'
        thetal_var[:]                = case_data['theta_il'][:]
    
    if (case_data['t'][0] != case_data['missing_value']):
        t_var                   = nc_file.createVariable('t', wp, ('t0','lev'))
        t_var.units             = 'K'
        t_var.standard_name     = 'absolute temperature'
        t_var[:]                = case_data['t'][:]
    
    #
    qt_var                       = nc_file.createVariable('qt', wp, ('t0','lev'))
    qt_var.units                 = 'kg kg-1'
    qt_var.standard_name         = 'mass_fraction_of_water_in_air'
    qt_var[:]                    = case_data['qt'][:]
    
    #
    u_var                        = nc_file.createVariable('ua', wp, ('t0','lev'))
    u_var.units                  = 'm s-1'
    u_var.standard_name          = 'eastward_wind'
    u_var[:]                     = case_data['u'][:]
    
    #
    v_var                        = nc_file.createVariable('va', wp, ('t0','lev'))
    v_var.units                  = 'm s-1'
    v_var.standard_name          = 'northward_wind'
    v_var[:]                     = case_data['v'][:]
    
    #
    p_var                        = nc_file.createVariable('pa', wp, ('t0','lev'))
    p_var.units                  = 'Pa'
    p_var.standard_name          = 'air_pressure'
    p_var[:]                     = case_data['levels'][:]
    
    #
    z_var                        = nc_file.createVariable('zh', wp, ('t0','lev'))
    z_var.units                  = 'm'
    z_var.standard_name          = 'height'
    z_var[:]                     = case_data['height'][:]
    
    #
    ps_var                       = nc_file.createVariable('ps', wp, ('t0'))
    ps_var.units                 = 'Pa'
    ps_var.standard_name         = 'surface_air_pressure'
    ps_var[:]                    = case_data['p_surf'][0]
    
    #
    ql_var                       = nc_file.createVariable('ql', wp, ('t0','lev'))
    ql_var.units                 = 'kg kg-1'
    ql_var.standard_name         = 'mass_fraction_of_cloud_liquid_water_in_air'
    ql_var[:]                    = case_data['ql'][:]
    
    #
    qi_var                       = nc_file.createVariable('qi', wp, ('t0','lev'))
    qi_var.units                 = 'kg kg-1'
    qi_var.standard_name         = 'mass_fraction_of_cloud_ice_water_in_air'
    qi_var[:]                    = case_data['qi'][:]
    
    #
    tke_var                      = nc_file.createVariable('tke', wp, ('t0','lev'))
    tke_var.units                = 'm2 s-2'
    tke_var.standard_name        = 'specific_turbulent_kinetic_energy'
    tke_var[:]                   = case_data['tke'][:]
    
    #
    ozone_var                    = nc_file.createVariable('o3', wp, ('t0','lev'))
    ozone_var.units              = 'kg kg-1'
    ozone_var.standard_name      = 'mole_fraction_of_ozone_in_air'
    ozone_var[:]                 = case_data['ozone'][:]
    
    if (case_data['stc'][0] != case_data['missing_value']):
        stc_var               = nc_file.createVariable('stc', wp, ('t0','nsoil'))
        stc_var.units         = 'K'
        stc_var.standard_name = 'initial profile of soil temperature'
        stc_var[:]            = case_data['stc'][:]
    
    if (case_data['smc'][0] != case_data['missing_value']):
        smc_var               = nc_file.createVariable('smc', wp, ('t0','nsoil'))
        smc_var.units         = 'kg'
        smc_var.standard_name = 'initial profile of soil moisture'
        smc_var[:]            = case_data['smc'][:]
        
    if (case_data['slc'][0] != case_data['missing_value']):
        slc_var               = nc_file.createVariable('slc', wp, ('t0','nsoil'))
        slc_var.units         = 'kg'
        slc_var.standard_name = 'initial profile of soil liquid moisture'
        slc_var[:]            = case_data['slc'][:]
    
    if (case_data['snicexy'][0] != case_data['missing_value']):
        snicexy_var               = nc_file.createVariable('snicexy', wp, ('t0','nsnow'))
        snicexy_var.units         = 'mm'
        snicexy_var.standard_name = 'initial profile of snow layer ice'
        snicexy_var[:]            = case_data['snicexy'][:]
    
    if (case_data['snliqxy'][0] != case_data['missing_value']):
        snliqxy_var               = nc_file.createVariable('snliqxy', wp, ('t0','nsnow'))
        snliqxy_var.units         = 'mm'
        snliqxy_var.standard_name = 'initial profile of snow layer liquid'
        snliqxy_var[:]            = case_data['snliqxy'][:]

    if (case_data['tsnoxy'][0] != case_data['missing_value']):
        tsnoxy_var               = nc_file.createVariable('tsnoxy', wp, ('t0','nsnow'))
        tsnoxy_var.units         = 'K'
        tsnoxy_var.standard_name = 'initial profile of snow layer temperature'
        tsnoxy_var[:]            = case_data['tsnoxy'][:]
    
    if (case_data['smoiseq'][0] != case_data['missing_value']):
        smoiseq_var               = nc_file.createVariable('smoiseq', wp, ('t0','nsoil'))
        smoiseq_var.units         = 'm3 m-3'
        smoiseq_var.standard_name = 'initial profile of equilibrium soil water content'
        smoiseq_var[:]            = case_data['smoiseq'][:]
    
    if (case_data['zsnsoxy'][0] != case_data['missing_value']):
        zsnsoxy_var               = nc_file.createVariable('zsnxosy', wp, ('t0','nsoil_plus_nsnow'))
        zsnsoxy_var.units         = 'm'
        zsnsoxy_var.standard_name = 'layer bottom depth from snow surface'
        zsnsoxy_var[:]            = case_data['zsnsoxy'][:]
    
    if (case_data['tiice'][0] != case_data['missing_value']):
        tiice_var               = nc_file.createVariable('tiice', wp, ('t0','nice'))
        tiice_var.units         = 'K'
        tiice_var.standard_name = 'sea ice internal temperature'
        tiice_var[:]            = case_data['tiice'][:]
    
    if (case_data['tslb'][0] != case_data['missing_value']):
        tslb_var               = nc_file.createVariable('tslb', wp, ('t0','nsoil'))
        tslb_var.units         = 'K'
        tslb_var.standard_name = 'soil temperature for RUC LSM'
        tslb_var[:]            = case_data['tslb'][:]
    
    if (case_data['smois'][0] != case_data['missing_value']):
        smois_var               = nc_file.createVariable('smois', wp, ('t0','nsoil'))
        smois_var.units         = 'none'
        smois_var.standard_name = 'volume fraction of soil moisture for RUC LSM'
        smois_var[:]            = case_data['smois'][:]
    
    if (case_data['sh2o'][0] != case_data['missing_value']):
        sh2o_var               = nc_file.createVariable('sh2o', wp, ('t0','nsoil'))
        sh2o_var.units         = 'none'
        sh2o_var.standard_name = 'volume fraction of unfrozen soil moisture for RUC LSM'
        sh2o_var[:]            = case_data['sh2o'][:]
    
    if (case_data['smfr'][0] != case_data['missing_value']):
        smfr_var               = nc_file.createVariable('smfr', wp, ('t0','nsoil'))
        smfr_var.units         = 'none'
        smfr_var.standard_name = 'volume fraction of frozen soil moisture for RUC LSM'
        smfr_var[:]            = case_data['smfr'][:]
    
    if (case_data['flfr'][0] != case_data['missing_value']):
        flfr_var               = nc_file.createVariable('flfr', wp, ('t0','nsoil'))
        flfr_var.units         = 'none'
        flfr_var.standard_name = 'flag for frozen soil physics for RUC LSM'
        flfr_var[:]            = case_data['flfr'][:]
        
    ps_forc_var                  = nc_file.createVariable('ps_forc', wp, ('time'))
    ps_forc_var.units            = 'Pa'
    ps_forc_var.standard_name    = 'forcing_surface_air_pressure'
    ps_forc_var[:]               = case_data['p_surf'][:]
    
    pa_forc_var                  = nc_file.createVariable('pa_forc', wp, ('time','lev'))
    pa_forc_var.units            = 'Pa'
    pa_forc_var.standard_name    = 'air_pressure_forcing'
    pa_forc_var[:]               = case_data['levels'][:]
    
    zh_forc_var                  = nc_file.createVariable('zh_forc', wp, ('time','lev'))
    zh_forc_var.units            = 'm'
    zh_forc_var.standard_name    = 'height_forcing'
    zh_forc_var[:]               = case_data['height'][:]
    
    if (nc_file.adv_ta == forcing_on):
        message = 'adv_ta is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnta_adv_var                    = nc_file.createVariable('tnta_adv', wp, ('time','lev'))
        # tnta_adv_var.units              = 'K s-1'
        # tnta_adv_var.standard_name      = 'tendency_of_air_temperature_due_to_advection'
        # tnta_adv_var[:]                 = np.swapaxes(case_data['tnta_adv'][:],0,1)
        
    if (nc_file.adv_qv == forcing_on):
        message = 'adv_qv is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnqv_adv_var                    = nc_file.createVariable('tnqv_adv', wp, ('time','lev'))
        # tnqv_adv_var.units              = 'kg kg-1 s-1'
        # tnqv_adv_var.standard_name      = 'tendency_of_specific_humidity_due_to_advection'
        # tnqv_adv_var[:]                 = np.swapaxes(case_data['tnqv_adv'][:],0,1)
    
    if (nc_file.adv_ua == forcing_on):
        message = 'adv_ua is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnua_adv_var                    = nc_file.createVariable('tnua_adv', wp, ('time','lev'))
        # tnua_adv_var.units              = 'm s-2'
        # tnua_adv_var.standard_name      = 'tendency_of_eastward_wind_due_to_advection'
        # tnua_adv_var[:]                 = np.swapaxes(case_data['tnua_adv'][:],0,1)
    
    if (nc_file.adv_va == forcing_on):
        message = 'adv_va is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnva_adv_var                    = nc_file.createVariable('tnva_adv', wp, ('time','lev'))
        # tnva_adv_var.units              = 'm s-2'
        # tnva_adv_var.standard_name      = 'tendency_of_northward_wind_due_to_advection'
        # tnva_adv_var[:]                 = np.swapaxes(case_data['tnva_adv'][:],0,1)
    
    if (nc_file.adv_theta == forcing_on):
        message = 'adv_theta is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tntheta_adv_var                    = nc_file.createVariable('tntheta_adv', wp, ('time','lev'))
        # tntheta_adv_var.units              = 'K s-1'
        # tntheta_adv_var.standard_name      = 'tendency_of_air_potential_temperature_due_to_advection'
        # tntheta_adv_var[:]                 = np.swapaxes(case_data['tntheta_adv'][:],0,1)
    
    if (nc_file.adv_thetal == forcing_on):
        tnthetal_adv_var                    = nc_file.createVariable('tnthetal_adv', wp, ('time','lev'))
        tnthetal_adv_var.units              = 'K s-1'
        tnthetal_adv_var.standard_name      = 'tendency_of_air_liquid_potential_temperature_due_to_advection'
        if (nc_file.forc_wap == forcing_on or nc_file.forc_wa == forcing_on):
            tnthetal_adv_var[:]                 = np.swapaxes(case_data['h_advec_thil'][:],0,1)
        else:
            tnthetal_adv_var[:]                 = np.swapaxes(case_data['h_advec_thil'][:] + case_data['v_advec_thil'][:],0,1)
    
    if (nc_file.adv_qt == forcing_on):
        tnqt_adv_var                    = nc_file.createVariable('tnqt_adv', wp, ('time','lev'))
        tnqt_adv_var.units              = 'kg kg-1 s-1'
        tnqt_adv_var.standard_name      = 'tendency_of_mass_fraction_of_water_in_air_due_to_advection'
        if (nc_file.forc_wap == forcing_on or nc_file.forc_wa == forcing_on):
            tnqt_adv_var[:]                 = np.swapaxes(case_data['h_advec_qt'][:],0,1)
        else:
            tnqt_adv_var[:]                 = np.swapaxes(case_data['h_advec_qt'][:] + case_data['v_advec_qt'][:],0,1)
    
    if (nc_file.adv_rv == forcing_on):
        message = 'adv_rv is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnrv_adv_var                    = nc_file.createVariable('tnrv_adv', wp, ('time','lev'))
        # tnrv_adv_var.units              = 'kg kg-1 s-1'
        # tnrv_adv_var.standard_name      = 'tendency_of_humidity_mixing_ratio_due_to_advection'
        # tnrv_adv_var[:]                 = np.swapaxes(case_data['tnrv_adv'][:],0,1)
    
    if (nc_file.adv_rt == forcing_on):
        message = 'adv_rt is turned on, but is not implemented in the proprietery CCPP SCM case format and cannot be used.'
        logging.critical(message)
        raise Exception(message)
        # tnrt_adv_var                    = nc_file.createVariable('tnrt_adv', wp, ('time','lev'))
        # tnrt_adv_var.units              = 'kg kg-1 s-1'
        # tnrt_adv_var.standard_name      = 'tendency_of_water_mixing_ratio_due_to_advection'
        # tnrt_adv_var[:]                 = np.swapaxes(case_data['tnrt_adv'][:],0,1)
    
    if (nc_file.forc_wap == forcing_on):
        wap_var                    = nc_file.createVariable('wap', wp, ('time','lev'))
        wap_var.units              = 'Pa s-1'
        wap_var.standard_name      = 'lagrangian_tendency_of_air_pressure'
        wap_var[:]                 = np.swapaxes(case_data['omega'][:],0,1)
    elif (nc_file.forc_wa == forcing_on):
        wa_var                    = nc_file.createVariable('wa', wp, ('time','lev'))
        wa_var.units              = 'm s-1'
        wa_var.standard_name      = 'upward_air_velocity'
        wa_var[:]                 = np.swapaxes(case_data['w_ls'][:],0,1)
        
    if (nc_file.forc_geo == forcing_on):
        ug_var                    = nc_file.createVariable('ug', wp, ('time','lev'))
        ug_var.units              = 'm s-1'
        ug_var.standard_name      = 'geostrophic_eastward_wind'
        ug_var[:]                 = np.swapaxes(case_data['u_g'][:],0,1)
        
        vg_var                    = nc_file.createVariable('vg', wp, ('time','lev'))
        vg_var.units              = 'm s-1'
        vg_var.standard_name      = 'geostrophic_northward_wind'
        vg_var[:]                 = np.swapaxes(case_data['v_g'][:],0,1)
        
    if (nc_file.nudging_ua != forcing_off):
        ua_nud_var                    = nc_file.createVariable('ua_nud', wp, ('time','lev'))
        ua_nud_var.units              = 'm s-1'
        ua_nud_var.standard_name      = 'nudging_eastward_wind'
        ua_nud_var[:]                 = np.swapaxes(case_data['u_nudge'][:],0,1)
    
    if (nc_file.nudging_va != forcing_off):
        va_nud_var                    = nc_file.createVariable('va_nud', wp, ('time','lev'))
        va_nud_var.units              = 'm s-1'
        va_nud_var.standard_name      = 'nudging_northward_wind'
        va_nud_var[:]                 = np.swapaxes(case_data['v_nudge'][:],0,1)
    
    if (nc_file.nudging_ta != forcing_off):
        ta_nud_var                    = nc_file.createVariable('ta_nud', wp, ('time','lev'))
        ta_nud_var.units              = 'K'
        ta_nud_var.standard_name      = 'nudging_air_temperature'
        ta_nud_var[:]                 = np.swapaxes(case_data['T_nudge'][:],0,1)
    
    if (nc_file.nudging_qt != forcing_off):
        qt_nud_var                    = nc_file.createVariable('qt_nud', wp, ('time','lev'))
        qt_nud_var.units              = 'kg kg-1'
        qt_nud_var.standard_name      = 'nudging_mass_fraction_of_water_in_air'
        qt_nud_var[:]                 = np.swapaxes(case_data['qt_nudge'][:],0,1)
    
    if (case_nml['case_config']['sfc_flux_spec']):
        nc_file.surface_forcing_temp = 'kinematic'
        nc_file.surface_forcing_moisture = 'kinematic'
        nc_file.surface_forcing_wind = 'z0'
        
        wpthetap_s_var                  = nc_file.createVariable('wpthetap_s', wp, ('time'))
        wpthetap_s_var.units            = 'K m s-1'
        wpthetap_s_var.standard_name    = 'surface_upward_potential_temperature_flux'
        wpthetap_s_var[:]               = case_data['sh_flux_sfc'][:]
        
        wpqtp_s_var                  = nc_file.createVariable('wpqvp_s', wp, ('time'))
        wpqtp_s_var.units            = 'kg kg-1 m s-1'
        wpqtp_s_var.standard_name    = 'surface_upward_specific_humidity_flux'
        wpqtp_s_var[:]               = case_data['lh_flux_sfc'][:]
        
        z0_var                  = nc_file.createVariable('z0', wp, ('time'))
        z0_var.units            = 'm'
        z0_var.standard_name    = 'surface_roughness_length_for_momentum_in_air'
        z0_var[:]               = case_nml['case_config']['sfc_roughness_length_cm']*1.0E-2
        
        
        ts_var                  = nc_file.createVariable('ts_forc', wp, ('time'))
        ts_var.units            = 'K'
        ts_var.standard_name    = 'forcing_surface_temperature'
        if np.any(case_data['T_surf'][:]):
            ts_var[:]               = case_data['T_surf'][:]
        else:
            ts_var[:]               = case_data['missing_value']
    else:
        nc_file.surface_forcing_temp = 'ts'
        nc_file.surface_forcing_wind = 'z0'
            
        z0_var                  = nc_file.createVariable('z0', wp, ('time'))
        z0_var.units            = 'm'
        z0_var.standard_name    = 'surface_roughness_length_for_momentum_in_air'
        z0_var[:]               = case_nml['case_config']['sfc_roughness_length_cm']*1.0E-2
        
        ts_var                  = nc_file.createVariable('ts_forc', wp, ('time'))
        ts_var.units            = 'K'
        ts_var.standard_name    = 'forcing_surface_temperature'
        ts_var[:]               = case_data['T_surf'][:]
        
    
    nc_file.close()
    
    return(fileOUT)

def write_SCM_nml_file(case_nml):
    filename = os.path.join(CASE_NML_DIR, case_nml['case_config']['case_name'] + '_dephy.nml')
    
    #Go through existing case namelist and only add necessary items to new DEPHY-based namelist
    
    #add _dephy to case (temporary - to differentiate from old format case)
    int_dict = {'case_name':case_nml['case_config']['case_name']+'_dephy',
                'input_type':1}
    
    nml_keys = case_nml['case_config'].todict().keys()
    if ('npz_type' in nml_keys):
        int_dict['npz_type'] = case_nml['case_config']['npz_type']
        if int_dict['npz_type'] == 'input' and 'vert_coord_file' in nml_keys:
            int_dict['vert_coord_file'] = case_nml['case_config']['vert_coord_file']
    
    if ('dt' in nml_keys):
        int_dict['dt'] = case_nml['case_config']['dt']
    
    #runtime is in netCDF file
    
    if ('output_dir' in nml_keys):
        int_dict['output_dir'] = case_nml['case_config']['output_dir']
        
    if ('model_ics' in nml_keys):
        int_dict['model_ics'] = case_nml['case_config']['model_ics']
    
    if ('lsm_ics' in nml_keys):
        int_dict['lsm_ics'] = case_nml['case_config']['lsm_ics']
    
    if ('do_spinup' in nml_keys):
        int_dict['do_spinup'] = case_nml['case_config']['do_spinup']
    
    if ('spinup_timesteps' in nml_keys):
        int_dict['spinup_timesteps'] = case_nml['case_config']['spinup_timesteps']
    
    if ('C_RES' in nml_keys):
        int_dict['C_RES'] = case_nml['case_config']['C_RES']
    
    #relax_time is in netCDF file
    
    #sfc_type is in netCDF file
    
    #sfc_flux_spec is in netCDF file
    
    if ('sfc_roughness_length_cm' in nml_keys):
        int_dict['sfc_roughness_length_cm'] = case_nml['case_config']['sfc_roughness_length_cm']
    
    if ('reference_profile_choice' in nml_keys):
        int_dict['reference_profile_choice'] = case_nml['case_config']['reference_profile_choice']
    
    if ('column_area' in nml_keys):
        int_dict['column_area'] = case_nml['case_config']['column_area']
    
    nml_dict = {'case_config':int_dict}
    
    nml = f90nml.namelist.Namelist(nml_dict)
    
    #print(nml)
    nml.write(filename, force=True)
    
    return(filename)

########################################################################################
#
########################################################################################    
def main():
    #read in arguments
    (case_name, use_area, debug) = parse_arguments()

    setup_logging(debug)

    (case_nml, error) = get_case_nml(case_name)
    if (error):
        message = 'The directory {0} does not contain a config file for case {1}'.format(CASE_NML_DIR, case_name)
        logging.critical(message)
        raise Exception(message)
    else:
        logging.info(case_nml)
        
    (case_data, error) = get_case_data(case_name)
    if (error):
        message = 'The directory {0} does not contain a data file for case {1}'.format(PROCESSED_CASE_DIR, case_name)
        logging.critical(message)
        raise Exception(message)
    else:
        logging.debug(case_data)
    
    fileOUT = write_SCM_case_file(case_nml, case_data, use_area)
    logging.debug("Created {}".format(fileOUT))
    
    write_SCM_nml_file(case_nml)
    
if __name__ == '__main__':
    main()
