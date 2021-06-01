# This file contain python codes to compute globally and yearly averaged datasets
# for 2D and 3D data fields
# The computations of feedback contributions have been made following software
# that accompanies the CAM5 kernels made by Pendergrass et al (2018).

##############################
### IMPORT Python packages ###
##############################

import numpy as np
import xarray as xr
import intake



##################################
## IMPORT CMIP DATA FROM GOOGLE ##
##################################

# Convert data catalog into a dictionary of xarray datasets using the intake-esm package
def drop_time_bounds(ds):
    if 'time_bounds' in ds.coords:
        ds = ds.drop('time_bounds')
    elif 'time_bnds' in ds.coords:
        ds = ds.drop('time_bnds')
    
    # Rename the spatial dimensions if necessary (not needed for CESM2)
    if ('longitude' in ds.dims) and ('latitude' in ds.dims):
        ds = ds.rename({'longitude':'lon', 'latitude':'lat'})
    return ds

def import_data(col, model_name, var_name):
    cat = col.search(
        source_id = model_name,
        experiment_id = ['abrupt-4xCO2', 'piControl'],
        table_id = 'Amon',
        variable_id = var_name,
        member_id = 'r1i1p1f1' # The only one for CESM2
        )
        
    ds_dict = cat.to_dataset_dict(preprocess = drop_time_bounds, zarr_kwargs={'consolidated': True, 'decode_times': False})
    
    return ds_dict
    
    
    
#################################
##### IMPORT RELEVANT FIELDS ####
#################################

def import_var(col, model_name, var_name):

    # Search for relevant datasets and load them in
    ds_dict = import_data(col, model_name, var_name)
    
    # Get rid of unnecessary coordinates and split in abrupt and control datasets
    for name, ds in ds_dict.items():
        
        ds = xr.decode_cf(ds)
        
        for coord in ds.coords:
            if coord not in ['lat', 'lon', 'plev', 'time']:
                ds = ds.drop(coord)    
    
        if 'abrupt' in name:
            ds_abr = ds
        elif 'Control' in name:
            ds_ctrl = ds
            
    # Compute the monthly averages for control experiment and track changes in the abrupt experiment compared to those
    ctrl_monthly_av = ds_ctrl[var_name].groupby('time.month').mean(dim='time')
    dVAR = ds_abr[var_name].groupby('time.month') - ctrl_monthly_av
            
    return dVAR.squeeze()
    
def import_salb(col, model_name):
    ds_rsds_dict = import_data(col, model_name, 'rsds')
    ds_rsus_dict = import_data(col, model_name, 'rsus')
    
    for name, ds_rsds in ds_rsds_dict.items():
        ds_rsus = ds_rsus_dict[name]
    
        ds_rsds = xr.decode_cf(ds_rsds)
        ds_rsus = xr.decode_cf(ds_rsus)
        
        for coord in ds_rsus.coords:
            if coord not in ['lat', 'lon', 'plev', 'time']:
                ds_rsus = ds_rsus.drop(coord)
                
        for coord in ds_rsds.coords:
            if coord not in ['lat', 'lon', 'plev', 'time']:
                ds_rsds = ds_rsds.drop(coord)
                
        if 'abrupt' in name:
            salb_abr = ( ds_rsus['rsus'] / ds_rsds['rsds'] ) * 100
        elif 'Control' in name:
            salb_ctrl = ( ds_rsus['rsus'] / ds_rsds['rsds'] ) * 100
        
    ctrl_monthly_av = salb_ctrl.groupby('time.month').mean(dim='time')
    dsalb = salb_abr.groupby('time.month') - ctrl_monthly_av
    return dsalb.squeeze()

    
     
    
##################################
###### 2D and 3D averages ########
##################################

def compute_global2D(Y):
    # Computes globally averaged yearly values for dataarray Y
    
    
    # Global mean average is computed via the surface integral
    # iint f(lat,lon) cos(lat) dS / iint cos(lat) dS
    # So, first: construct help variable containing cosine effective weight values
    weight = np.cos(np.deg2rad(Y['lat'])) * xr.ones_like(Y['lon'])
    
    # Then compute averages (converting integrals to sums)
    Y_globalMean = (
        ( Y * weight ).sum(dim=['lat', 'lon']) /
        weight.sum(dim=['lat', 'lon'])
    ).squeeze()
    
    return Y_globalMean
    

def compute_global3D(Y, p_thickness):
    # Computes globally averaged yearly values for dataarray Y
    
    
    # For 3D variables, the tropopause needs to be masked out
    # p_tropopause = 30000 - 20000 cos(lat)
    p_tropopause = 30000 - 20000 * np.cos(np.deg2rad(Y['lat']))
    Y_masked = ( ( Y['plev'] > p_tropopause ) * Y )
    
    
    # Global mean average is computed via the surface integral
    # iint f(lat,lon) cos(lat) dS / iint cos(lat) dS
    # So, first: construct help variable containing cosine effective weight values
    latlon_weight = np.cos(np.deg2rad(Y['lat'])) * xr.ones_like(Y['lon'])
    p_weight = p_thickness
    weight = latlon_weight * p_weight
    
    # Then compute averages (converting integrals to sums)
    Y_globalMean = (
        ( Y_masked * weight ).sum(dim=['lat', 'lon', 'plev']) /
        weight.sum(dim=['lat', 'lon', 'plev'])
    ).squeeze()
    
    return Y_globalMean   
    
    
def compute_global3D_fb(Y, p_thickness, kernel):
    # Computes globally averaged yearly values for dataarray Y
        
    
    # Global mean average is computed via the surface integral
    # iint f(lat,lon) cos(lat) dS / iint cos(lat) dS
    # So, first: construct help variable containing cosine effective weight values
    weight = np.cos(np.deg2rad(Y['lat'])) * xr.ones_like(Y['lon'])
    
    # Then compute averages (converting integrals to sums)
    Y_globalMean = (
        ( (p_thickness * kernel * Y.groupby('time.month')).sum(dim='plev', skipna = True) * weight ).sum(dim=['lat', 'lon']) /
        weight.sum(dim=['lat', 'lon'])
        ).groupby('time.year').mean(dim='time').squeeze()
    
    return Y_globalMean   
    
def compute_2D_field_from_3D(Y, p_thickness, kernel):
    Y_field = (p_thickness * kernel * Y.groupby('time.month')).sum(dim='plev', skipna = True)
    Y_field_yearly = Y_field.groupby('time.year').mean(dim='time')
    Y_field_yearly = Y_field_yearly.compute()
    return Y_field_yearly
    
    
    
############################################
# COMPUTE 2D/3D FEEDBACKS AND GLOBAL MEANS #
############################################

def compute_feedback2D(kernel, dvar, fb_name, var_name):

    # Compute feedback 2Dfield per year
    try:
        fb_field = xr.open_dataarray("Data/fields/" + fb_name + "_field.nc")
        print("Found dataset " + fb_name + "_field -- skipping recomputation")
    except:
        fb_field = (kernel * dvar.groupby('time.month')).groupby('time.year').mean('time').compute()
        fb_field = fb_field.to_netcdf("Data/fields/" + fb_name + "_field.nc")
        fb_field = xr.open_dataarray("Data/fields/" + fb_name + "_field.nc")
        print("Computed " + fb_name + "_field")
    
    # Compute feedback global mean per yearly
    try:
        fb = xr.open_dataarray("Data/global/" + fb_name + ".nc")
        print("Found dataset " + fb_name + " -- skipping recomputation")
    except:
        fb = compute_global2D(fb_field)
        fb = fb.compute()
        fb.to_netcdf("Data/global/" + fb_name + ".nc")
        print("Computed " + fb_name)
        
    # Compute global mean for variable
    try:
        dvar_GM = xr.open_dataarray("Data/global/" + var_name + ".nc")
        print("Found dataset " + var_name + " -- skipping recomputation")
    except:
        dvar_GM = compute_global2D(dvar.groupby('time.year').mean('time'))
        dvar_GM = dvar_GM.compute()
        dvar_GM.to_netcdf("Data/global/" + var_name + ".nc")
        print("Computed " + var_name)
        
    return fb, dvar_GM
    
    

def compute_feedback3D(kernel, pthick, dvar, fb_name, var_name):

    # Make height integration (so a field is obtained)
    try:
        fb_field = xr.open_dataarray("Data/fields/" + fb_name + "_field.nc")
        print("Found dataset " + fb_name + "_field -- skipping recomputation")  
    except:
        fb_field = compute_2D_field_from_3D(dvar, pthick['dp'], kernel)
        fb_field.to_netcdf("Data/fields/" + fb_name + "_field.nc")
        fb_field = xr.open_dataarray("Data/fields/" + fb_name + "_field.nc")
        print("Computed " + fb_name + "_field")
    
    try:
        fb = xr.open_dataarray("Data/global/" + fb_name + ".nc")
        print("Found dataset " + fb_name + " -- skipping recomputation")
    except:
        fb = compute_global2D(fb_field)
        fb = fb.compute()
        fb.to_netcdf("Data/global/" + fb_name + ".nc")
        print("Computed " + fb_name)
    
    return fb 



##################################
### NORMALIZATION OF Q KERNELS ###
##################################

## Calculate the change in moisture per degree warming at constant relative humidity.
def calcsatspechum(ta, plev):
    # We only care about the monthly average
    t = ta.squeeze().groupby('time.month').mean(dim='time')
    p = plev.squeeze() * xr.ones_like(t)
    p = p/100 # needs to be in hPa
    
    # formulae from Buck (1982)
    es = (1.0007+(3.46e-6*p)) * 6.1121 * np.exp(17.502*(t-273.15) / (240.97+(t-273.15)));
    wsl = .622*es /(p-es) # Saturation mixing ratio wrt liquid water
    es = (1.0003+(4.18e-6*p)) * 6.1115 * np.exp(22.452*(t-273.15) / (272.55+(t-273.15)));
    wsi = .622*es/(p-es) # Saturation mixing ratio wrt ice
    
    # Below freezing we only care about ice, above only about liquid water.
    ws = wsl * (t >= 273.15) + wsi * (t < 273.15)
    
    qs = ws / (1 + ws) # Saturation specific humidity
    
    return qs

def comp_moisture_change(col, model_name, dta):
    # Obtain initial q1
    ds_dict = import_data(col, model_name, 'hus')
    for name, ds in ds_dict.items():
        if 'Control' in name:
            ds = xr.decode_cf(ds)
            for coord in ds.coords:
                if coord not in ['lat', 'lon', 'plev', 'time']:
                    ds = ds.drop(coord)
            q1 = ds['hus'].groupby('time.month').mean(dim='time')
            ds_q_ctrl = ds
    
    ds_dict_ta = import_data(col, model_name, 'ta')
    for name, ds in ds_dict_ta.items():
        ds = xr.decode_cf(ds)
        for coord in ds.coords:
            if coord not in ['lat', 'lon', 'plev', 'time']:
                ds = ds.drop(coord)
        if 'abrupt' in name:
            ds_ta_abr = ds
        elif 'Control' in name:
            ds_ta_ctrl = ds
    
    qs1 = ( calcsatspechum(ds_ta_ctrl['ta'].isel(plev=slice(0,17)), ds_q_ctrl['plev']) )
    qs2 = ( calcsatspechum(ds_ta_abr['ta'].isel(plev=slice(0,17)), ds_q_ctrl['plev']) )
    
    # Compute the change of qs2 and qs1
    dqsdt = (qs2 - qs1) / (dta.squeeze().groupby('time.month').mean(dim='time'))
    # Constant relative humidity
    rh = q1 / qs1

    dqdt = rh * dqsdt
    
    return dqdt
            
    
##############################    
### COMBINATION FUNCTION  ####
##############################

def compute_feedback_timeseries(col, model_name):

    fb = {}
    #####
    # 1 # TEMPERATURE FEEDBACK
    #####
    
    
    ###
    ### Contribution of skin temperature
    ###
    
    dts = import_var(col,model_name, 'ts')
    
    # Read-in the kernel related to ts:
    ts_kernel = xr.open_dataset("kernels/ts.kernel.nc")
    ts_kernel = ts_kernel.drop('time')
    ts_kernel = ts_kernel.assign_coords({'month': ts_kernel.time + 1}).swap_dims({'time':'month'}) # Convert 'time' dim to a 'month' coordinate
    
    dLW_ts, dts_GM = compute_feedback2D(-ts_kernel['FLNT'], dts, 'dLW_ts', 'dts') # Minus sign since feedback contribution of LW has minus sign
    dLW_CS_ts, dts_GM = compute_feedback2D(-ts_kernel['FLNTC'], dts, 'dLW_CS_ts', 'dts')
    
    ###
    ### Contribution of atmospheric temperature
    ###
    
    dta = import_var(col, model_name, 'ta')
    
    # Read-in the kernel related to ta
    ta_kernel = xr.open_dataset("kernels/t.kernel.plev.nc")
    # Only keep the 'FLNT' variable
    for var in ta_kernel.data_vars:
        if var not in ['FLNT', 'FLNTC']:
            ta_kernel = ta_kernel.drop(var)            

    ta_kernel = ta_kernel.drop('time').drop('lev_p') # Only keep the relevant coordinates (time will be coverted one line below)
    ta_kernel = ta_kernel.assign_coords({'month': ta_kernel.time + 1}).swap_dims({'time':'month'}) # Covert 'time' dims to a 'month' coordinate
    ta_kernel = ta_kernel.swap_dims({'ncl1' : 'plev'})  
    
    
    # Read in some pressure weight computed along the kernels (functions as a weight to accurately average over elevation/pressure level)
    pdiff = xr.open_dataset("kernels/dp_plev.nc")/100
    pdiff = pdiff.drop('time')
    pdiff = pdiff.assign_coords({'month': pdiff.time + 1}).swap_dims({'time':'month'})
   
    ## Filter out the stratosphere
    # Filtering is done with a crude approximation of the (pressure)height of the tropopause:
    # p_tropopause = 30000 - 20000 cos(lat)
    p_tropopause = 30000 - 20000 * np.cos(np.deg2rad(dta['lat']))
    # Only use the plev values that are relevant (the first 17);
    # then mask out the values above the tropopause
    dta_masked = ( dta['plev'].isel(plev=slice(0,17)) > p_tropopause ) * dta
    
    dLW_ta = compute_feedback3D(-ta_kernel['FLNT'], pdiff, dta_masked, 'dLW_ta', 'dta')# Minus sign since feedback contribution of LW has minus sign
    dLW_CS_ta = compute_feedback3D(-ta_kernel['FLNTC'], pdiff, dta_masked, 'dLW_CS_ta', 'dta')
    
    fb['temp'] = (dLW_ts + dLW_ta)
    fb['temp-clearsky'] = (dLW_CS_ts + dLW_CS_ta)
    
    
    #####
    # 1A# PLANCK FEEDBACK CONTRIBUTION TO TEMPERATURE FEEDBACK
    #####
    
    # To compute the planck feedback we project the surface temperature (ts) into the height;
    # This mimics the values for atmospheric temperature (ta) that would happen if there were no lapse-rate
    dts3d = dts * xr.ones_like(dta)

    # then mask out the values above the tropopause
    dts3d_masked = ( dts3d['plev'].isel(plev=slice(0,17)) > p_tropopause ) * dts3d
    
    dLW_planck_atm = compute_feedback3D(-ta_kernel['FLNT'],pdiff, dts3d_masked, 'dLW_planck_atm', 'dts3D')
    dLW_planck = (dLW_planck_atm + dLW_ts)
    dLW_planck.to_netcdf("Data/global/dLW_planck.nc")
    
    dLW_CS_planck_atm = compute_feedback3D(-ta_kernel['FLNTC'],pdiff, dts3d_masked, 'dLW_CS_planck_atm', 'dts3D')
    dLW_CS_planck = (dLW_CS_planck_atm + dLW_CS_ts)
    dLW_CS_planck.to_netcdf("Data/global/dLW_CS_planck.nc")
    
    
    fb['planck'] = dLW_planck
    fb['planck-clearsky'] = dLW_CS_planck
    
    #####
    # 1B# LAPSE-RATE FEEDBACK CONTRIBUTION TO TEMPERATURE FEEDBACK
    #####
    
    dlr = dta - dts3d
    dlr_masked = ( dlr['plev'].isel(plev=slice(0,17)) > p_tropopause ) * dlr
    
    dLW_lr = compute_feedback3D(-ta_kernel['FLNT'], pdiff, dlr_masked, 'dLW_lr', 'dlr') # Minus sign since feedback contribution of LW has minus sign
    dLW_CS_lr = compute_feedback3D(-ta_kernel['FLNTC'], pdiff, dlr_masked, 'dLW_CS_lr', 'dlr')
    
    fb['LR'] = dLW_lr
    fb['LR-clearksy'] = dLW_CS_lr
    
    
    
    #####
    # 2 # ALBEDO FEEDBACK
    #####
    
    dsalb = import_salb(col, model_name)
    
    # Read-in the kernel related to albedo
    alb_kernel = xr.open_dataset("kernels/alb.kernel.nc")
    alb_kernel = alb_kernel.drop('time')
    alb_kernel = alb_kernel.assign_coords({'month': alb_kernel.time + 1}).swap_dims({'time':'month'}) # Convert 'time' dim to a 'month' coordinate
    
    dSW_alb, dsalb_GM = compute_feedback2D(alb_kernel['FSNT'], dsalb, 'dSW_alb', 'dsalb')
    dSW_CS_alb, dsalb_GM = compute_feedback2D(alb_kernel['FSNTC'], dsalb, 'dSW_CS_alb', 'dsalb')
        
    fb['albedo'] = dSW_alb 
    fb['albedo-clearsky'] = dSW_CS_alb
    
    #####
    # 3 # WATER VAPOUR FEEDBACK
    #####
    
    dq = import_var(col, model_name, 'hus')
    
    # Load kernels
    q_kernel = xr.open_dataset("kernels/q.kernel.plev.nc")
    for var in q_kernel.data_vars:
        if var not in ['FLNT', 'FSNT', 'FLNTC', 'FSNTC']:
            q_kernel = q_kernel.drop(var)
        
    q_kernel = q_kernel.drop('time').drop('lev_p') # Only keep the relevant coordinates (time will be coverted one line below)
    q_kernel = q_kernel.assign_coords({'month': q_kernel.time + 1}).swap_dims({'time':'month'}) # Covert 'time' dims to a 'month' coordinate

    # Then compute the change in moisture at constant relative humidity (used for normalizing the kernel)
    # Here we only use the monthly averages, as it is too computationally heavy otherwise;
    # this shouldn't be a problem, as we are only interested in a good estimation for the change in moisture,
    # for which this forms an approximation anyways.
    try:
        q_LW_kernel = xr.open_dataarray("kernels/CESM2/q_LW_kernel.nc")
        q_SW_kernel = xr.open_dataarray("kernels/CESM2/q_SW_kernel.nc")
        q_LW_CS_kernel = xr.open_dataarray("kernels/CESM2/q_LW_CS_kernel.nc")
        q_SW_CS_kernel = xr.open_dataarray("kernels/CESM2/q_SW_CS_kernel.nc")
    except:
        try:
            dqdt = xr.open_dataarray("kernels/CESM2/dqdt.nc")
        except:
            dqdt = comp_moisture_change(col, model_name, dta)
            dqdt = dqdt.compute()
            dqdt.to_netcdf("kernels/CESM2/dqdt.nc")
            
        q_LW_kernel = q_kernel['FLNT']/dqdt
        q_SW_kernel = q_kernel['FSNT']/dqdt
        q_LW_CS_kernel = q_kernel['FLNTC']/dqdt
        q_SW_CS_kernel = q_kernel['FSNTC']/dqdt
        
        q_LW_kernel = q_LW_kernel.compute()
        q_SW_kernel = q_SW_kernel.compute()
        q_LW_CS_kernel = q_LW_CS_kernel.compute()
        q_SW_CS_kernel = q_SW_CS_kernel.compute()
        
        q_LW_kernel.to_netcdf("kernels/CESM2/q_LW_kernel.nc")
        q_SW_kernel.to_netcdf("kernels/CESM2/q_SW_kernel.nc")
        q_LW_CS_kernel.to_netcdf("kernels/CESM2/q_LW_CS_kernel.nc")
        q_SW_CS_kernel.to_netcdf("kernels/CESM2/q_SW_CS_kernel.nc")
    
    
    dq_masked = ( dq['plev'].isel(plev=slice(0,17)) > p_tropopause ) * dq
    dLW_q = compute_feedback3D(-q_LW_kernel, pdiff, dq_masked, 'dLW_q', 'dq') # Minus sign since feedback contribution of LW has minus sign
    dSW_q = compute_feedback3D(q_SW_kernel, pdiff, dq_masked, 'dSW_q', 'dq')
    dLW_CS_q = compute_feedback3D(-q_LW_CS_kernel, pdiff, dq_masked, 'dLW_CS_q', 'dq')
    dSW_CS_q = compute_feedback3D(q_SW_CS_kernel, pdiff, dq_masked, 'dSW_CS_q', 'dq')
    
    fb['WV-LW'] = dLW_q
    fb['WV-SW'] = dSW_q
    fb['WV-LW-clearsky'] = dLW_CS_q
    fb['WV-SW-clearsky'] = dSW_CS_q
    
    return fb
    
    
def compute_GMST_imbalance(col, model_name):

    dtas = import_var(col, model_name, 'tas')
    drsdt = import_var(col, model_name, 'rsdt')
    drsut = import_var(col, model_name, 'rsut')
    drlut = import_var(col, model_name, 'rlut')
    drsut_cs = import_var(col, model_name, 'rsutcs')
    drlut_cs = import_var(col, model_name, 'rlutcs')
    
    dIMB = drsdt - drsut - drlut
    dIMB_CS = drsdt - drsut_cs - drlut_cs
    
    ###
    # near-surface atmosphere temperature 'tas'
    ###
    
    # Compute 2Dfield per year
    try:
        dtas_field = xr.open_dataarray("Data/fields/dtas_field.nc")
        print("Found dataset dtas_field -- skipping recomputation")
    except:
        dtas_field = dtas.groupby('time.year').mean('time').compute()
        dtas_field.to_netcdf("Data/fields/dtas_field.nc")
        dtas_field = xr.open_dataarray("Data/fields/dtas_field.nc")
        print("Computed dtas_field")
    
    # Compute global mean per yearl
    try:
        dGMST = xr.open_dataarray("Data/global/dGMST.nc")
        print("Found dataset dGMST -- skipping recomputation")
    except:
        dGMST = compute_global2D(dtas_field)
        dGMST = dGMST.compute()
        dGMST.to_netcdf("Data/global/dGMST.nc")
        print("Computed dGMST")
        
    ###
    # Radiative imbalance (full sky)
    ###
        
    try:
        dIMB_field = xr.open_dataarray("Data/fields/dIMB_field.nc")
        print("Found dataset dIMB_field -- skipping recomputation")
    except:
        dIMB_field = dIMB.groupby('time.year').mean('time').compute()
        dIMB_field.to_netcdf("Data/fields/dIMB_field.nc")
        dIMB_field = xr.open_dataarray("Data/fields/dIMB_field.nc")
        print("Computed dIMB_field")
    
    try:
        dIMB = xr.open_dataarray("Data/global/dIMB.nc")
        print("Found dataset dIMB -- skipping recomputation")
    except:
        dIMB = compute_global2D(dIMB_field)
        dIMB = dIMB.compute()
        dIMB.to_netcdf("Data/global/dIMB.nc")
        print("Computed dIMB")
    
    ###
    # Radiative imbalance (clear sky)
    ###
    
    try:
        dIMB_CS_field = xr.open_dataarray("Data/fields/dIMB_CS_field.nc")
        print("Found dataset dIMB_CS_field -- skipping recomputation")
    except:
        dIMB_CS_field = dIMB_CS.groupby('time.year').mean('time').compute()
        dIMB_CS_field.to_netcdf("Data/fields/dIMB_CS_field.nc")
        dIMB_CS_field = xr.open_dataarray("Data/fields/dIMB_CS_field.nc")
        print("Computed dIMB_CS_field")
    
    try:
        dIMB_CS = xr.open_dataarray("Data/global/dIMB_CS.nc")
        print("Found dataset dIMB_CS -- skipping recomputation")
    except:
        dIMB_CS = compute_global2D(dIMB_CS_field)
        dIMB_CS = dIMB_CS.compute()
        dIMB_CS.to_netcdf("Data/global/dIMB_CS.nc")
        print("Computed dIMB_CS")
    
    
    
    return dGMST, dIMB, dIMB_CS