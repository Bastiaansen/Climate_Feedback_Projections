# This file contain python codes to compute globally and yearly averaged datasets
# for 2D and 3D data fields

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
        fb = xr.open_dataarray("Data/global/" + fb_name + " .nc")
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


def open_compute_2D(var, file_name, global_file_name):
    
    # Compute 2D field time series
    try:
        var_field = xr.open_dataarray("Data/fields/" + file_name)
        print("Found dataset " + file_name + " -- skipping recomputation")
    except:
        var_field = var.groupby('time.year').mean('time').compute()
        var_field.to_netcdf("Data/fields/" + file_name)
        var_field = xr.open_dataarray("Data/fields/" + file_name)
        print("Computed " + file_name)
        
     # Compute global average
    try:
        var_global = xr.open_dataarray("Data/global/" + global_file_name)
        print("Found dataset " + global_file_name + " -- skipping recomputation")
    except:
        var_global = compute_global2D(var_field)
        var_global = var_global.compute()
        var_global.to_netcdf("Data/global/" + global_file_name)
        print("Computed " + global_file_name)
        
    return var_field, var_global


    
def compute_2D_field_from_3D_obs(Y, p_thickness):
    Y_field = (p_thickness * Y.groupby('time.month')).sum(dim='plev', skipna = True)
    Y_field_yearly = Y_field.groupby('time.year').mean(dim='time')
    Y_field_yearly = Y_field_yearly.compute()
    return Y_field_yearly

def open_compute_3D(var, file_name, global_file_name):

    pdiff = xr.open_dataset("kernels/dp_plev.nc")/100
    pdiff = pdiff.drop('time')
    pdiff = pdiff.assign_coords({'month': pdiff.time + 1}).swap_dims({'time':'month'})
    pthick = pdiff
    
    # Filter out the stratosphere
    p_tropopause = 30000 - 20000 * np.cos(np.deg2rad(var['lat']))
    var_masked = ( var['plev'].isel(plev=slice(0,17)) > p_tropopause ) * var
    
    # Make height integration (so a 2D field is obtained)
    try:
        var_field = xr.open_dataarray("Data/fields/" + file_name)
        print("Found dataset " + file_name + " -- skipping recomputation")
    except:
        var_field = compute_2D_field_from_3D_obs(var, pthick['dp'])
        var_field = var_field.compute()
        var_field.to_netcdf("Data/fields/" + file_name)
        var_field = xr.open_dataarray("Data/fields/" + file_name)
        print("Computed " + file_name)
        
     # Compute global average
    try:
        var_global = xr.open_dataarray("Data/global/" + global_file_name)
        print("Found dataset " + global_file_name + " -- skipping recomputation")
    except:
        var_global = compute_global2D(var_field)
        var_global = var_global.compute()
        var_global.to_netcdf("Data/global/" + global_file_name)
        print("Computed " + global_file_name)
        
    return var_field, var_global


def compute_observable_field(col, model_name):

    dts = import_var(col, model_name, 'ts')
    dta = import_var(col, model_name, 'ta')
    dts3d = dts * xr.ones_like(dta)
    dsalb = import_salb(col, model_name)
    dq = import_var(col, model_name, 'hus')
    
    dlr = dta - dts3d
    
    
    drsdt = import_var(col, model_name, 'rsdt')
    drsut = import_var(col, model_name, 'rsut')
    drlut = import_var(col, model_name, 'rlut')
    
    ## 2D vars
    drsdt_field, dsrdt = open_compute_2D(drsdt, 'drsdt_field.nc', 'drsdt.nc')
    drsut_field, drsut = open_compute_2D(drsut, 'drsut_field.nc', 'drsut.nc')
    drlut_field, drlut = open_compute_2D(drlut, 'drlut_field.nc', 'drlut.nc')
    dsalb_field, dsalb = open_compute_2D(dsalb, 'dsalb_field.nc', 'dsalb.nc')
    dts_field, dts = open_compute_2D(dts, 'dts_field.nc', 'dts.nc')
    
    ## 3D vars
    dlr_field, dlr = open_compute_3D(dlr, 'dlr_field.nc', 'dlr.nc')
    dq_field, dq = open_compute_3D(dq, 'dq_field.nc', 'dq.nc')
    
    return dlr