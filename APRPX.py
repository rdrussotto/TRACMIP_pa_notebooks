# New version of the APRP code designed to work with XArray, 
# in a situation where we have multiple models and runs
# in one XArray dataset, as it is in the TRACMIP climatology.
# (don't load directly from NetCDF4)
#
# Should do all models and months simultaneously.
#
# Rick Russotto
# Started 18 March 2019

import numpy as np
import xarray as xr

#Issues: 
#XArray doesn't support masked arrays...
#So will have to either extract regular NumPy arrays
#or 
#But original functions not designed to work with multiple models,
#so will have to write functions that work with extracted
#NumPy arrays for each model.
#Write wrapper to loop for all models. 
#Need to fill in dimensions of time, lat, lon


#Main function to run the other functions
#Loop through input list of models, run APRP for each one, 
#and return xarray dataset containing the relevant results
#
#Input arguments:
#ds: XArray dataset to run APRP on
#exp1: experiment considered the base or control state
#exp2: experiment representing the perturbed state
#models: list of models to run it on
#
#Return: XArray dataset containing changes in TOA net SW radiation
#due to changes in cloud properties, non-cloud atmosphere, and 
#surface albedo ("cloud", "noncloud" and "surface" variables)--
#add more later if needed
def aprp_main_TRACMIP(ds, exp1, exp2, models):
    
    #Pre-initialize 4D arrays for each variable
    cloud = np.empty([12, len(ds.lat), len(ds.lon), len(models)])
    noncloud = np.empty([12, len(ds.lat), len(ds.lon), len(models)])
    surface = np.empty([12, len(ds.lat), len(ds.lon), len(models)])
    
    #Run APRP for each model and put relevant variables into the 
    #arrays
    i = 0 #model counter
    for model in models:
        
        print('Running APRP calculations for model: '+model)
        
        #Calculations for basic state
        dict1A = setup(ds, exp1, model)
        dict1B = parameters(dict1A)
        #Calculations for perturbed state
        dict2A = setup(ds, exp2, model)
        dict2B = parameters(dict2A)
        #Calculations for change between the states
        dictC = d_albedo(dict1A, dict1B, dict2A, dict2B)
        
        #Convert masked arrays to arrays with nans 
        cloud_nan = dictC['cloud'].data
        cloud_nan[np.where(dictC['cloud'].mask == True)] = np.nan
        noncloud_nan = dictC['noncloud'].data
        noncloud_nan[np.where(dictC['noncloud'].mask == True)] = np.nan
        surface_nan = dictC['surface'].data
        surface_nan[np.where(dictC['surface'].mask == True)] = np.nan
        
        #Populate the empty arrays for the variables of interest
        cloud[:,:,:,i] = dictC['cloud']
        noncloud[:,:,:,i] = dictC['noncloud']
        surface[:,:,:,i] = dictC['surface']
                
        i = i + 1
        
        
    #Construct the dataset with the results and return
    results = xr.Dataset({'cloud': (['time', 'lat', 'lon', 'model'], cloud), 
                          'noncloud': (['time', 'lat', 'lon', 'model'], noncloud), 
                          'surface': (['time', 'lat', 'lon', 'model'], surface)},
                         coords={'time': ds.time.values, 
                                 'lat': ds.lat.values,
                                 'lon': ds.lon.values, 
                                 'model': models})
    
    return results

#Wrapper for 2 different datasets, for just one model and month
def aprp_main_TRACMIP_2Datasets(ds_clim, ds_perturb):
    
    
    #Run APRP and put relevant variables into the arrays
    #Wait: the "setup" routine subsetted by model--need to fix this up
    #Make model and experiment arguments optional
    
    #Calculations for basic state
    dict1A = setup(ds_clim)
    dict1B = parameters(dict1A)
    #Calculations for perturbed state
    dict2A = setup(ds_perturb)
    dict2B = parameters(dict2A)
    #Calculations for change between the states
    dictC = d_albedo(dict1A, dict1B, dict2A, dict2B)
    
    #Convert masked arrays to arrays with nans 
    cloud_nan = dictC['cloud'].data
    cloud_nan[np.where(dictC['cloud'].mask == True)] = np.nan
    noncloud_nan = dictC['noncloud'].data
    noncloud_nan[np.where(dictC['noncloud'].mask == True)] = np.nan
    surface_nan = dictC['surface'].data
    surface_nan[np.where(dictC['surface'].mask == True)] = np.nan
    
    #Construct the dataset with the results and return
    results = xr.Dataset({'cloud': (['lat', 'lon'], dictC['cloud']), 
                          'noncloud': (['lat', 'lon'], dictC['noncloud']), 
                          'surface': (['lat', 'lon'], dictC['surface'])},
                         coords={'lat': ds_clim.lat.values,
                                 'lon': ds_clim.lon.values})
    
    
    return results


#Setup function replacing "loadNetCDF"
#In the original code, took multi-annual means in this function; 
#don't need to do that for TRACMIP (already done),
#might need to for other applications.
#
#Wait: I ran this twice, once for each file, 
#along with the "parameters" function.
#Can combine "loadNetCDF" and "parameters" 
#into one function for each experiment. 
#
#Input arguments: 
#"ds" (an XArray dataset)
#"exp": Experiment to subset to
#"model": 
def setup(ds, exp='none', model='none'):
    if exp == 'none' and model == 'none':
        ds1 = ds.copy()
    elif exp != 'none' and model == 'none':
        ds1 = ds.sel(exp=exp)
    elif exp == 'none' and model != 'none': 
        ds1 = ds.sel(model=model)
    else:
        ds1 = ds.sel(exp=exp, model=model)

    #Extract lat and lon as NumPy arrays
    lat = ds.lat.values
    lon = ds.lon.values
    
    #Extract each of the necessary
    #monthly mean variables as 3D NumPy arrays.
    #Have tested, should automatically result in 
    #time, lat, lon dimensions.
    #(Skip step of taking monthly means--already done)
    
    m_rsds = ds1['rsds'].values 
    m_rsus = ds1['rsus'].values 
    m_rsut = ds1['rsut'].values 
    m_rsdt = ds1['rsdt'].values 
    m_rsutcs = ds1['rsutcs'].values 
    m_rsdscs = ds1['rsdscs'].values 
    m_rsuscs = ds1['rsuscs'].values 
    m_clt = ds1['clt'].values 
    
    #Mask out any values greater than 10^10
    m_rsds = np.ma.masked_greater(m_rsds,1.e10)
    m_rsus = np.ma.masked_greater(m_rsus,1.e10)
    m_rsut = np.ma.masked_greater(m_rsut,1.e10)
    m_rsdt = np.ma.masked_greater(m_rsdt,1.e10)
    m_rsutcs = np.ma.masked_greater(m_rsutcs,1.e10)
    m_rsdscs = np.ma.masked_greater(m_rsdscs,1.e10)
    m_rsuscs = np.ma.masked_greater(m_rsuscs,1.e10)
    m_clt = np.ma.masked_greater(m_clt,1.e10)
    
    #Calculate overcast versions of rsds, rsus, rsut
    
    #Mask zero values of cloud fraction so you 
    #don't calculate overcast values in clear-sky pixels
    m_clt = np.ma.masked_values(m_clt, 0)
    #c = m_clt/100. #convert from percentage to fraction
    c = m_clt #For TRACMIP, clt is already a fraction
    m_rsdsoc = (m_rsds-(1.-c)*(m_rsdscs))/c  #Can derive this algebraically from Taylor et al., 2007, Eq. 3
    m_rsusoc = (m_rsus-(1.-c)*(m_rsuscs))/c
    m_rsutoc = (m_rsut-(1.-c)*(m_rsutcs))/c
    
    #Mask zero values of the downward SW radiation (polar night)
    m_rsds = np.ma.masked_values(m_rsds, 0)
    m_rsdscs = np.ma.masked_values(m_rsdscs, 0)
    m_rsdsoc = np.ma.masked_values(m_rsdsoc, 0)
    m_rsdt = np.ma.masked_values(m_rsdt, 0)
    
    #Return dictionary with all the variables calculated here 
    #(called "dictA" because calculated in first function called)
    dictA = dict()
    dictA['rsds'] = m_rsds
    dictA['rsus'] = m_rsus
    dictA['rsut'] = m_rsut
    dictA['rsdt'] = m_rsdt
    dictA['rsutcs'] = m_rsutcs
    dictA['rsdscs'] = m_rsdscs
    dictA['rsuscs'] = m_rsuscs
    dictA['clt'] = m_clt
    dictA['lat'] = lat
    dictA['lon'] = lon
    dictA['rsdsoc'] = m_rsdsoc
    dictA['rsusoc'] = m_rsusoc
    dictA['rsutoc'] = m_rsutoc
    dictA['c'] = c #Cloud fraction as fraction, not %
    
    return dictA


###  ----- Same as original APRP code below this line -----  ###


#Calculate the tuning parameters for the idealized single-layer radiative transfer model
#for the individual time period (i.e. control or warmed)
#See Figure 1 of Taylor et al., 2007, and other parts of that paper. Equations referenced are from there.
#
#Based on Ting's "parameters.m".
#
#Inputs: the dictionary output by loadNetCDF
#Outputs: a dictionary of additional outputs
def parameters(dictA):
    #Clear-sky parameters
    a_clr = dictA['rsuscs']/dictA['rsdscs'] #Surface albedo   
    Q = dictA['rsdscs']/dictA['rsdt'] #Ratio of incident surface flux to insolation
    mu_clr = dictA['rsutcs']/dictA['rsdt']+Q*(1.-a_clr) #Atmospheric transmittance (Eq. 9)  #"Invalid value in divide"
    ga_clr = (mu_clr-Q)/(mu_clr-a_clr*Q) #Atmospheric scattering coefficient (Eq. 10)

    #Overcast parameters
    a_oc = dictA['rsusoc']/dictA['rsdsoc'] #Surface albedo
    Q = dictA['rsdsoc']/dictA['rsdt'] #Ratio of incident surface flux to insolation
    mu_oc = dictA['rsutoc']/dictA['rsdt']+Q*(1.-a_oc) #Atmospheric transmittance (Eq. 9)
    ga_oc = (mu_oc-Q)/(mu_oc-a_oc*Q) #Atmospheric scattering coefficient (Eq. 10)   
    
    #Calculating cloudy parameters based on clear-sky and overcast ones 
    #Difference between _cld and _oc: _cld is due to the cloud itself, as opposed to 
    #scattering and absorption from all constituents including clouds in overcast skies.
    mu_cld = mu_oc / mu_clr            #Eq. 14
    ga_cld = (ga_oc-1.)/(1.-ga_clr)+1. #Eq. 13  
    
    #Save the relevant variables to a dictionary for later use
    dictB = dict()
    dictB['a_clr'] = a_clr
    dictB['a_oc'] = a_oc
    dictB['mu_clr'] = mu_clr
    dictB['mu_cld'] = mu_cld
    dictB['ga_clr'] = ga_clr
    dictB['ga_cld'] = ga_cld
    
    #Ting saved a cloud fraction variable here--I did this in earlier function instead.
    
    return dictB





#Calculations for the differences between time periods
def d_albedo(dict1A, dict1B, dict2A, dict2B):
    
    #First, Ting set cloud values that were masked in one time period 
    #equal to the value in the other time period, assuming no cloud changes.
    #I'll take these variables out of the dictionary before modifying them. 
    a_oc1 = dict1B['a_oc']
    a_oc2 = dict2B['a_oc']
    a_oc2[a_oc2.mask == True] = a_oc1[a_oc2.mask == True]
    a_oc1[a_oc1.mask == True] = a_oc2[a_oc1.mask == True]
    
    mu_cld1 = dict1B['mu_cld']
    mu_cld2 = dict2B['mu_cld']
    mu_cld2[mu_cld2.mask == True] = mu_cld1[mu_cld2.mask == True]
    mu_cld1[mu_cld1.mask == True] = mu_cld2[mu_cld1.mask == True]
    
    ga_cld1 = dict1B['ga_cld']
    ga_cld2 = dict2B['ga_cld']
    ga_cld2[ga_cld2.mask == True] = ga_cld1[ga_cld2.mask == True]
    ga_cld1[ga_cld1.mask == True] = ga_cld2[ga_cld1.mask == True]
    
    #Now a bunch of calls to the "albedo" function to see how the albedo changes as a result of 
    #...the changes to each of the radiative components. 
    
    #Retrieve other variables from dictionaries to make calls to albedo shorter/more readable
    c1 = dict1A['c']
    c2 = dict2A['c']
    a_clr1 = dict1B['a_clr']
    a_clr2 = dict2B['a_clr']
    mu_clr1 = dict1B['mu_clr']
    mu_clr2 = dict2B['mu_clr']
    ga_clr1 = dict1B['ga_clr']
    ga_clr2 = dict2B['ga_clr']
    
    #Base state albedo
    A1 = albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)
    A2 = albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2)
    
    #Change in albedo due to each component (Taylor et al., 2007, Eq. 12b)
    dA_c =      .5*(albedo(c2, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c1, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2))    
    
    dA_a_clr =  .5*(albedo(c1, a_clr2, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr1, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld2))
                
    dA_a_oc =   .5*(albedo(c1, a_clr1, a_oc2, mu_clr1, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc1, mu_clr2, mu_cld2, ga_clr2, ga_cld2))
               
    dA_mu_clr = .5*(albedo(c1, a_clr1, a_oc1, mu_clr2, mu_cld1, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr1, mu_cld2, ga_clr2, ga_cld2))               
                 
    dA_mu_cld = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld2, ga_clr1, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld1, ga_clr2, ga_cld2))                 
                 
    dA_ga_clr = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr2, ga_cld1)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr1, ga_cld2))
                 
    dA_ga_cld = .5*(albedo(c1, a_clr1, a_oc1, mu_clr1, mu_cld1, ga_clr1, ga_cld2)-A1)+.5*(
                 A2-albedo(c2, a_clr2, a_oc2, mu_clr2, mu_cld2, ga_clr2, ga_cld1))
                 
    #Set changes due to overcast or cloudy sky parameters, or changes to clouds themselves, to zero
    #...if cloud fraction is less than 3% in either time period
    dA_a_oc[dict1A['c'] < .03] = 0
    dA_a_oc[dict2A['c'] < .03] = 0
    dA_mu_cld[dict1A['c'] < .03] = 0
    dA_mu_cld[dict2A['c'] < .03] = 0
    dA_ga_cld[dict1A['c'] < .03] = 0
    dA_ga_cld[dict2A['c'] < .03] = 0
    dA_c[dict1A['c'] < .03] = 0
    dA_c[dict2A['c'] < .03] = 0
    
    #Combine different components into changes due to surface albedo, atmospheric clear-sky and atmospheric cloudy-sky
    dA_a = dA_a_clr + dA_a_oc                #Eq. 16a
    dA_cld = dA_mu_cld + dA_ga_cld + dA_c    #Eq. 16b
    dA_clr = dA_mu_clr + dA_ga_clr           #Eq. 16c
    
    #Set all planetary albedo changes = zero when incoming solar radaition is zero    
    #(This will replace NaNs with zeros in the polar night--affects annual means)
    dA_a[dict2A['rsdt']<0.1] = 0
    dA_clr[dict2A['rsdt']<0.1] = 0
    dA_cld[dict2A['rsdt']<0.1] = 0
    dA_a_clr[dict2A['rsdt']<0.1] = 0
    dA_a_oc[dict2A['rsdt']<0.1] = 0
    dA_mu_cld[dict2A['rsdt']<0.1] = 0
    dA_ga_cld[dict2A['rsdt']<0.1] = 0
    dA_c[dict2A['rsdt']<0.1] = 0
    dA_mu_clr[dict2A['rsdt']<0.1] = 0
    dA_ga_clr[dict2A['rsdt']<0.1] = 0
    
    
    #Calculate radiative effects in W/m^2 by multiplying negative of planetary albedo changes by downward SW radation
    #(This means positive changes mean more downward SW absorbed)
    surface = -dA_a*dict2A['rsdt']   #Radiative effect of surface albedo changes
    surface[dict2A['rsdt']<0.1] = 0
    surface = np.ma.masked_outside(surface, -100, 100) # Ting called this "boundary for strange output"
    
    cloud = -dA_cld*dict2A['rsdt']   #Radiative effect of cloud changes
    cloud[dict2A['rsdt']<0.1] = 0
    cloud = np.ma.masked_outside(cloud, -100, 100) # Ting called this "boundary for strange output"
    
    noncloud = -dA_clr*dict2A['rsdt'] #Radiative effect of non-cloud SW changes (e.g. SW absorption)
    noncloud[dict2A['rsdt']<0.1] = 0
    
    #Broken down further into the individual terms in Eq. 16
    surface_clr = -dA_a_clr*dict2A['rsdt']    #Effects of surface albedo in clear-sky conditions
    surface_clr[dict2A['rsdt']<0.1] = 0
    
    surface_oc = -dA_a_oc*dict2A['rsdt']      #Effects of surface albedo in overcast conditions
    surface_oc[dict2A['rsdt']<0.1] = 0
    
    cloud_c = -dA_c*dict2A['rsdt']            #Effects of changes in cloud fraction
    cloud_c[dict2A['rsdt']<0.1] = 0
    
    cloud_ga = -dA_ga_cld*dict2A['rsdt']      #Effects of atmospheric scattering in cloudy conditions
    cloud_ga[dict2A['rsdt']<0.1] = 0
    
    cloud_mu = -dA_mu_cld*dict2A['rsdt']      #Effects of atmospheric absorption in cloudy conditions
    cloud_mu[dict2A['rsdt']<0.1] = 0
    
    noncloud_ga = -dA_ga_clr*dict2A['rsdt']   #Effects of atmospheric scattering in clear-sky conditions
    noncloud_ga[dict2A['rsdt']<0.1] = 0
    
    noncloud_mu = -dA_mu_clr*dict2A['rsdt']   #Effects of atmospheric absorption in clear-sky conditions
    noncloud_mu[dict2A['rsdt']<0.1] = 0
    
    #Calculate more useful radiation output
    CRF = dict1A['rsut'] - dict1A['rsutcs'] - dict2A['rsut'] + dict2A['rsutcs'] #Change in cloud radiative effect
    cs = dict1A['rsutcs'] - dict2A['rsutcs']  #Change in clear-sky upward SW flux at TOA
    
    #Define a dictionary to return all the variables calculated here
    dictC = dict()
    dictC['A1'] = A1
    dictC['A2'] = A2
    dictC['dA_c'] = dA_c
    dictC['dA_a_clr'] = dA_a_clr
    dictC['dA_a_oc'] = dA_a_oc
    dictC['dA_mu_clr'] = dA_mu_clr
    dictC['dA_mu_cld'] = dA_mu_cld
    dictC['dA_ga_clr'] = dA_ga_clr
    dictC['dA_ga_cld'] = dA_ga_cld
    dictC['dA_a'] = dA_a
    dictC['dA_cld'] = dA_cld
    dictC['dA_clr'] = dA_clr
    dictC['surface'] = surface
    dictC['cloud'] = cloud
    dictC['noncloud'] = noncloud
    dictC['surface_clr'] = surface_clr
    dictC['surface_oc'] = surface_oc
    dictC['cloud_c'] = cloud_c
    dictC['cloud_ga'] = cloud_ga
    dictC['cloud_mu'] = cloud_mu
    dictC['noncloud_ga'] = noncloud_ga
    dictC['noncloud_mu'] = noncloud_mu
    dictC['CRF'] = CRF
    dictC['cs'] = cs
    
    return dictC





    
    
#Function to calculate the planetary albedo, A.
#Inputs: (see Fig. 1 of Taylor et al., 2007)
#   c: fraction of the region occupied by clouds
#   a_clr: clear sky surface albedo (SW flux up / SW flux down)
#   a_oc: overcast surface albedo
#   mu_clr: clear-sky transmittance of SW radiation
#   mu_cld: cloudy-sky transmittance of SW radiation
#   ga_clr: clear-sky atmospheric scattering coefficient
#   ga_cld: cloudy-sky atmospheric scattering coefficient
def albedo(c, a_clr, a_oc, mu_clr, mu_cld, ga_clr, ga_cld): #Labeled with equation numbers from Taylor et al. 2007
    mu_oc = mu_clr*mu_cld                                                            #Eq. 14
    ga_oc =  1. - (1.-ga_clr)*(1.-ga_cld)                                            #Eq. 13
    A_clr = mu_clr*ga_clr + mu_clr*a_clr*(1.-ga_clr)*(1.-ga_clr)/(1.-a_clr*ga_clr)   #Eq. 7 (clear-sky)
    A_oc = mu_oc*ga_oc + mu_oc*a_oc*(1.-ga_oc)*(1.-ga_oc)/(1.-a_oc*ga_oc)            #Eq. 7 (overcast sky)
    A = (1-c)*A_clr + c*A_oc                                                         #Eq. 15
    return A
    