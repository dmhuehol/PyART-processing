# -*- coding: utf-8 -*-
"""
STARTING POINT FOR PYART PROCESSING

DESCRIPTION: Use this script to easily process and plot radar data in PyART.

MANUAL VARIABLES (These must be changed by hand):
    inpath: File path containing the radar data to be processed.
    outpath: File path where images and/or processed data should be saved
    wildcard: String found in all of the data files. Used to automatically determine the radar type.
    scan_strat: String describing the scan strategy ('PPI', 'RHI', or 'Sector').
    
    snow_rate_bool: True/False on whether to calculate the Rasmussen snow rate (rescaled reflectivity).
    vdiv_bool: True/False on whether to calculate the vertical divergence of horizontal velocity (RHI only)
    
    contour_bool: True/False whether to overlay contours of a secondary field over a plot of a base field.
    base_field: Field to be plotted normally.
    contour_field: Field to be contoured over the base field.
    contour_levels: Array of values at which to draw contours, e.g. [10,20]
    
    azi_overlay: Dictionary containing settings for drawing RHI azimuths on a PPI plot.
    
    dealias_bool: True/False on whether to dealias velocity data or leave folded.
    save_cfradial_bool: True/False on whether to save a CF/Radial data files containing dealiased velocity data.
    
SEMIAUTOMATIC VARIABLES (Take care of themselves for KASPR, CHILL, and NEXRAD, but can be controlled manually):
    radar_type: Used throughout the toolkit to discriminate between the different radars.
    fields: Data types observed.
    ranges: Ranges corresponding to the fields.
    plot_bool: True/False on whether or not to make plots. Leave at True unless only using PyART to save off data.
    x_lim: x-limit for image output.
    y_lim: y-limit for image output.
    colorbar_labels: Labels for the colorbars corresponding to the fields and ranges.
    name2dealias: Name of the folded velocity field.
    new_name: String that PyART will name the dealiased velocity field.
    nyquist_vel: Nyquist velocity used by dealiaser. MUST BE SET MANUALLY FOR NEXRAD.
    cmaps: Colormaps to use when plotting the fields.
    rhoHV_mask: Dictionary containing settings for rhoHV filter. Enabled and set to retain data between 0.45 and 1.2 by default.
    NCP_mask: Dictionary containing settings for NCP filter. Disabled by default.
    SNR_mask: Dictionary containing settings for SNR filter. Disabled by default.
    Zdr_offset: Dictionary containing settings for accounting for Zdr offset on some radars. Disabled for CHILL and NEXRAD. Enabled and accounts for 1.2 dB offset on KASPR.
    mountain_clutter_bool: True/False to run/not run the mountain removal code. MUST BE SET MANUALLY FOR NEXRAD.
    
    
    

Created on Thu Aug 06 14:01:55 2015

@author: thecakeisalie

Version date: 11/17/2019
Daniel Hueholt
North Carolina State University
Undergraduate Research Assistant at Environment Analytics

"""
# Load necessary packages (set appropriate working directory!)
import os
os.chdir('C:\\Users\\lmtomkin\\Documents\\winter_storms\\PyART-processing')
import colormap
import gen_fun
import run_fun
import gc
import time
import numpy as np
from multiprocessing import Process
import sys

######### Define Variables #############
### Path Variables
#inpath = 'C:\\Users\\lmtomkin\\Documents\\winter_storms\\data\\KASPR\\'+sys.argv[1]+'\\'+sys.argv[2]
#outpath = 'C:\\Users\\lmtomkin\\Documents\\winter_storms\\images\\KASPR\\'+sys.argv[1]+'\\'+sys.argv[2]
inpath = 'C:\\Users\\lmtomkin\\Documents\\winter_storms\\data\\KASPR\\20191201\\RHI'
outpath = 'C:\\Users\\lmtomkin\\Documents\\winter_storms\\images\\KASPR\\20191201\\RHI'

### File and Data Variables ###
wildcard = 'KASPR' #Common wildcards are below
#CHL: CSU-CHILL S-band
#CHX: CSU-CHILL X-band
#KASPR: SBU Ka-band
#KCYS, KFTG: NEXRAD radars in FRONT network (other NEXRAD data works too)
# KASPR, CHILL, and NEXRAD are the only fully-implemented radars for now
#KXA: Dallas HF-S
#KXAS: NBC5 StormRanger

scan_strat = sys.argv[2] #Possible entries are below
#PPI: Plan view at a specific tilt angle
#RHI: Cross section along a specific azimuth
#Sector: Plan view at a specific tilt angle, with a confined set of azimuths. ONLY PARTIALLY IMPLEMENTED
# PyART can be used with other scan strategies but these are not yet supported in this toolkit.

# Determine radar type
if wildcard == 'CHX':
    radar_type = 'CHILL'
elif wildcard == 'CHL':
    radar_type = 'CHILL'
elif wildcard == 'KASPR':
    radar_type = 'KASPR'
else:
    if len(wildcard)==4:
        radar_type = 'NEXRAD'
    else:
        print("WARNING: Unknown radar type! Processing with NEXRAD settings to give best chance of success.")
        radar_type = 'NEXRAD'
    
# List of strings (field names)
#HF-S: ['DBZH','DBZV','ZDR','RHOHV','PHIDP','SNRHC','SNRVC','dealiasVELH']
#KASPR: ['correlation_coefficient','differential_phase','differential_reflectivity','linear_depolarization_ratio','linear_depolarization_ratio_hc_hh','mean_doppler_velocity','mean_doppler_velocity_folded','reflectivity','reflectivity_xpol_htx','reflectivity_xpol_vtx','snr','snr_xpol_htx','snr_xpol_vtx','spectrum_width']
#StormRanger:
#
#fields: strings, refer to the data types observed
#ranges: numerical tuples, the ranges corresponding to the fields. ORDER OF THE RANGES MUST MATCH ORDER OF FIELDS
if radar_type=='CHILL':
    fields = ['reflectivity','dealiased_velocity','corrected_differential_reflectivity','spectrum_width','cross_correlation_ratio','normalized_coherent_power','one_way_differential_phase','two_way_differential_phase']
    ranges = [(-5,25),(-20,20),(-1,2),(0,4),(0.4,1),(0,1),(-60,-120),(-1.5,1.5)] #CSU-CHILL (winter) X-band
elif radar_type=='KASPR':
    fields = ['correlation_coefficient','differential_reflectivity','PyART_dealiased_velocity','reflectivity','spectrum_width','linear_depolarization_ratio','snr']
    ranges = [(0.5,1),(-2,2),(-35,35),(-5,40),(0,3),(-40,-20),(0,100)] #KASPR (commonly-used) (winter)
elif radar_type=='NEXRAD':
    fields = ['reflectivity','dealiased_velocity','spectrum_width','cross_correlation_ratio','differential_reflectivity']
    ranges = [(-5,50),(-50,50),(0,12),(0,1),(-4,4)]

# Other common ranges
#ranges = [(-5,25),(-40,40),(-1,2),(0,8),(0.4,1),(0,1),(-30,-60),(-0.5,0.5)] #CSU-CHILL (winter) S-band
#ranges = [(-5,65),(-40,40),(-3,5),(0,8),(0.5,1),(0,1),(-5,5)] #CSU-CHILL (summer)
#ranges = #StormRanger
#ranges = [(0,60),(0,60),(-20,0),(0,1),(0,180),(0,1),(0,1),(-40,40)] #HF-S, with broken Zdr

### Plotting Variables ###
plot_bool = True #Determines whether to make plots or not

# Set x and y limits for radar plot
if scan_strat != 'RHI':
    x_lim = [-60,60]
    y_lim = [-60,60]
    if radar_type == 'NEXRAD':
        x_lim = [-175,175]
        y_lim = [-175,175]
if scan_strat == 'RHI':
    x_lim = [0,60]
    y_lim = [0,9]
    if radar_type == 'KASPR':
        x_lim = [-40,40]
        y_lim = [0,12]

# Other useful x-limits
#x_lim = [-375,375] #HF-S PPI
#x_lim =  #StormRanger
#x_lim = [0,75] #CSU-CHILL X-band max RHI range
# Other useful y-limits
#y_lim = [-375,375] #HF-S PPI
#y_lim = #StormRanger
#y_lim = [0,16] #CSU-CHILL RHI (summer)

# List of strings (colorbar labels)
#colorbar_labels = ['DBZ (dBZ)','DBZ (dBZ)','ZDR (DB)','rhoHV','PhiDP','SNR','SNR','V (m/s)'] #HF-S
#colorbar_labels = #StormRanger
#colorbar_labels = ['rhoHV','PhiDP','ZDR (DB)','V (m/s)','DBZ (dBZ)','Width (m/s)'] #KASPR few
if radar_type=='CHILL':
    colorbar_labels = ['DBZ (dBZ)','V (m/s)','ZDR (dB)','Width (m/s)','rhoHV','NCP','PhiDP (deg)','KDP (deg/km)'] #CSU-CHILL
elif radar_type=='KASPR':
    colorbar_labels = ['rhoHV','Zdr (dB)','V (m/s)','Z (dBZ)','Spectral Width (m/s)','LDR (dB)','SNR'] #KASPR common
elif radar_type=='NEXRAD':
    colorbar_labels = ['DBZ (dBZ)','V (m/s)','Width (m/s)','rhoHV','Zdr (dB)']

### Dealiasing Variables ###
dealias_bool = True
save_cfradial_bool = False #Save the radar data with dealiased velocity in a CF/Radial file

if radar_type=='CHILL':
    name2dealias = 'corrected_velocity' #CSU-CHILL
    new_name = 'dealiased_velocity'
    if wildcard=='CHL':
        nyquist_vel = 27.5039 #S-band CHILL
    elif wildcard=='CHX':
        nyquist_vel = 25.893 #X-band CHILL
elif radar_type=='KASPR':
    name2dealias = 'mean_doppler_velocity_folded' #KASPR
    new_name = 'PyART_dealiased_velocity' #KASPR
    nyquist_vel= 9.999 #KASPR
elif radar_type=='NEXRAD':
    name2dealias = 'velocity'
    new_name = 'dealiased_velocity'
    nyquist_vel = 26.389 #KCYS 2018 02 01, this can vary and will need to be detected automatically
#name2dealias = 'VELH' #HF-S
#new_name = 'dealiasVELH' #HF-S
#nyquist_vel = 8.5048 #HF-S
#name2dealias = #StormRanger
#new_name = #StormRanger

#######################################

# Adject inpath and outpath for easier writing
inpath = inpath + '\\'
outpath = outpath + '\\'

print(inpath)
print(outpath)

# Load color maps
# L_range is the luminance range. If the colors come out too dark then raise the minimum luminance value.
max_luminance = 100
min_luminance = 0
LCH = colormap.LCH_Spiral()[0]
LCH_zdr = colormap.LCH_Spiral(nc = 100, np = .3, offset = 0, reverse = 1, L_range = [max_luminance, min_luminance], name = 'LCH_zdr')[0]
LCH_wid = colormap.LCH_Spiral(nc = 100, np = .3, offset = 45, reverse = 0, L_range = [max_luminance, min_luminance], name = 'LCH_wid')[0]
Int = colormap.PID_Integer()
IntCHILL = colormap.PID_Integer_CHILL()
cuckooPalette = colormap.cuckoo()

# Set colormaps
#CAUTION: DO NOT USE 'seismic' FOR VELOCITY DATA. 'seismic' is asymmetric around its zero value, which skews the apparent magnitude
#of negative values relative to positive values. 'RdBu_r' is preferred.
#cmaps = [LCH,LCH,LCH_zdr,'bone_r','cividis','copper','copper','RdBu_r'] #HF-S
#cmaps = #StormRanger
if radar_type=='CHILL':
    cmaps = [LCH,'RdBu_r',LCH_zdr,LCH_wid,'bone_r','copper','magma',cuckooPalette] #CSU-CHILL
elif radar_type=='KASPR':
    cmaps = ['bone_r',LCH_zdr,'RdBu_r',LCH,LCH_wid,'inferno','copper'] #KASPR common
elif radar_type=='NEXRAD':
    cmaps = [LCH,'RdBu_r',LCH_wid,'bone_r',LCH_zdr]

# Make filelist
filelist = gen_fun.get_filelist(inpath, wildcard, False)

# Plausibility check on inputs
if len(cmaps)!=len(colorbar_labels):
    raise IndexError('Check number of colormaps and colorbars!')


# Data quality
#   Remove values outside of a given range for some variable.
#   Common Z ranges: (-5, 50) (KASPR), (-10, 80) (HF-S)
#   Common PhiDP ranges: (0, 180) (HF-S)
#   Common rhoHV ranges: (0.45, 1.2) (good for most cases)
#   Common NCP ranges: (0.05, 1.2) (good for most cases)
#   Most common approach is to mask on SNR if available, or rhoHV if not. Sometimes an NCP mask is also needed.
#   Z, Zdr, and PhiDP masks risk removing real echo if used incautiously.

Z_mask = {
        "bool": False,
        "range": (-5,45)
        }
Zdr_mask = {
        "bool": False,
        "range": (-2,2)
        }
PhiDP_mask = {
        "bool": False,
        "range": (-5, 5)
        }
rhoHV_mask = {
        "bool": True,
        "range": (0.45, 1.2)
        }
NCP_mask = {
        "bool": False,
        "range": (0.15, 1.2)
        }
SNR_mask = {
        "bool": False,
        "range": (5, 100)
        }

#   Account for Zdr offset on radars such as KASPR or HF-S
if radar_type=='CHILL':
    Zdr_offset = {
            "bool": False,
            "offset": 0.0
            }
elif radar_type=='KASPR':
    Zdr_offset = {
            "bool": True,
            "offset": 1.2}
else:
    Zdr_offset = {
            "bool": False,
            "offset": 0.0}
       
#   Logicals for derived data
snow_rate_bool = True #Derive the Rasmussen snow rate from reflectivity
vdiv_bool = True #Calculate the vertical divergence of horizontal dealiased velocity (RHI only)

#   Other logicals
if radar_type=='CHILL':
    mountain_clutter_bool = True #Remove mountain clutter from CHILL data automatically
elif radar_type=='KASPR':
    mountain_clutter_bool = False #No mountains in New York City!
elif radar_type=='NEXRAD':
    mountain_clutter_bool = True #Needs to be set manually

#   Contour overlay settings
contour_bool = False
base_field = 'dealiased_velocity' #Background field
contour_field = 'reflectivity' #Field that will be contoured over the base field
contour_levels = [15] #Values for contours

#   Settings controlling RHI azimuth overlay on PPI plot
if scan_strat == 'PPI':
    azi_overlay = {
            "bool": False,
            "azi_lines": [134,224],
            "color": "#105456",
            "linewidth": 4.2
            }
elif scan_strat == 'RHI':
    azi_overlay = {
        "bool": False,
        "azi_lines": [0],
        "color": "#105456",
        "linewidth": 4.2
        }
        

#   Parse through filelist
# Processing is conducted within individual processes that are started and ended with each file.
# Otherwise, Python's memory use will grow until the computer's memory is expended,
# crashing the system. This is not a memory leak in PyART! It's just a consequence of how Python's memory handling
# deals with long-running processes.

if __name__== '__main__':
        length_filelist = np.size(filelist)
        print("Processing in progress!")
        for item in range(0,length_filelist):
            filelist_ind=filelist[item]
            # Instantiate a new process.
            p = Process(target=run_fun.parse_filelist,args=(filelist_ind, inpath, outpath, radar_type, fields, ranges, plot_bool, 
                                                                cmaps, colorbar_labels, x_lim, y_lim, scan_strat, 
                                                                dealias_bool, save_cfradial_bool, name2dealias, new_name, nyquist_vel, Z_mask, Zdr_mask, PhiDP_mask,
                                                                rhoHV_mask, NCP_mask, SNR_mask, Zdr_offset, snow_rate_bool, vdiv_bool, mountain_clutter_bool,
                                                                contour_bool, base_field, contour_field, contour_levels, azi_overlay))
            p.start()
            #print(p.is_alive())
            p.join()
            p.terminate() #End the process
            numleft = (length_filelist - item -1)
            print(numleft) #Displays how many files remain to be processed in current job
            
        print("Completed!")

gc.collect()
