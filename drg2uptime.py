#!/usr/bin/env python3
# AKIMOTO
# 2023-07-05
# 2023-07-06 update1 (adding the polar graph)
# 2023-07-07 update2 (calculating LST, adding xtick/ytick-right/left of rcParams and ProgressBar)
# 2023-09-24 update3 (add sckedule list)
# 2024-03-07
# 2024-09-23

import os
import sys
import datetime
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_sun
from astropy.time import Time

plt.rcParams["xtick.direction"]     = "in"       
plt.rcParams["ytick.direction"]     = "in"       
plt.rcParams["xtick.minor.visible"] = True       
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["xtick.top"]           = True
plt.rcParams["xtick.bottom"]        = True
plt.rcParams["ytick.left"]          = True
plt.rcParams["ytick.right"]         = True 
plt.rcParams["xtick.major.size"]    = 5          
plt.rcParams["ytick.major.size"]    = 5          
plt.rcParams["xtick.minor.size"]    = 3          
plt.rcParams["ytick.minor.size"]    = 3          
plt.rcParams["axes.grid"]           = True
plt.rcParams["grid.color"]          = "lightgray"
plt.rcParams["axes.labelsize"]      = 12
plt.rcParams["font.size"]           = 12

def RaDec2AltAz_LST_Sun(object_ra: float, object_dec: float, observation_time: float, delta_time: float, latitude: float, longitude: float, height: float) -> float :
    
    delta_time_sample = 0
    if 10 <= delta_time < 100 :
        delta_time_sample = int(delta_time / 10)
    elif 100 <= delta_time < 1000 :
        delta_time_sample = int(delta_time / 20)
    else :
        delta_time_sample = int(delta_time / 100)
    
    
    # calculating AZ-EL
    location_geocentrice = EarthLocation.from_geocentric(latitude, longitude, height, unit=u.m)
    location_geodetic    = EarthLocation.to_geodetic(location_geocentrice)
    location_lon_lat     = EarthLocation(lon=location_geodetic.lon, lat=location_geodetic.lat, height=location_geodetic.height)
    
    delta_time           = np.linspace(0, delta_time, delta_time_sample +1) * u.second
    obstime              = Time(observation_time, scale="utc") + delta_time
    object_ra_dec        = SkyCoord(ra=object_ra * u.deg, dec=object_dec * u.deg)
    AltAz_coord_azel     = AltAz(location=location_lon_lat, obstime=obstime)
    object_altaz         = object_ra_dec.transform_to(AltAz_coord_azel)

    # calculating LST
    one_day_hour         = np.linspace(0, 24, 24*10) * u.hour
    object_ra_dec        = SkyCoord(ra=object_ra * u.deg, dec=object_dec * u.deg)
    obstime_lst          = Time(observation_time, scale="utc" ,location=location_lon_lat) + one_day_hour
    AltAz_coord_lst      = AltAz(location=location_lon_lat, obstime=obstime_lst)
    lst_altaz            = object_ra_dec.transform_to(AltAz_coord_lst)
    lst                  = obstime_lst.sidereal_time('apparent')
    
    
    return obstime.datetime, object_altaz.az.deg, object_altaz.alt.deg, lst.hour, lst_altaz.az.deg, lst_altaz.alt.deg, obstime, AltAz_coord_azel

def target_RaDec_rectime(drg) :
    
    target_ra_dec  = {}
    target_rectime = defaultdict(list)
    d = 0
    target_skd_list = []
    drg_open = open(drg, "r").readlines()
    for drg_line in drg_open :
        if "2000.0" in drg_line :
            
            drg_line_split = drg_line.split()[:8]
            
            target_ra_dec[drg_line_split[0]] = SkyCoord(["%s %s" % (":".join(drg_line_split[2:5]), ":".join(drg_line_split[5:8]))], unit=(u.hourangle, u.deg)).to_string()[0]

        if "PREOB" in drg_line :

            drg_line_split = drg_line.split()
            del drg_line_split[1:4], drg_line_split[3:]

            conv2datetime_start = datetime.datetime.strptime(drg_line_split[1], "%y%j%H%M%S")
            conv2datetime_end   = conv2datetime_start + datetime.timedelta(seconds=int(drg_line_split[2]))
            target_skd_list.append([drg_line_split[0], conv2datetime_start, conv2datetime_end, drg_line_split[2]])
            
            drg_line_split[1] = conv2datetime_start.strftime("%Y-%m-%dT%H:%M:%S")
            
            target_rectime[drg_line_split[0]].append(" ".join(drg_line_split[1:3]))
            
            if d == 0 :
                observation_start_date = drg_line_split[1]
            observation_end_date = drg_line_split[1]
            d += 1

    return target_ra_dec, target_rectime, observation_start_date, observation_end_date, target_skd_list


#
# antenna coordinate in xml file
#
antenna_position_x, antenna_position_y, antenna_position_z, antenna_name = -3502544.5870, +3950966.2350, +3566381.1920, "Yamaguchi-32m antenna"
#antenna_position_x, antenna_position_y, antenna_position_z, antenna_name = -3961788.9740, +3243597.4920, +3790597.6920, "Hitachi-32m antenna"

drg_input = sys.argv[1]

print(f"Extracting target coordinates and an observation schedule from \"{drg_input}\"")
target_ra_dec_out, target_rectime_out, observation_start_date, observation_end_date, target_skd_list  = target_RaDec_rectime(drg_input)

print("Plot...")
fig_azel , axs_azel  = plt.subplots(2, 1, figsize=(12,9), sharex=True, tight_layout=True)
fig_polar, axs_polar = plt.subplots(1, 1, figsize=(8, 9)             , tight_layout=True, subplot_kw={'projection': 'polar'})
fig_lst  , axs_lst   = plt.subplots(2, 1, figsize=(12,9), sharex=True, tight_layout=True)
for i, target in enumerate(target_ra_dec_out.keys()) :
    
    label_check = False
    
    target_ra, target_dec = target_ra_dec_out[target].split()
    
    cm = plt.colormaps.get_cmap("tab20")
    rgb = cm.colors[i]
    
    for target_rectime_length in target_rectime_out[target] :

        target_rectime_start, target_length = target_rectime_length.split()

        # Az & EL
        target_datetime, target_az, target_el, target_lst, _, _, _, _ = RaDec2AltAz_LST_Sun(float(target_ra), float(target_dec), target_rectime_start, float(target_length), antenna_position_x, antenna_position_y, antenna_position_z)

        axs_azel[0].plot(target_datetime, target_az, label=f"{target}" if label_check == False else "", color=rgb) # az
        axs_azel[1].plot(target_datetime, target_el, label=f"{target}" if label_check == False else "", color=rgb) # el
        
        axs_polar.plot(target_az * np.pi / 180.0, target_el, label=f"{target}" if label_check == False else "", color=rgb)
        
        label_check = True
    
    # LST
    _, _, _, target_lst, target_lst_az, target_lst_el, sun_one_day_hour, sun_coordinate = RaDec2AltAz_LST_Sun(float(target_ra), float(target_dec), observation_start_date, (Time(observation_end_date)-Time(observation_start_date)).datetime.total_seconds(), antenna_position_x, antenna_position_y, antenna_position_z)
    
    target_lst_az_el_zip = list(map(list, zip(target_lst, target_lst_az, target_lst_el)))
    target_lst_az_el_zip.sort()
    target_lst_az_el_zip = np.array(target_lst_az_el_zip)
    
    axs_lst[0].plot(target_lst_az_el_zip[:,0], target_lst_az_el_zip[:,1], label=f"{target}", color=rgb) # az
    axs_lst[1].plot(target_lst_az_el_zip[:,0], target_lst_az_el_zip[:,2], label=f"{target}", color=rgb) # el
    
# sun
sun = get_sun(sun_one_day_hour).transform_to(sun_coordinate)
axs_azel[0].plot(sun_one_day_hour.datetime, sun.az, ":", label="sun", color="k") # sun az
axs_azel[1].plot(sun_one_day_hour.datetime, sun.alt, ":", label="sun", color="k") # sun el


formatter = mdates.DateFormatter("%H:%M")
axs_azel[0].xaxis.set_major_formatter(formatter)
axs_azel[0].set_ylim(0,360)
axs_azel[0].set_yticks(np.linspace(0,360,9))
axs_azel[0].set_ylabel("AZ (deg)")
axs_azel[0].legend(ncols=7, bbox_to_anchor=(0, 1.15), loc='upper left', borderaxespad=0, fontsize=10)
axs_azel[1].xaxis.set_major_formatter(formatter)
axs_azel[1].set_ylim(0,90)
axs_azel[1].set_yticks(np.linspace(0,90,10))
axs_azel[1].set_xlabel(f"The observation time: {observation_start_date} to {observation_end_date} UT in {antenna_name}")
axs_azel[1].set_ylabel("EL (deg)")

axs_polar.set_rticks(np.linspace(0,90,10))
axs_polar.set_rmax(0)
axs_polar.set_rmin(90)
axs_polar.set_theta_direction(-1)
axs_polar.set_theta_offset(np.pi/2)
axs_polar.legend(ncols=5, bbox_to_anchor=(0, 1.15), loc='upper left', borderaxespad=0, fontsize=10)

axs_lst[0].set_xlim(0,24)
axs_lst[0].set_ylim(0,360)
axs_lst[0].set_xticks(np.linspace(0,24,25))
axs_lst[0].set_yticks(np.linspace(0,360,9))
axs_lst[0].set_ylabel("AZ (deg)")
axs_lst[0].legend(ncols=7, bbox_to_anchor=(0, 1.15), loc='upper left', borderaxespad=0, fontsize=10)
axs_lst[1].set_xlim(0,24)
axs_lst[1].set_ylim(0,90)
axs_lst[1].set_xticks(np.linspace(0,24,25))
axs_lst[1].set_yticks(np.linspace(0,90,10))
axs_lst[1].set_xlabel(f"LST in {antenna_name}")
axs_lst[1].set_ylabel("EL (deg)")


skd_list = "   target     start (UT)               end (UT)      integ (sec)  slew (sec)\n"
for i in range(len(target_skd_list)) :
    if i < len(target_skd_list)-1 :
        slew_time = (target_skd_list[i+1][1] - target_skd_list[i][2]).total_seconds()
    scan_target = target_skd_list[i][0]
    scan_start  = target_skd_list[i][1].strftime("%Y-%m-%dT%H:%M:%S")
    scan_end    = target_skd_list[i][2].strftime("%Y-%m-%dT%H:%M:%S")
    scan_integ  = target_skd_list[i][3]
    if i == len(target_skd_list)-1 : slew_time = 0; print(slew_time)
    skd_list += "%10s   %20s   %20s   %+5s   %5d\n"  % (scan_target, scan_start, scan_end, scan_integ, slew_time)
print(skd_list)
fig, axs = plt.subplots(1, 1, figsize=(5,12), tight_layout=True)
plt.text(0, 0, skd_list, size="x-small")
plt.gca().axis('off')
        


plt.show()



