"""
scanner.py 

Author: Thomas Plunkett

Date: 24/03/24

Purpose:

Defines functions needed to create the schedule dashboard.

"""
# Import necessary packages 
import pandas as pd
from datetime import datetime, date, timedelta
import numpy as np
from dateutil import tz
import os
import time
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.coordinates import Angle
from astroplan import Observer, FixedTarget, ObservingBlock
from astroplan.constraints import AtNightConstraint, AirmassConstraint, MoonSeparationConstraint, AltitudeConstraint
from astroplan.scheduling import Transitioner, Schedule, SequentialScheduler, PriorityScheduler
from astroplan.plots import plot_sky
import matplotlib.pyplot as plt  
import plotly.express as px
from repeat_constraint import *
import gspread

# Set up the location
date_str = str(datetime.now().date())
bisdee_tier = EarthLocation.from_geodetic(lon = 147.2878*u.deg, lat = -42.4311*u.deg, height = 646*u.m)
GHO = Observer(location=bisdee_tier, name='Greenhill Observatory', timezone='Australia/Tasmania')

# Define necessary functions 
def convert_RA(ra_deg):
    ang = Angle(ra_deg, unit = u.deg)
    sex_ang = ang.to_string('hourangle', sep = ':')
    return sex_ang

def convert_dec(dec_deg):
    ang = Angle(dec_deg, unit = u.deg)
    sex_ang = ang.to_string(unit = 'degree', sep = ':')
    return sex_ang

def read_targets(date=date_str):
    """
    Read the targets from a Google spreadsheet, labelled by date
    
    params:
    date - Today's date (str)
    
    return:
    target_df - A dataframe containing the list of targets
    """
    key_df = pd.read_csv('key.csv')
    jsonfile = key_df['json'][0]
    keystr = key_df['key'][0]
    
    gc = gspread.service_account(filename=jsonfile)
    targets = gc.open_by_key(keystr)
    sheet_name = date
    try:
        target_df = pd.DataFrame(targets.worksheet(sheet_name).get_all_records())
    except:
        target_df = pd.DataFrame()
        
    return target_df

def construct_schedule(ms_df):
    """
    Make a dataframe containing the night's schedule, ordered on priority and compatibility
    with observing constraints (moon sep, airmass, at night)
    
    params:
    ms_df - The target dataframe
    
    return:
    sched - The schedule dataframe 
    """
    # Define constants relating to the H50
    read_out = 10 * u.second
    slew_rate = 8.0*u.deg/u.second # For 50cm
    transitioner = Transitioner(slew_rate, {'filter':{'default': 30*u.second}})
    sched = pd.DataFrame()
    
    # Calculate the global constraints on the observing session
    start, end = GHO.tonight()
    global_constraints = [MoonSeparationConstraint(min=30*u.deg), AltitudeConstraint(min = 20*u.deg, boolean_constraint = True), AirmassConstraint(max = 3, boolean_constraint = False), AtNightConstraint.twilight_astronomical()]
    
    # Loop over the dataframe and find requested targets, add to observing block
    
    # We start by building subsets based on the requested no. of repeats
    for k in range(0, ms_df['REPEATS'].max()+1):
        blocks = []
        target_coord = []
        target = []
        filter_list = []
        df_sub = ms_df[ms_df.REPEATS >= k]
        df_sub = df_sub.reset_index(drop=True)
       
        # Loop through this subset and create FixedTarget objects to perform calculations on
        for i in range(0, len(df_sub)):
            target_coord += [SkyCoord(ra = df_sub['RA'][i], dec = df_sub['DEC'][i], unit = (u.hourangle, u.deg))]
            target += [FixedTarget(coord = target_coord[i], name = df_sub['ID'][i])]
            
            # Check the no. of images requested is greater than zero
            if df_sub['IMAGES'][i] > 0:
                fltr_list = df_sub['FILTERS'][i].split(',')
                
                # Iterate through filters
                for fltr in fltr_list:
                    # If repeats are requested, add constraint of not reobserving for at least 30 mins
                    if k > 0:
                        repeat_sub = sched[sched.target == str(df_sub['ID'][i])].reset_index(drop=True)
                        
                        # Find the time when the last observation ended to use in constraint
                        if len(repeat_sub) != 0:
                            end_time = Time(repeat_sub['end time (UTC)'][k*len(fltr_list)-1], format='iso', scale = 'utc')
                            b = ObservingBlock.from_exposures(target[i], df_sub['PRIORITY'][i], df_sub['EXP. TIME [s]'][i]*u.second, df_sub['IMAGES'][i], read_out, configuration = {'filter': fltr}, constraints = [RepeatConstraint(t_end = end_time)])
                            blocks.append(b)
                    
                    else:
                        b = ObservingBlock.from_exposures(target[i], df_sub['PRIORITY'][i], df_sub['EXP. TIME [s]'][i]*u.second, df_sub['IMAGES'][i], read_out, configuration = {'filter': fltr})
                        blocks.append(b)
                
        # Initialize the sequential scheduler with the constraints and transitioner
        seq_scheduler = SequentialScheduler(constraints = global_constraints, observer = GHO, transitioner = transitioner)
        #pri_scheduler = PriorityScheduler(constraints = global_constraints, observer = GHO, transitioner = transitioner)

        # Initialize a Schedule object, to contain the new schedule
        sequential_schedule = Schedule(start, end)
        #priority_schedule = Schedule(start, end)

        # Call the schedule with the observing blocks and schedule to schedule the blocks
        seq_scheduler(blocks, sequential_schedule)
        #pri_scheduler(blocks, priority_schedule)

        # Convert to a schedule table
        sched_temp = seq_scheduler.schedule.to_table(show_transitions=False).to_pandas()
        #sched = pri_scheduler.schedule.to_table(show_transitions=False).to_pandas()

        # Clean up the table
        for p in range(0, len(sched_temp)): 
            filter_list += [sched_temp['configuration'][p]['filter']]

        sched_temp['Filters'] = filter_list

        del sched_temp['configuration']
        
        sched = pd.concat([sched, sched_temp])
        sched = sched.sort_values('start time (UTC)')
        sched = sched.reset_index(drop=True)
                
    return sched.sort_values('start time (UTC)')
                          
def construct_sky(sched):
    """
    Create a polar plot of the sky, with zenith angle as radius and azimuth as angle. 
    
    params:
    sched -  The schedule dataframe 
    
    return:
    fig - The plotly figure 
    
    """
    # Create a list of unique targets in the schedule
    target_list = sched.target.unique()
    sky_df = pd.DataFrame()
    
    # Iterate through each target and create a new dataframe with alt/az
    for ob_id in target_list:
        sched_sub = sched[sched.target == ob_id].reset_index(drop=True)
        target_coord = SkyCoord(ra = sched_sub.ra[0], dec = sched_sub.dec[0], unit = (u.deg, u.deg))
        target = FixedTarget(coord = target_coord, name = ob_id)
        
        times = Time(sched_sub['start time (UTC)'].to_list(), format='iso', scale='utc')
        
        altitude = GHO.altaz(times, target).alt * (1/u.deg)
        azimuth = GHO.altaz(times, target).az * (1/u.deg) 
        
        df_temp = pd.DataFrame({'ID':[ob_id]*len(altitude), 'Z': 90-altitude, 'Az': azimuth, 'Time': sched_sub['start time (UTC)'].to_list()})
        
        sky_df = pd.concat([sky_df, df_temp])
           
    fig = px.scatter_polar(sky_df, r='Z', theta='Az', color = 'ID', template='plotly_dark')
    fig.update_polars(radialaxis={'range':[0,90]})
        
    return fig

    