"""
scanner.py 

Author: Thomas Plunkett

Date: 29/06/23

Purpose:

Defines functions needed to create schedule dashboard, including functions to scan 
MOA and GAIA alerts pages

"""

# Import necessary packages 
import pandas as pd
from datetime import datetime, date, timedelta
import numpy as np
from dateutil import tz
import os
import time
from astropy.coordinates import Angle
from astroplan import Observer, FixedTarget, ObservingBlock
from astroplan.scheduling import Transitioner
from astroplan.constraints import TimeConstraint, AtNightConstraint, AirmassConstraint, MoonSeparationConstraint
from astroplan.scheduling import SequentialScheduler, PriorityScheduler
from astroplan.scheduling import Schedule
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, SkyCoord

bisdee_tier = EarthLocation.from_geodetic(lon = 147.3*u.deg, lat = -42.43*u.deg, height = 646*u.m)
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

def get_tdelta(path):
    """
    A function to calculate the time in seconds since a file was last created/modified
    
    params: 
    path - The filename/path
    
    return:
    t_delta - The time difference in seconds
    """
    if os.path.isfile(path):
        timestamp = os.path.getmtime(path)
        datestamp = datetime.fromtimestamp(timestamp)
        t_delta = Time(datestamp, format='datetime', scale='local') - Time(datetime.now(), format='datetime', scale='local')
    
    return abs(t_delta.sec)
    

def get_MOA():
    """
    Function to retrieve dataframes from MOA alerts page for all, high mag. and anomoly events
    by scraping the hmtl tables using Pandas. Updates the tables every 2 hours, but saves for faster running.
    
    return:
    moa_df - The full alerts page table for all events of the year
    moa_df2 - The high magnification alerts for the year
    moa_df3 - The anomoly alerts tables for the year
    """
    try:
        if os.path.isfile('MOA_ALL.csv') and get_tdelta('MOA_ALL.csv') < (2*3600):
            moa_df = pd.read_csv('MOA_ALL.csv')
        else:
            moa_df = pd.read_html('https://www.massey.ac.nz/~iabond/moa/alert2023/alert.php')[0] # All alerts for the year
            moa_df.to_csv('MOA_ALL.csv')
        
        moa_df2 = pd.read_html('https://www.massey.ac.nz/~iabond/moa/alert2023/himag.php') # Himag alerts
        #moa_df2.to_csv('MOA_HIMAG.csv')
            
        if os.path.isfile('MOA_ANOM.csv') and get_tdelta('MOA_ANOM.csv') < (2*3600):
            moa_df3 = pd.read_csv('MOA_ANOM.csv')
        else:
            moa_df3 = pd.read_html('http://iral2.ess.sci.osaka-u.ac.jp/~moa/anomaly/2023/index.html')[0] # Anomolies
            moa_df3.to_csv('MOA_ANOM.csv')
    except:
        raise
        if os.path.isfile('MOA_ALL.csv') and os.path.isfile('MOA_ANOM.csv'):
            moa_df = pd.read_csv('MOA_ALL.csv')
            moa_df2 = pd.read_html('https://www.massey.ac.nz/~iabond/moa/alert2023/himag.php')
            moa_df3 = pd.read_csv('MOA_ANOM.csv')
        else:
            print('Unable to access MOA alerts pages and no previous files could be found! Check your internet connection...')
            moa_df = pd.DataFrame()
            moa_df2 = pd.DataFrame()
            moa_df3 = pd.DataFrame()

    return moa_df, moa_df2, moa_df3

def construct_MOA():
    """
    Function to construct the desired dataframe for displaying interesting MOA events
    
    return:
    final_df - A combined dataframe with the desired anomolies and high mag. events from MOA
    """
    himag_df = pd.DataFrame()
    anom_df = pd.DataFrame()
    str_list = []
    moa_df, moa_df2, moa_df3 = get_MOA()
    
    try:
        # Make the dataframe for high mag. events
        for i in range(0, len(moa_df2)):
            himag_df = pd.concat([himag_df, moa_df[moa_df['tE(days)'] == float(moa_df2[i][2][1].split(' ')[0])]])  
        himag_sub = himag_df.copy()
        himag_sub = himag_sub[himag_sub['Ibase'] < 19.0]
        himag_sub['Comment'] = ['High Magnification']*len(himag_sub)
        
        # Make the dataframe for the anomolies
        for i in range(0, len(moa_df3)):
            anom_df = pd.concat([anom_df, moa_df[moa_df['ID'] == moa_df3['MOA ID'][i].replace('MOA-','')]])
            str_list += [str(moa_df3['comment'][i]) + ' , mass ratio: ' + str(moa_df3['q'][i])]
        
        anom_df['Comment'] = str_list
        
        anom_sub = anom_df.copy()
        anom_sub = anom_sub[anom_sub['Ibase']<20.0]
        
        # Finally construct the final MOA dataframe by combining previous
        final_df = pd.concat([anom_sub, himag_sub])
        del final_df['tE(days)']
        del final_df['Assessment']
        del final_df['Amax']
        final_df = final_df.reset_index(drop=True)
        for j in range(0, len(final_df)):
            t0_date = datetime.strptime(final_df['t0'][j][:-3].replace('-',' '), "%Y %b %d").date()
            if t0_date < (date.today() - timedelta(days = 30)):
                final_df = final_df.drop(j)
        final_df = final_df.reset_index(drop=True)
        #final_df.sort_values('t0', ascending = False)
        final_df['ID'] = 'MOA-'+final_df['ID']
    
    except:
        raise 
        print('Unable to construct the MOA dataframe... Check internet connection...')
        final_df = pd.DataFrame()
   
    return final_df 

def get_gaia():
    """
    A function to either download or read the full GAIA alerts csv into a pandas datafra,e
    """
    try:
        # Read in the events csv from GAIA alerts website
        if os.path.isfile('GAIA_ALL.csv') and get_tdelta('GAIA_ALL.csv') < (24*3600):
            gaia_df = pd.read_csv('GAIA_ALL.csv')
        else:
            gaia_df = pd.read_csv('http://gsaweb.ast.cam.ac.uk/alerts/alerts.csv')
            gaia_df.to_csv('GAIA_ALL.csv')
            
    except:
        if os.path.isfile('GAIA_ALL.csv'):
            gaia_df = pd.read_csv('GAIA_ALL.csv')
        else:
            print('Unable to construct the Gaia dataframe... Check internet connection')
            gaia_df = pd.DataFrame()
    
    return gaia_df
            

def construct_gaia():
    """
    A function to construct a dataframe with interesting gaia microlensing
    candidates, in the same format as the moa dataframes. Updates the table every 24 hours.
    
    Return:
    gaia_final - The gaia dataframe
    """
    gaia_df = get_gaia() 
    ind_list = []
    
    if len(gaia_df) != 0:
        # Sort so only visible, candidate microlensing events appear
        gaia_df = gaia_df[gaia_df[' DecDeg'] < -(90-42.5)]
        gaia_df = gaia_df.reset_index(drop=True)

        # Find a subset of event in 2023
        for i in range(0, len(gaia_df)):
            if type(gaia_df[' Comment'][i])==str:
                if 'microlensing' in gaia_df[' Comment'][i] and '2023' in gaia_df[' Date'][i]:
                    ind_list += [i]

        # Make a subset that drops faint targets
        gaia_sub = gaia_df.iloc[ind_list]
        gaia_sub = gaia_sub[gaia_sub[' HistoricMag'] < 20.0]
        gaia_sub = gaia_sub.sort_values(' HistoricMag')
        gaia_sub = gaia_sub.reset_index(drop=True)

        # If the event is older than ~ 3 months, ignore it
        for k in range(0, len(gaia_sub)):
            g_date = datetime.strptime(gaia_sub[' Date'][k].replace('-',' '), "%Y %m %d %H:%M:%S").date()
            if g_date < (date.today() - timedelta(days = 90)):
                gaia_sub = gaia_sub.drop(k)

        # Match the dataframe format to the MOA events
        gaia_final = pd.DataFrame()
        gaia_final['ID'] = gaia_sub['#Name']
        gaia_final['RA (J2000.0)'] = convert_RA(gaia_sub[' RaDeg']*u.deg)
        gaia_final['Dec (J2000.0)'] = convert_dec(gaia_sub[' DecDeg']*u.deg)
        gaia_final['t0'] = [np.nan]*len(gaia_sub)
        gaia_final['Ibase'] = gaia_sub[' HistoricMag']
        gaia_final['Comment'] = gaia_sub[' Comment']
    
    else:
        gaia_final = pd.DataFrame(columns = gaia_df.columns)
    
    return gaia_final

def get_GRB():
    """
    A function to retrieve a list of GRBs from the IceCube website and read into a
    Pandas Dataframe
    """
    cols = ['GRB_name','T0','ra','decl','pos_error','T90', 'MJD']
    try:
        if os.path.isfile('GRB_ALL.csv') and get_tdelta('GAIA_ALL.csv') < (2*3600):
            grb_df = pd.read_csv('GRB_ALL.csv')
        else:
            grb_raw = pd.read_csv('https://user-web.icecube.wisc.edu/~grbweb_public/Summary_table.txt', sep=' ', skipfooter=1, skiprows=3, skipinitialspace=True, header=None, engine='python')
            grb_df = grb_raw[[0,2,3,4,5,6,14]]
            grb_df.columns = cols
            grb_df.to_csv('GRB_ALL.csv')
    except:
        if os.path.isfile('GRB_ALL.csv'):
            grb_df = pd.read_csv('GRB_ALL.csv')
        else: 
            print('Unable to construct the GRB dataframe... Check internet connection!')
            grb_df = pd.DataFrame(columns = cols)
            
    return grb_df

def construct_GRB():
        """
        A function to construct the final GRB dataframe ready to join with the microlensing events
        """
        grb_df = get_GRB()
        current_mjd = Time(datetime.now(), format='datetime', scale='local').mjd
        
        if len(grb_df) != 0:
            # Sort so only visible, candidate microlensing events appear
            grb_sub = grb_df[grb_df['decl'] < -(90-42.5)]
            grb_sub = grb_sub.reset_index(drop=True)
            
            # If the event is older than ~ 3 months, ignore it
            for k in range(0, len(grb_sub)):
                if grb_sub['MJD'][k] < (current_mjd - 90):
                    grb_sub = grb_sub.drop(k)
                    
            # Match the dataframe format to the MOA events
            grb_final = pd.DataFrame()
            grb_final['ID'] = grb_sub['GRB_name']
            grb_final['RA (J2000.0)'] = convert_RA(grb_sub['ra']*u.deg)
            grb_final['Dec (J2000.0)'] = convert_dec(grb_sub['decl']*u.deg)
            grb_final['t0'] = ['']*len(grb_sub)
            grb_final['Ibase'] = ['']*len(grb_sub)
            grb_final['Comment'] = ['']*len(grb_sub)
            
            for t in range(0, len(grb_sub)):
                grb_final['t0'][t] = str(grb_sub['GRB_name'][t][3:9]) + ' ' + str(grb_sub['T0'][t])
                grb_final['Comment'][t] = 'T90: ' + str(grb_sub['T90'][t]) + ' seconds, Pos. Error: ' + str(grb_sub['pos_error'][t]) + ' degrees'
            
        else:
            grb_final = pd.DataFrame(columns = grb_df.columns)
            
        return grb_final

def construct_final():
    """
    A function to create a master dataframe with all interesting microlensing events from
    MOA and Gaia alerts pages and GRBs from IceCube website
    
    Return:
    ms_df - The master dataframe for microlensing 
    """
    moa_final = construct_MOA()
    gaia_final = construct_gaia()
    grb_final = construct_GRB()
    
    ms_df = pd.concat([moa_final, gaia_final, grb_final])
    ms_df = ms_df.reset_index(drop=True)
    ms_df['Priority'] = [0.5]*len(ms_df)
    ms_df['Exp. Time (s)'] = [30.0]*len(ms_df)
    ms_df['No. of exposures'] = [0]*len(ms_df)
    ms_df['No. of repeats'] = [0]*len(ms_df)
    ms_df['Filters'] = ["r'"]*len(ms_df)
    del ms_df[ms_df.columns[0]]
    
    return ms_df

def construct_schedule(ms_df):
    """
    Make a dataframe containing the night's schedule, ordered on priority and compatibility
    with observing constraints (moon sep, airmass, at night)
    
    params:
    ms_df - The microlensing dataframe
    
    return:
    sched - The schedule dataframe 
    """
    target_coord = []
    target = []
    blocks = []
    read_out = 10 * u.second
    slew_rate = 8.0*u.deg/u.second # For 50cm
    transitioner = Transitioner(slew_rate, {'filter':{'default': 30*u.second}})
    
    # Create FixedTarget objects for all events in the dataframe
    for i in range(0, len(ms_df)):
        target_coord += [SkyCoord(ra = ms_df['RA (J2000.0)'][i], dec = ms_df['Dec (J2000.0)'][i], unit = (u.hourangle, u.deg))]
        target += [FixedTarget(coord = target_coord[i], name = ms_df['ID'][i])]
    
    # Calculate the constraints on the observing session
    start, end = GHO.tonight()
    global_constraints = [MoonSeparationConstraint(min=30*u.deg), AirmassConstraint(max = 3, boolean_constraint = False), AtNightConstraint.twilight_astronomical()]
    
    # Loop over the dataframe and find requested targets, add to observing block
    for i in range(0, len(ms_df)):
        if ms_df['No. of exposures'][i] != 0:
            for fltr in ms_df['Filters'][i].split(','):
                for k in range(0, int(ms_df['No. of repeats'][i])+1):
                    b = ObservingBlock.from_exposures(target[i], ms_df['Priority'][i], ms_df['Exp. Time (s)'][i]*u.second, ms_df['No. of exposures'][i], read_out, configuration = {'filter': fltr})
                    blocks.append(b)
                
    # Initialize the sequential scheduler with the constraints and transitioner
    seq_scheduler = SequentialScheduler(constraints = global_constraints, observer = GHO, transitioner = transitioner)
    
    # Initialize a Schedule object, to contain the new schedule
    sequential_schedule = Schedule(start, end)
    
    # Call the schedule with the observing blocks and schedule to schedule the blocks
    seq_scheduler(blocks, sequential_schedule)
    
    # Convert to a schedule table
    sched = seq_scheduler.schedule.to_table(show_transitions=False).to_pandas()
    sched['Filters'] = ['']*len(sched)
    
    # Clean up the table
    for p in range(0, len(sched)):
        sched['Filters'][p] = sched['configuration'][p]['filter']
    del sched['configuration']
    
    return sched
                          
    
    