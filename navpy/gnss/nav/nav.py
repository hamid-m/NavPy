"""
Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""

import numpy as np
import scipy.linalg as la
import navpy as navpy
from ..satorbit import satorbit

def code_phase_LS(raw_meas, gps_ephem, lat=0.0, lon=0.0, alt=0.0, rxclk=0.0):
    """
    Calculate code phase least square solution
    
    Parameters
    ----------
    raw_meas: prn_class object, containing the pseudorange measurements
    gps_ephem: ephem_class object, containing the satellite ephemeris
    
    lat: (optional, degrees), initial guess of latitude
    lon: (optional, degrees), initial guess of longitude
    alt: (optional, m), initial guess of altitude
    rxclk: (optional, seconds), initial guess of Receiver clock bias
    
    Returns
    -------
    lat: (deg) Latitude
    lon: (deg) Longitude
    alt: (m) Altitude
    rxclk: (sec) Receiver Clock Bias
    
    """
    # If there are less than 3 satellites, don't do anything
    SV_avbl = np.nonzero(raw_meas.is_dataValid(range(32)))[0]
    
    if(len(SV_avbl) < 3):
        return lat, lon, alt, rxclk
    
    # Begin computing position using Least Squares
    t_tx = raw_meas.get_TOW()*np.ones(len(SV_avbl))
    delta_time = np.zeros(len(SV_avbl))
    
    # Iterate because time at transmission is not known ...
    for k in xrange(5):
        ecef = navpy.lla2ecef(lat,lon,alt)
        t_tx = raw_meas.get_TOW()*np.ones(len(SV_avbl)) - delta_time
        
        # Build satellite information
        clk = satfn.compute_sat_clk_bias(gps_ephem,np.vstack((SV_avbl,t_tx)).T)
        x,y,z = satfn.compute_sat_pos(gps_ephem,np.vstack((SV_avbl,t_tx)).T)
        
        # Build Geometry Matrix H = [LOS 1]
        rho = np.sqrt((x-ecef[0])**2 + (y-ecef[1])**2 + (z-ecef[2])**2)  #LOS Magnitude
        H = np.vstack( ( -(x-ecef[0])/rho, -(y-ecef[1])/rho, -(z-ecef[2])/rho, np.ones(len(SV_avbl)) ) ).T
        
        # Calculate estimated pseudorange
        PR_hat = np.sqrt((x-ecef[0])**2 + (y-ecef[1])**2 + (z-ecef[2])**2) + rxclk*np.ones(len(SV_avbl))
        
        # Innovation: Difference between the measurement (corrected for satellite clock bias)
        #             and PR_hat
        dy = (raw_meas.get_pseudorange(SV_avbl) + clk*2.99792458e8) - PR_hat
        
        # Measurement Covariance
        RR = raw_meas.get_PR_cov(SV_avbl)

        # Least Square Solution
        dx = la.inv(H.T.dot(RR.dot(H))).dot(H.T.dot(RR)).dot(dy)
        # ... Correction
        ecef += dx[0:3]
        rxclk += dx[3]
        # ... Reclaculate Lat, Lon, Alt
        lat, lon, alt = navpy.ecef2lla(ecef)
        
        # Recalculate delta_time to correct for time at transmission
        x,y,z = compute_sat_pos(gps_ephem,np.vstack((SV_avbl,t_tx)).T)
        delta_time = (np.sqrt((x-ecef[0])**2 + (y-ecef[1])**2 + (z-ecef[2])**2))/2.99792458e8
        # ... Let's go to the next iteration ...
     
    return lat, lon, alt, rxclk