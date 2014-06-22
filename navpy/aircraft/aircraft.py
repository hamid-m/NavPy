"""
Aircraft Class and Functions

Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""

import numpy as np
import navpy.gnss as gnss

class ac_class:
    """
    Aircraft Class
    """
    
    def __init__(self,fname=''):
        self._file = fname
        self.imu = imu_class()
        self.gps = gnss.rx_class()
        self.airdata = airdata_class()
        
        
class imu_class:
    """
    IMU Class
    """
    def __init__(self):
        self._p = np.nan
        self._q = np.nan
        self._r = np.nan
        self._ax = np.nan
        self._ay = np.nan
        self._az = np.nan
        self._hx = np.nan
        self._hy = np.nan
        self._hz = np.nan
    
    def set_imu(self,accel,gyro,mag):
        self._ax, self._ay, self._az = accel
        self._p, self._q, self._r = gyro
        self._hx, self._hy, self._hz = mag
    
    def get_accel(self):
        return np.array([self._ax, self._ay, self._az])
    
    def get_gyro(self):
        return np.array([self._p, self._q, self._r])
        
    def get_mag(self):
        return np.array([self._hx, self._hy, self._hz])

class airdata_class:
    """
    Air Data Class
    """
    def __init__(self):
        self._ias = np.nan
        self._alpha = np.nan
        self._beta = np.nan
    
    def set_imu(self,v,alpha,beta):
        self._ias, self._alpha, self._beta = v, alpha, beta
    
    def get_ias(self):
        return self._ias
    
    def get_alpha(self):
        return self._alpha
        
    def get_beta(self):
        return self._beta
    