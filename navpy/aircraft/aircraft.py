"""
Aircraft Class and Functions

Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""

import numpy as np
import navpy.gnss as gnss
import navpy.utils as _utils

class ac_class:
    """
    Aircraft Class
    """
    
    def __init__(self,log_file='',model_file=''):
        self._logfile = log_file
        self._model = model_file
        self.h = 0.0
        self.V = 0.0
        self.alpha = 0.0
        self.beta = 0.0
        self.p = 0.0
        self.q = 0.0
        self.r = 0.0
        self.phi = 0.0
        self.theta = 0.0
        self.psi = 0.0
        self.ax = 0.0
        self.ay = 0.0
        self.az = 0.0
        self.lat = 0.0
        self.lon = 0.0
        self.alt = 0.0

class fcs_class:
    """
    Flight Control System Class
    """        
    
    def __init__(self):
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
    
    # ====================== API's ==========================
    def set_imu(self,accel,gyro,mag):
        accel,N1 = _utils.input_check_Nx3(accel)
        gyro,N2 = _utils.input_check_Nx3(gyro)
        mag,N3 = _utils.input_check_Nx3(mag)
        
        if( (N1 != N2) or (N1 != N3) ):
            raise ValueError("Wrong Input Dimension")
        
        if(N1 > 1):
            self._ax = accel.reshape(N1,3)[:,0]
            self._ay = accel.reshape(N1,3)[:,1]
            self._az = accel.reshape(N1,3)[:,2]
            self._p = gyro.reshape(N2,3)[:,0]
            self._q = gyro.reshape(N2,3)[:,1]
            self._r = gyro.reshape(N2,3)[:,2]
            self._hx = mag.reshape(N3,3)[:,0]
            self._hy = mag.reshape(N3,3)[:,1]
            self._hz = mag.reshape(N3,3)[:,2]
        else:
            self._ax, self._ay, self._az = accel
            self._p, self._q, self._r = gyro
            self._hx, self._hy, self._hz = mag
            
    def get_accel(self):
        return np.array([self._ax, self._ay, self._az]).T
    
    def get_gyro(self):
        return np.array([self._p, self._q, self._r]).T
        
    def get_mag(self):
        return np.array([self._hx, self._hy, self._hz]).T

class airdata_class:
    """
    Air Data Class
    """
    def __init__(self):
        self._ias = np.nan
        self._alpha = np.nan
        self._beta = np.nan
    
    def set_airdata(self,v,alpha,beta):
        v,N1 = _utils.input_check_Nx1(v)
        alpha,N2 = _utils.input_check_Nx1(alpha)
        beta,N3 = _utils.input_check_Nx1(beta)
        
        if( (N1 != N2) or (N1 != N3) ):
            raise ValueError("Wrong Input Dimension")
            
        self._ias = v
        self._alpha = alpha
        self._beta = beta
        
    def get_ias(self):
        return self._ias
    
    def get_alpha(self):
        return self._alpha
        
    def get_beta(self):
        return self._beta
    