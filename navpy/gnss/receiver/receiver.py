"""
Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""

import numpy as np

class rx_class:
    """
    GNSS Receiver Class
    This class contains the GNSS Solution for the receiver
    and the raw measurement object (from prn_class)
    
    Adhika Lie, 06/21/2014
    """
    def __init__(self):
        self.TOW = 0
        self.lat = 0
        self.lon = 0
        self.alt = 0
        self.clkbias = 0
        self.sig_N = 0
        self.sig_E = 0
        self.sig_D = 0
        
        self.rawdata = prn_class()
        
        self.INIT = False
        
class prn_class:
    """
    PRN Class
    This class stores raw measurement data to each GPS Satellite (PRN)
    Fields inside this class are set and accessed via APIs
    Available methods / APIs and their prototypes are:
    1. Pseudorange:
        set_pseudorange(PR,PR_std,sv)
        get_pseudorange(sv)
        get_PR_cov(sv)
    2. Accumulated Delta Range or Carrier Phase:
        set_carrierphase(ADR,ADR_std,sv)
        get_carrierphase(sv)
        get_ADR_cov(sv)
    3. Doppler:
        set_doppler(DR,sv)
        * NEED TO ADD GET *
    4. CNo:
        set_CNo(CNo,sv)
        * NEED TO ADD GET *
    5. Flags (not listed here)
    
    Adhika Lie, 06/21/2014
    """
    def __init__(self):
        self._PR = np.nan*np.ones(32)    # Pseudorange
        self._ADR = np.nan*np.ones(32)   # Carrier Phase
        self._DR = np.nan*np.ones(32)    # Doppler
        self._PR_std = np.nan*np.ones(32)
        self._ADR_std = np.nan*np.ones(32)
        self._locktime = np.nan*np.ones(32)
        self._CNo = np.nan*np.ones(32)
        self._dataValid = [False for i in xrange(32)]
        
    ######################################################################
    #                            prn_class APIs
    ######################################################################
    
    # 1. PSEUDORANGE
    # ====================================================================
    def set_pseudorange(self,PR,PR_std,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        PR,N2 = _utils.input_check_Nx1(PR)
        PR_std,N3 = _utils.input_check_Nx1(PR_std)

        if(np.any(sv>=32)):
            raise TypeError('sv > 32')

        if( (N1!=N2) or (N1!=N3) ):
            raise TypeError('Incompatible size')
            
        if(N1==1):
            sv=[sv]
            PR = [PR]
            PR_std = [PR_std]

        for i in xrange(N1):
            self._PR[sv[i]] = PR[i]
            self._PR_std[sv[i]] = PR_std[i]
    
    def get_pseudorange(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]
        return np.array([self._PR[prn] for prn in sv])
    
    def get_PR_cov(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]
        
        return np.diag([1./self._PR_std[prn]**2 for prn in sv])
    
    # 2. CARRIER PHASE / ACCUMULATED DOPPLER RANGE 
    # ====================================================================  
    def set_carrierphase(self,ADR,ADR_std,locktime,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        ADR,N2 = _utils.input_check_Nx1(ADR)
        ADR_std,N3 = _utils.input_check_Nx1(ADR_std)
        locktime,N4 = _utils.input_check_Nx1(locktime)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        if( (N1!=N2) or (N1!=N3) or (N1!=N4)):
            raise TypeError('Incompatible size')
        
        if(N1==1):
            sv=[sv]
            ADR = [ADR]
            ADR_std = [ADR_std]
            locktime = [locktime]
        
        for i in xrange(N1):
            self._ADR[sv[i]] = ADR[i]
            self._ADR_std[sv[i]] = ADR_std[i]
            self._locktime[sv[i]] = locktime[i]

    def get_carrierphase(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]
        return np.array([self._ADR[prn] for prn in sv])
    
    def get_ADR_cov(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]
        
        return np.diag([1./self._ADR_std[prn]**2 for prn in sv])
    
    # 3. DOPPLER      
    # ====================================================================
    def set_doppler(self,DR,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        DR,N2 = _utils.input_check_Nx1(DR)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        if( N1!=N2 ):
            raise TypeError('Incompatible Size')
        
        if(N1==1):
            sv=[sv]
            DR = [DR]
            
        for i in xrange(N1):
            self._DR[sv[i]] = DR[i]
 
    # 4. CARRIER-TO-NOISE 
    # ====================================================================
    def set_CNo(self,CNo,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        CNo,N2 = _utils.input_check_Nx1(CNo)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        if( N1!=N2 ):
            raise TypeError('Incompatible Size')
            
        if(N1==1):
            sv=[sv]
            CNo = [CNo]
            
        for i in xrange(N1):
            self._CNo[sv[i]] = CNo[i]

    # 5. FLAGS           
    # ====================================================================
    def _set_dataValid(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')

        if(N1==1):
            sv=[sv]
            
        for i in xrange(N1):
            self._dataValid[sv[i]] = True
    
    def _reset_dataValid(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]
        for i in xrange(N1):
            self._dataValid[sv[i]] = False
    
    def check_dataValid(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]    
            
        for prn in sv:
            self._dataValid[prn] = ~( np.any(np.isnan(self._PR[prn])) or np.any(np.isnan(self._ADR[prn])) ) 
        
    def is_dataValid(self,sv):
        sv,N1 = _utils.input_check_Nx1(sv)
        if(np.any(sv>=32)):
            raise TypeError('sv > 32')
        
        if(N1==1):
            sv=[sv]
        return [self._dataValid[prn] for prn in sv]