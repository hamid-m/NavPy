import numpy as np
import navpy as _navpy
import navpy.utils as _utils
import datetime

class ephem_class:
    """
    Satellite Ephemeris Class
    This class contains the orbital parameters for GPS PRN 0 to 31
    The current method is to read from RINEX file.
    Can be expanded to more standard ephemeris file
    
    Adhika Lie, 06/21/2014
    """
    def __init__(self):
        self.toc = [[] for i in xrange(0,32)]
        self.af0 = [[] for i in xrange(0,32)]
        self.af1 = [[] for i in xrange(0,32)]
        self.af2 = [[] for i in xrange(0,32)]
        
        self.iode = [[] for i in xrange(0,32)]
        self.Crs = [[] for i in xrange(0,32)]
        self.delta_n = [[] for i in xrange(0,32)]
        self.M0 = [[] for i in xrange(0,32)]
        
        self.Cuc = [[] for i in xrange(0,32)]
        self.e = [[] for i in xrange(0,32)]
        self.Cus = [[] for i in xrange(0,32)]
        self.sqrt_a = [[] for i in xrange(0,32)]
        
        self.toe = [[] for i in xrange(0,32)]
        self.Cic = [[] for i in xrange(0,32)]
        self.Loa = [[] for i in xrange(0,32)]
        self.Cis = [[] for i in xrange(0,32)]
        
        self.i = [[] for i in xrange(0,32)]
        self.Crc = [[] for i in xrange(0,32)]
        self.perigee = [[] for i in xrange(0,32)]
        self.ra_rate = [[] for i in xrange(0,32)]
        
        self.i_rate = [[] for i in xrange(0,32)]
        self.GPS_week = [[] for i in xrange(0,32)]
        
        self.iodc = [[] for i in xrange(0,32)]
        self.URA = [[] for i in xrange(0,32)]
        
    def read_RINEX(self,ephem_file):
        """
        Read Ephemeris File from RINEX file
        
        Usage Example:
        --------------
        ephem_file = 'examples/brdc1680.13n'
        gps_ephem = gnss.ephem_class()
        gps_ephe.read_RINEX(ephem_file)
        
        Returns:
        ---------
        None
        
        Dependencies
        -------------
        ephem_class._RINEX2array()
        
        Adhika Lie, 06/21/2014
        """
        ephem_array = _RINEX2array(ephem_file)
        ephem_array = ephem_array[ephem_array[:,0].argsort()]
        for i in xrange(0,ephem_array.shape[0]):
            prn = int(ephem_array[i,0]-1)
            
            self.toc[prn].append(ephem_array[i,19])
            self.af0[prn].append(ephem_array[i,20])
            self.af1[prn].append(ephem_array[i,21])
            self.af2[prn].append(ephem_array[i,22])
            
            self.iode[prn].append(ephem_array[i,17])
            self.Crs[prn].append(ephem_array[i,13])
            self.delta_n[prn].append(ephem_array[i,2])
            self.M0[prn].append(ephem_array[i,1])
            
            self.Cuc[prn].append(ephem_array[i,10])
            self.e[prn].append(ephem_array[i,3])
            self.Cus[prn].append(ephem_array[i,11])
            self.sqrt_a[prn].append(ephem_array[i,4])
            
            self.toe[prn].append(ephem_array[i,16])
            self.Cic[prn].append(ephem_array[i,14])
            self.Loa[prn].append(ephem_array[i,5])
            self.Cis[prn].append(ephem_array[i,15])
            
            self.i[prn].append(ephem_array[i,6])
            self.Crc[prn].append(ephem_array[i,12])
            self.perigee[prn].append(ephem_array[i,7])
            self.ra_rate[prn].append(ephem_array[i,8])
            
            self.i_rate[prn].append(ephem_array[i,9])
            self.GPS_week[prn].append(ephem_array[i,18])
            
            self.iodc[prn].append(ephem_array[i,23])

    def _RINEX2array(ephem_file):
        """
        Read Ephemeris File from RINEX file and convert it to numpy array
        """
        fid = open(ephem_file,'r')
    
        # First 8 lines are irrelevant
        while(1):
            tline = fid.readline().lstrip().rstrip(' \n')
            if(tline == 'END OF HEADER'):
                break   
    
        satNo = np.zeros(32)
    
        while(1):
            tline = fid.readline()
            if(tline==''):
                break
            #print tline
            # Parsing first line
            prn = int(tline[0:2])
            year = int(tline[2:5].replace('D','E'))
            month = int(tline[5:8].replace('D','E'))
            day = int(tline[8:11].replace('D','E'))
            hour = int(tline[11:14].replace('D','E'))
            minute = int(tline[14:17].replace('D','E'))
            sec = float(tline[17:22].replace('D','E'))
        
            Af0 = float(tline[22:19*1+22].replace('D','E'))
            Af1 = float(tline[19*1+22:19*2+22].replace('D','E'))
            Af2 = float(tline[19*2+22:19*3+22].replace('D','E'))
        
            # Parsing second line
            tline = fid.readline()
            IODE = float(tline[3:19*1+3].replace('D','E'))
            Crs = float(tline[19*1+3:19*2+3].replace('D','E'))
            delta_n = float(tline[19*2+3:19*3+3].replace('D','E'))
            M0 = float(tline[19*3+3:19*4+3].replace('D','E'))
        
            # Parsing third line
            tline = fid.readline()
            Cuc = float(tline[3:19*1+3].replace('D','E'))
            e = float(tline[19*1+3:19*2+3].replace('D','E'))
            Cus = float(tline[19*2+3:19*3+3].replace('D','E'))
            sqrt_a = float(tline[19*3+3:19*4+3].replace('D','E'))
        
            # Parsing fourth line
            tline = fid.readline()
            Toe = float(tline[3:19*1+3].replace('D','E'))
            Cic = float(tline[19*1+3:19*2+3].replace('D','E'))
            Loa = float(tline[19*2+3:19*3+3].replace('D','E'))
            Cis = float(tline[19*3+3:19*4+3].replace('D','E'))
        
            # Parsing fifth line
            tline = fid.readline()
            i = float(tline[3:19*1+3].replace('D','E'))
            Crc = float(tline[19*1+3:19*2+3].replace('D','E'))
            perigee = float(tline[19*2+3:19*3+3].replace('D','E'))
            ra_rate = float(tline[19*3+3:19*4+3].replace('D','E'))
        
            # Parsing sixth line
            tline = fid.readline()
            i_rate = float(tline[3:19*1+3].replace('D','E'))
            GPS_week = float(tline[19*3+3:19*4+3].replace('D','E'))
        
            # Parsing seventh and eighth line
            tline = fid.readline()
            IODC = float(tline[19*3+3:19*4+3].replace('D','E'))
            tline = fid.readline()
        
            # Calculate Toc
            #weekNo, Toc = calendar2gpstime(year,month,day,hour,minute,sec)
        
            date_diff = datetime.date(year,month,day) - datetime.date(1980,1,6)
            Nsecs = date_diff.total_seconds()
        
            weekNo, rem = divmod(Nsecs,604800) # There are 604800 seconds/week
            # rem is the number of seconds into the week
        
            # Add to rem, the number of seconds into the day
            dsec = hour*3600.0 + minute*60.0 + sec
            Toc = rem + dsec
        
            new_line = [prn, M0, delta_n, e, sqrt_a, Loa, i, perigee, ra_rate,\
                        i_rate, Cuc, Cus, Crc, Crs, Cic, Cis, Toe, IODE, GPS_week,\
                        Toc, Af0, Af1, Af2, IODC]
            try:
                gps_ephem = np.vstack((gps_ephem,new_line))
            except NameError:
                gps_ephem = np.array(new_line)
        
        return gps_ephem

def compute_sat_pos(ephem,sv_t):
    """ 
    Compute Satellite Position given ephemeris class
    Based on Table 8 GPS Blue Book Vol 1 (pp. 138)
    
    Parameters
    ----------
    ephem: ephem_class object 
    sv_t: Nx2 iterable object
        Column 0: Satellite PRN Number [0-31]
        Column 1: Time at which the position of PRN on Column 0 is desired
    
    Returns
    -------
    x, y, z: (m) ECEF Position of each satellite specified on Column 0 of sv_t
            at the time specified on Column 1 of sv_t
    
    Adhika Lie, 06/21/2014
    """
    sv_t, N = _utils.input_check_Nx2(sv_t)
    
    if(N==1):
        sv = [int(sv_t[0])]
        t_input = [sv_t[1]]
    else:
        sv = sv_t[:,0].astype(int)
        t_input = sv_t[:,1]
    
    # Pre-allocate
    x = np.nan*np.ones(N)
    y = np.nan*np.ones(N)
    z = np.nan*np.ones(N)

    k = 0
    for prn in sv:
        # Find closest time of ephemeris
        closest_dTOE = np.min ( np.abs(np.array(ephem.toe[prn]) - t_input[k]) )
        idx = np.nonzero( ( np.abs(np.array(ephem.toe[prn]) - t_input[k]) ) == closest_dTOE)[0]
        
        ephemeris_epoch = ephem.toe[prn][idx]
        #print ("PRN = %d, ephemeris_epoch = %f" % (prn+1, ephemeris_epoch))
        
        a = ephem.sqrt_a[prn][idx]**2       # Compute a (meters)
        n0 = np.sqrt(_navpy.wgs84.GM_GPS /a**3)
        
        t = t_input[k] - ephemeris_epoch
        
        n = n0 + ephem.delta_n[prn][idx]    # Corrected mean motion
        M = ephem.M0[prn][idx] + n*t       # Mean anomaly
        
        # Load eccentricity
        e = ephem.e[prn][idx]
        Ek = _calcEA(M,e)
        
        # Compute true anomaly
        nu = np.arctan2( np.sqrt(1-e**2) * np.sin(Ek) / (1-e*np.cos(Ek)),
                         (np.cos(Ek)-e) / (1-e*np.cos(Ek)) )
        
        #Eccentric Anomaly
        Ek = np.arccos((e + np.cos(nu))/(1+e*np.cos(nu)))

        # Argument of perigee
        omega = ephem.perigee[prn][idx] 
        # Compute the argument of latitude
        phi = nu + omega
        
        # Compute sinusoidal terms for the ephemeris data
        del_u = ephem.Cus[prn][idx] * np.sin(2*phi) + ephem.Cuc[prn][idx] * np.cos(2*phi)
        del_r = ephem.Crs[prn][idx] * np.sin(2*phi) + ephem.Crc[prn][idx] * np.cos(2*phi)
        del_i = ephem.Cis[prn][idx] * np.sin(2*phi) + ephem.Cic[prn][idx] * np.cos(2*phi)
        
        u = phi+del_u    # Argument of latitude correction
        r = a*(1-e*np.cos(Ek)) + del_r # Radius correction
        i = ephem.i[prn][idx] + del_i + ephem.i_rate[prn][idx]*t
        
        x0 = r*np.cos(u)   # Satellite x position in orbital plane
        y0 = r*np.sin(u)   # Satellite y position in orbital plane
        
        #print ("PRN = %d, x0 = %f, y0 = %f" % (prn+1, x0, y0))
        
        # Corrected longitude of ascending node for node rate and
        # Earth's rotation
        node = ephem.Loa[prn][idx] + (ephem.ra_rate[prn][idx] - _navpy.wgs84.omega_E_GPS)*t -\
                _navpy.wgs84.omega_E_GPS*ephem.toe[prn][idx]
                
        # X, Y, Z ECEF location of the satellite
        
        x[k] = x0*np.cos(node)-y0*np.cos(i)*np.sin(node)
        y[k] = x0*np.sin(node)+y0*np.cos(i)*np.cos(node)
        z[k] = y0*np.sin(i)
            
        #print ("PRN = %d,x = %f" % (prn+1, x[k,prn]))
        #print ("PRN = %d,y = %f" % (prn+1, y[k,prn]))
        #print ("PRN = %d,z = %f" % (prn+1, z[k,prn]))
        k += 1
            
    if(N==1):
        x = x.reshape(len(sv))
        y = y.reshape(len(sv))
        z = z.reshape(len(sv))
    return x, y, z

def compute_sat_clk_bias(ephem,sv_t):
    """ 
    Compute Satellite Clock Correction Term
    GPS Blue Book Vol 1 pp. 133
    
    Parameters
    ----------
    ephem: ephem_class object 
    sv_t: Nx2 iterable object
        Column 0: Satellite PRN Number [0-31]
        Column 1: Time at which the position of PRN on Column 0 is desired
    
    Returns
    -------
    clkbias: (sec) Satellite clock bias estimate  of each satellite specified
            on Column 0 of sv_t, at the time specified on Column 1 of sv_t
    
    Adhika Lie, 06/21/2014
    """
    # Constant for Relativistic Time Correction
    F = -2*np.sqrt(_navpy.wgs84.GM_GPS)/2.99792458e8**2
    
    sv_t, N = _utils.input_check_Nx2(sv_t)
    
    if(N==1):
        sv = [int(sv_t[0])]
        t_input = [sv_t[1]]
    else:
        sv = sv_t[:,0].astype(int)
        t_input = sv_t[:,1]
    
    # Pre-allocate
    clk_bias = np.nan*np.ones(N)
    k = 0
    for prn in sv:
        # Find closest time of ephemeris
        closest_dTOE = np.min ( np.abs(np.array(ephem.toe[prn]) - t_input[k]) )
        idx = np.nonzero( ( np.abs(np.array(ephem.toe[prn]) - t_input[k]) ) == closest_dTOE)[0]

        ephemeris_epoch = ephem.toe[prn][idx]
        
        n0 = np.sqrt(_navpy.wgs84.GM_GPS /ephem.sqrt_a[prn][idx]**6)
        
        t = t_input[k] - ephemeris_epoch
        
        n = n0 + ephem.delta_n[prn][idx]    # Corrected mean motion
        M = ephem.M0[prn][idx] + n*t       # Mean anomaly
        
        # Load eccentricity
        e = ephem.e[prn][idx]
        Ek = _calcEA(M,e)
        
        dtr = F*ephem.e[prn][idx]*ephem.sqrt_a[prn][idx]*np.sin(Ek)
        
        af0 = ephem.af0[prn][idx]
        af1 = ephem.af1[prn][idx]
        af2 = ephem.af2[prn][idx]
        
        clk_bias[k] = af0 + af1*t + af2*t**2 + dtr

        k+=1
        
    if(N==1):
        clk_bias = clk_bias.reshape(len(sv))
    
    return clk_bias
    

def _calcEA(M,e):
    """
    Calculate Mean Anomaly, used for satellite position computation
    Parameters
    ----------
    M: Mean anomaly
    e: eccentricity
    
    Return
    ------
    E: Eccentric Anomaly
    
    Adhika Lie, 09/21/2014
    """
    Etemp = M
    ratio = 1
    while(abs(ratio)>1e-8):
        f_E = Etemp - e*np.sin(Etemp) - M
        f_Eprime = 1-e*np.cos(Etemp)
        ratio = f_E/f_Eprime
        
        if(abs(ratio)>1e-8):
            Etemp = Etemp - ratio
        else:
            E = Etemp
    return E

def calendar2gpstime(year,month,day,hour,minute,sec):
    """
    Calculate GPS Time (Week and TOW) from a Julian calendar epoch
    
    Parameters
    ----------
    year: Julian calendar year
    month: Month 1 to 12
    day: Date
    hour: Time
    minute: Time
    sec: Time
    
    Return
    ------
    weekNo, TOW: GPS Week number and time of week in seconds
    
    Adhika Lie, 06/21/2014
    """
    date_diff = datetime.date(year,month,day) - datetime.date(1980,1,6)
    Nsecs = date_diff.total_seconds()
    
    weekNo, rem = divmod(Nsecs,604800) # There are 604800 seconds/week
    # rem is the number of seconds into the week
    
    # Add to rem, the number of seconds into the day
    dsec = hour*3600.0 + minute*60.0 + sec
    TOW = rem + dsec
    
    return weekNo,TOW