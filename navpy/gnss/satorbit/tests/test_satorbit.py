"""
Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""
import navpy.gnss.satorbit as satorbit
import unittest
import numpy as np

assert_almost_equal = np.testing.assert_almost_equal

class TestSatOrbitClass(unittest.TestCase):
    def test_sat_pos_computation(self):
        ephem_file = 'test_data/brdc1680.13n'
        TOW = 87300.00
        
        # Reference File
        ref_fname = 'test_data/igs_2013_06_17_00_15_00.sp3'
        ref_file  = open(ref_fname)
        
        ref_sat_pos = np.nan*np.zeros((32,3))
        for line in ref_file:
            line = line.strip()
            if (line.startswith('//') or line.startswith('#')):
                continue
            if (line.startswith('PG')):
                prn = int(line[2:4])-1
                data = line[4:].strip().split()
                ref_sat_pos[prn,:] = [float(data[0]), float(data[1]), float(data[2])]

        # Calculate
        gps_ephem = satorbit.ephem_class()
        gps_ephem.read_RINEX(ephem_file)

        sv_t = np.vstack((range(0,32),TOW*np.ones(32))).T

        # Compare with igs_2013_06_17_00_15.sp3
        x,y,z = satorbit.compute_sat_pos(gps_ephem,sv_t) 
        sat_pos = np.vstack((x,y,z)).T/1000.0
        
        for i in xrange(0,32):
            for e1, e2 in zip(sat_pos[i,:], ref_sat_pos[i,:]):
                if(np.isnan(e2)):
                    continue
                self.assertAlmostEqual(e1,e2,places=0)
    
    def test_clk_bias_computation(self):
        ephem_file = 'test_data/brdc1680.13n'
        TOW = 87300.00
        
        # Reference File
        ref_fname = 'test_data/igs_clk_2013_06_17_00_15_00.clk'
        ref_file  = open(ref_fname)
        
        ref_clk_bias = np.nan*np.zeros(32)
        for line in ref_file:
            line = line.strip()
            if (line.startswith('//') or line.startswith('#')):
                continue
            if (line.startswith('AS G')):
                prn = int(line[4:6])-1
                data = line[6:].strip().split()
                ref_clk_bias[prn] = float(data[7])

        # Calculate
        gps_ephem = satorbit.ephem_class()
        gps_ephem.read_RINEX(ephem_file)

        sv_t = np.vstack((range(0,32),TOW*np.ones(32))).T

        # Compare with igs_2013_06_17_00_15.sp3
        clk = satorbit.compute_sat_clk_bias(gps_ephem,sv_t) 
        
        for e1, e2 in zip(clk, ref_clk_bias):
            if (np.isnan(e2)):
                continue
            self.assertAlmostEqual(e1,e2,places=1)

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSatOrbitClass)
    unittest.TextTestRunner(verbosity=2).run(suite)    
