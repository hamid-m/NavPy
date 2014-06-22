"""
Copyright (c) 2014 NavPy Developers. All rights reserved.
Use of this source code is governed by a BSD-style license that can be found in
LICENSE.txt
"""

# So that nosetests would be able to find test data
import os
import sys

orig_wd = os.getcwd()

def setUp():
    """
    test package setup:  change working directory to root of test package, so that 
    relative path to test data will work.
    """
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

def tearDown():
    global orig_wd
    os.chdir(orig_wd)