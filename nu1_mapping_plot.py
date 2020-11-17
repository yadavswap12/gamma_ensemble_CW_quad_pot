# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 16:28:17 2018

@author: yadav
"""

# code to generate nu1 mapping
import math
import cmath
import numpy
import matplotlib    #pylab is submodule in matplotlib
data = numpy.loadtxt("input_output_files/rho_2.0/mapping/mapping_output_nu1_c1=2.0_c0=1.0_10000points_iter1.0.txt",float)
x = data[:,0]
y = data[:,1]
matplotlib.pylab.plot(x,y)
#matplotlib.pyplot.plot(x,y)
#pylab.plot(x,y)    # Use this step to run in anaconda