# -*- coding: utf-8 -*-
"""
Created on Thu Nov 15 16:32:54 2018

@author: yadav
"""

# code to generate nu2 mapping
import math
import cmath
import numpy
import matplotlib    #pylab is submodule in matplotlib
data = numpy.loadtxt("input_output_files/rho_2.0/mapping/mapping_output_nu2_c1=1_c0=0.5_10000points.txt",float)
x = data[:,0]
y = data[:,1]
matplotlib.pylab.plot(x,y)
#matplotlib.pyplot.plot(x,y)
#pylab.plot(x,y)    # Use this step to run in anaconda