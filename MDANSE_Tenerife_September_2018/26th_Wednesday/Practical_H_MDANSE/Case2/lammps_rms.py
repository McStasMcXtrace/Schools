"""
Script to analyze the data block extracted from a LAMMPS log file
and compute the average and RMS deviations of the data.

Usage: python lammps_rms.py input_file 
"""

import sys
import numpy as np

filename = str(sys.argv[1])

f = open(filename, 'r')

lines = f.readlines()
nlines = len(lines)

# Identify data block
for i in range(nlines):
    if lines[i][0:4] == 'Step':
        lini = i
    if lines[i][0:4] == 'Loop':
        lend = i

ndata = lend - lini - 1

properties = lines[lini].split()[1:]

# Create array to store all dat
t = np.zeros((ndata, len(properties)))


# Get data
for i in range(ndata):
    j = lini+1+i
    data = lines[j].split()[1:]
    t[i,:] = data

averages = np.sum(t,axis=0) / float(ndata)

print " "
print "Number of data points = ", ndata
print "Property, Average, R.M.S., R.M.S./average:"
print " "
for i in range(len(properties)):
    delta = (t[:,i] - averages[i])**2
    rms = np.sqrt ( np.sum(delta) / float(ndata) )
    if np.abs(averages[i]) > 0:
        fluct = rms/np.abs(averages[i])
    else:
        fluct = 0.0
    print ' ', properties[i], ' = ', averages[i], "   ,   ", rms, "   ,   ", fluct
    

