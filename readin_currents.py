#! usr/bin/python2
import os
import numpy as np

path = '.'
omegar = np.array([])
hr = np.array([])
hrstd = np.array([])

# for num, itemName in enumerate(os.listdir(path)):
for itemName in os.listdir(path):
    if itemName.startswith("diatomic"):
        with open(itemName,'r') as f:
          for line in f:
              a = line.split()     
              for j, ele in enumerate(a):
                  if ele=='omega_r':
                      omegar = np.append(omegar, float(a[j+2][:-1]))
                  if ele=='JR':
                      # print float(a[j+2][:-1]) # the last character to be left
                      # out is comma
                      hr = np.append(hr, float(a[j+2][:-1]))
                      hrstd = np.append(hrstd, float(a[j+5]))

filename = 'omegar_white_hr-1.txt'
np.savetxt(filename, np.c_[omegar, hr, hrstd])
