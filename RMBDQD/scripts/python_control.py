#!/usr/bin/python3.6
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import subprocess

np = 4
main_f90 = 'tp'
input_file = 'poissonModal.dat'
mpirun = 'mpirun -n ' + str(np) +' ' + ' '+  main_f90  + ' '  + input_file

subprocess.run(mpirun, shell=True)
print('Done')
