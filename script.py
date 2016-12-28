#!/usr/bin/python

import subprocess
import sys

# Set parameter
L = 1.
Nx = 200
T = 1.
Nt = 1000
alpha = 1.
theta = 0.5

params = "%s %s %s %s %s %s" % (L,Nx,T,Nt,alpha,theta)
subprocess.call('./main %s' % params,shell=True)
