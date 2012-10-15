#!/usr/bin/env python

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
# A simple code for visualizing the phase variation of the non-coplanar 
# baselines term in radio interferometry.
# Author: N. Iniyan
# e-mail: iniyan@ast.uct.ac.za
# Last modified: Oct. 2012

# Version 0.2 - Thanks to Thuso Simon for showing me how to better use numpy arrays and get rid of loops.
# If you feel this code is fast, it's due to him!


import os
import sys
import numpy as np
import pylab as pl
import sphfn # shared library sphfn.so made using f2py.

# Nominal choices for weightparm and cutoffparm

weightparm = 3          # value of 'weighting exponent selector'
                        # parameter fed to sphfn

cutoffparm = 4          # could be user-specified as 3 or 4.

#===========================================
# The python version of this function for calculating the PSWF 
# along one dimension was written by Tony Willis.

def compute_pswf(image_dimension):

 # set up gridding constants
 iflag = 0
 isupp  = 2 * cutoffparm
 ialpha = weightparm

 # compute array of PSWF function values along one dimension.
 nr = int(image_dimension / 2)  + 1
 conv_arr = np.zeros((image_dimension,),np.float32)
 eta = 0.0
 j = image_dimension-1
 norm_factor = 1.0 / float(nr-1)
 for i in xrange(1,nr+1):
  eta = float(nr-i) * norm_factor
  fx, ierr = sphfn.sphfn(ialpha,isupp,iflag,eta)
  conv_arr[i-1] = fx
  if i > 1:
   conv_arr[j] = conv_arr[i-1]
   j = j - 1
 return conv_arr

#===========================================
def wtermcalc(lmdim,delta_lm,maxw,color):
 lmdim = int(lmdim);
 delta_lm = float(delta_lm);
 maxw = float(maxw);

 #Compute spheroidal function
 conv_array_m = conv_array_l = compute_pswf(lmdim);

 #make 2-D spheroidal matrix
 conv_array = np.zeros((lmdim,lmdim));
 for m in xrange(0,lmdim):
  conv_array[m,:] = conv_array_m[m] * conv_array_l;

 # Make L, M matrices of the required dimensions that contain the 
 # delta_l and delta_m values from the image centre as elements.
 L,M = np.indices((lmdim,lmdim));
 L = L * (- 1) + lmdim / 2
 L = L.astype(np.float)*delta_lm
 M = M - lmdim / 2
 M = M.astype(np.float)*delta_lm

 # Calculate (w-term phases*sphfn) and their Fourier transforms for w=0, maxw/2 and maxw.
 norm_factor = np.sqrt(1.0-L**2-M**2); 
 
 w_factor = np.exp(-2 * np.pi * 1j * ( maxw * (norm_factor-1.0) ))/norm_factor
 wph = np.angle(w_factor) * conv_array;
 ftwph = np.fft.fft2(wph);
 ftwph = np.sqrt(ftwph.real**2+ftwph.imag**2)
 ftwph = np.fft.fftshift(ftwph);

 w_factor = np.exp(-2 * np.pi * 1j * ( maxw/2.0 * (norm_factor-1.0) ))/norm_factor
 wph2 = np.angle(w_factor) * conv_array;
 ftwph2 = np.fft.fft2(wph2);
 ftwph2 = np.sqrt(ftwph2.real**2+ftwph2.imag**2)
 ftwph2 = np.fft.fftshift(ftwph2);

 #for w=0, the tapering PSWF is all there is.
 wph3 = conv_array;
 ftwph3 = np.fft.fft2(wph3)
 ftwph3 = np.sqrt(ftwph3.real**2+ftwph3.imag**2)
 ftwph3 = np.fft.fftshift(ftwph3);
  
 #Plotting...
 fig = pl.figure();
 fig.add_subplot(2,3,1,title='w=Max_w');
 pl.imshow(wph,cmap=color);
 fig.add_subplot(2,3,2,title='w=Max_w/2');
 pl.imshow(wph2,cmap=color);
 fig.add_subplot(2,3,3,title='w=0 (just the tapering PSWF)');
 pl.imshow(wph3,cmap=color);

 fig.add_subplot(2,3,4,title='Fourier Transform');
 pl.imshow(ftwph,cmap=color);
 fig.add_subplot(2,3,5,title='Fourier Transform');
 pl.imshow(ftwph2,cmap=color);
 fig.add_subplot(2,3,6,title='Fourier Transform');
 pl.imshow(ftwph3,cmap=color);

 pl.show();
 
#=============================
def main( argv ):
  
  if len(argv) == 4:
    wtermcalc(argv[1],argv[2],argv[3],'jet')
  elif len(argv) == 5:
    wtermcalc(argv[1],argv[2],argv[3],argv[4])
  else:
    print """Usage: python wterm.py <image-dim> <delta_lm> <w-planes> [cmap (def='jet')]
	       Eg.: python wterm.py 1024 0.0007 128 (or)
                    python wterm.py 1024 0.0007 128 gray 
    """
  return
#=============================
if __name__ == "__main__":
  main(sys.argv)
