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
# Modified from wterm.py to illustrate zero-padding and frequency interpolation.
# Author: N. Iniyan
# e-mail: iniyan@ast.uct.ac.za
# Last modified: Oct. 2012


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
def wtermcalc(lmdim,delta_lm,maxw,zeropad,color):
 lmdim = int(lmdim);
 delta_lm = float(delta_lm);
 maxw = float(maxw);
 zeropad = int(zeropad);

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

 if maxw!=0: 
  w_factor = np.exp(-2 * np.pi * 1j * ( maxw * (norm_factor-1.0) ))/norm_factor
  wph = np.angle(w_factor) * conv_array;
 else:
  wph = conv_array; #for w=0, the tapering PSWF is all there is.

 # zero-padding input
 padlm = zeropad * lmdim;
 padarr = np.zeros((padlm,padlm));
 padarr[padlm/2-lmdim/2:padlm/2+lmdim/2,padlm/2-lmdim/2:padlm/2+lmdim/2] = wph;

 ftwph = np.fft.fft2(wph);
 ftwph = np.sqrt(ftwph.real**2+ftwph.imag**2)
 ftwph = np.fft.fftshift(ftwph);
  
 ftpad = np.fft.fft2(padarr);
 ftpad = np.sqrt(ftwph.real**2+ftwph.imag**2)
 # No need to shift here because of the zero-padding???

 #Plotting...
 fig = pl.figure();
 fig.add_subplot(2,2,1,title='Input signal');
 pl.imshow(wph,cmap=color);
 fig.add_subplot(2,2,2,title='Zero-padded input');
 pl.imshow(padarr,cmap=color);

 pltft1=fig.add_subplot(2,2,3,title='Fourier Transform');
 pl.imshow(ftwph,cmap=color);
 pltft2=fig.add_subplot(2,2,4,title='Fourier Transform');
 pl.imshow(ftpad,cmap=color);

 pl.show();
 
#=============================
def main( argv ):
  
  if len(argv) == 5:
    wtermcalc(argv[1],argv[2],argv[3],argv[4],'jet')
  elif len(argv) == 6:
    wtermcalc(argv[1],argv[2],argv[3],argv[4],argv[5])
  else:
    print """Usage: python zeropad.py <image-dim> <delta_lm> <w-planes> <oversampling factor> [cmap (def='jet')]
	       Eg.: python zeropad.py 1024 0.0007 128 2 (or)
                    python zeropad.py 1024 0.0007 128 2 gray 
    """
  return
#=============================
if __name__ == "__main__":
  main(sys.argv)
