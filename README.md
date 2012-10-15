wtermcalc
=========

Python code to visualize the phase variation of the non-coplanar baselines term in radio interferometry.

This code computes the phase variation of the non-coplanar baselines term (or the w-term) for w-values of w, w/2 and 0
and plots all three along with their corresponding Fourier transforms.

The w-term image is tapered with a prolate spheroidal wave function (PSWF) before its F.T is taken. For w=0, the PSWF
is all there is.

For more information, refer Frater and Doherty 1980 and Cornwell et al. 2008.

Libraries used - Matplotlib, Numpy.

Included files:
===============

sphfn.so - Shared library that computes the PSWF; generated using f2py.

wterm.py - Computes and plots the phase variation in the w-term as a function of distance from the phase centre.

zeropad.py - Modified version of wterm.py that illustrates the effects of zero-padding in the spatial domain.