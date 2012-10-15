wtermcalc
=========

Python code to visualize the phase variation of the non-coplanar baselines term in radio interferometry.

This code computes the phase variation of the non-coplanar baselines term (or the w-term) for w-values of w, w/2 and 0
and plots all three along with their corresponding Fourier transforms.

The w-term image is tapered with a prolate spheroidal wave function (PSWF) before its F.T is taken. For w=0, the PSWF
is all there is.

The source includes the shared library sphfn.so generated using f2py.

For more information, refer Frater and Doherty 1980 and Cornwell et al. 2008.

Libraries used - Matplotlib, Numpy.