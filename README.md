# 2D-Fourier-imaging
This is the code package for magnetic Fourier imaging with NV ensemble in 2D case. It entails the uniform analysis of the reconstruction method and the non-uniform case as well. This code is wirtrn by Ke-Hang Zhu. And the idea is under the discussion with Ruolan Xue, Mark Ku and Professor Amir Yacoby in Harvard.

The code for generating 2D geometry lattice is taken from Matlab's reservatory: geom2d. Here I borrow 2 parts from that package:triangleGrid.m and hexagonalGrid.m. 
In order to perform non-uniform Fourier transformation, we nned multi-variable Lomb-Scargle algorithm, this part is transferred into R studio and the part for manipulating R code is included in the manual.
