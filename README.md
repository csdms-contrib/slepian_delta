# slepian_delta
Spatiospectral localization of data from the GRACE/GRACE-FO missions using Slepian functions as originally developed by
<a href="https://polarice.geo.arizona.edu/">C. Harig</a> &amp; <a href="http://www.frederik.net">F. J. Simons</a>.
This package is hosted by CSDMS at
http://csdms.colorado.edu/wiki/Model:SLEPIAN_Delta





## Citation information
Please cite our work as appropriate if you find this package useful.

**To cite the method:**<br>
Harig, Christopher and Frederik J. Simons. 
Mapping Greeenland's mass loss in space and time.
<i>Proc. Natl. Acad. Sc.</i>, 109(49), 19934-19937, 2012.
<a href="http://dx.doi.org/10.1073/pnas.1206785109">doi:10.1073/pnas.1206785109</a>

**To cite the code itself:**<br>
Persistent DOI of the latest release: [![DOI](https://zenodo.org/badge/6548/csdms-contrib/slepian_delta.svg)](https://zenodo.org/badge/latestdoi/6548/csdms-contrib/slepian_delta)

**To cite the application to Greenland:**<br>
Harig, Christopher and Frederik J. Simons.  
Ice mass loss in Greenland, the Gulf of Alaska, and the Canadian Archipelago: Seasonal cycles and decadal trends.  
<i>Geophys. Res. Let.</i>, 43 (7), 3150-3159, 2016.
<a href="http://dx.doi.org/10.1002/2016GL067759">http://dx.doi.org/10.1002/2016GL067759</a>

Harig, Christopher and Frederik J. Simons.  
Mapping Greeenland's mass loss in space and time.  
<i>Proc. Natl. Acad. Sc.</i>, 109(49), 19934-19937, 2012.
<a href="http://dx.doi.org/10.1073/pnas.1206785109">doi:10.1073/pnas.1206785109</a>


**To cite the application to Antarctica:**<br>
Harig, Christopher and Frederik J. Simons.  
Accelerated West Antarctic ice mass loss continues to outpace East Antarctica gains. 
<i>Earth Planet. Sci. Let.</i>, 415, 134-141, 2015.
<a href="http://dx.doi.org/10.1016/j.epsl.2015.01.029">http://dx.doi.org/10.1016/j.epsl.2015.01.029</a>

Latest release hosted on [![View slepian_delta on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/81106-slepian_delta)



## Where to start? Slepian localization and GRACE: a recipe

Here is a general outline to study regions and signals in GRACE data.

1. Set up your data and directory structure, and then with GRACE2PLMT 
read in the GRACE data files into a matrix for use in Matlab. This 
function does corrections for C2,0 from the SLR values of Cheng and 
Tapley, [2004] and degree 1 values from Swenson, [2008].

2. Decide on your choice of basis, depending on your region of interest 
and the bandwidth you want. Using GRACE2SLEPT, project the results of 
GRACE2PLMT into your chosen basis. We recommend that your basis is 
chosen based on a set of synthetic experiments which estimate the 
leakage/recovery tradeoffs.

3. Next run SLEPT2RESID to fit a choice of functions (e.g. lines, 
quadratics, etc.) to the Slepian coefficients. Note: if you want to 
remove a model of GIA then you should do that before this step, using 
a function like CORRECT4GIA.

4. If you want to examine the total mass trend in a region, this 
information is calculated and returned as a result from SLEPT2RESID. 
To summarize, each Slepian coefficient up to the Shannon truncation is 
multiplied by the integral (found using INTEGRATEBASIS) of the corresponding 
function over the region. This total data is then fit with TIMESERIESFIT 
which allows fitting with a custom data variance.

5. If you want the total map of mass change, multiply the difference 
between estimated signal coefficients by the corresponding Slepian 
function, and add them up. Remember the appropriate units. This sum can 
then be expanded to space using PLM2XYZ By limiting the set of coefficients 
you use, you can instead make this for any date span, such as for several 
years or a single year.



Note: The demos above and the package in general expects a couple of
environmental variables to be set (e.g. your storage location:
$IFILES) and a directory tree for data storage to already be created
(i.e. subdirectories in $IFILES). The programs will initially display
errors instructing the user to set up these shell variables and
folders. Also, if you have an open pool of Matlab workers, a few
functions such as KERNELCP will take advantage of the parallel
computing resources and run much faster than otherwise.




## Function List:

### Administrative
| | |
| --- | --- |
| LICENSE.txt | The license file |
| README.md | The masthead for the package |
| README_Delta.md | This file! |

### Geographic
| | |
| --- | --- |
| baffing.m | An outline for a region encompassing the glaciers around Baffin Island |
| ellesmereg.m | An outline for a region encompassing the glaciers around Ellesmere Island |
| getglaciers.m | Return the XY coordinates for glaciers from the Randolph Glacier Inventory |

### Analysis
| | |
| --- | --- |
| Clmlmp2Cab.m | Given a spectral covariance matrix, turns it into a Slepian covariance matrix |
| Clmlmp2Crrp.m | Given a spectral covariance matrix, evaluates it in space |
| coarsenmask.m | Given a spatial mask, coarsen its resolution by several methods |
| cov2plm.m | Given a spectral covariance matrix, generates spherical harmonics realizations |
| geopotential.m | Plots geopotential fields for Earth, Moon and Mars |
| getglaciers.m | Get XY coordinates for glaciers in RGI database |
| gldas2TWSt.m | Read the GLDAS products and calculate Terrestrial Water Storage (TWS) fields |
| grace2plmt.m | Turns monthly GRACE data files into a single matrix for time-dependent analysis |
| grace2slept.m | Transform the result of GRACE2PLMT into a Slepian basis
| grs.m | Computes parameters for a certain geodetic reference system
| hash.m | Makes a hash (NOTE: This is one of the few files not written by us! It was made by Michael Kleder) |
| integratebasis.m | Integrates Slepian eigenfunctions given as spherical harmonics expansions |
| lovenums.m | Returns elastic Love numbers for a certain Earth model |
| maskfromgmt.m | Given a grid and polgon, generate a gridded mask using Generic Mapping Tools |
| net2mat.m | Saves specified variables from netcdf files as MAT files |
| netvarread.m | Read a netcef file and get the variable list in the file |
| periodfit.m | Find and fit periodic cycles through a data set
| plm2avg.m | Integrates and averages spherical harmonic expansions 
| plm2pot.m | Reads in and scales geopotential coeffficients
| plmresid2cov.m | Turns GRACE residual time series into a spherical-harmonic spectral covariance matrix 
| plmt2diff.m | Turns monthly GRACE data matrix into a month-to-month difference map
| plmt2resid.m | Turns monthly GRACE data matrix into residuals after fitting analysis in the spherical harmonic basis
| POMME4.m | Reads in and scales geomagnetic coeffficients of the POMME-4 model
| readgldas.m | Read in the GLDAS products and save them in Matlab friendly formats |
| resid2plot.m | Plots GRACE residual time series
| slepresid2cov.m | Turns GRACE residual time series into a Slepian covariance matrix 
| slept2resid.m | Turns monthly GRACE data matrix into residuals after fitting in the Slepian basis 
| timeseriesfit.m | Fits polynomial functions to time series with an F-test criterion 
| topography.m | Plots topography fields for Earth, Moon and Mars
| TWSt2slept.m | Project the TWS fields you previously made into a Slepian basis |






## Other important stuff

This software is distributed under the GNU Public License v2, which can be
found at http://geoweb.princeton.edu/people/simons/license.html  and also
copied below.

Copyright (C) 2014. Developer can be contacted by email. 

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version. 

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
details. 

You can receive a copy of the GNU General Public License by writing to the
Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA. 
