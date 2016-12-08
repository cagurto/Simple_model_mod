
# coding: utf-8

# In[ ]:

import matplotlib.pylab as plb
import radmc3dPy
import numpy as np
import os
from radmc3dPy.natconst import *
import sys
# Import the simple models:
import SimpleDisk

# go to work folder, create model folders
# cd /alma/home/cagurto/radmc3d/radmc-3d/version_0.40/SIMPLE/SimpleDisk
cd /alma/home/cagurto/radmc3d/radmc-3d/version_0.40/SIMPLE/SimpleDisk/SimpleDisk-mod
mkdir run1_noExt
mkdir run1_Ext
mkdir run1_Ext_noStar


# In[ ]:

# Create an envelope model without external radiation:
cd run1_noExt
cp ../problem_params.inp .
cp ../dustkappa_osshenn_thinextra.inp .
SimpleDisk.simpleRadmc3Dmodel(disk=False)
# run the model
os.system('radmc3d mctherm setthreads 40')
os.system('radmc3d sed incl 67 phi 30 setthreads 40')


# In[ ]:

# Read the model density and temperature
data = radmc3dPy.analyze.readData(dtemp=True,ddens=True,binary=False)

# Set plotting ranges:
Tmin = 0.
Tmax = 50.

# Plot the density
PltDens = plb.contourf(data.grid.x/au, np.pi/2.-data.grid.y, np.log10(data.rhodust[:,:,0,0].T), 60)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbD = plb.colorbar(PltDens)
cbD.set_label(r'$\log_{10}{\rho}$')
plb.show()

# Plot the temperature
Tdust = data.dusttemp[:,:,0,0].T.clip(Tmin,Tmax)
PltTemp = plb.contourf(data.grid.x/au, np.pi/2.-data.grid.y,Tdust, 30)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbT = plb.colorbar(PltTemp)
cbT.set_label('T [K]')
plb.show()


# In[ ]:

# Create an envelope model *WITH* external radiation:
cd ../run1_Ext
cp ../problem_params.inp .
cp ../dustkappa_osshenn_thinextra.inp .
SimpleDisk.simpleRadmc3Dmodel(disk=False,isrf=True,shisrf=True)
os.system('radmc3d mctherm setthreads 40')


# In[ ]:

# Read the model density and temperature
data2 = radmc3dPy.analyze.readData(dtemp=True,ddens=True,binary=False)

# Set plotting ranges:
Tmin = 0.
Tmax = 50.

# Plot the density
PltDens = plb.contourf(data2.grid.x/au, np.pi/2.-data.grid.y, np.log10(data.rhodust[:,:,0,0].T), 60)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbD = plb.colorbar(PltDens)
cbD.set_label(r'$\log_{10}{\rho}$')
plb.show()

# Plot the temperature
Tdust = data.dusttemp[:,:,0,0].T.clip(Tmin,Tmax)
PltTemp = plb.contourf(data2.grid.x/au, np.pi/2.-data.grid.y,Tdust, 30)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbT = plb.colorbar(PltTemp)
cbT.set_label('T [K]')
plb.show()


# In[ ]:




# In[ ]:

# Create an envelope model *WITH* external radiation but *NO STAR*:
cd ../run1_Ext_noStar
cp ../problem_params.inp .
cp ../dustkappa_osshenn_thinextra.inp .
SimpleDisk.simpleRadmc3Dmodel(disk=False,isrf=True,shisrf=False)
rm stars.inp
os.system('radmc3d mctherm setthreads 40')


# In[ ]:

# Read the model density and temperature
data = radmc3dPy.analyze.readData(dtemp=True,ddens=True,binary=False)

# Set plotting ranges:
Tmin = 0.
Tmax = 50.

# Plot the density
PltDens = plb.contourf(data.grid.x/au, np.pi/2.-data.grid.y, np.log10(data.rhodust[:,:,0,0].T), 60)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbD = plb.colorbar(PltDens)
cbD.set_label(r'$\log_{10}{\rho}$')
plb.show()

# Plot the temperature
Tdust = data.dusttemp[:,:,0,0].T.clip(Tmin,Tmax)
PltTemp = plb.contourf(data.grid.x/au, np.pi/2.-data.grid.y,Tdust, 30)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbT = plb.colorbar(PltTemp)
cbT.set_label('T [K]')
plb.show()


# In[ ]:

# Create an envelope model *WITH* external radiation but *NO STAR* 1e6 photons:
import matplotlib.pylab as plb
import radmc3dPy
import numpy as np
import os
from radmc3dPy.natconst import *
import sys
# Import the simple models:
import SimpleDisk
mkdir run2_Ext_noStar
cd run2_Ext_noStar
cp ../problem_params.inp .
cp ../dustkappa_osshenn_thinextra.inp .
SimpleDisk.simpleRadmc3Dmodel(disk=False,isrf=True,shisrf=False)
rm stars.inp
os.system('radmc3d mctherm setthreads 40')
os.system('radmc3d sed incl 67 phi 30 setthreads 40')


# In[ ]:

# Read the model density and temperature
data = radmc3dPy.analyze.readData(dtemp=True,ddens=True,binary=False)

# Set plotting ranges:
Tmin = 0.
Tmax = 14.

# Plot the density
PltDens = plb.contourf(data.grid.x/au, np.pi/2.-data.grid.y, np.log10(data.rhodust[:,:,0,0].T), 60)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbD = plb.colorbar(PltDens)
cbD.set_label(r'$\log_{10}{\rho}$')
plb.show()

# Plot the temperature
Tdust = data.dusttemp[:,:,0,0].T.clip(Tmin,Tmax)
PltTemp = plb.contourf(data.grid.x/au, np.pi/2.-data.grid.y,Tdust, 30)
plb.xlabel('r [AU]')
plb.ylabel(r'$\pi/2-\theta$')
plb.xscale('log')
cbT = plb.colorbar(PltTemp)
cbT.set_label('T [K]')
plb.show()


# In[ ]:

#Make an image
radmc3dPy.image.makeImage(npix=1000,sizeau=20000.,wav=3000.,incl=67,posang=80)
imag=radmc3dPy.image.readImage()
radmc3dPy.image.plotImage(imag,au=True,dpc=240.,log=True,maxlog=None,cmap='CMRmap',bunit='snu')
radmc3dPy.image.plotImage(imag,au=True,dpc=240.,log=True,vmin=-4,vmax=-3,cmap='CMRmap',bunit='snu')

