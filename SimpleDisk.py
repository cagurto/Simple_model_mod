# Code to set up a simple radmc-3d disk models with envelope

# Import the radmc3dPy module
import matplotlib.pylab as plb
import numpy as np
import os

# import the radmc3dPy library
import radmc3dPy
from radmc3dPy.natconst import *
import radmc3dPy.analyze as analyze

def ISradField(G=None,grid=None,ppar=None,show=False,write=False):
# Determine the G scale parameter of the FUV field
   G0 = 1.7
   if G==None and ppar.has_key('G'):
      G = ppar['G']
   elif G==None:
      G = G0

   erg2eV = 6.242e11  # conversion from electronvolt to erg
   
# Get the wavelength and frequency ranges of the model from grip object
   wav = grid.wav
   nwav = grid.nwav
   freq = grid.freq
   nfreq = grid.nfreq
   eV = hh * freq * erg2eV   # erg -> eV


# Create the black body components of the Black (1994) radiation field:
# see: https://ui.adsabs.harvard.edu/#abs/1994ASPC...58..355B/abstract
   p = [0., 0., 0., 0., -1.65, 0.]
   T = [7500.,4000.,3000.,250.,23.3,2.728]
   W = [1e-14, 1e-13, 4e-13, 3.4e-09, 2e-4, 1.]
   lp = [0.4,0.75,1.,1.,140.,1060.]
   
   IBlack = np.zeros(nwav)

   for i in range(len(p)):
       IBlack = IBlack + ( 2.*hh*freq**3. / cc**2. * (wav / lp[i])**p[i] * W[i] / (np.exp(hh * freq / kk / T[i])-1))
# don't worry about the warnings about division by zero, it comes from the (np.exp(hh * freq / kk / T[i])-1) part.
   
# Create the Draine (1978) radiation field using his original formula, 
# see eq. 11 in https://ui.adsabs.harvard.edu/#abs/1978ApJS...36..595D/abstract
   IDraineEv = (1.658E6 * eV) - (2.152e5 * eV**2) + (6.919E3 * eV**3)   # in photons cm^-2 s^-1 sr^-1 eV^-1
   IDraineErg = IDraineEv * hh**2 * freq * erg2eV                        # in erg cm^-2 s^-1 sr^-1 Hz^-1
   
   IDraine = IDraineErg * G/G0  # scale the FUV Draine field
# The formula was designed on the 5 eV (0.24 micron) to 13.6 eV (0.09 micron) range,
# limit the calculated intensities
   IDraine[wav < 0.09117381] = 0.0
   IDraine[wav > 0.24799276] = 0.0
   IBlack[wav < 0.24799276] = 0.0           # limit the Black field as well below 0.24 micron

# Combine the expressions for the different wavelength ranges:
   Iisrf = IBlack+IDraine
#
# Plot the results if asked
#
# Measurements to overplot from Allen book ISRF, unit: erg s-1 cm-2 mum-1 
   if show:
      wavObs = np.array([0.091, 0.1, 0.11, 0.13, 0.143, 0.18, 0.2, 0.21, 0.216, 0.23, 0.25, \
              0.346, 0.435, 0.55, 0.7, 0.9, 1.2, 1.8, 2.2, 2.4, 3.4, 4, 5, 12,  \
              25, 60, 100, 200, 300, 400, 600, 1000])
      freqObs_hz = cc / (wavObs / 1e4)
      FlamObs  = np.array([1.07e-2, 1.47e-2, 2.04e-2, 2.05e-2, 1.82e-2, 1.24e-2, 1.04e-2,   \
              9.61e-3, 9.17e-3, 8.25e-3, 7.27e-3, 1.3e-2, 1.5e-2, 1.57e-2,      \
              1.53e-2, 1.32e-2, 9.26e-3, 4.06e-3, 2.41e-3, 1.89e-3, 6.49e-4,    \
              3.79e-4, 1.76e-4, 1.7e-4, 6.0e-5, 4.6e-5, 7.3e-5, 2.6e-5, 5.4e-6, \
              1.72e-6, 3.22e-6, 7.89e-6])
      InuObs = (wavObs / 1e4)**2. / cc * FlamObs / (4.*3.14) * 1e4
# Plot   
      plb.xscale('log')
      plb.yscale('log')
      plb.plot(wavObs, InuObs * freqObs_hz, 'ro')
      plb.plot(wav, Iisrf * freq, color='black')
#
# Write to file is asked
#
   if write:   
      fname = 'external_source.inp'
      print 'Writing '+fname
      wfile = open(fname, 'w')
      wfile.write('%d\n'%2)     # this is the format, should be 2 always
      wfile.write('%d\n'%nwav)  # number of wavelength
      for ilam in range(nwav):  # first wavelength then Inu ISRF
          wfile.write('%.9e\n'%wav[ilam])
      for ilam in range(nwav):
          wfile.write('%.9e\n'%Iisrf[ilam])
      wfile.close()
            
def simpleRadmc3Dmodel(envelope=True, disk=True, isrf=False, shisrf=False, G=1.7):
# Read the parameters from the problem_params.inp file 
   modpar = analyze.readParams()

# Make a local copy of the ppar dictionary
   ppar = modpar.ppar

# Write out the important used parameters:

# Stellar properties:
   print "\nStellar properties:\n"
   print "T_star = ", ppar['tstar'][0], "K"
   print "R_star = ", ppar['rstar'][0] / rs, "Rsol"
   print "L_star = ", (ppar['rstar'][0]/rs)**2 * (ppar['tstar'][0]/5772.)**4, "L_sol"
# Grid parameters
   print "\nGrid parameters:\n"
   print "coordinate system:", ppar['crd_sys']
   print "x coordinate boundaries:", np.array(ppar['xbound']) / au, "in au"
   print "nx = ", ppar['nx']
   print "y coordinate boundaries:", ppar['ybound']
   print "ny = ", ppar['ny']
   if ppar.has_key('zbound'):
      print "z coordinate boundaries:", ppar['zbound']
      print "nz = ", ppar['nz']
   else:
      print "z coordinate is not activated (2D model)!"
   if ppar.has_key('xres_nlev'):
      print "Refinement along x axis:"
      print "in ", ppar['xres_nstep'], "steps"
# Wavelength grid
   print "\nWavelength grid:\n"
   print "Wavelength ranges", ppar['wbound'], "in micron"
   print "Bin number in range", ppar['nw']
# Envelope parameters:
   print "\nEnvelope parameters:\n"
   if envelope == True:
      print "Envelope is included in the model!"
      print "Density at 1 AU = ", ppar['rho0Env']
      print "Density power law index = ", ppar['prhoEnv']
      print "Opening angle [deg] = ", ppar['thetac_deg']
      print "Truncation radius [au] = ", ppar['rTrunEnv'] / au
   else:
      print "*NO* envelope is included!"
# Disk parameters:
   print "\nDisk parameters:\n"
   if disk == True:
      print "Disk is included in the model!"
      print "Mdisk (dust+gas) [MS] = ", ppar['mdisk'] / ms
      print "Rin [au] = ", ppar['rin'] / au
      print "Rdisk [au] = ", ppar['rdisk'] / au
      print "H(rdisk) = ", ppar['hrdisk']
      print "Power law of H = ", ppar['plh']
      print "Power law of surface density = ", ppar['plsig1']
      print "Dust-to-gas = ", ppar['dusttogas']
   else:
      print "*NO* disk is included!"

# --------------------------------------------------------------------------------------------
# Create the grid
# --------------------------------------------------------------------------------------------

# create the radmc3dGrid object:
   grid = analyze.radmc3dGrid()
# create the wavelength grid
   grid.makeWavelengthGrid(ppar=ppar)
# create the spatial grid
   grid.makeSpatialGrid(ppar=ppar)

# --------------------------------------------------------------------------------------------
# Create the input stellar radiation field
# --------------------------------------------------------------------------------------------

   radSources = analyze.radmc3dRadSources(ppar=ppar, grid=grid)
   radSources.getStarSpectrum(tstar=ppar['tstar'], rstar=ppar['rstar'])

# --------------------------------------------------------------------------------------------
# Create the dust density distribution 
# --------------------------------------------------------------------------------------------

# Create a radmc3dData object, this will contain the density
   data = analyze.radmc3dData(grid)

# Creat a grid for the 
   xx,yy = np.meshgrid(grid.x, np.pi/2.-grid.y)
   xx = xx.swapaxes(0,1)
   yy = yy.swapaxes(0,1)
    
#======================== Backgroud density =========================#

   rho_bg = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) + (ppar['bgdens'] * ppar['dusttogas'])

#============ Envelope density + Cavity + Cut off radius ============#

   rho_env = np.zeros([grid.nx, grid.ny, grid.nz, 1], dtype=np.float64) #array for the envelope density
   if envelope == True:
      thetac  = np.deg2rad(90) - np.deg2rad(ppar['thetac_deg']) #convert from degrees to radian
      for iz in range(grid.nz):
          for iy in range(grid.ny):
              for ix in range(grid.nx):
                  if np.abs(yy[ix,iy]) <= np.abs(thetac) and xx[ix,iy] >= ppar['rTrunEnv']:
                     rho_env[ix,iy,iz,0] = ppar['rho0Env'] * 1./(1. + (xx[ix,iy]/ppar['rTrunEnv'])**ppar['prhoEnv'])

#========================= Flarring disk ============================#

# set up the arrays containing the disk density
   rho_disk_tot  = np.zeros([grid.nx, grid.ny, grid.nz,1], dtype=np.float64)  # gas + dust
   rho_disk_dust = np.zeros([grid.nx, grid.ny, grid.nz,1], dtype=np.float64)  # only dust

   if disk == True:
# 1) first make the gas and dust disk (i.e. total mass)
      rr, th = np.meshgrid(grid.x, grid.y)
      z0 = np.zeros([grid.nx, grid.nz, grid.ny], dtype=np.float64)
      zz   = rr * np.cos(th)
      rcyl = rr * np.sin(th)

# Calculate the pressure scale height as a function of r, phi
      hp = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
      dum = ppar['hrdisk'] * (rcyl/ppar['rdisk'])**ppar['plh'] * rcyl
      dum = dum.swapaxes(0,1)
      for iz in range(grid.nz):  # copy the (x,y) grid to each z coordinates
          hp[:,:,iz] = dum
       
# We assume that sig0 is 1 and calculate the surface density profile:
      sigma = np.zeros([grid.nx, grid.ny, grid.nz], dtype=np.float64)
      dum1 = 1.0 * (rcyl/ppar['rdisk'])**ppar['plsig1']

   # Adding the smoothed inner rim
      if (ppar.has_key('srim_rout') & ppar.has_key('srim_plsig')):
         sig_srim = 1.0 * (ppar['srim_rout']*ppar['rin'] / ppar['rdisk'])**ppar['plsig1']
         dum2     = sig_srim * (rcyl / (ppar['srim_rout']*ppar['rin']))**ppar['srim_plsig']
         p    = -5.0
         dum  = (dum1**p + dum2**p)**(1./p)
      else:
         dum = dum1
      dum = dum.swapaxes(0,1)
      for iz in range(grid.nz):  # copy the (x,y) grid to each z coordinates
         sigma[:,:,iz] = dum

   # setting the disk surface density to 0 outside of the rin and rdisk range:
      for iy in range(grid.ny):
          ii = (rcyl[iy,:]<ppar['rin'])|(rcyl[iy,:]>ppar['rdisk'])
          sigma[ii,iy,:] = 0.0

# We calculate the disk density using the above surface density profile:
      for iz in range(grid.nz):
          for iy in range(grid.ny):
              rho_disk_tot[:,iy,iz,0] = sigma[:,iy,iz] /          \
                   (hp[:,iy,iz] * np.sqrt(2.0*np.pi)) *         \
                   np.exp(-0.5 * ((zz[iy,:])-z0[:,iz,iy]) *     \
                   ((zz[iy,:])-z0[:,iz,iy]) /                   \
                   (hp[:,iy,iz]*hp[:,iy,iz]))
                    
# Now we calculate the mass in rho_disk_tot and scale the density to get back the
# desired disk mass (['mdisk'] parameter):
      if ppar.has_key('mdisk') and ppar['mdisk']!=0.:
   # Calculate the volume of each grid cell
         vol  = grid.getCellVolume()
         mass = (rho_disk_tot[:,:,:,0]*vol).sum(0).sum(0).sum(0)
         rho_disk_tot = rho_disk_tot * (ppar['mdisk']/mass)
                    
# 2) Now calculate the dust density by scaling the disk mass with the dust/gas ratio:
      rho_disk_dust = np.array(rho_disk_tot) * ppar['dusttogas']
                    
#============== Adding up the density contributions =================#

   data.rhodust = rho_env + rho_disk_dust + rho_bg


# --------------------------------------------------------------------------------------------
# Now write out everything 
# --------------------------------------------------------------------------------------------
   print "\n Writing the input files:\n"
#Frequency grid
   grid.writeWavelengthGrid(old=False)
#Spatial grid
   grid.writeSpatialGrid(old=False)
#Input radiation field
   radSources.writeStarsinp(ppar=ppar, old=False)
# Write the external radiation field input file if needed
   ISradField(G,grid=grid,ppar=ppar,show=shisrf,write=isrf)
#Dust density distribution
   data.writeDustDens(binary=False, old=False)
#radmc3d.inp
   radmc3dPy.setup.writeRadmc3dInp(modpar=modpar)
#Master dust opacity file
   opac=analyze.radmc3dDustOpac()
   opac.writeMasterOpac(ext=ppar['dustkappa_ext'], scattering_mode_max=ppar['scattering_mode_max'], old=False)
