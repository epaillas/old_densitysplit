import sys
import numpy as np
from scipy.io import FortranFile
from astropy.io import fits
import argparse
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import eval_legendre
from scipy.integrate import simps

def fits_to_ascii(
  input_filename, output_filename, cosmology,
  is_random=False, equal_weights=False, zrange=None
):
  # open fits file
  with fits.open(input_filename) as hdul:
    cat = hdul[1].data

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat['Z'] > zmin) & (cat['Z'] < zmax)
    cat = cat[ind]

  # convert redshifts to distances
  dist = cosmology.ComovingDistance(cat['Z'])
  x = dist * np.cos(cat['DEC'] * np.pi / 180) * np.cos(cat['RA'] * np.pi / 180)
  y = dist * np.cos(cat['DEC'] * np.pi / 180) * np.sin(cat['RA'] * np.pi / 180)
  z = dist * np.sin(cat['DEC'] * np.pi / 180)

  if not equal_weights:
    weight = cat['WEIGHT_FKP']
    if not is_random:
      weight *= cat['WEIGHT_SYSTOT'] * (cat['WEIGHT_CP'] + cat['WEIGHT_NOZ'] - 1)
  else:
    weight = np.ones(len(cat))

  #write result to output file
  cout = np.c_[x, y, z, weight]
  np.savetxt(output_filename, cout)

def fits_to_unformatted(
  input_filename, output_filename, cosmology,
  is_random=False, equal_weights=False, zrange=None
):
  # open fits file
  with fits.open(input_filename) as hdul:
    cat = hdul[1].data

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat['Z'] > zmin) & (cat['Z'] < zmax)
    cat = cat[ind]

  # convert redshifts to distances
  dist = cosmology.ComovingDistance(cat['Z'])
  x = dist * np.cos(cat['DEC'] * np.pi / 180) * np.cos(cat['RA'] * np.pi / 180)
  y = dist * np.cos(cat['DEC'] * np.pi / 180) * np.sin(cat['RA'] * np.pi / 180)
  z = dist * np.sin(cat['DEC'] * np.pi / 180)

  if not equal_weights:
    weight = cat['WEIGHT_FKP']
    if not is_random:
      weight *= cat['WEIGHT_SYSTOT'] * (cat['WEIGHT_CP'] + cat['WEIGHT_NOZ'] - 1)
  else:
    weight = np.ones(len(cat))

  #write result to output file
  cout = np.c_[x, y, z, weight]
  nrows, ncols = np.shape(cout)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()

def fits_to_cute(
  input_filename, output_filename, is_random=False,
  equal_weights=False, zrange=None
):
  # open fits file
  with fits.open(input_filename) as hdul:
    cat = hdul[1].data

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat['Z'] > zmin) & (cat['Z'] < zmax)
    cat = cat[ind] 

  ra = cat['RA']
  dec = cat['DEC']
  redshift = cat['Z']

  if not equal_weights:
    if is_random:
      weight = cat['WEIGHT_FKP']
    else:
      weight = cat['WEIGHT_FKP'] * cat['WEIGHT_SYSTOT']
  else:
    weight = np.ones(len(cat))

  #write result to output file
  cout = np.c_[ra, dec, redshift, weight]
  np.savetxt(output_filename, cout)


def revolver_to_unformatted(
  input_filename, output_filename, cosmology,
  equal_weights=False, zrange=None
):
  # open numpy file
  cat = np.load(input_filename)

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat[:,2] > zmin) & (cat[:,2] < zmax)
    cat = cat[ind] 

  # convert redshifts to distances
  dist = cosmology.ComovingDistance(cat[:,2])
  x = dist * np.cos(cat[:,1] * np.pi / 180) * np.cos(cat[:,0] * np.pi / 180)
  y = dist * np.cos(cat[:,1] * np.pi / 180) * np.sin(cat[:,0] * np.pi / 180)
  z = dist * np.sin(cat[:,1] * np.pi / 180)

  if not equal_weights:
    weight = cat[:,3]
  else:
    weight = np.ones(len(cat))

  #write result to output file
  cout = np.c_[x, y, z, weight]
  nrows, ncols = np.shape(cout)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()


def patchy_to_unformatted(
  input_filename, output_filename, cosmology,
  is_random=False, equal_weights=False, zrange=None
):
  # open text file
  cat = np.genfromtxt(input_filename)

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat[:,2] > zmin) & (cat[:,2] < zmax)
    cat = cat[ind]

  # convert redshifts to distances
  dist = cosmology.ComovingDistance(cat[:,2])
  x = dist * np.cos(cat[:,1] * np.pi / 180) * np.cos(cat[:,0] * np.pi / 180)
  y = dist * np.cos(cat[:,1] * np.pi / 180) * np.sin(cat[:,0] * np.pi / 180)
  z = dist * np.sin(cat[:,1] * np.pi / 180)

  if not equal_weights:
    if is_random:
      weight = (cat[:,5] * cat[:,6]) / (1 + 10000 * cat[:,3])
    else:
      weight = (cat[:,6] * cat[:,7]) / (1 + 10000 * cat[:,4])
  else:
    weight = np.ones(len(cat))

  #write result to output file
  cout = np.c_[x, y, z, weight]
  nrows, ncols = np.shape(cout)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()


def patchy_to_revolver(
  input_filename, output_filename,
  is_random=False, zrange=None
):
  # open text file
  cat = np.genfromtxt(input_filename)

  if zrange is not None:
    zmin, zmax = zrange
    ind = (cat[:,2] > zmin) & (cat[:,2] < zmax)
    cat = cat[ind]

  ra = cat[:,0]
  dec = cat[:,1]
  z = cat[:,2]

  if is_random:
    fkp = 1 / (1 + 10000 * cat[:,3])
    veto = cat[:,5]
    cp = cat[:,6]
  else:
    fkp = 1 / (1 + 10000 * cat[:,4])
    veto = cat[:,6]
    cp = cat[:,7]

  #write result to output file
  cout = np.c_[ra, dec, z, fkp, cp, veto]
  np.savetxt(output_filename, cout)

def read_array_2d(input_filename):
  data = np.genfromtxt(input_filename)
  dim1 = np.unique(data[:,0])
  dim2 = np.unique(data[:,1])

  vary_dim2 = False
  if data[0,0] == data[1,0]:
      vary_dim2 = True

  result = np.zeros([len(dim1), len(dim2)])
  counter = 0
  if vary_dim2:
      for i in range(len(dim1)):
          for j in range(len(dim2)):
              result[i, j] = data[counter, 2]
              counter += 1
  else:
      for i in range(len(dim2)):
          for j in range(len(dim1)):
              result[j, i] = data[counter, 2]
              counter += 1
  return dim1, dim2, result

def get_multipole(ell, s, mu, xi_smu):
  multipole = np.zeros(xi_smu.shape[0])
  if mu.min() < 0:
      factor = 2
      mumin = -1
  else:
      factor = 1
      mumin=0
  for i in range(xi_smu.shape[0]):
      mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3, ext=3)
      xaxis = np.linspace(mumin, 1, 1000)
      lmu = eval_legendre(ell, xaxis)
      yaxis = mufunc(xaxis) * (2 * ell + 1) / factor * lmu
      multipole[i] = simps(yaxis, xaxis)
  return s, multipole