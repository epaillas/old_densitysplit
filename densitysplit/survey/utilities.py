import sys
import numpy as np
from scipy.io import FortranFile
from astropy.io import fits
import argparse

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
    if is_random:
      weight = cat['WEIGHT_FKP']
    else:
      weight = cat['WEIGHT_FKP'] * cat['WEIGHT_SYSTOT']
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
