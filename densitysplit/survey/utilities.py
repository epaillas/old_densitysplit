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
  x = dist * np.sin(cat['DEC'] * np.pi / 180) * np.cos(cat['RA'] * np.pi / 180)
  y = dist * np.sin(cat['DEC'] * np.pi / 180) * np.sin(cat['RA'] * np.pi / 180)
  z = dist * np.cos(cat['DEC'] * np.pi / 180)

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
  x = dist * np.sin(cat[:,1] * np.pi / 180) * np.cos(cat[:,0] * np.pi / 180)
  y = dist * np.sin(cat[:,1] * np.pi / 180) * np.sin(cat[:,0] * np.pi / 180)
  z = dist * np.cos(cat[:,1] * np.pi / 180)

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
