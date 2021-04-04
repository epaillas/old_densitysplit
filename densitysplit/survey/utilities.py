import sys
import numpy as np
from scipy.io import FortranFile
from astropy.io import fits
import argparse

def fits_to_unformatted(
  input_filename, output_filename, cosmology
):
  # open fits file
  with fits.open(input_filename) as hdul:
    cat = hdul[1].data

  # convert redshifts to distances
  dist = cosmology.ComovingDistance(cat['Z'])
  x = dist * np.sin(cat['DEC'] * np.pi / 180) * np.cos(cat['RA'] * np.pi / 180)
  y = dist * np.sin(cat['DEC'] * np.pi / 180) * np.sin(cat['RA'] * np.pi / 180)
  z = dist * np.cos(cat['DEC'] * np.pi / 180)

  # write result to output file
  cout = np.c_[x, y, z]
  nrows, ncols = np.shape(cout)
  f = FortranFile(output_filename, 'w')
  nrows, ncols = np.shape(cout)
  f.write_record(nrows)
  f.write_record(ncols)
  f.write_record(cout)
  f.close()