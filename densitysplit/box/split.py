import sys
from os import path
import numpy as np
import subprocess
from scipy.io import FortranFile

def generate_centres(
  ncentres, output_filename, sampling='uniform',
  sampling_filename=None, xmin=0, xmax=0,
  ymin=0, ymax=0, zmin=0,
  zmax=0, output_format='unformatted'
):
  np.random.seed(0)

  if sampling == 'subsampling':
    # check if file exists
    if not path.isfile(sampling_filename):
      raise FileNotFoundError(f'{sampling_filename} does not exist.')
    # check if this is a numpy file
    if '.npy' in sampling_filename:
      centres = np.load(sampling_filename)
    else:
      # if not, check if it is a text file
      try:
        centres = np.genfromtxt(sampling_filename)
      except:
        # else, check if it is an unformatted file
        try:
          fin = FortranFile(sampling_filename, 'r')
          nrows = fin.read_ints()[0]
          ncols = fin.read_ints()[0]
          pos = fin.read_reals(dtype=np.float64).reshape(nrows, ncols)
          idx = np.random.choice(nrows, size=ncentres, replace=False)
          centres = pos[idx]
        except:
          sys.exit('Format of sampling file not recognized.')

  elif sampling == 'uniform':
      x = np.random.uniform(xmin, xmax, ncentres)
      y = np.random.uniform(ymin, ymax, ncentres)
      z = np.random.uniform(zmin, zmax, ncentres)
      centres = np.c_[x, y, z]
  else:
      sys.exit('Sampling type not recognized')

  centres = centres.astype('float64')
  if output_format == 'unformatted':
    f = FortranFile(output_filename, 'w')
    nrows, ncols = np.shape(centres)
    f.write_record(nrows)
    f.write_record(ncols)
    f.write_record(centres)
    f.close()
  elif output_format == 'ascii':
    np.savetxt(output_filename, centres)
  elif output_format == 'numpy':
    np.save(output_filename, centres)
  else:
    sys.exit('Output format not recognized.')

  return centres
  

def filtered_density(
  centres_filename, tracers_filename, output_filename,
  filter_type, filter_size, ngrid, box_size,
  nthreads=1, dim1_min=0, dim1_max=None,
  output_format='unformatted'
):

  # check if files exist
  if not path.isfile(centres_filename):
    raise FileNotFoundError(f'{centres_filename} does not exist.')

  if not path.isfile(tracers_filename):
    raise FileNotFoundError(f'{tracers_filename} does not exist.')


  if dim1_max == None:
    if filter_type == 'tophat':
            dim1_max = filter_size
    elif filter_type == 'gaussian':
            dim1_max = 5 * filter_size

  binpath = path.join(path.dirname(__file__),
    'bin', '{}_filter.exe'.format(filter_type))

  cmd = [
    binpath, centres_filename, tracers_filename,
    output_filename, str(box_size), str(dim1_min),
    str(dim1_max), str(filter_size), str(ngrid),
    str(nthreads)
  ]

  subprocess.call(cmd)

  # open filter file
  f = FortranFile(output_filename, 'r')
  smoothed_delta = f.read_ints()[0]
  smoothed_delta = f.read_reals(dtype=np.float64)
  f.close()

  if output_format != 'unformatted':
    if output_format == 'npy':
      subprocess.call(['rm', output_filename])
      np.save(output_filename, smoothed_delta)
    elif output_format == 'ascii':
      np.savetxt(output_filename, smoothed_delta)
    else:
      print('Output format not recognized. Using unformatted F90 file.')
  
  return smoothed_delta


def split_centres(
  centres_filename, filter_filename, quantiles,
  handle=None, output_format='unformatted'
):

  # read centres
  # first check if this is a numpy file
  if '.npy' in centres_filename:
    centres = np.load(centres_filename)
  else:
    # if not, check if it is a text file
    try:
      centres = np.genfromtxt(centres_filename)
    except:
      # else, check if it is an unformatted file
      try:
        fin = FortranFile(centres_filename, 'r')
        nrows = fin.read_ints()[0]
        ncols = fin.read_ints()[0]
        centres = fin.read_reals(dtype=np.float64).reshape(nrows, ncols)
      except:
        sys.exit('Format of centres file not recognized.')

  # read smoothed densities
  # first check if this is a numpy file
  if '.npy' in filter_filename:
    smoothed_delta = np.load(filter_filename)
    ncentres = len(smoothed_delta)
  else:
    # if not, check if it is a text file
    try:
      smoothed_delta = np.genfromtxt(filter_filename)
      ncentres = len(smoothed_delta)
    except:
      # else, check if it is an unformatted file
      try:
        f = FortranFile(filter_filename, 'r')
        ncentres = f.read_ints()[0]
        smoothed_delta = f.read_reals(dtype=np.float64)
        f.close()
      except:
        sys.exit('Format of filter file not recognized.')
  idx = np.argsort(smoothed_delta)

  # sort centres using smoothed densities
  sorted_centres = centres[idx]

  # generate quantiles
  binned_centres = {}
  for i in range(1, quantiles + 1):
      binned_centres['DS{}'.format(i)] = sorted_centres[int((i-1)*ncentres/quantiles):int(i*ncentres/quantiles)]
      cout = binned_centres['DS{}'.format(i)]

      if handle != None:
            output_filename = handle + '_DS{}'.format(i)
      else:
          output_filename = centres_filename.split('.unf')[0] + '_DS{}'.format(i)
      
      if output_format == 'unformatted':
        output_filename += '.unf'
        f = FortranFile(output_filename, 'w')
        f.write_record(np.shape(cout)[0])
        f.write_record(np.shape(cout)[1])
        f.write_record(cout)
        f.close()
      
      elif output_format == 'ascii':
        output_filename += '.dat'
        np.savetxt(output_filename, cout)

      elif output_format == 'npy':
        np.save(output_filename, cout)
      
      else:
        sys.exit('Output format not recognized.')

  return binned_centres
