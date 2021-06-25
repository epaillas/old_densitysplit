import numpy as np
from scipy.io import FortranFile
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import eval_legendre
from scipy.integrate import simps


def save_as_unformatted(data, filename):
    '''
    Saves a numpy array as an unformatted
    Fortran 90 file that can be handled by
    this package's numerical routines.

    Parameters:  data: ND array_like
                 Array to be saved.

                 filename: str
                 Name of the output file.
    '''
    data = np.asarray(data)

    nrows, ncols = np.shape(data)
    f = FortranFile(filename, 'w')
    nrows, ncols = np.shape(data)
    f.write_record(nrows)
    f.write_record(ncols)
    f.write_record(data)
    f.close()


def read_array_2d(filename):
    '''
    Read a two dimensional from an ascii
    file.

    Parameters:  filename: str
                 Name of the ascii file containing
                 the array
    '''
    data = np.genfromtxt(filename)
    dim1 = np.unique(data[:, 0])
    dim2 = np.unique(data[:, 1])

    vary_dim2 = False
    if data[0, 0] == data[1, 0]:
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
    '''
    Calculate the ell multipole moment
    from the redshift-space correlation
    function xi_smu.

    Parameters:  ell: int
                 Multipole moment to calculate.

                 s: 1D array_like
                 Radial bins of the correlation function.

                 mu: 1D array_like
                 Mu bins of the correlation function, where
                 mu is the cosine of the angle with respect
                 to the line of sight.

                 xi_smu: 2D array_like
                 Correlation function binned in s and mu.
    '''
    multipole = np.zeros(xi_smu.shape[0])
    if mu.min() < 0:
        factor = 2
        mumin = -1
    else:
        factor = 1
        mumin = 0
    for i in range(xi_smu.shape[0]):
        mufunc = InterpolatedUnivariateSpline(mu, xi_smu[i, :], k=3, ext=3)
        xaxis = np.linspace(mumin, 1, 1000)
        lmu = eval_legendre(ell, xaxis)
        yaxis = mufunc(xaxis) * (2 * ell + 1) / factor * lmu
        multipole[i] = simps(yaxis, xaxis)
    return s, multipole
