# Load Packages

import numpy
#from numba import njit
from matplotlib import pyplot
#get_ipython().magic('matplotlib inline')
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq
import scipy.fftpack
from matplotlib import rcParams
import matplotlib.ticker


from class_definition import *

#Load Data

#####################################################################################################################
#@njit
def load_data(path,file_name,headers,data_format):
    """
    Load data from an external file in nek5000 format.
    
    Arguments:
    ----------
    file_name (str)             : name of the data file to load
    headers (int)               : row where data starts (skips headers).
    data_format (str)           : FFormat of the file name, could be Nek5000 "ASCII" or "csv"
    
    Returns:
    --------
    volume_aux (volume)         : data loaded
                                      - volume_aux
                                          |-> .time []           : current time for the data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : x-velocity
                                                |-> .V []        : y-velocity
                                                |-> .W []        : z-velocity
                                                |-> .P []        : Pressure
                                                '-> .T []        : Temperature

    """

    volume_aux = volume()
    
    volume_aux.history = file_name

    if (data_format == "ASCII"):
        data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=headers)
        volume_aux.mesh.x = data[:,0]
        volume_aux.mesh.y = data[:,1]
        volume_aux.mesh.z = data[:,2]
        volume_aux.field.U = data[:,3]
        volume_aux.field.V = data[:,4]
        volume_aux.field.W = data[:,5]
        volume_aux.field.P = data[:,6]
        volume_aux.field.T = data[:,7]
    elif (data_format == "csv"):
        data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=headers,delimiter=',')
        volume_aux.mesh.x = data[:,9]
        volume_aux.mesh.y = data[:,10]
        volume_aux.mesh.z = data[:,11]
        volume_aux.field.U = data[:,2]
        volume_aux.field.V = data[:,3]
        volume_aux.field.W = data[:,4]
        volume_aux.field.P = data[:,0]
        volume_aux.field.T = data[:,1]
    elif (data_format == "2D_csv"):
        data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=headers,delimiter=',')
        volume_aux.mesh.x = data[:,8]
        volume_aux.mesh.y = data[:,9]
        volume_aux.mesh.z = data[:,10]
        volume_aux.field.U = data[:,5]
        volume_aux.field.V = data[:,6]
        volume_aux.field.W = data[:,7]
        volume_aux.field.P = data[:,0]
        volume_aux.field.T = data[:,1]
    elif (data_format == "oF_2D_csv"):
        data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=headers,delimiter=',')
        volume_aux.mesh.x = data[:,5]
        volume_aux.mesh.y = data[:,6]
        volume_aux.mesh.z = data[:,7]
        volume_aux.field.U = data[:,1]
        volume_aux.field.V = data[:,2]
        volume_aux.field.W = data[:,3]
        volume_aux.field.P = data[:,4]
        volume_aux.field.T = data[:,0]
    elif (data_format == "oF_csv"):
        data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=headers,delimiter=',')
        volume_aux.mesh.x = data[:,5]
        volume_aux.mesh.y = data[:,6]
        volume_aux.mesh.z = data[:,7]
        volume_aux.field.U = data[:,1]
        volume_aux.field.V = data[:,2]
        volume_aux.field.W = data[:,3]
        volume_aux.field.P = data[:,4]
        volume_aux.field.T = data[:,0]
    elif (data_format == "test_csv"):
        data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=headers,delimiter=',')
        volume_aux.mesh.x = data[:,10]
        volume_aux.mesh.y = data[:,11]
        volume_aux.mesh.z = data[:,12]
        volume_aux.field.U = data[:,0]
        volume_aux.field.V = data[:,1]
        volume_aux.field.W = data[:,2]
        volume_aux.field.P = data[:,3]
        volume_aux.field.T = data[:,4]
        volume_aux.field.uu = data[:,6]
        volume_aux.field.vv = data[:,5]
        volume_aux.field.ww = data[:,7]
        volume_aux.field.PP = data[:,8]
        volume_aux.field.TT = data[:,9]

    return volume_aux

#####################################################################################################################
#@njit
def load_monitoring_points_data(path,file_name):
    """
    load monitoring points data from an external file.
    
    Arguments:
    ----------
    path (str)                  : File's absolute path
    file_name (str)             : Name of the data file to load
    ti (float)                  : Initial time to be load
    
    Returns:
    --------
    mp_aux (volume)             : volumes coordinates of the mesh:
                                      - mp_aux
                                          |-> .time []           : current time for the data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : x-velocity
                                                |-> .V []        : y-velocity
                                                |-> .W []        : z-velocity
                                                |-> .P []        : Pressure
                                                '-> .T []        : Temperature

    """

    mp_aux = volume()
    
    data = numpy.loadtxt(fname=path + 'Data/' + file_name,skiprows=1)

    mp_aux.history = file_name

    mp_aux.time = data[:,0]
    mp_aux.field.U = data[:,1]
    mp_aux.field.V = data[:,2]
    mp_aux.field.W = data[:,3]
    mp_aux.field.P = data[:,4]
    mp_aux.field.T = data[:,5]
    
    
    return mp_aux
