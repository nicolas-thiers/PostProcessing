# Load Packages

import numpy
#from numba import njit
from matplotlib import pyplot
#get_ipython().magic('matplotlib inline')
from scipy.interpolate import interp1d
from scipy.fftpack import fft, fftfreq
import scipy.fftpack
from matplotlib import rcParams



from class_definition import *

#Load Data

#####################################################################################################################
#@njit
def load_data(file_name,headers):
    """
    Load data from an external file in nek5000 format.
    
    Arguments:
    ----------
    file_name (str)             : name of the data file to load
    headers (int)               : row where data starts (skips headers).
    
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
    
    data = numpy.loadtxt(fname=file_name,skiprows=headers)
    
    volume_aux.mesh.x = data[:,0]
    volume_aux.mesh.y = data[:,1]
    volume_aux.mesh.z = data[:,2]
    volume_aux.field.U = data[:,3]
    volume_aux.field.V = data[:,4]
    volume_aux.field.W = data[:,5]
    volume_aux.field.P = data[:,6]
    volume_aux.field.T = data[:,7]
    
    return volume_aux

#####################################################################################################################
#@njit
def load_monitoring_points_data(file_name):
    """
    load monitoring points data from an external file.
    
    Arguments:
    ----------
    file_name (str)             : name of the data file to load
    
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
    
    data = numpy.loadtxt(fname=file_name,skiprows=1)
    
    mp_aux.time = data[:,0]
    mp_aux.field.U = data[:,1]
    mp_aux.field.V = data[:,2]
    mp_aux.field.W = data[:,3]
    mp_aux.field.P = data[:,4]
    mp_aux.field.T = data[:,5]
    
    
    return mp_aux


# Processing data


#####################################################################################################################
#@njit
def slice_plane(data,a,b,c,d):
    """
    select data over the plane:
        pi : a*x + b*y + c*z = d.
    ** Current version only supports points from the mesh located over the mesh, dosnt interpolate data.
    
    Arguments:
    ----------
    data (volume)               : data to be sliced over the plane.
                                      - data
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


    a (real)                    : x-component of the plane's normal vector.
    b (real)                    : y-component of the plane's normal vector.
    c (real)                    : z-component of the plane's normal vector.
    d (real)                    : position parameter
    
    Returns:
    --------
    volume_aux (volume)         : filtered data:
                                      - data
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
    
    volume_aux.field.U = data.field.U[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    volume_aux.field.V = data.field.V[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    volume_aux.field.W = data.field.W[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    volume_aux.field.P = data.field.P[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    volume_aux.field.T = data.field.T[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]

    
    volume_aux.mesh.x = data.mesh.x[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    volume_aux.mesh.y = data.mesh.y[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    volume_aux.mesh.z = data.mesh.z[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]

    return volume_aux


#####################################################################################################################
#@njit
def compute_gradients(data):
    """
    Compute gradient of unstructured data.
    
    Arguments:
    ----------
    avg_data (volume)         : average field data: <U>
                                      - avg_data
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    rms_data (volume)         : average square field  data: <UU>
                                      - rms_data 
                                          |-> .time []           : current time for the average data square
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity square
                                                |-> .V []        : Mean y-velocity square
                                                |-> .W []        : Mean z-velocity square
                                                |-> .P []        : Mean Pressure square
                                                '-> .T []        : Mean Temperature square
    
    Returns:
    --------
    rms_fluctuation_out (volume)  : rms fluctuation field data: <u'u'>
                                      - rms_fluctuation_out
                                          |-> .time []           : current time for the rms fluctuation data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : rms fluctuation of x-velocity
                                                |-> .V []        : rms fluctuation of y-velocity
                                                |-> .W []        : rms fluctuation of z-velocity
                                                |-> .P []        : rms fluctuation of pressure
                                                '-> .T []        : rms fluctuation of temperature

    """
    # work in progess
    
    return gradient_out

#####################################################################################################################
#@njit
def monitoring_point(mp_data,mp_,nmp_,t0_):
    """
    Filter data with monitoring points info and output the desired monitor point history from the desired initial 
    time.
    
    Arguments:
    ----------
    mp_data (volume)            : Data loaded from hpts.out file:
                                      - mp_data
                                          |-> .time []           : Current time 
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
                                                
    mp_ (integer)               : Id number from the desire monitoring point.
    nmp_ (integer)              : Total number of monitoring point.
    t0_ (integer)               : Initial time to clip data.
    
    Returns:
    --------
    point_out (volume)          : Data from monitoring point:
                                      - point_out
                                          |-> .time []           : Current time 
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
    
    point_aux = volume()
    
    point_aux.time = mp_data.time[mp_::nmp_][t0_:]
    
    N = 1000
    window = numpy.ones(N)/N
    
    point_aux.field.U = mp_data.field.U[mp_::nmp_][t0_:]
    point_aux.field.V = mp_data.field.V[mp_::nmp_][t0_:]
    point_aux.field.W = mp_data.field.W[mp_::nmp_][t0_:]
    point_aux.field.P = mp_data.field.P[mp_::nmp_][t0_:]
    point_aux.field.T = mp_data.field.T[mp_::nmp_][t0_:]
    
    #point_aux.field.U = mp_data.field.U[mp_::nmp_][t0_:] - numpy.convolve(mp_data.field.U[mp_::nmp_][t0_:],window,'same')
    #point_aux.field.V = mp_data.field.V[mp_::nmp_][t0_:] - numpy.convolve(mp_data.field.V[mp_::nmp_][t0_:],window,'same')
    #point_aux.field.W = mp_data.field.W[mp_::nmp_][t0_:] - numpy.convolve(mp_data.field.W[mp_::nmp_][t0_:],window,'same')
    #point_aux.field.P = mp_data.field.P[mp_::nmp_][t0_:] - numpy.convolve(mp_data.field.P[mp_::nmp_][t0_:],window,'same')
    #point_aux.field.T = mp_data.field.T[mp_::nmp_][t0_:] - numpy.convolve(mp_data.field.T[mp_::nmp_][t0_:],window,'same')
    
    #point_aux.field.U = mp_data.field.U[mp_::nmp_][t0_:]- numpy.mean(mp_data.field.U[mp_::nmp_][t0_:])
    #point_aux.field.V = mp_data.field.V[mp_::nmp_][t0_:]- numpy.mean(mp_data.field.V[mp_::nmp_][t0_:])
    #point_aux.field.W = mp_data.field.W[mp_::nmp_][t0_:]- numpy.mean(mp_data.field.W[mp_::nmp_][t0_:])
    #point_aux.field.P = mp_data.field.P[mp_::nmp_][t0_:]- numpy.mean(mp_data.field.P[mp_::nmp_][t0_:])
    #point_aux.field.T = mp_data.field.T[mp_::nmp_][t0_:]- numpy.mean(mp_data.field.T[mp_::nmp_][t0_:]) 
    
    point_out = point_aux

    return point_out

#####################################################################################################################
#@njit
def smooth_data(data,N):
    """
    Filter data with monitoring points info and output the desired monitor point history from the desired initial 
    time.
    
    Arguments:
    ----------
    data (volume)               : Data to be loaded:
                                      - mp_data
                                          |-> .time []           : Current time 
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
    N (integer)                 : Size of the window used to smooth data.                                            
    
    Returns:
    --------
    smooth_data_out (volume)    : Smoothed data:
                                      - point_out
                                          |-> .time []           : Current time 
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
    
    window = numpy.ones(N)/N
    
    smooth_data_out = volume()
    
    smooth_data_out.time = data.time[N:][:-N]
    smooth_data_out.mesh.x = data.mesh.x
    smooth_data_out.mesh.y = data.mesh.y
    smooth_data_out.mesh.z = data.mesh.z
    
    smooth_data_out.field.U = data.field.U[N:][:-N] - numpy.convolve(data.field.U,window,'same')[N:][:-N]
    smooth_data_out.field.V = data.field.V[N:][:-N] - numpy.convolve(data.field.V,window,'same')[N:][:-N]
    smooth_data_out.field.W = data.field.W[N:][:-N] - numpy.convolve(data.field.W,window,'same')[N:][:-N]
    smooth_data_out.field.P = data.field.P[N:][:-N] - numpy.convolve(data.field.P,window,'same')[N:][:-N]
    smooth_data_out.field.T = data.field.T[N:][:-N] - numpy.convolve(data.field.T,window,'same')[N:][:-N]

    return smooth_data_out

#####################################################################################################################
#@njit
def dft_monitoring_point(point_):
    """
    Compute Discrete Fourier Transform for all fields in target point.
    
    Arguments:
    ----------
    point_ (volume)             : Data loaded from hpts.out file:
                                      - point_out
                                          |-> .time []           : Current time 
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
    
    Returns:
    --------
    dft_point_out (volume)      : Discrete Fourier Transform from the monitoring point's signal :
                                      - dft_point_out
                                          |-> .time []           : range of frecuencys 
                                          |-> .mesh
                                          |     |-> .x []        : -- not used --
                                          |     |-> .y []        : -- not used --
                                          |     '-> .z []        : -- not used --
                                          '-> .field
                                                |-> .U []        : Frecuency spectrum of x-velocity
                                                |-> .V []        : Frecuency spectrum of y-velocity
                                                |-> .W []        : Frecuency spectrum of z-velocity
                                                |-> .P []        : Frecuency spectrum of Pressure
                                                '-> .T []        : Frecuency spectrum of Temperature

    """
    
    point_aux = volume() 
    
    x = numpy.linspace(min(point_.time), max(point_.time) , len(point_.time))
    U_interp = interp1d(x, point_.field.U ) 
    V_interp = interp1d(x, point_.field.V )
    W_interp = interp1d(x, point_.field.W )
    P_interp = interp1d(x, point_.field.P )
    T_interp = interp1d(x, point_.field.T )
    
    x = numpy.linspace(0, len(point_.time)/(max(point_.time)-min(point_.time)), len(scipy.fftpack.fft(U_interp.y)))
    
    point_aux.time = x
    point_aux.field.U = scipy.fftpack.fft(U_interp.y)
    point_aux.field.V = scipy.fftpack.fft(V_interp.y)
    point_aux.field.W = scipy.fftpack.fft(W_interp.y)
    point_aux.field.P = scipy.fftpack.fft(P_interp.y)
    point_aux.field.T = scipy.fftpack.fft(T_interp.y)
    
    point_out = point_aux

    return point_out


######################
###   Calculator   ###
######################

#######################################################################################################
def calc_sum(data_1,data_2):
    """
    Addition operation of two volume class.
        A = B + C
    
    Arguments:
    ----------
    data_1 (volume)             : Volume 1
                                      - data_1
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    data_2 (volume)             : Volume 2
                                      - data_2
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    
    Returns:
    --------
    data_aux (volume)             : Sum volume 
                                      - data_aux
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature

    """
    
    data_aux = volume()
    
    data_aux.time = data_1.time
    
    data_aux.mesh.x = data_1.mesh.x
    data_aux.mesh.y = data_1.mesh.y
    data_aux.mesh.z = data_1.mesh.z
    
    data_aux.field.U = data_1.field.U + data_2.field.U
    data_aux.field.V = data_1.field.V + data_2.field.V
    data_aux.field.W = data_1.field.W + data_2.field.W
    data_aux.field.P = data_1.field.P + data_2.field.P
    data_aux.field.T = data_1.field.T + data_2.field.T
    
    return data_aux

#######################################################################################################
def calc_res(data_1,data_2):
    """
    Subtraction operation of two volume class.
        A = B - C
    
    Arguments:
    ----------
    data_1 (volume)             : Volume 1
                                      - data_1
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    data_2 (volume)             : Volume 2
                                      - data_2
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    
    Returns:
    --------
    data_aux (volume)           : subtract volume 
                                      - data_aux
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature

    """
    
    data_aux = volume()
    
    data_aux.time = data_1.time
    
    data_aux.mesh.x = data_1.mesh.x
    data_aux.mesh.y = data_1.mesh.y
    data_aux.mesh.z = data_1.mesh.z
    
    data_aux.field.U = data_1.field.U - data_2.field.U
    data_aux.field.V = data_1.field.V - data_2.field.V
    data_aux.field.W = data_1.field.W - data_2.field.W
    data_aux.field.P = data_1.field.P - data_2.field.P
    data_aux.field.T = data_1.field.T - data_2.field.T
    
    return data_aux

#######################################################################################################
def calc_mul(data_1,data_2):
    """
    Subtraction operation of two volume class.
        A = B * C
    
    Arguments:
    ----------
    data_1 (volume)             : Volume 1
                                      - data_1
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    data_2 (volume)             : Volume 2
                                      - data_2
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    
    Returns:
    --------
    data_aux (volume)           : Product volume 
                                      - data_aux
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature

    """
    
    data_aux = volume()
    
    data_aux.time = data_1.time
    
    data_aux.mesh.x = data_1.mesh.x
    data_aux.mesh.y = data_1.mesh.y
    data_aux.mesh.z = data_1.mesh.z
    
    data_aux.field.U = data_1.field.U * data_2.field.U
    data_aux.field.V = data_1.field.V * data_2.field.V
    data_aux.field.W = data_1.field.W * data_2.field.W
    data_aux.field.P = data_1.field.P * data_2.field.P
    data_aux.field.T = data_1.field.T * data_2.field.T
    
    return data_aux

#######################################################################################################
def calc_div(data_1,data_2):
    """
    Subtraction operation of two volume class.
        A = B / C
    
    Arguments:
    ----------
    data_1 (volume)             : Volume 1
                                      - data_1
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    data_2 (volume)             : Volume 2
                                      - data_2
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    
    Returns:
    --------
    data_aux (volume)           : subtract volume 
                                      - data_aux
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature

    """
    
    data_aux = volume()
    
    data_aux.time = data_1.time
    
    data_aux.mesh.x = data_1.mesh.x
    data_aux.mesh.y = data_1.mesh.y
    data_aux.mesh.z = data_1.mesh.z
    
    data_aux.field.U = data_1.field.U / data_2.field.U
    data_aux.field.V = data_1.field.V / data_2.field.V
    data_aux.field.W = data_1.field.W / data_2.field.W
    data_aux.field.P = data_1.field.P / data_2.field.P
    data_aux.field.T = data_1.field.T / data_2.field.T
    
    return data_aux

#######################################################################################################
def calc_power(data_1,N):
    """
    Power operation of two volume class.
        A = B ** N
    
    Arguments:
    ----------
    data_1 (volume)             : Volume 1
                                      - data_1
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature
    N (real)                    : Exponent
    
    Returns:
    --------
    data_aux (volume)           : Power value 
                                      - data_aux
                                          |-> .time []           : current time for the average data
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center
                                          |     |-> .y []        : y-coordinate of the volume's center
                                          |     '-> .z []        : z-coordinate of the volume's center
                                          '-> .field
                                                |-> .U []        : Mean x-velocity
                                                |-> .V []        : Mean y-velocity
                                                |-> .W []        : Mean z-velocity
                                                |-> .P []        : Mean Pressure
                                                '-> .T []        : Mean Temperature

    """
    
    data_aux = volume()
    
    data_aux.time = data_1.time
    
    data_aux.mesh.x = data_1.mesh.x
    data_aux.mesh.y = data_1.mesh.y
    data_aux.mesh.z = data_1.mesh.z
    
    data_aux.field.U = data_1.field.U ** N
    data_aux.field.V = data_1.field.V ** N
    data_aux.field.W = data_1.field.W ** N
    data_aux.field.P = data_1.field.P ** N
    data_aux.field.T = data_1.field.T ** N
    
    return data_aux


# Ploting Data

#####################################################################################################################
#@njit
def plot_3d(abs_,ord_,field_,data,contour_,number_of_level,field_min_,field_max_,scale_factor,out_path):
    """
    Plot data in two dimensional heat map.
    
    Arguments:
    ----------
    abs_ (str)                  : Space dimension to be plotted over x-axis.
    ord_ (str)                  : Space dimension to be plotted over y-axis.
    field_ (str)                : Field to be plotted in the heat map.
    data (volume)               : Data to be plotted:
                                      - data
                                          |-> .time []           : Current time 
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center.
                                          |     |-> .y []        : y-coordinate of the volume's center.
                                          |     '-> .z []        : z-coordinate of the volume's center.
                                          '-> .field
                                                |-> .U []        : x-velocity.
                                                |-> .V []        : y-velocity.
                                                |-> .W []        : z-velocity.
                                                |-> .P []        : Pressure.
                                                '-> .T []        : Temperature.
    contour_ (str)              : Add contour lines to he plot.
                                      - possible values:
                                          |-> "yes"              : iso-lines enable.
                                          '-> "no"               : iso-lines disable.
    number_of_level (int)       : Number of iso-lines when contour lines are enable.
    field_min_ ("str" or float) : minimal value in the field data to be plotted when float type is used.
                                      - possible string values:
                                          '-> "min"              : auto-adjust to the minimum value in the field 
                                                                   data.
    field_max_ ("str" or float) : maximal value in the field data to be plotted when float type is used.
                                      - possible string values:
                                          '-> "max"              : auto-adjust to the maximum value in the field 
    scale_factor (int)          : scale factor of the output plot.                                    
    out_path (str)              : output path where store graph.
    
    Returns:
    --------
    

    """
    
    x = getattr(getattr(data, "mesh"), abs_)
    y = getattr(getattr(data, "mesh"), ord_)
    field = getattr(getattr(data, "field"), field_)
    
    #pyplot.figure(figsize=(max(x)*scale_factor,max(y)*scale_factor))
    
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 14
    
    delta=1e-10
    
    pyplot.figure(figsize=(max(x)*scale_factor,max(y)*scale_factor))
    x_ticks = numpy.arange(min(x), max(x)+delta, (max(x)-min(x))/5)
    y_ticks = numpy.arange(min(y), max(y)+delta, (max(y)-min(y))/10)
    
    if field_min_ == "min" :
        min_ = min(field)
    else :
        min_ = field_min
    if field_max_ == "max" :
        max_ = max(field) + delta
    else :
        max_ = field_max_ + delta
        
    resolution = numpy.arange(min_, max_, (max_ - min_)/501)
    iso_levels = numpy.arange(min_, max_, (max_ - min_)/number_of_level)
    
    #resolution = numpy.arange(min(field), max(field)+delta, (max(field)-min(field))/100)
    #iso_levels = numpy.arange(min(field), max(field)+delta, (max(field)-min(field))/number_of_level)
    
    pyplot.title('Reduced Temperature at mid-large plane. \n')
    pyplot.xlim([min(x)-delta,max(x)+delta])
    pyplot.ylim([min(y)-delta,max(y)+delta])
    pyplot.xticks(x_ticks, rotation=75)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.tricontourf(x,y,field,
                       levels = resolution,
                       cmap='coolwarm')    
    pyplot.colorbar()
    
    if contour_ == "yes" :
        pyplot.tricontour(x,y,field,
                          levels = iso_levels,
                          ls='-.',
                          lw=0.5,
                          colors='black');
    pyplot.savefig(out_path+abs_+ord_+field_+'.eps')
    return

#####################################################################################################################
#@njit
def plot_over_line(abs_,field_,data,scale_factor,out_path):
    """
    Plot data previosuly sliced in a x-y dispersion graph.
    
    Arguments:
    ----------
    abs_ (str)                  : Space or time dimension to be plotted over x-axis.
    field_ (str)                : Field to be plotted in the graph.
    data (volume)               : Data to be plotted:
                                      - data
                                          |-> .time []           : Current time 
                                          |-> .mesh
                                          |     |-> .x []        : x-coordinate of the volume's center.
                                          |     |-> .y []        : y-coordinate of the volume's center.
                                          |     '-> .z []        : z-coordinate of the volume's center.
                                          '-> .field
                                                |-> .U []        : x-velocity.
                                                |-> .V []        : y-velocity.
                                                |-> .W []        : z-velocity.
                                                |-> .P []        : Pressure.
                                                '-> .T []        : Temperature.
    scale_factor (int)          : scale factor of the output plot.                                    
    out_path (str)              : output path where store graph.                                      
    Returns:
    --------
    

    """
    pyplot.figure(figsize=(2*scale_factor,1*scale_factor));
 
    if abs_ == "time" :
        x = getattr(data, abs_)
    else :
        x = getattr(getattr(data, "mesh"), abs_)
    y = getattr(getattr(data, "field"), field_) 


    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 14
    
    delta=1e-12
    
    x_ticks = numpy.arange(min(x), max(x)+delta, (max(x)-min(x))/10)
    y_ticks = numpy.arange(min(y), max(y)+delta, (max(y)-min(y))/10)
    
    x,y = zip(*sorted(zip(x,y)))
    
    delta_y = (max(y) - min(y))
    pyplot.xlim([min(x),max(x)])
    pyplot.ylim([min(y)-delta_y*0.1,max(y)+delta_y*0.1])
    pyplot.xticks(x_ticks, rotation=75)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.plot(x,y,
                color='#2929e3',
                ls='-',
                lw=1)
    pyplot.savefig(out_path+abs_+field_+'.eps');
    return

#####################################################################################################################
#@njit
def plot_dft(dft_point_,field_,scale_factor,out_path):
    
    pyplot.figure(figsize=(2*scale_factor,1*scale_factor));
    
    if field_ == "U":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.U[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude X-Velocity [-]')
    elif field_ == "V":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.V[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude Y-Velocity [-]')
    elif field_ == "W":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.W[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude Z-Velocity [-]')
    elif field_ == "P":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.P[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude Pressure [-]')        
    else:
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.T[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude Temperature [-]')
    
    x = dft_point_.time[0:len(dft_point_.time)][1:]
    pyplot.xlabel('Dimensionless Frecuency [-]')

    
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 14
    
    delta=1e-15
    
    #x_ticks = numpy.arange(min(x), max(x)/2, (max(x)/2-min(x))/10)
    x_ticks = numpy.arange(0, 2, 0.1)
    y_ticks = numpy.arange(min(y), max(y)+delta, (max(y)-min(y))/10)
    
    #x,y = zip(*sorted(zip(x,y)))   

    #pyplot.xlim([min(x),max(x)/2])
    pyplot.xlim(0.05, 2)
    pyplot.ylim([min(y)*1.1,max(y)*1.1])
    pyplot.xticks(x_ticks, rotation=75)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.semilogy(x, y,
                color='#2929e3',
                ls='-',
                lw=1)    
    pyplot.savefig(out_path+'dft_mp'+field_+'.eps')
    return
