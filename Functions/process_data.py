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
   
    volume_aux.history = data.history + '_slicePlane_pi:' +str(a)+'x'+str(b)+'y'+str(c)+'z'+str(d) 

    # Interpolated Mesh
    #int_size = 1000
    #volume_aux.mesh.x = numpy.linspace(min(data.mesh.x),max(data.mesh.x),num=int_size)
    #volume_aux.mesh.y = numpy.linspace(min(data.mesh.y),max(data.mesh.y),num=int_size)
    #for i in range(0,int_size):
	#    volume_aux.mesh.z[i] = (d - a*volume_aux.mesh.x[i] - b*volume_aux.mesh.y[i]) / c

    

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
def slice_line(data,point_1,point_2):
    """
    select data over the line which connects point_1 with point_2:
        l:      (x1,y1,z1) + t * (x2,y2,z2).
    ** Current version only supports points from the mesh located over the mesh, dosnt interpolate data.
    
    Arguments:
    ----------
    data (volume)               : data to be sliced over the line.
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
   
    volume_aux.history = data.history + '_sliceLine_l:' + point_1 + 'to' + point_2

    ### work in progress

    #volume_aux.field.U = data.field.U[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]

    #volume_aux.field.U = data.field.U[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    #volume_aux.field.V = data.field.V[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    #volume_aux.field.W = data.field.W[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    #volume_aux.field.P = data.field.P[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    #volume_aux.field.T = data.field.T[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]

    
    #volume_aux.mesh.x = data.mesh.x[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    #volume_aux.mesh.y = data.mesh.y[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]
    #volume_aux.mesh.z = data.mesh.z[a*data.mesh.x[:]+b*data.mesh.y[:]+c*data.mesh.z[:] == d]

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
    
    point_aux.history = mp_data.history + '_sensor_' + str(mp_) #+ 't0_' + str(t0_)

    i = 0
    while(mp_data.time[mp_::nmp_][i] < t0_ and i <= len(mp_data.time[mp_::nmp_])-2):
        i = i+1
    
    point_aux.time = mp_data.time[mp_::nmp_][i:]
    
    
    point_aux.field.u = mp_data.field.u[mp_::nmp_][i:]
    point_aux.field.v = mp_data.field.v[mp_::nmp_][i:]
    point_aux.field.w = mp_data.field.w[mp_::nmp_][i:]
    point_aux.field.P = mp_data.field.P[mp_::nmp_][i:]
    point_aux.field.T = mp_data.field.T[mp_::nmp_][i:]
    
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
    
    #window = numpy.ones(N)/N
    
    smooth_data_out = volume()
    
    smooth_data_out.history = data.history + "_smooth"

    #smooth_data_out.time = data.time[N:][:-N]
    smooth_data_out.time = data.time
    smooth_data_out.mesh.x = data.mesh.x
    smooth_data_out.mesh.y = data.mesh.y
    smooth_data_out.mesh.z = data.mesh.z
    
    #smooth_data_out.field.U = data.field.U[N:][:-N] - numpy.convolve(data.field.U,window,'same')[N:][:-N]
    #smooth_data_out.field.V = data.field.V[N:][:-N] - numpy.convolve(data.field.V,window,'same')[N:][:-N]
    #smooth_data_out.field.W = data.field.W[N:][:-N] - numpy.convolve(data.field.W,window,'same')[N:][:-N]
    #smooth_data_out.field.P = data.field.P[N:][:-N] - numpy.convolve(data.field.P,window,'same')[N:][:-N]
    #smooth_data_out.field.T = data.field.T[N:][:-N] - numpy.convolve(data.field.T,window,'same')[N:][:-N]

    smooth_data_out.field.u = data.field.u - numpy.average(data.field.u)
    smooth_data_out.field.v = data.field.v - numpy.average(data.field.v)
    smooth_data_out.field.w = data.field.w - numpy.average(data.field.w)
    smooth_data_out.field.P = data.field.P - numpy.average(data.field.P)
    smooth_data_out.field.T = data.field.T - numpy.average(data.field.T)

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
    U_interp = interp1d(x, point_.field.u ) 
    V_interp = interp1d(x, point_.field.v )
    W_interp = interp1d(x, point_.field.w )
    P_interp = interp1d(x, point_.field.P )
    T_interp = interp1d(x, point_.field.T )
    
    x = numpy.linspace(0, len(point_.time)/(max(point_.time)-min(point_.time)), len(scipy.fftpack.fft(U_interp.y)))
    
    point_aux.history = point_.history + "_dft"

    point_aux.time = x
    point_aux.field.u = scipy.fftpack.fft(U_interp.y)
    point_aux.field.v = scipy.fftpack.fft(V_interp.y)
    point_aux.field.w = scipy.fftpack.fft(W_interp.y)
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

