# Load Packages

import numpy
import csv 
#from matplotlib import pyplot
#from scipy.interpolate import interp1d
#from scipy.fftpack import fft, fftfreq
#import scipy.fftpack
#from matplotlib import rcParams
#import matplotlib.ticker

from class_definition import *
from create_class import *

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

    fieldsDict = {}
    csvFile = csv.DictReader(open( path + file_name ), delimiter=',')   
    for i in csvFile.fieldnames[:-3]:
        fieldsDict[i] = []

    data = numpy.loadtxt(fname=path + file_name,skiprows=headers,delimiter=',')

	create_classes(volume_=volume_aux , classDict_ = fieldsDict)	

	for i in range(0,len(csvFile.fieldnames[:-3])):
	    setattr(field,list(fieldsDict)[i],data[:,i])

	print("Done loading Data file")

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
