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
    pyplot.savefig(out_path+"Images/"+abs_+ord_+field_+'.eps')
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
    pyplot.savefig(out_path+"Images/"+abs_+field_+'.eps');
    return

#####################################################################################################################
#@njit
def plot_dft(dft_point_,field_,scale_factor,out_path):
    
    pyplot.figure(figsize=(2*scale_factor,1.2*scale_factor));
    
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
    x_ticks = numpy.arange(0 , 2 , 0.2)
    y_ticks = numpy.arange(min(y), max(y)+delta, (max(y)-min(y))/10)
    
    #x,y = zip(*sorted(zip(x,y)))   

    pyplot.xlim(0,2)
    #pyplot.xlim([min(x),max(x)/2])
    pyplot.ylim([min(y)*1.1,max(y)*1.1])
    pyplot.xticks(x_ticks, rotation=0)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.semilogy(x, y,
                color='#2929e3',
                ls='-',
                lw=1)    
    pyplot.savefig(out_path+'Images/dft_mp_'+field_+'.'+'png')
    return
