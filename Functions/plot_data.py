# Load Packages

import numpy
#from numba import njit
from matplotlib import pyplot
#get_ipython().magic('matplotlib inline')
from scipy.interpolate import interp1d,griddata
from scipy.fftpack import fft, fftfreq
import scipy.fftpack
from matplotlib import rcParams
from matplotlib.colors import LogNorm

import colormaps as cmaps
pyplot.register_cmap(name='viridis', cmap=cmaps.viridis)
pyplot.set_cmap(cmaps.viridis)
pyplot.register_cmap(name='magma', cmap=cmaps.magma)
pyplot.set_cmap(cmaps.magma)

from class_definition import *

# Ploting Data

#####################################################################################################################
#@njit
def plot_3d(abs_,ord_,field_,data,contour_,number_of_level,field_min_,field_max_,scale_factor,out_path,plot_format = "png",normalized=True,cmaps='viridis'):
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
 
    if (normalized == True):
        x = (x-min(x))/(max(x)-min(x))
        y = (y-min(y))/(max(y)-min(y))
   
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    rcParams.update({'figure.autolayout': True})

    delta=1e-10

    if field_min_ == "min" :
        min_ = min(field) - delta
    else :
        min_ = field_min_ - delta
    if field_max_ == "max" :
        max_ = max(field) + delta
    else :
        max_ = field_max_ + delta

    pyplot.figure(figsize=(max(x)*1.5*scale_factor,max(y)*scale_factor))
    x_ticks = numpy.linspace(min(x),max(x),num=5) 
    y_ticks = numpy.linspace(min(y),max(y),num=5) 

    resolution = numpy.arange(min_-delta, max_+delta, (max_ - min_)/number_of_level)
    iso_levels = numpy.arange(min_, max_, (max_ - min_)/number_of_level)
    
    #resolution = numpy.arange(min(field), max(field)+delta, (max(field)-min(field))/100)
    #iso_levels = numpy.arange(min(field), max(field)+delta, (max(field)-min(field))/number_of_level)
    
    #pyplot.title('Reduced Temperature at mid-large plane. \n')
    pyplot.xlim([min(x)-delta,max(x)+delta])
    pyplot.ylim([min(y)-delta,max(y)+delta])
    pyplot.xticks(x_ticks, rotation=75)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.tricontourf(x,y,field,
                       levels = resolution,
                       cmap=cmaps)    
    pyplot.colorbar()
    
    if contour_ == "yes" :
        pyplot.tricontour(x,y,field,
                          levels = iso_levels,
                          ls='-.',
                          lw=0.5,
                          colors='black');
    pyplot.savefig(out_path+"Images/Slices/" + data.history[20:-4] + "_" + field_ + '.' + plot_format )
    pyplot.close("all")
    return

#####################################################################################################################
#@njit
def plot_streamlines(abs_,ord_,vectorX_,vectorY_,field_,data_,outPath_,plotFormat_ = "png",normalized=True,cmap_='viridis',numberOfLevel_=11,fieldMin_="min",fieldMax_="max",scaleFactor_=10,xrange_=['min','max'],yrange_ = ['min','max']):
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
    
    x1 = getattr(getattr(data_, "mesh"), abs_)
    y1 = getattr(getattr(data_, "mesh"), ord_)
    
    fieldX1 = getattr(getattr(data_, "field"), vectorX_)
    fieldY1 = getattr(getattr(data_, "field"), vectorY_)
    if (field_ == "magnitud"):
        field1 = numpy.sqrt(fieldX1**2+fieldY1**2)
    else:
        field1 = getattr(getattr(data_, "field"), field_)
        print(field1) 

    #####    Interpolation to regular grid    #####
    x = numpy.linspace(x1.min(),x1.max(),100)
    y = numpy.linspace(y1.min(),y1.max(),100)
    xi, yi = numpy.meshgrid(x,y,sparse=True) 
    px = x1.flatten()
    py = y1.flatten()
    pu = fieldX1.flatten()
    pv = fieldY1.flatten()
    pfield = field1.flatten()
    fieldX = griddata(zip(px,py), pu, (xi,yi))
    fieldY = griddata(zip(px,py), pv, (xi,yi))
    field = griddata(zip(px,py), pfield, (xi,yi))

    if (normalized == True):
        x = (x-min(x))/(max(x)-min(x))
        y = (y-min(y))/(max(y)-min(y))

   
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    #rcParams.update({'figure.autolayout': True})

    delta=1e-10
    if fieldMin_ == "min" :
        fieldMin = numpy.amin(field) - delta
    else :
        fieldMin = fieldMin_ - delta
    if fieldMax_ == "max" :
        fieldMax = numpy.amax(field) + delta
    else :
        fieldMax = fieldMax_ + delta
    resolution = numpy.arange(fieldMin - delta, fieldMax + delta, (fieldMax - fieldMin)/numberOfLevel_)


    
    lw = 5*field / numpy.amax(field1)

    pyplot.figure(figsize=(max(x)*scaleFactor_,max(y)*scaleFactor_))

    if (xrange_ == ['min','max']):
        pyplot.xlim([min(x)-delta,max(x)+delta])
        x_ticks = numpy.linspace(min(x),max(x),num=5) 
    else :
        pyplot.xlim([xrange_[0]-delta,xrange_[1]+delta])
        x_ticks = numpy.linspace(xrange_[0],xrange_[1],num=5) 
    if (yrange_ == ['min','max']):
        pyplot.ylim([min(y)-delta,max(y)+delta])
        y_ticks = numpy.linspace(min(y),max(y),num=5) 
    else :
        pyplot.ylim([yrange_[0]-delta,yrange_[1]+delta])
        y_ticks = numpy.linspace(yrange_[0],yrange_[1],num=5) 

    pyplot.xticks(x_ticks, rotation=75)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.contourf(x,y,field,alpha=0.5,cmap=cmap_)
    pyplot.colorbar()
    pyplot.streamplot(x, y, fieldX, fieldY, density=1,color='k',linewidth=lw)
    '''
    pyplot.tricontourf(x,y,field,
                       levels = resolution,
                       cmap=cmaps)    
    '''
    pyplot.savefig(outPath_+"Images/streamLines/" + data_.history[20:-4] + "_" + field_+'.'+plotFormat_,format=plotFormat_,bbox_inches='tight',rasterized=True,dpi=300)
    pyplot.close("all")
    return


#####################################################################################################################
#@njit
def plot_PDF(abs_,ord_,field_,data,x_range,y_range,scale_factor,out_path,plot_format = "png"):
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
    pyplot.figure(figsize=(scale_factor,scale_factor))
    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    rcParams.update({'figure.autolayout': True})

    x = getattr(getattr(data, "field"), abs_) 
    y = getattr(getattr(data, "field"), ord_)
    xedges = numpy.linspace(x_range[0],x_range[1],1000)
    yedges = numpy.linspace(y_range[0],y_range[1],1000)
    
    print(xedges)
    print(yedges)

    H, xedges, yedges = numpy.histogram2d(x, y, bins=(xedges,yedges))
    H = H.T
    X,Y = numpy.meshgrid(xedges,yedges)
    pyplot.contourf(X[1:,1:],Y[1:,1:],H,cmap = "viridis",norm = LogNorm())
    pyplot.colorbar()
    pyplot.contour(X[1:,1:],Y[1:,1:],H,colors = "black",norm = LogNorm())
    pyplot.savefig(out_path+"Images/" + data.history[20:] + "_" + field_ + '.' + plot_format )
    pyplot.close("all")
    return

#####################################################################################################################
#@njit
def plot_over_line(abs_,field_,data,scale_factor,out_path,plot_format = "png",labels = {'x_label':'','y_label':'','title':''}):
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
    pyplot.figure(figsize=(2*scale_factor,1.2*scale_factor));
 
    if abs_ == "arc_length":
        x = getattr(getattr(data, "field"), abs_)
    elif abs_ == "time" :
        x = getattr(data, abs_)
    else :
        x = getattr(getattr(data, "mesh"), abs_)
    y = getattr(getattr(data, "field"), field_) 

    '''
    if field_ == "U":
        pyplot.ylabel('X-Velocity [-]')
    elif field_ == "V":
        pyplot.ylabel('Y-Velocity [-]')
    elif field_ == "W":
        pyplot.ylabel('Z-Velocity [-]')
    elif field_ == "P":
        pyplot.ylabel('Pressure [-]')        
    else:
        pyplot.ylabel('Reduced Temperature [-]')

    if abs_ == "arc_length":
        pyplot.xlabel('')
    elif abs_ == "x":
        pyplot.xlabel('x')
    elif abs_ == "y":
        pyplot.xlabel('y')
    elif abs_ == "z":
        pyplot.xlabel('z')
    elif abs_ == "time":
        pyplot.xlabel('time')        
    '''
    pyplot.xlabel(labels['x_label'])
    pyplot.ylabel(labels['y_label'])
    pyplot.title(labels['title'])

    rcParams['font.family'] = 'serif'
    rcParams['font.size'] = 12
    
    delta=1e-12

    x_ticks = numpy.arange( min(x) , max(x)+delta , (max(x)-min(x)+delta)/10)
    y_ticks = numpy.arange( min(y) , max(y)+delta , (max(y)-min(y)+delta)/10)
    
    x,y = zip(*sorted(zip(x,y)))
    
    #delta_y = (max(y) - min(y))
    pyplot.xlim([min(x),max(x)])
    pyplot.ylim([min(y),max(y)])
    pyplot.xticks(x_ticks, rotation=0 )
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.plot(x,y,
                color='black',
                ls='-',
                lw=1)
    if (abs_ == 'time'):
        pyplot.savefig(out_path+"Images/monitoringPoints/" + data.history[20:] + "_"  + abs_ + field_ + '.' + plot_format);
    else:
        pyplot.savefig(out_path+"Images/Lines/" + data.history[20:-4] + "_"  + abs_ + field_ + '.' + plot_format);
    pyplot.close("all") 
    return

#####################################################################################################################
#@njit
def plot_dft(dft_point_,field_,scale_factor,out_path,plot_format = "png"):

    pyplot.figure(figsize=(2*scale_factor,1.2*scale_factor));
    
    if field_ == "u":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.u[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude X-Velocity [-]')
    elif field_ == "v":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.v[0:len(dft_point_.time)][1:])
        pyplot.ylabel('Amplitude Y-Velocity [-]')
    elif field_ == "w":
        y = 2.0/len(dft_point_.time) * numpy.abs(dft_point_.field.w[0:len(dft_point_.time)][1:])
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
    rcParams['font.size'] = 12
    
    delta=1e-15
    
    x_ticks = numpy.arange(min(x), max(x)/2, (max(x)/2-min(x))/10)
    x_ticks = numpy.arange(0 , 4 , 0.2)
    y_ticks = numpy.arange(min(y), max(y)+delta, (max(y)-min(y)+delta)/10)
    
    #x,y = zip(*sorted(zip(x,y)))   

    pyplot.xlim(0,4)
    #pyplot.xlim([min(x),max(x)/2])
    pyplot.ylim([min(y)*1.1,max(y)*1.1])
    pyplot.xticks(x_ticks, rotation=0)
    pyplot.yticks(y_ticks)
    pyplot.grid()
    pyplot.plot(x, y,
                color='black',
                ls='-',
                lw=1)    
    pyplot.savefig(out_path+'Images/monitoringPoints/'+ dft_point_.history[20:] + field_ +'.' + plot_format)
    pyplot.close("all") 
    return
