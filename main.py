import sys 
import os
path = "/home/nicolas/Dropbox/Doctorado/PostProcessing"
sys.path.append(os.path.abspath(path))

# Load Functions

from class_definition import *
from functions import *

# Main program


def main():
    
    #path to data, must be in absolute form
    path = '/home/nicolas/Dropbox/Doctorado/Simulaciones/3D/base_case/'
    
    # Determinar la fila donde comienza la data
    #!head -n 411 high_ray.fld02 | tail -n 4
    
    #data = load_data(path + 'high_ray.fld13',409)
    mp_data = load_monitoring_points_data(path + 'hpts.out')
    #data_rms = load_data(path + 'rmshigh_ray.fld03',409)
    #data_avg = load_data(path + 'avghigh_ray.fld03',409)
    '''
    plane_1 = slice_plane(data_avg,0,1,0,0)
    
    plot_3d("x","z","T",
            plane_1,
            "yes",10,"min","max",
            30,path)
    plot_3d("x","z","W",
            plane_1,
            "no",10,"min","max",
            30,path)
    
    line_1 = slice_plane(slice_plane(data_avg,0,1,0,0),0,0,1,0)
    line_2 = slice_plane(slice_plane(data_avg,0,1,0,0),1,0,0,0)    

    
    plot_over_line("x","T",line_1,5,path)
    plot_over_line("x","W",line_1,5,path)
    plot_over_line("z","T",line_2,5,path)
    '''
    point_1 = monitoring_point(mp_data,1,5,30000)
    smooth_point_1 = smooth_data(point_1,100)

    plot_over_line("time","U",point_1,5,path)
    plot_over_line("time","V",point_1,5,path)
    plot_over_line("time","W",point_1,5,path)
    plot_over_line("time","P",point_1,5,path)
    plot_over_line("time","T",point_1,5,path)

    dft_point_1 = dft_monitoring_point(smooth_point_1)
    
    plot_dft(dft_point_1,"U",5,path)
    plot_dft(dft_point_1,"V",5,path)
    plot_dft(dft_point_1,"W",5,path)
    plot_dft(dft_point_1,"P",5,path)
    plot_dft(dft_point_1,"T",5,path)
    
    '''
    data_prime = calc_res(data_rms,calc_power(data_avg,2))
    
    plane_1 = slice_plane(data_prime,0,1,0,0)
    
    plot_3d("x","z","U",
            plane_1,
            "yes",10,"min","max",
            30,path)
    plot_3d("x","z","V",
            plane_1,
            "no",10,"min","max",
            30,path)
    plot_3d("x","z","W",
            plane_1,
            "yes",10,"min","max",
            30,path)
    plot_3d("x","z","P",
            plane_1,
            "yes",10,"min","max",
            30,path)
    plot_3d("x","z","T",
            plane_1,
            "yes",10,"min","max",
            30,path)
    
    return data,mp_data,data_rms,data_avg
    '''
    return mp_data


main()
'''
path = '/home/nicolas/Dropbox/Doctorado/Simulaciones/3D/thermal_actuator/fase/phi_0.66/'
mp_data = load_monitoring_points_data(path + 'hpts.out')

point_1 = monitoring_point(mp_data,4,7,2000)
point_2 = monitoring_point(mp_data,6,7,2000)

x1 = point_1.time #- 0.44
y1 = (point_1.field.T - numpy.mean(point_1.field.T))
x2 = point_2.time
y2 = (point_2.field.T - numpy.mean(point_2.field.T))*13



pyplot.xlim(100, 110)
pyplot.ylim(-0.07,0.07)
pyplot.xlabel("Dimensionless time [-]")
pyplot.ylabel("Reduced Temperature [-]")
pyplot.plot(x1,y1,
               color = 'b',
               alpha = 0.5,
               label = 'Actuator');
pyplot.plot(x2,y2,
               color = 'r',
               alpha = 1,
               label = 'BLC sensor');
pyplot.legend(loc=0);


path1 = '/home/nicolas/Dropbox/Doctorado/Simulaciones/3D/thermal_actuator/fase/phi_1/'
mp_data1 = load_monitoring_points_data(path1 + 'hpts.out')

point_1 = monitoring_point(mp_data1,4,7,2000)
point_2 = monitoring_point(mp_data1,6,7,2000)

x1 = point_1.time - 0.3
y1 = (point_1.field.T - numpy.mean(point_1.field.T))
x2 = point_2.time
y2 = (point_2.field.T - numpy.mean(point_2.field.T))*13



pyplot.xlim(120, 130)
pyplot.ylim(-0.07,0.07)
pyplot.xlabel("Dimensionless time [-]")
pyplot.ylabel("Reduced Temperature [-]")
pyplot.plot(x1,y1,
               color = 'b',
               alpha = 0.5,
               label = 'Actuator');
pyplot.plot(x2,y2,
               color = 'r',
               alpha = 1,
               label = 'BLC sensor');
pyplot.legend(loc=0);
'''
