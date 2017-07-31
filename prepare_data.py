#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'

# create a new 'VisItNek5000Reader'
avg_data = VisItNek5000Reader(FileName= 'Data/avghigh_ray.nek5000')
avg_data.Meshes = ['mesh']
avg_data.Materials = []
avg_data.CellArrays = []
avg_data.PointArrays = ['pressure', 'temperature', 'x_velocity', 'y_velocity', 'z_velocity']

rms_data = VisItNek5000Reader(FileName= 'Data/rmshigh_ray.nek5000')
rms_data.Meshes = ['mesh']
rms_data.Materials = []
rms_data.CellArrays = []
rms_data.PointArrays = ['pressure', 'temperature', 'x_velocity', 'y_velocity', 'z_velocity']

# create a new 'Temporal Statistics'
avg_temporalStatistics = TemporalStatistics(Input=avg_data)
avg_temporalStatistics.ComputeAverage = 1
avg_temporalStatistics.ComputeMinimum = 0
avg_temporalStatistics.ComputeMaximum = 0
avg_temporalStatistics.ComputeStandardDeviation = 0

rms_temporalStatistics = TemporalStatistics(Input=rms_data)
rms_temporalStatistics.ComputeAverage = 1
rms_temporalStatistics.ComputeMinimum = 0
rms_temporalStatistics.ComputeMaximum = 0
rms_temporalStatistics.ComputeStandardDeviation = 0

# create a new 'Programmable Filter'
programmableFilter1 = ProgrammableFilter(Input=[avg_data, rms_data])
programmableFilter1.Script ='\
	meanU = inputs[0].PointData["x_velocity"]\n\
	meanUU = inputs[1].PointData["x_velocity"]\n\
	meanV = inputs[0].PointData["y_velocity"]\n\
	meanVV = inputs[1].PointData["y_velocity"]\n\
	meanW = inputs[0].PointData["z_velocity"]\n\
	meanWW = inputs[1].PointData["z_velocity"]\n\
	meanP = inputs[0].PointData["pressure"]\n\
	meanPP = inputs[1].PointData["pressure"]\n\
	meanT = inputs[0].PointData["temperature"]\n\
	meanTT = inputs[1].PointData["temperature"]\n\
	output.PointData.append(meanU,"<u>")\n\
	output.PointData.append(meanV,"<v>")\n\
	output.PointData.append(meanW,"<w>")\n\
	output.PointData.append(meanP,"<P>")\n\
	output.PointData.append(meanP,"<T>")\n\
	output.PointData.append(meanVV-meanV*meanV,"<v\'v\'>")\n\
	output.PointData.append(meanUU-meanU*meanU,"<u\'u\'>")\n\
	output.PointData.append(meanWW-meanW*meanW,"<w\'w\'>")\n\
	output.PointData.append(meanPP-meanP*meanP,"<P\'>")\n\
	output.PointData.append(meanTT-meanT*meanT,"<T\'>")\n'

# save data
#SaveData('Data/avg_data.csv', proxy=avg_temporalStatistics, Precision=5,
#    UseScientificNotation=1,
#    WriteAllTimeSteps=0,
#    FieldAssociation='Points')
#
#SaveData('Data/rms_data.csv', proxy=rms_temporalStatistics, Precision=5,
#    UseScientificNotation=1,
#    WriteAllTimeSteps=0,
#    FieldAssociation='Points')

SaveData('Data/data.csv', proxy=programableFilter1, Precision=5,
    UseScientificNotation=1,
    WriteAllTimeSteps=0,
    FieldAssociation='Points')
