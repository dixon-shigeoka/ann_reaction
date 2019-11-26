import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
import time
from sklearn import preprocessing
#from vtk import *
import vtk
from vtk.util.numpy_support import vtk_to_numpy

#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_length = abspath+ '/learning_data/data_length.npy'
abspath_model = abspath + '/learned_model/stanford_model.json'
abspath_eval = abspath + '/output/x.csv'
abspath_answer = abspath + '/output/y.csv'
abspath_randt = abspath + '/output/randt.csv'
abspath_flow = abspath + '/data/#flow.0103107.vts'
print(abspath_x)

state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])
select  = np.array([0,5])

# load a vtk file as input
reader = vtk.vtkXMLStructuredGridReader()
reader.SetFileName(abspath_flow)
reader.Update()

# Get the coordinates of nodes in the mesh
nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

#The "Temperature" field is the third scalar in my vtk file
temperature_vtk_array = reader.GetOutput().GetPointData().GetArray('Temperature [K]')
pressure_vtk_array    = reader.GetOutput().GetPointData().GetArray('Pressure [MPa]')
H2_vtk_array    = reader.GetOutput().GetPointData().GetArray('  H2')

#Get the coordinates of the nodes and their temperatures
nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

temperature_np = vtk_to_numpy(temperature_vtk_array)
pressure_np    = vtk_to_numpy(pressure_vtk_array)
H2_np    = vtk_to_numpy(pressure_vtk_array)

#make histogram data
temp_min = np.amin(temperature_np)
temp_max = np.amax(temperature_np)
pres_min = np.amin(pressure_np)
pres_max = np.amax(pressure_np)
print(temp_min, temp_max, pres_min, pres_max)
hist_temp = np.histogram(temperature_np,pressure_np,bins=1000,range=[])

#プロット生成

plt.hist(train_x[:,1],bins=200)
plt.savefig("O.png")
plt.close()
plt.hist(train_x[:,0],bins=200)
plt.savefig("temp.png")
plt.close()

