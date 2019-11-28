import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.colors
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
abspath_flow = abspath + '/data/#flow.0124310.vts'
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
print('get scalar data')
temperature_vtk_array = reader.GetOutput().GetPointData().GetArray('Temperature [K]')
pressure_vtk_array    = reader.GetOutput().GetPointData().GetArray('Pressure [MPa]')
O_vtk_array    = reader.GetOutput().GetPointData().GetArray('   O')

#Get the coordinates of the nodes and their temperatures
print('change to numpy from vtk')
nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

temperature_np = vtk_to_numpy(temperature_vtk_array)
pressure_np    = vtk_to_numpy(pressure_vtk_array)
O_np    = vtk_to_numpy(O_vtk_array)
print(O_np.size)

#Data scanning and delete columun
print('data scanning')
counter_list = []

for i in range(O_np.size):
    if(O_np[i] < 0.01 and O_np[i] > 0 and temperature_np[i] > 400):
        counter_list.append(i)

filt_temp = np.zeros([1,len(counter_list)])
filt_pres = np.zeros([1,len(counter_list)])

for i in range(len(counter_list)):
    filt_temp[0,i] = temperature_np[counter_list[i]]
    filt_pres[0,i] = pressure_np[counter_list[i]]


#make histogram data
fig = plt.figure()
ax = fig.add_subplot(111)
temp_min = np.amin(filt_temp)
temp_max = np.amax(filt_temp)
pres_min = np.amin(filt_pres)
pres_max = np.amax(filt_pres)
print(temp_min, temp_max, pres_min, pres_max)

x = np.ravel(filt_temp)
y = np.ravel(filt_pres)
hist = ax.hist2d(x,y,bins=100,norm=matplotlib.colors.LogNorm(),cmap=cm.jet)
xlabel = 'Temperature [K]'
ylabel = 'Pressure [MPa]'
xticklabels = ax.get_xticklabels()
yticklabels = ax.get_yticklabels()
#ax.set_xticklabels(xticklabels,fontsize=9,fontname='Times New Roman')
#ax.set_yticklabels(yticklabels,fontsize=9,fontname='Times New Roman')
ax.set_xlabel(xlabel,fontsize=10)
ax.set_ylabel(ylabel,fontsize=10)
fig.colorbar(hist[3],ax=ax)
#ax.set_title(title,fontsize=30,fontname='Times New Roman')

print(hist)

plt.savefig("histogram_plot.png")

#プロット生成

#plt.hist(train_x[:,1],bins=200)
#plt.savefig("O.png")
#plt.close()
#plt.hist(train_x[:,0],bins=200)
#plt.savefig("temp.png")
#plt.close()

