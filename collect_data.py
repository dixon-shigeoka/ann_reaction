import numpy as np
import ctypes as c
import os
import time
import subprocess
import glob
import vtk
from vtk.util.numpy_support import vtk_to_numpy


#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_length = abspath + '/learning_data/data_length.npy'
abspath_library = abspath + '/library'

seed = 7
np.random.seed(seed)

#ctypesの引数設定
fn = np.ctypeslib.load_library("library/mtsdriver.so",".")
fn.imtss_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
]
fn.imtss_.restype = c.c_void_p

#fn2 = np.ctypeslib.load_library("1dsteady.so",".")
#fn2.steadystanford_.argtypes = [
#    c.POINTER(c.c_double),
#    c.POINTER(c.c_double),
#    np.ctypeslib.ndpointer(dtype=np.float64),
#    c.POINTER(c.c_double),
#]
#fn2.steadystanford_.restype = c.c_void_p

fn3 = np.ctypeslib.load_library("library/pidriver.so",".")
fn3.pointimplicit_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
]
fn3.pointimplicit_.restype = c.c_void_p


def Combustion():
    # set the condition
    aemn = 1.e-20
    delt_base = 1.e-9 # base timeの時間刻み
    delt_mts  = 1.e-9 # mts loopの時間刻み
    range_mts = 1     # mts loopを回す回数
    list_num  = 9     # 要素数
    tcounter  = 0     # 温度を変化させるループのカウンター
    #total_num = 1   # 合計値
    train_x = np.zeros([1,12],dtype=float)
    train_y = np.zeros([1,11],dtype=float)
    omega_ave_data = np.zeros([1,9],dtype=float)
    omega_ave = np.zeros([1,9],dtype=float)
    train_x_zeros = np.zeros([1,1])
    train_y_zeros = np.zeros([1,1])
    data_length = np.zeros([1,1])
    atmp = 0
    counter = 0     # データ数カウンター

    for i in range(1):  #質量分率，温度，密度を変化させるループ(base time)

        data_length = np.append(data_length,counter)
        t1 = time.time()
        #atmp = 3000 + float(tcounter*25)
        #atmp = 1663 + float(tcounter*10)
        atmp = 1763
        dprs = 3343232
        #dprs = 2000000
        #dprs = 1e5 + float(tcounter*25)
        #dprs =1.01325e5 * 1.5
        aMi = np.array([2.016,32.000,1.008,16.000,17.008,18.016,33.008,34.016,28.016])
        aMolarRatio = np.array([2,1,0,0,0,0,0,0,0])
        aMolarRatio = np.where(aMolarRatio == 0,aemn,aMolarRatio)
        aMolarRatio[8] = 0
        aTotalMass = np.dot(aMi,aMolarRatio)
        aMixMolarRatio = np.multiply(aMi,aMolarRatio)
        aYi = aMixMolarRatio/aTotalMass
        aYi  = aYi.reshape((1,9))
        equiv_error = 1.e-1
        tcounter = tcounter + 1
        print('data number is ', counter)
        print('tcounter is ', tcounter)

        print(atmp,dprs,aYi)

        while (equiv_error > 1.E-4 or atmp < 2000) :     #時間方向にmts計算を行い訓練データを格納するループ(base timeからの変化)
            #while (aYi[0,2] < 2e-5) :

            counter = counter + 1
            delt_mem = 0.e0

            train_x_append = np.append(train_x_zeros,atmp)
            train_x_append = np.append(train_x_append,dprs)
            train_x_append = np.append(train_x_append,aYi)


            for ii in range(range_mts):    #deltをパラメータにする分だけmtsを回す

                #delt_rand = (np.random.rand(1)*0.5 + 0.5)*delt_mts
                delt_rand = delt_mts
                delt_mem = delt_mem + delt_rand #このループ内で回った時間

                #入力データを格納
                temp_before = atmp
                train_x_append_in_mtsloop = np.append(train_x_append,delt_mem)
                train_x_moment = np.delete(train_x_append_in_mtsloop,0,0)
                train_x_moment = train_x_moment.reshape((1,12))
                train_x = np.concatenate([train_x,train_x_moment],axis=0)

                atmp = c.c_double(atmp)     #ctypes形式でラップ
                dprs = c.c_double(dprs)
                delt = c.c_double(delt_mem)

                #byrefでポインタにして渡す，mtsの実行
                fn.imtss_(c.byref(atmp),c.byref(dprs),aYi,c.byref(delt))
                #fn3.pointimplicit_(c.byref(atmp),c.byref(dprs),aYi,c.byref(delt))
                atmp = atmp.value
                dprs = dprs.value
                delt = delt.value

                #MTS結果を教師データとして格納
                train_y_append = np.append(train_y_zeros,atmp)
                train_y_append = np.append(train_y_append,dprs)
                train_y_append = np.append(train_y_append,aYi)
                train_y_moment = np.delete(train_y_append,0,0)
                train_y_moment = train_y_moment.reshape((1,11))
                train_y = np.concatenate([train_y,train_y_moment],axis=0)
                #omega_ave_data = np.vstack((omega_ave_data,omega_ave))
                #print(train_x_append,atmp)


                atmp = train_x_append[1]
                dprs = train_x_append[2]
                aYi  = train_x_append[3:12]
                aYi  = aYi.reshape((1,9))
                atmp = c.c_double(atmp)     #ctypes形式でラップ
                dprs = c.c_double(dprs)
                delt = c.c_double(delt_base)
                #H2_before = aYi[0,0]

                #byrefでポインタにして渡す，mtsの実行
                fn.imtss_(c.byref(atmp),c.byref(dprs),aYi,c.byref(delt))
                #fn3.pointimplicit_(c.byref(atmp),c.byref(dprs),aYi,c.byref(delt))
                atmp = atmp.value
                dprs = dprs.value
                delt = delt.value

                #equiv_error = abs(aYi[0,0] - H2_before)
                equiv_error = abs(atmp - temp_before)


    t2 = time.time()
    print('mts loop time is ', t2-t1)
    counter = counter + 1
    data_length = np.append(data_length,counter)


    # OUTPUT
    data_length = np.delete(data_length,0,0)
    print(data_length)
    train_x = np.delete(train_x,0,0)
    train_y = np.delete(train_y,0,0)
    train_x = np.delete(train_x,10,1) # delete N2 mass fraction
    train_x = np.delete(train_x,10,1) # delete delt
    train_y = np.delete(train_y,10,1) # delete N2 mass fraction
    train_y = np.delete(train_y,0,1) # delete Temperature
    train_y = np.delete(train_y,0,1) # delete Pressure
    print(train_x.shape, train_y.shape)
    #omega_ave = np.delete(omega_ave,0,0)
    np.save(abspath_x,train_x)
    np.save(abspath_y,train_y)
    np.save(abspath_length,data_length)
    abspath_temp = abspath + '/learning_data/temp_%s.csv' % delt_mts
    np.savetxt(abspath_temp,train_x[:,0],delimiter=',')


def SteadyDetonation():

    train_x = np.zeros([1,12],dtype=float)
    train_y = np.zeros([1,11],dtype=float)
    train_x_zeros = np.zeros([1,1])
    train_y_zeros = np.zeros([1,1])
    data_length = np.zeros([1,1])
    atmp = 0
    shock_speed = 2836
    with open('com','w') as f:
        print(shock_speed,file=f)
    subprocess.call(['/share/library/1dsteady.out'])
    test_x = np.loadtxt('/share/1dsteady.csv',delimiter=',')
    print(test_x.shape)

def FlowFileConverter():

    flow_file_list = sorted(glob.glob('/share/data/Planar_flow/vts/#flow*'))
    flow_file_list = sorted(glob.glob('/share/data/Planar_flow/vts/#flow.0150206.vts'))
    before_file_list = sorted(glob.glob('/share/data/Planar_flow/vts/#before*'))
    before_file_list = sorted(glob.glob('/share/data/Planar_flow/vts/#before.0150206.vts'))
    print(flow_file_list)
    print(before_file_list)
    train_x    = np.zeros([1,10])
    train_y    = np.zeros([1,10])

    for i_file in range(len(before_file_list)):
        # load a vtk file as input
        print('reading : ', before_file_list[i_file])
        reader = vtk.vtkXMLStructuredGridReader()
        reader.SetFileName(before_file_list[i_file])
        reader.Update()

        # get the coordinates of nodes in the mesh
        nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

        #the "temperature" field is the third scalar in my vtk file
        print('get scalar data')
        temperature_vtk_array = reader.GetOutput().GetPointData().GetArray('Temperature [K]')
        pressure_vtk_array    = reader.GetOutput().GetPointData().GetArray('Pressure [MPa]')
        h2_vtk_array          = reader.GetOutput().GetPointData().GetArray('  H2')
        o2_vtk_array          = reader.GetOutput().GetPointData().GetArray('  O2')
        h_vtk_array           = reader.GetOutput().GetPointData().GetArray('   H')
        o_vtk_array           = reader.GetOutput().GetPointData().GetArray('   O')
        oh_vtk_array          = reader.GetOutput().GetPointData().GetArray('  OH')
        h2o_vtk_array         = reader.GetOutput().GetPointData().GetArray(' H2O')
        ho2_vtk_array         = reader.GetOutput().GetPointData().GetArray(' HO2')
        h2o2_vtk_array        = reader.GetOutput().GetPointData().GetArray('H2O2')

        #get the coordinates of the nodes and their temperatures
        print('change to numpy from vtk')
        nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
        x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

        temperature_np = vtk_to_numpy(temperature_vtk_array)
        pressure_np    = vtk_to_numpy(pressure_vtk_array)
        h2_np          = vtk_to_numpy(h2_vtk_array)
        o2_np          = vtk_to_numpy(o2_vtk_array)
        h_np           = vtk_to_numpy(h_vtk_array)
        o_np           = vtk_to_numpy(o_vtk_array)
        oh_np          = vtk_to_numpy(oh_vtk_array)
        h2o_np         = vtk_to_numpy(h2o_vtk_array)
        ho2_np         = vtk_to_numpy(ho2_vtk_array)
        h2o2_np        = vtk_to_numpy(h2o2_vtk_array)

        temperature_np = temperature_np.reshape(len(temperature_np),1)
        pressure_np    = pressure_np.reshape(len(temperature_np),1)
        h2_np          = h2_np.reshape(len(temperature_np),1)
        o2_np          = o2_np.reshape(len(temperature_np),1)
        h_np           = h_np.reshape(len(temperature_np),1)
        o_np           = o_np.reshape(len(temperature_np),1)
        oh_np          = oh_np.reshape(len(temperature_np),1)
        h2o_np         = h2o_np.reshape(len(temperature_np),1)
        ho2_np         = ho2_np.reshape(len(temperature_np),1)
        h2o2_np        = h2o2_np.reshape(len(temperature_np),1)
        #x              = x.reshape(len(temperature_np),1)

        train_x_tmp = np.concatenate([temperature_np,pressure_np ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,h2_np      ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,o2_np      ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,h_np       ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,o_np       ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,oh_np      ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,h2o_np     ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,ho2_np     ],axis=1)
        train_x_tmp = np.concatenate([train_x_tmp    ,h2o2_np    ],axis=1)
        #train_x_tmp = np.concatenate([train_x_tmp    ,x          ],axis=1)

        train_x  = np.vstack([train_x,train_x_tmp])
        print(train_x.shape)


    for i_file in range(len(flow_file_list)):
        # load a vtk file as input
        print('reading : ', flow_file_list[i_file])
        reader = vtk.vtkXMLStructuredGridReader()
        reader.SetFileName(flow_file_list[i_file])
        reader.Update()

        # get the coordinates of nodes in the mesh
        nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

        #the "temperature" field is the third scalar in my vtk file
        print('get scalar data')
        temperature_vtk_array = reader.GetOutput().GetPointData().GetArray('Temperature [K]')
        pressure_vtk_array    = reader.GetOutput().GetPointData().GetArray('Pressure [MPa]')
        h2_vtk_array          = reader.GetOutput().GetPointData().GetArray('  H2')
        o2_vtk_array          = reader.GetOutput().GetPointData().GetArray('  O2')
        h_vtk_array           = reader.GetOutput().GetPointData().GetArray('   H')
        o_vtk_array           = reader.GetOutput().GetPointData().GetArray('   O')
        oh_vtk_array          = reader.GetOutput().GetPointData().GetArray('  OH')
        h2o_vtk_array         = reader.GetOutput().GetPointData().GetArray(' H2O')
        ho2_vtk_array         = reader.GetOutput().GetPointData().GetArray(' HO2')
        h2o2_vtk_array        = reader.GetOutput().GetPointData().GetArray('H2O2')

        #get the coordinates of the nodes and their temperatures
        print('change to numpy from vtk')
        nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
        x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]

        temperature_np = vtk_to_numpy(temperature_vtk_array)
        pressure_np    = vtk_to_numpy(pressure_vtk_array)
        h2_np          = vtk_to_numpy(h2_vtk_array)
        o2_np          = vtk_to_numpy(o2_vtk_array)
        h_np           = vtk_to_numpy(h_vtk_array)
        o_np           = vtk_to_numpy(o_vtk_array)
        oh_np          = vtk_to_numpy(oh_vtk_array)
        h2o_np         = vtk_to_numpy(h2o_vtk_array)
        ho2_np         = vtk_to_numpy(ho2_vtk_array)
        h2o2_np        = vtk_to_numpy(h2o2_vtk_array)

        temperature_np = temperature_np.reshape(len(temperature_np),1)
        pressure_np    = pressure_np.reshape(len(temperature_np),1)
        h2_np          = h2_np.reshape(len(temperature_np),1)
        o2_np          = o2_np.reshape(len(temperature_np),1)
        h_np           = h_np.reshape(len(temperature_np),1)
        o_np           = o_np.reshape(len(temperature_np),1)
        oh_np          = oh_np.reshape(len(temperature_np),1)
        h2o_np         = h2o_np.reshape(len(temperature_np),1)
        ho2_np         = ho2_np.reshape(len(temperature_np),1)
        h2o2_np        = h2o2_np.reshape(len(temperature_np),1)
        #x              = x.reshape(len(temperature_np),1)

        train_y_tmp = np.concatenate([temperature_np,pressure_np ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,h2_np      ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,o2_np      ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,h_np       ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,o_np       ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,oh_np      ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,h2o_np     ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,ho2_np     ],axis=1)
        train_y_tmp = np.concatenate([train_y_tmp    ,h2o2_np    ],axis=1)
        #train_y_tmp = np.concatenate([train_y_tmp    ,x      ],axis=1)

        train_y = np.vstack([train_y,train_y_tmp])
        print(train_y.shape)


    train_x = np.delete(train_y,0,axis=0)
    train_y = np.delete(train_y,0,axis=0)

    #------------------------
    #pick up for learning
    #------------------------

    pickup = 'off'

    if (pickup == 'on'):
        #Data scanning and delete columun
        print('data scanning')
        counter_list = []
        train_x = np.delete(train_x    , 0,axis=0)
        train_y = np.delete(train_y    , 0,axis=0)

        #exclude zero
        print('exclude zero')
        train_x_tmp = np.delete(train_x,np.where(train_x == 0), axis=0)
        train_y_tmp = np.delete(train_y,np.where(train_x == 0), axis=0)
        train_x     = train_x_tmp
        train_y     = train_y_tmp
        print(train_x_tmp.shape, train_y_tmp.shape)

        #separate equil. and reac. zone
        print('separate zone : equil. and reac.')
        train_y_eq = np.delete(train_y_tmp,np.where(train_x_tmp[:,2] > 0.025), axis=0)
        train_x_eq = np.delete(train_x_tmp,np.where(train_x_tmp[:,2] > 0.025), axis=0)
        train_y_re = np.delete(train_y_tmp,np.where(train_x_tmp[:,2] <=0.025), axis=0)
        train_x_re = np.delete(train_x_tmp,np.where(train_x_tmp[:,2] <=0.025), axis=0)
        train_y_re = np.delete(train_y_re ,np.where(train_x_re[:,0]  <   400), axis=0)
        train_x_re = np.delete(train_x_re ,np.where(train_x_re[:,0]  <   400), axis=0)
        print(train_x_eq.shape, train_x_re.shape)

        #extract random rows from equil. zone
        print('extract random rows from equil. zone')
        row_i  = np.random.choice(np.arange(train_x_eq.shape[0]),size=train_x_re.shape[0],replace=False)
        row_i.sort()
        train_x_tmp = train_x_eq[row_i,:]
        train_y_tmp = train_y_eq[row_i,:]
        print(train_x_tmp.shape)

        #concatenate equil. and reac. zone
        train_x = np.vstack([train_x_tmp,train_x_re])
        train_y = np.vstack([train_y_tmp,train_y_re])

        tmp_diff = np.zeros([train_x.shape[0],2])
        for i in range(train_x.shape[0]):
            tmp_diff[i,0] = i
            tmp_diff[i,1] = abs(train_x[i,0]-train_y[i,0])

            col_num = 1
            tmp_diff_sort = tmp_diff[np.argsort(tmp_diff[:,col_num])]

    #delete column
    #train_x = np.delete(train_x_tmp,-1,axis=1)
    train_y = np.delete(train_y,0,axis=1)
    train_y = np.delete(train_y,0,axis=1)
    #print(train_x.shape,train_y.shape)
    #print(train_x.dtype,train_y.dtype)
    print(train_x[1,:],train_y[1,:])

    print(train_x.shape,train_y.shape)
    np.save(abspath_x,train_x)
    np.save(abspath_y,train_y)
    #train_x = np.hstack([train_x,train_y])
    abspath_csv = abspath + '/learning_data/train_y.csv'
    np.savetxt(abspath_csv,train_y,delimiter=',')

if __name__ == '__main__':

    #choose your object
    #Combustion()
    #SteadyDetonation()
    FlowFileConverter()
