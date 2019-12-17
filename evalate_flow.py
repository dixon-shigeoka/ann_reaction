import numpy as np
import tensorflow as tf
from tensorflow.python.keras.models import load_model
import time
import ctypes as c
import os
from keras.models import Sequential, model_from_json
from keras.layers import Dense,Activation
from keras import backend as K
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import glob

def vtk_write(train_y,ann_data):
    flow_file_list = sorted(glob.glob('/share/data/Planar_flow/vts/#flow*'))
    # load a vtk file as input
    print('reading : ', flow_file_list[0])
    reader = vtk.vtkXMLStructuredGridReader()
    reader.SetFileName(flow_file_list[0])
    reader.Update()

    # get the coordinates of nodes in the mesh
    nodes_vtk_array= reader.GetOutput().GetPoints().GetData()

    #get the coordinates of the nodes
    print('change to numpy from vtk')
    nodes_nummpy_array = vtk_to_numpy(nodes_vtk_array)
    x,y,z= nodes_nummpy_array[:,0] , nodes_nummpy_array[:,1] , nodes_nummpy_array[:,2]
    gridToVTK("./structured",x,y,z,pointData={"temp":train_y[:,0],"ann_temp":ann_data[:,0]})


#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_length = abspath+ '/learning_data/data_length.npy'
abspath_model = abspath + '/learned_model/stanford_model.json'
abspath_eval = abspath + '/output/eval.csv'
abspath_answer = abspath + '/output/answer.csv'
abspath_hist = abspath + '/output/time_hist.csv'
abspath_randt = abspath + '/output/randt.csv'
abspath_stat = abspath + '/cpp/resource/statistics.csv'
print(abspath_x)


#ctypesの引数設定
fn = np.ctypeslib.load_library("library/mtsdriver.so",".")
fn.imtss_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
]
fn.imtss_.restype = c.c_void_p

fn1 = np.ctypeslib.load_library("library/make_constant.so",".")
fn1.make_constant_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
]
fn1.make_constant_.restype = c.c_void_p

fn2 = np.ctypeslib.load_library("library/make_properties.so",".")
fn2.make_properties_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
]
fn2.make_properties_.restype = c.c_void_p

fn3 = np.ctypeslib.load_library("library/pidriver.so",".")
fn3.pointimplicit_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
]
fn3.pointimplicit_.restype = c.c_void_p

# training dataの読み込み
#print('train_y',train_y)
train_x = np.load(abspath_x)
train_y = np.load(abspath_y)
train_x[:,1] = train_x[:,1] * 1e6
data_length_float = np.load(abspath_length)
data_length = data_length_float.astype(np.int64)
state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])



#######
#USER DIFINE
input_num     = 10
output_num    = 8
mts_loop      = 1  #base timeからのmts loopの回数
data_num      = 1  #train_x中で初期値とする状態量の番号
start         = 0  #検証を開始するstep数
########


#------------------------------------
#初期状態の決定と熱的状態量の算出
#------------------------------------
train_int_zeros = np.zeros([1,1])
#dtmp = train_x[start,0]
#dprs = train_x[start,1]
dtmp_init = 3646.3380148366131834337
dprs_init = 2029699.3854148103855550289
#dtmp_init =1355
#dprs_init = 1.01325e5 * 1.5

aMi = np.array([2.016,32.000,1.008,16.000,17.008,18.016,33.008,34.016,28.016])
aMolarRatio = np.array([2,1,0,0,0,0,0,0,0])
#aMolarRatio = np.where(aMolarRatio == 0,aemn,aMolarRatio)
#aMolarRatio[8] = 0
aTotalMass = np.dot(aMi,aMolarRatio)
aMixMolarRatio = np.multiply(aMi,aMolarRatio)
aYi = aMixMolarRatio/aTotalMass
aYi = np.array([0.0205572317184883613,0.0989824859011038255,0.0045745792559854330,0.0345493752659215073,0.1533132157188197564,
                0.6875744446241545127,0.0004020743187675550,0.0000465931967590170,0])
aYi_init  = aYi.reshape((1,9))
#aYi_init = train_x[start,2:11]

train_int_append = np.append(train_int_zeros,dtmp_init)
train_int_append = np.append(train_int_append,dprs_init)
train_int_append = np.append(train_int_append,aYi_init)
train_int_moment = np.delete(train_int_append,0,0)
train_int_moment = np.delete(train_int_moment,-1,0)
train_int_moment = train_int_moment.reshape((1,10))

#dtmp = c.c_double(dtmp_init)     #ctypes形式でラップ
#dprs = c.c_double(dprs_init)
#totaldens = c.c_double(totaldens)
#aeng = c.c_double(aeng)
#
#fn1.make_constant_(c.byref(dtmp),c.byref(dprs),aYi_init,c.byref(totaldens),c.byref(aeng))
#
#totaldens = totaldens.value
#aeng = aeng.value

#train_x = train_x[0:200,:]
#train_y = train_y[0:200,:]
#train_x = np.delete(train_x,10,1)
#train_x = np.delete(train_x,10,1)
#train_y = np.delete(train_y,10,1)
#train_y = np.delete(train_y,0,1)
#train_y = np.delete(train_y,0,1)

#print('test data')
#print(train_int_moment)

#------------------------------------

#------------------------------------
#標準化または対数正規化
#------------------------------------

#train_x_log = np.log(train_x)

#standardization
#mean_x = np.mean(train_x_log,axis=0)
#std_x = np.std(train_x_log,axis=0)
statistics = np.loadtxt(abspath_stat,dtype='float',delimiter=',')
mean_x = statistics[0,:]
std_x  = statistics[1,:]

#train_x = (train_int_moment - mean_x) / std_x

#------------------------------------

print(train_y.shape)
#train_x = train_x[7602,:]
#train_y = train_y[7602,:]
#train_x = train_x.reshape(1,10)
#train_y = train_y.reshape(1,10)

#train_y_mts = np.zeros([1,11],dtype=float)
#------------------------------------
#正当データの復元
#------------------------------------
#for i in range (train_x.shape[0]):
#    train_x_zeros = np.zeros([1,1])
#    train_y_zeros = np.zeros([1,1])
#    #train_x = np.zeros([1,12],dtype=float)
#    equiv_error = 0
#    counter = 0
#    #dtmp = dtmp_init
#    #dprs = dprs_init
#    #aYi  = aYi_init
#    dtmp = train_x[i,0]
#    dprs = train_x[i,1]
#    aYi  = train_x[i,2:10]
#    aYi  = np.append(aYi,0)
#    delt_mts = 2.e-10
#    time_hist_ode = np.zeros([1,1])
#    time_hist_pred = np.zeros([1,1])
#    t_diff        = np.zeros([1,1])
#
#    #train_x_append = np.append(train_x_zeros,dtmp)
#    #train_x_append = np.append(train_x_append,dprs)
#    #train_x_append = np.append(train_x_append,aYi)
#
#    #入力データを格納
#    temp_before = dtmp
#    #train_x_append_in_mtsloop = np.append(train_x_append,delt_mem)
#    #train_x_moment = np.delete(train_x_append_in_mtsloop,0,0)
#    #train_x_moment = train_x_moment.reshape((1,12))
#    #train_x = np.concatenate([train_x,train_x_moment],axis=0)
#
#    dtmp = c.c_double(dtmp)     #ctypes形式でラップ
#    dprs = c.c_double(dprs)
#    #delt = c.c_double(delt_mem)
#    delt = c.c_double(delt_mts)
#
#    #------------------------------------
#    # prediction
#    #------------------------------------
#    tbefore = time.time()
#
#    #byrefでポインタにして渡す，mtsの実行
#    fn.imtss_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(delt))
#    #fn3.pointimplicit_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(delt))
#
#    tafter = time.time()
#    t_diff[0,0] = tafter - tbefore
#    time_hist_ode = np.concatenate([time_hist_ode,t_diff],axis=0)
#    #------------------------------------
#
#    dtmp = dtmp.value
#    dprs = dprs.value
#    delt = delt.value
#
#    #MTS結果を教師データとして格納
#    train_y_append = np.append(train_y_zeros,dtmp)
#    train_y_append = np.append(train_y_append,dprs)
#    train_y_append = np.append(train_y_append,aYi)
#    train_y_moment = np.delete(train_y_append,0,0)
#    train_y_moment = train_y_moment.reshape((1,11))
#    train_y_mts = np.concatenate([train_y_mts,train_y_moment],axis=0)
#
##train_x = np.delete(train_x,0,0)
#train_y_mts = np.delete(train_y_mts,0,0)
#print('train_x')
#print(train_x)
#print('train_y')
#print(train_y)
#print('train_y_mts')
#print(train_y_mts)

#------------------------------------


#------------------------------------
#------------------------------------
#evaluation
#------------------------------------
#------------------------------------

#modelの読み込み
#abspath_weight = abspath + '/learned_model/stanford_weight.hdf5'
#model = model_from_json(open(abspath_model,'r').read())
#model.load_weights(abspath_weight)
model = load_model(abspath + '/learned_model/stanford_model.h5')

eval_data = np.zeros([1,input_num])
eval_zeros = np.zeros([1,1])
#train = train_x[(start+data_length[data_num-1]),:]

#standardization
train_int_moment = np.log(train_int_moment)
train = (train_int_moment - mean_x) / std_x

#print(train_x[(start+data_length[data_num-1]),:])
#print(train)

train_int = train.reshape(1,input_num)

print('PREDICTION')
# threshold for start chemical reaction
threshold    = (np.log(500) - mean_x[0]) / std_x[0]
aemn         = np.full(10,1e-20)

for i in range(train_x.shape[0]):

    #------------------------------------
    # make conserved properties
    #------------------------------------
    dtmp = c.c_double(train_x[i,0])     #ctypes形式でラップ
    dprs = c.c_double(train_x[i,1])
    aYi  = train_x[i,2:10]
    aYi  = np.append(aYi,0)
    totaldens = 0
    aeng = 0
    totaldens = c.c_double(totaldens)
    aeng = c.c_double(aeng)

    fn1.make_constant_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(totaldens),c.byref(aeng))

    totaldens = totaldens.value
    aeng = aeng.value

    #------------------------------------
    # prediction
    #------------------------------------
    # preprocess data
#    tbefore = time.time()

    train_int = train_x[i,:]
    train_int = train_int.reshape(1,input_num)

    if (train_int[0,0] > threshold) :
        #print('o', i, train_int)
        train_int[0,:] = (np.log(train_int[0,:] + aemn) - mean_x) / std_x
        eval_moment = model.predict(train_int)
    else :
        eval_moment = train_int[0,2:10]
        print('x')

#    tafter = time.time()
#    t_diff[0,0] = tafter - tbefore
    #------------------------------------

    #保存量とモル分率から温度、圧力の算出
    aYi = eval_moment/np.sum(eval_moment)

    #dtmp = c.c_double(dtmp)     #ctypes形式でラップ
    #dprs = c.c_double(dprs)
    totaldens = c.c_double(totaldens)
    aeng = c.c_double(aeng)
    #print('aYi',aYi, np.sum(aYi))

    #------------------------------------
    # make thermodynamic properties from conservations
    #------------------------------------
#    tbefore = time.time()

    fn2.make_properties_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(totaldens),c.byref(aeng))

#    tafter = time.time()
#    t_diff[0,0] = t_diff[0,0] + tafter - tbefore
#    time_hist_pred = np.concatenate([time_hist_pred,t_diff],axis=0)
#    t_diff[0,0] = 0
    #------------------------------------

    dtmp = dtmp.value
    dprs = dprs.value
    eval_append = np.append(eval_zeros,dtmp)
    eval_append = np.append(eval_append,dprs)
    eval_moment = np.append(eval_append,aYi)
    eval_moment = np.delete(eval_moment,0,0)
    eval_next = eval_moment.reshape((1,input_num))
    eval_data = np.concatenate([eval_data,eval_next],axis=0)

ann_data = np.delete(eval_data,0,axis=0)
vtk_write(train_y,ann_data)
print('ann_data')
print(ann_data[0,:])
#------------------------------------
