import numpy as np
import time
import ctypes as c
import os
from keras.models import Sequential, model_from_json
from keras.layers import Dense,Activation
from keras import backend as K


#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_length = abspath+ '/learning_data/data_length.npy'
abspath_model = abspath + '/learned_model/stanford_model.json'
abspath_eval = abspath + '/output/eval.csv'
abspath_answer = abspath + '/output/answer.csv'
abspath_randt = abspath + '/output/randt.csv'
print(abspath_x)


#ctypesの引数設定
fn1 = np.ctypeslib.load_library("make_constant.so",".")
fn1.make_constant_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
]
fn1.make_constant_.restype = c.c_void_p

fn2 = np.ctypeslib.load_library("make_properties.so",".")
fn2.make_properties_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
]
fn2.make_properties_.restype = c.c_void_p


# MTS training dataの読み込み
#print('train_y',train_y)
train_x = np.load(abspath_x)
train_y = np.load(abspath_y)
data_length_float = np.load(abspath_length)
data_length = data_length_float.astype(np.int64)
state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])



#######
#USER DIFINE
input_num     = 10
output_num    = 8
mts_loop      = 1  #base timeからのmts loopの回数
data_num      = 6  #train_x中で初期値とする状態量の番号
start         = 200  #検証を開始するstep数
########

#熱的状態量の算出
dtmp = train_x[start,0]
dprs = train_x[start,1]
aYi = train_x[start,2:10]
totaldens = 0
aeng = 0

dtmp = c.c_double(dtmp)     #ctypes形式でラップ
dprs = c.c_double(dprs)
totaldens = c.c_double(totaldens)
aeng = c.c_double(aeng)

fn1.make_constant_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(totaldens),c.byref(aeng))

totaldens = totaldens.value
aeng = aeng.value

#train_x = train_x[0:200,:]
#train_y = train_y[0:200,:]
train_x = np.delete(train_x,10,1)
train_x = np.delete(train_x,10,1)
train_y = np.delete(train_y,10,1)
train_y = np.delete(train_y,0,1)
train_y = np.delete(train_y,0,1)

print(data_length)
print(train_x[(start+data_length[data_num-1]),:])

#標準化または対数正規化
#Yiのみ対数変換
train_x = np.log(train_x)
#train_y = np.log(train_y)

#standardization
mean_x = np.mean(train_x,axis=0)
#mean_y = np.mean(train_y,axis=0)
std_x = np.std(train_x,axis=0)
#std_y = np.std(train_y,axis=0)
train_x = (train_x - mean_x) / std_x
#train_y = (train_y - mean_y) / std_y


#min-max normalize
#min_x = np.min(train_x,axis=0)
#min_y = np.min(train_y,axis=0)
#max_x = np.max(train_x,axis=0)
#max_y = np.max(train_y,axis=0)
#train_x = (train_x - min_x)/(max_x - min_x)
#train_y = (train_y - min_y)/(max_y - min_y)


t1 = time.time()


#modelの読み込み
abspath_weight = abspath + '/learned_model/stanford_weight.hdf5'
model = model_from_json(open(abspath_model,'r').read())
model.load_weights(abspath_weight)



#evaluation

eval_data = np.zeros([1,input_num])
eval_zeros = np.zeros([1,1])
train = train_x[(start+data_length[data_num-1]),:]

train_int = train.reshape(1,input_num)
eval_moment = model.predict(train_int)

data_range = data_length[data_num] - data_length[data_num-1] - 1
eval_range = data_range/mts_loop - start - 1 #startから最終ステップまでの長さ
eval_range = int(eval_range)


for ii in range(eval_range) :

    #保存量とモル分率から温度、圧力の算出
    dtmp = 0
    dprs = 0
    aYi = eval_moment/np.sum(eval_moment)

    dtmp = c.c_double(dtmp)     #ctypes形式でラップ
    dprs = c.c_double(dprs)
    totaldens = c.c_double(totaldens)
    aeng = c.c_double(aeng)
    #print('aYi',aYi, np.sum(aYi))

    fn2.make_properties_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(totaldens),c.byref(aeng))

    dtmp = dtmp.value
    dprs = dprs.value
    totaldens = totaldens.value
    aeng = aeng.value
    #print('temp',dtmp, 'dprs',dprs)

    eval_append = np.append(eval_zeros,dtmp)
    eval_append = np.append(eval_append,dprs)
    eval_moment = np.append(eval_append,aYi)
    eval_moment = np.delete(eval_moment,0,0)
    eval_next = eval_moment.reshape((1,input_num))
    eval_data = np.concatenate([eval_data,eval_next],axis=0)
    eval_next = np.log(eval_next)
    eval_next = (eval_next - mean_x) / std_x
    #eval_next = eval_moment.reshape((1,input_num))
    #eval_data = np.concatenate([eval_data,eval_next],axis=0)
    eval_moment = model.predict(eval_next)

t2 = time.time()
print('eval time is ',t2-t1)


#標準化または正規化の再変換
#min-max normalize
#train_y = train_y*(max_y - min_y) + min_y
#eval_data = eval_data*(max_y - min_y) + min_y
#train_x = train_x*(max_x - min_x) + min_x

#standardization
#eval_data = std_y * eval_data + mean_y
#train_y = std_y * train_y + mean_y
train_x = std_x * train_x + mean_x


#exp conversion
#train_y = np.exp(train_y)
#eval_data = np.exp(eval_data)
train_x = np.exp(train_x)

#train = train_x[(data_length[data_num]),:]
#print(train)
print(train_x[start,:])
print(train_y[start,:])
#print(eval_data)

#train_data = train_x[start*mts_loop::mts_loop,]
#answer_data = train_y[start*mts_loop::mts_loop,]
#train_data = train_x[start*mts_loop::mts_loop,:]
answer_moment = train_y[start+data_length[data_num-1]::mts_loop,:]
answer_data = answer_moment[0:data_range-start+1,:]
ann_data = np.delete(eval_data,0,axis=0)

np.savetxt(abspath_eval,ann_data,delimiter=',')
#np.savetxt(abspath_eval,train_data,delimiter=',')
np.savetxt(abspath_answer,answer_data,delimiter=',')
#np.savetxt(abspath_randt,train_x,delimiter=',')
