import numpy as np
import time
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
abspath_weight = abspath + '/learned_model/stanford_weight.hdf5'
abspath_eval = abspath + '/output/eval.csv'
abspath_answer = abspath + '/output/answer.csv'
abspath_randt = abspath + '/output/randt.csv'
print(abspath_x)


# MTS training dataの読み込み
#print('train_y',train_y)
train_x = np.load(abspath_x)
train_y = np.load(abspath_y)
data_length_float = np.load(abspath_length)
data_length = data_length_float.astype(np.int64)
#train_x = train_x[:,5] #H mass fraction
#train_y = train_y[:,5] #H mass fraction
train_x = train_x[:,np.array([0,5])] #H mass fraction
train_y = train_y[:,np.array([0,5])] #H mass fraction
#train_x = np.delete(train_x,10,1) # delete N2 mass fraction
#train_x = np.delete(train_x,10,1) # delete delt
#train_y = np.delete(train_y,10,1) # delete N2 mass fraction
#train_y = np.delete(train_y,0,1) # delete Temperature
#train_y = np.delete(train_y,0,1) # delete Pressure
#print(train_x)
#print(train_y)

# for validation (input : time)
#train_x = np.array(range(3145))
#train_x = np.random.rand(3145)*1.0 + 0.5
#rand_t = train_x*3145/np.sum(train_x)
#for i in range(3145) :
#    train_x[i] = np.sum(rand_t[0:i])
#print(train_x)


#標準化または対数正規化
#log conversion
train_x = np.log(train_x)
train_y = np.log(train_y)

#standardization
#mean_x = np.mean(train_x,axis=0)
#mean_y = np.mean(train_y,axis=0)
#std_x = np.std(train_x,axis=0)
#std_y = np.std(train_y,axis=0)
#train_x = (train_x - mean_x) / std_x
#train_y = (train_y - mean_y) / std_y

#min-max normalize
min_x = np.min(train_x,axis=0)
min_y = np.min(train_y,axis=0)
max_x = np.max(train_x,axis=0)
max_y = np.max(train_y,axis=0)
train_x = (train_x - min_x)/(max_x - min_x)
train_y = (train_y - min_y)/(max_y - min_y)


t1 = time.time()

#modelの読み込み
model = model_from_json(open(abspath_model,'r').read())
model.load_weights(abspath_weight)

#for validation (input:time)
#eval_data = model.predict(train_x)



#evaluation

#######
#USER DIFINE
mts_loop = 1 #base timeからのmts loopの回数
data_num = 1 #number of temperature
#data_ini = data_num
start = 200  #検証を開始するstep数
########

#eval_data = np.zeros([1,10])
eval_data = np.zeros([1,2])
#train = train_x[(start+data_length[data_num]),:]
train = train_x[start,:] # if only use 1 property
print(data_length)
#train = train.reshape(1,11) #delt有りの場合
train = train.reshape(1,2) #deltなしの場合
eval_moment = model.predict(train)
data_range = data_length[data_num] - data_length[data_num-1] + 1
eval_range = data_range/mts_loop - start #startから最終ステップまでの長さ
eval_range = int(eval_range)
for ii in range(eval_range) :
    eval_data = np.concatenate([eval_data,eval_moment],axis=0)
    #delt有りの場合
    #eval_next = np.append(eval_moment,1.e-8)
    #eval_next = eval_next.reshape(1,11)
    #deltなしの場合
    eval_next = eval_moment.reshape(1,2)
    #eval_next = eval_moment.reshape(1,1)
    eval_moment = model.predict(eval_next)

t2 = time.time()
print('eval time is ',t2-t1)

#標準化または正規化の再変換
#standardization
#eval_data = std_y * eval_data + mean_y
#train_y = std_y * train_y + mean_y
#train_x = std_x * train_x + mean_x

#min-max normalize
train_y = train_y*(max_y - min_y) + min_y
eval_data = eval_data*(max_y - min_y) + min_y

#exp conversion
eval_data = np.exp(eval_data)
train_y = np.exp(train_y)
train_x = np.exp(train_x)

#train = train_x[(data_length[data_num]),:]
#print(train)
print(eval_data)

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
