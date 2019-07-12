import numpy as np
import time
import os
from keras.models import Sequential
from keras.layers import Dense,BatchNormalization,Activation, RepeatVector
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import TimeDistributed
from keras import backend as K
from keras.initializers import he_normal
from sklearn.model_selection import KFold, train_test_split


### add for TensorBoard
import keras.callbacks
import keras.backend.tensorflow_backend as KTF
import tensorflow as tf

#old_session = KTF.get_session()
#
#session = tf.Session('')
#KTF.set_session(session)
#KTF.set_learning_phase(1)
###


# fix random seed for reproducibility
seed = 7
np.random.seed(seed)


#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_model = abspath + '/learned_model/stanford_model.json'
abspath_model_temp = abspath + '/learned_model/temp_model.json'
abspath_weight_temp = abspath + '/learned_model/temp_weight.hdf5'
abspath_model_mass = abspath + '/learned_model/mass_model.json'
abspath_weight_mass = abspath + '/learned_model/mass_weight.hdf5'
abspath_tflog = abspath + '/tflog/'
print(abspath_x)

# MTS training dataの読み込み
#print('train_y',train_y)
train_x = np.load(abspath_x)
train_y = np.load(abspath_y)
state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])
#train_y = train_y[:,np.array([0,5])] #H mass fraction
#train_x = np.delete(train_x,10,1) # delete N2 mass fraction
#train_x = np.delete(train_x,10,1) # delete delt
#train_y = np.delete(train_y,10,1) # delete N2 mass fraction
#train_y = np.delete(train_y,0,1) # delete Temperature
#train_y = np.delete(train_y,0,1) # delete Pressure
print(train_x.shape,train_y.shape)


#for validation, Input: t,  Output:Yi
#train_x = np.array(range(3145))


#単精度と倍精度の変更
#K.set_floatx('float32')
#print('K.float is ', K.floatx())
#train_x = K.cast_to_floatx(train_x)
#train_y = K.cast_to_floatx(train_y)
#print('train_x shape is ',train_x.dtype)

#標準化または対数正規化
#log conversion
#train_x_thermo = train_x[:,0:2]
#train_x_thermo = train_x_thermo
#train_x_Yi = np.log(train_x[:,2:11])
#train_x_delt = train_x[:,11].reshape(train_x.shape[0],1)
#train_y_thermo = train_y[:,0:2]
#train_y_thermo = train_y_thermo
#train_y_Yi = np.log(train_y[:,2:11])
#train_x = np.concatenate([train_x_thermo,train_x_Yi],axis=1)
#train_x = np.concatenate([train_x,train_x_delt],axis=1)
#train_y = np.concatenate([train_y_thermo,train_y_Yi],axis=1)
train_x = np.log(train_x)
train_y = np.log(train_y)

#train_x = train_x[20:200,:]
#train_y = train_y[20:200,:]

print(train_x[0:10,0])
print(train_y[0:10,0])

#standardization
mean_x = np.mean(train_x,axis=0)
mean_y = np.mean(train_y,axis=0)
std_x = np.std(train_x,axis=0)
std_y = np.std(train_y,axis=0)
train_x = (train_x - mean_x) / std_x
train_y = (train_y - mean_y) / std_y

#min-max normalize
#min_x = np.min(train_x,axis=0)
#min_y = np.min(train_y,axis=0)
#max_x = np.max(train_x,axis=0)
#max_y = np.max(train_y,axis=0)
#train_x = (train_x - min_x)/(max_x - min_x)
#train_y = (train_y - min_y)/(max_y - min_y)



##自作cost functionの実装
def loss_logadd(y_true,y_pred):
    first_log = K.log(K.abs(y_pred))
    second_log = K.log(K.abs(y_true))
    return  K.mean(K.square(K.abs(y_pred - y_true) + first_log - second_log),axis=-1)

def mix_mse_mape(y_true,y_pred):
    first_log  = K.abs(y_pred)
    second_log = K.abs(y_true)
    return  K.mean(K.square(K.abs(y_pred - y_true) + K.abs(y_true - y_pred)/y_true),axis=-1)

#time record
t1 = time.time()


#####################################
#Hold-Out method for validation
#モデル作成宣言

input_num       = 1
output_num      = 1
state_num       = int(input_num/output_num)
select          = np.array([5])
seq_in_length   = 40
seq_out_length  = 1
n_hidden        = 32
batch_size      = 16
epochs          = 2000

time_length = train_x.shape[0]-seq_in_length+1
X = np.zeros([time_length,seq_in_length,train_x.shape[1]])
Y = np.zeros([time_length,seq_out_length,train_y.shape[1]])
#Y = np.zeros([train_y.shape[0]-seq_length+1,train_y.shape[1]])

print(X.shape,Y.shape)

for i in range(time_length) :
    X[i,:seq_in_length,:] = train_x[i:i+seq_in_length,:]
    Y[i,:seq_out_length,:] = train_y[i:i+seq_out_length,:]

train_x = X[:,:,select] #H mass fraction
state_x = state_x[select]

print(train_x[0:2,:,0])
print(train_y[0:2,0])
print(train_x.shape,Y.shape)

#各状態量に対しての個別のモデルを作成するためのループ
for i in range(state_num) :
    #(train_x,test_x,train_y,test_y) = train_test_split(X,Y,test_size=0.1,shuffle=True,random_state=seed)
    if (state_num == 1) :
        train_y_state = Y[:,:,select]
    else :
        train_y_state = Y[:,:,select[i]] #H mass fraction

    train_y_state = train_y_state.reshape([time_length,1])
    print(train_y_state.shape)


    # define X-fold cross validation
    model = Sequential()

    #パラメータ設定
    #model.add(RepeatVector(seq_inout_length,input_dim=input_num))

    #stateful=Trueではbatch_input_shapeで三次元配列を与える必要あり(batch-size込)
    model.add(LSTM(units=n_hidden,
                   return_sequences=False,
                   stateful=False,
                   batch_input_shape=(None,seq_in_length,input_num)))

    #model.add(TimeDistributed(Dense(units=output_num)))
    model.add(Activation('tanh'))
    model.add(Dense(output_num))
    model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mae'])
    #model.compile(optimizer='adam',loss='mape',metrics=['acc'])
    #model.compile(optimizer='adam',loss=mix_mse_mape,metrics=['acc'])
    model.summary()

    ### make callbacks
    ### add for TensorBoard
    tb_cb = keras.callbacks.TensorBoard(log_dir=abspath_tflog,
                                        #histogram_freq=1,
                                        write_grads=True,
                                        write_images=True
    )
    abspath_weight = abspath + '/learned_model/weight_history/stanford_weights_{epoch:02d}_{loss:.2E}_{val_loss:.2E}.hdf5'
    #cp_cb = keras.callbacks.ModelCheckpoint(filepath=abspath_weight, monitor="val_loss", verbose=0,
    #                                        save_best_only=True
    #)
    es_cb = keras.callbacks.EarlyStopping(monitor='val_loss', patience=100, verbose=1,mode='min')
    #cbks = [tb_cb,cp_cb]
    cbks = [tb_cb]
    ###


    #学習実行
    verbose = 2
    #history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,
    #                    callbacks=cbks,validation_data=(test_x,test_y))
    history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,callbacks=cbks,shuffle=False)

    #評価
    #score = model.evaluate(test_x, test_y, verbose=1)
    score = model.evaluate(train_x, train_y_state, verbose=1)
    print(model.metrics_names[1], score[1])
    #cvscores.append(score[1])
    #score = model.evaluate(train_x,train_y_state,verbose=1)
    print('Train loss:', score[0])

    #modelの保存
    if (state_num != 1) :
        abspath_weight_state = abspath + '/learned_model/%s_weight.hdf5' % state_x[i]
        abspath_model_state  = abspath + '/learned_model/%s_model.json'  % state_x[i]
        open(abspath_model_state,"w").write(model.to_json())
        model.save_weights(abspath_weight_state)
    print('Train accuracy:', score[1])
    #####################################




#####################################
##k-fold cross validation for validation
#
## define K-fold cross validation
#kfold = KFold(n_splits=5, shuffle=True, random_state=seed)
#cvscores = []
#
#for train, test in kfold.split(X, Y):
#
#    #モデル作成宣言
#    model = Sequential()
#
#    #パラメータ設定
#    model.add(Dense(36,input_dim=2))
#    model.add(BatchNormalization())
#    model.add(Activation('relu'))
#    model.add(Dense(36))
#    model.add(BatchNormalization())
#    model.add(Activation('relu'))
#    model.add(Dense(2))
#    model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mae'])
#    #model.compile(optimizer='adam',loss=loss_logadd,metrics=['mae'])
#    model.summary()
#
#
#    ### add for TensorBoard
#    tb_cb = keras.callbacks.TensorBoard(log_dir=abspath_tflog,
#                                        #                                    histogram_freq=1,
#                                        write_grads=True,
#                                        write_images=True
#    )
#    abspath_weight = abspath + '/learned_model/stanford_weights_{epoch:02d}_{loss:.2E}_{val_loss:.2E}.hdf5'
#    cp_cb = keras.callbacks.ModelCheckpoint(filepath=abspath_weight, monitor="val_loss", verbose=1,
#                                            save_best_only=True
#    )
#    cbks = [tb_cb,cp_cb]
#    ###
#
#
#    #学習実行
#    epochs = 10
#    batch_size = 64
#    verbose = 1
#    history = model.fit(X[train],Y[train],batch_size,epochs,verbose,
#                        callbacks=cbks,validation_data=(X[test],Y[test]))
#
#
#    #評価
#    score = model.evaluate(X[test], Y[test], verbose=1)
#    print(model.metrics_names[1], score[1])
#    cvscores.append(score[1])
#####################################


t2 = time.time()
print('errapsed time is ',t2-t1)
#print((np.mean(cvscores), np.std(cvscores)))

#modelの保存
abspath_weight = abspath + '/learned_model/stanford_weight.hdf5'
open(abspath_model,"w").write(model.to_json())
model.save_weights(abspath_weight)

### add for TensorBoard
#KTF.set_session(old_session)
###
