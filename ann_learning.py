import numpy as np
import time
import os
import shutil
import scipy.io
from keras.models import Sequential
from keras.layers import Dense,BatchNormalization,Activation,Dropout
from keras import backend as K
from keras.initializers import he_normal
from tensorflow.python.framework.graph_util import convert_variables_to_constants
import keras.backend.tensorflow_backend as KTF
import tensorflow as tf
from sklearn.model_selection import KFold, train_test_split


### add for TensorBoard
import keras.callbacks
import keras.backend.tensorflow_backend as KTF
import tensorflow as tf

old_session = KTF.get_session()
sess = tf.Session('')
KTF.set_session(sess)
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
abspath_weight = abspath + '/learned_model/weight_history/'
abspath_tflog = abspath + '/tflog/'
print(abspath_weight)

#履歴の削除
shutil.rmtree(abspath_weight)
os.mkdir(abspath_weight)

# MTS training dataの読み込み
#print('train_y',train_y)
train_x = np.load(abspath_x)
train_y = np.load(abspath_y)
state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])
#train_y = train_y[:,np.array([0,5])] #H mass fraction
train_x = np.delete(train_x,10,1) # delete N2 mass fraction
train_x = np.delete(train_x,10,1) # delete delt
train_y = np.delete(train_y,10,1) # delete N2 mass fraction
train_y = np.delete(train_y,0,1) # delete Temperature
train_y = np.delete(train_y,0,1) # delete Pressure
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
#train_y = np.log(train_y)

#train_x = train_x[0:200,:]
#train_y = train_y[0:200,:]

#standardization
mean_x = np.mean(train_x,axis=0)
#mean_y = np.mean(train_y,axis=0)
std_x = np.std(train_x,axis=0)
#std_y = np.std(train_y,axis=0)
train_x = (train_x - mean_x) / std_x
#train_y = (train_y - mean_y) / std_y
mean = mean_x.reshape(1,10)
std = std_x.reshape(1,10)
statistics_x = np.concatenate([mean,std],axis=0)
print(statistics_x)

#saving statistic data for c++
abspath_stat = abspath + '/cpp/resource/statistics.csv'
abspath_mean = abspath + '/cpp/resource/mean.csv'
abspath_std = abspath + '/cpp/resource/std.csv'
np.savetxt(abspath_stat,statistics_x,delimiter=',')
np.savetxt(abspath_mean,mean_x,delimiter=',')
np.savetxt(abspath_std,std_x,delimiter=',')

#min-max normalize
#min_x = np.min(train_x,axis=0)
#min_y = np.min(train_y,axis=0)
#max_x = np.max(train_x,axis=0)
#max_y = np.max(train_y,axis=0)
#train_x = (train_x - min_x)/(max_x - min_x)
#train_y = (train_y - min_y)/(max_y - min_y)

X = train_x
Y = train_y

##自作cost functionの実装
def loss_log(y_true,y_pred):
    first_log = K.log(K.abs(y_pred))
    second_log = K.log(K.abs(y_true))
    #return  K.mean(K.square(K.abs(first_log - second_log)),axis=-1)
    return  K.mean(K.square(first_log - second_log),axis=-1)

def loss_logaddmse(y_true,y_pred):
    first_log = K.log(K.abs(y_pred))
    second_log = K.log(K.abs(y_true))
    #return  K.mean(K.square(K.abs(y_pred - y_true) + K.log(y_pred/y_true)),axis=0)
    #return  K.mean(K.square(K.abs(y_pred - y_true) + tf.div(y_pred,y_true)),axis=-1)
    return  K.mean(K.square(K.abs((y_pred - y_true) + first_log - second_log)),axis=-1)

def mix_mse_mape(y_true,y_pred):
    first_log  = K.abs(y_pred)
    second_log = K.abs(y_true)
    return  K.mean(K.square(K.abs(y_pred - y_true) + K.abs(y_pred - y_true)/y_true),axis=-1)

#time record
t1 = time.time()


#####################################
#Hold-Out method for validation
#モデル作成宣言

input_num = 10
output_num = 8
state_num = int(input_num/output_num)
#select = np.array([5])

#train_x = train_x[:,select] #H mass fraction
#state_x = state_x[select]


#各状態量に対しての個別のモデルを作成するためのループ
for i in range(state_num) :
    #(train_x,test_x,train_y,test_y) = train_test_split(X,Y,test_size=0.1,shuffle=True,random_state=seed)
    if (output_num == 1) :
        if (state_num == 1) :
            train_y_state = train_y[:,select]
        else :
            train_y_state = train_y[:,select[i]] #H mass fraction
    else :
        train_y_state = train_y

    #print(train_y_state)

    model = Sequential()

    #パラメータ設定
    #model.add(InputLayer(input_shape=(input_num,)))
    model.add(Dense(30,input_dim=input_num))
    #model.add(BatchNormalization())
    model.add(Activation('relu'))
    #model.add(Dropout(0.5))
    model.add(Dense(30))
    #model.add(BatchNormalization())
    model.add(Activation('relu'))
    #model.add(Dropout(0.5))
    #model.add(Dense(30))
    #model.add(BatchNormalization())
    #model.add(Activation('relu'))
    #model.add(Dropout(0.5))
    model.add(Dense(output_num))
    model.add(Activation('sigmoid'))
    #model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mae'])
    #model.compile(optimizer='adam',loss='mean_squared_logarithmic_error',metrics=['mae'])
    model.compile(optimizer='adam',loss=loss_log,metrics=['mae'])
    #model.compile(optimizer='adam',loss=loss_logaddmse,metrics=['mae'])
    #model.compile(optimizer='adam',loss=mix_mse_mape,metrics=['mae'])
    model.summary()

    for layer in model.layers:
        print(layer.output_shape)

    ### make callbacks
    ### add for TensorBoard
    tb_cb = keras.callbacks.TensorBoard(log_dir=abspath_tflog,
                                        #histogram_freq=1,
                                        write_grads=True,
                                        write_images=True
    )
    abspath_weight = abspath + '/learned_model/weight_history/stanford_weights_{epoch:02d}_{loss:.2e}_{val_loss:.2e}.hdf5'
    cp_cb = keras.callbacks.ModelCheckpoint(filepath=abspath_weight, monitor="val_loss", verbose=0,
                                            save_best_only=True
    )
    #es_cb = keras.callbacks.EarlyStopping(monitor='val_loss', patience=200, verbose=1,mode='min')
    es_cb = keras.callbacks.EarlyStopping(monitor='val_loss', patience=2000, verbose=1,mode='min')
    cbks = [tb_cb,cp_cb,es_cb]
    #cbks = [tb_cb]
    ###


    #学習実行
    epochs = 100000
    batch_size = 2048
    verbose = 2
    #history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,callbacks=cbks)
    history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,callbacks=cbks,validation_split=0.1)

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
#    model.add(Dense(30,input_dim=input_num))
#    #model.add(BatchNormalization())
#    model.add(Activation('relu'))
#    #model.add(Dropout(0.1))
#    model.add(Dense(30))
#    #model.add(BatchNormalization())
#    model.add(Activation('relu'))
#    #model.add(Dropout(0.1))
#    model.add(Dense(output_num))
#    model.add(Activation('sigmoid'))
#    #model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mae'])
#    model.compile(optimizer='adam',loss=loss_logadd2,metrics=['mae'])
#    #model.compile(optimizer='adam',loss=mix_mse_mape,metrics=['mae'])
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
#    epochs = 20000
#    batch_size = 4096
#    verbose = 1
#    history = model.fit(X[train],Y[train],batch_size,epochs,verbose,
#                        callbacks=cbks,validation_data=(X[test],Y[test])
#
#
#    #評価
#    score = model.evaluate(X[test], Y[test], verbose=1)
#    print(model.metrics_names[1], score[1])
#    cvscores.append(score[1])
######################################


t2 = time.time()
print('errapsed time is ',t2-t1)
#print((np.mean(cvscores), np.std(cvscores)))

#saving model for python
abspath_weight = abspath + '/learned_model/stanford_weight.hdf5'
open(abspath_model,"w").write(model.to_json())
model.save_weights(abspath_weight)

#saving model for python
abspath_model = abspath + '/learned_model/stanford_model.h5'
model.save(abspath_model, include_optimizer=False)

#saving model for c++
abspath_model_for_c = abspath + '/learned_model/model_for_c++/'
saver = tf.compat.v1.train.Saver()
saver.save(sess, "./learned_model/model_for_c++/model.ckpt")
saver.export_meta_graph("./learned_model/model_for_c++/model.meta")
tf.io.write_graph(sess.graph.as_graph_def(), "./learned_model/model_for_c++/" ,"graph.pb")


### add for TensorBoard
#KTF.set_session(old_session)
###
