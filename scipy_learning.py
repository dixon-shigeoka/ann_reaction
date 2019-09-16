import numpy as np
import time
import os
import scipy.io
from keras.models import Sequential
from keras.layers import Dense,BatchNormalization,Activation,Dropout
from keras import backend as K
from keras.initializers import he_normal
from sklearn.model_selection import KFold, train_test_split

from scipy.optimize import minimize
from tensorflow import keras
from tensorflow.keras import backend as K  # pylint: disable=import-error
from tensorflow.python.keras.callbacks import BaseLogger, CallbackList, History  # pylint: disable=no-name-in-module
from tensorflow.python.keras.optimizer_v2 import optimizer_v2   # pylint: disable=no-name-in-module
from tqdm import trange, tqdm_notebook


class GradientObserver(optimizer_v2.OptimizerV2):
    """
    Implements the Keras Optimizer interface in order to accumulate gradients for
    each mini batch. Gradients are then read at the end of the epoch by the ScipyOptimizer. 
    """

    def __init__(self):
        super(GradientObserver, self).__init__('GradientObserver')
        self._vars = []
        self._grads = {}

    def get_updates(self, loss, params):
        """
        Build the graph nodes that accumulate gradients.
        """
        self.updates = []
        grads = self.get_gradients(loss, params)
        for param, grad in zip(params, grads):
            shape = K.int_shape(param)
            var = K.zeros(shape)
            self._vars.append(var)
            self.updates.append(K.update_add(var, grad))
        return self.updates

    def get_gradient_values(self):
        """
        Read gradient values (at epoch end).
        """
        values = []
        for g in self._grads.values():
            values.append(K.eval(g))
        for v in self._vars:
            values.append(K.eval(v))
        return values

    def clear(self):
        """
        Clear gradient values (used at epoch start)
        """
        self._grads = {}
        for var in self._vars:
            K.set_value(var, np.zeros(var.shape))

    def _create_slots(self, var_list):
        pass

    def _resource_apply_dense(self, grad, var, **apply_kwargs):
        self._grads[var.name] = grad

    def _resource_apply_sparse(self, grad, handle, indices, apply_state):
        if handle.name in self._grads:
            dense_grad = self._grads[handle.name]
        else:
            dense_grad = np.zeros(handle.shape.as_list())
            self._grads[handle.name] = dense_grad

        for i, idx in enumerate(indices.numpy()):
            dense_grad[idx] += grad[i]

    def get_config(self):
        config = super(GradientObserver, self).get_config()
        return config


class GeneratorWrapper(keras.utils.Sequence):
    """
    Converts fit() into fit_generator() interface.
    """

    def __init__(self, inputs, outputs):
        self._inputs = inputs
        self._outputs = outputs

    def __len__(self):
        return 1

    def __getitem__(self, index):
        assert index == 0
        return self._inputs, self._outputs


class ScipyOptimizer(object):
    """
    Invokes the underlying model in order to obtain the cost and gradients for the function
    being optimized.
    """

    def __init__(self, model):
        self._model = model
        self._layers = [layer for layer in model._layers if layer.weights]
        self._weights_size = 0
        for layer in self._layers:
            for w in layer.weights:
                if not w.trainable:
                    continue
                self._weights_size += np.prod(w.shape)

    def _update_weights(self, x):
        x_offset = 0
        for layer in self._layers:
            w_list = []
            w_trainable = [w.trainable for w in layer.weights]
            batch_update = False not in w_trainable
            for w in layer.weights:
                if not w.trainable:
                    continue
                shape = w.get_shape()
                w_size = np.prod(shape)
                value = np.array(x[x_offset:x_offset+w_size]).reshape(shape)
                if batch_update:
                    w_list.append(value)
                else:
                    K.set_value(w, value)
                x_offset += w_size
            if batch_update:
                layer.set_weights(w_list)
        assert x_offset == self._weights_size

    def _collect_weights(self):
        x_values = np.empty(self._weights_size)
        x_offset = 0
        for layer in self._layers:
            w_trainable = [w.trainable for w in layer.weights]
            for var, trainable in zip(layer.get_weights(), w_trainable):
                if not trainable:
                    continue
                w_size = var.size
                x_values[x_offset:x_offset+w_size] = var.reshape(-1)
                x_offset += w_size
        assert x_offset == self._weights_size
        return x_values

    def _fun_generator(self, x, generator, state):
        self._model.optimizer.clear()
        self._update_weights(x)
        callbacks = state['callbacks']

        if not state['in_epoch']:
            callbacks.on_epoch_begin(state['epoch'])
            state['in_epoch'] = True

        cost_sum = 0
        iterator = trange(len(generator)) if state['verbose'] else range(
            len(generator))

        state['epoch_logs'] = {}
        epoch_logs = state['epoch_logs']

        for batch_index in iterator:
            inputs, outputs = generator[batch_index]
            if isinstance(inputs, list):
                isize = inputs[0].shape[0]
            else:
                isize = inputs.shape[0]
            batch_logs = {'batch': batch_index, 'size': isize}
            callbacks.on_batch_begin(batch_index, batch_logs)
            outs = self._model.train_on_batch(inputs, outputs)
            if not isinstance(outs, list):
                outs = [outs]
            for lbl, v in zip(self._model.metrics_names, outs):
                batch_logs[lbl] = v
                epoch_logs[lbl] = epoch_logs.get(lbl, 0.0) + v
            callbacks.on_batch_end(batch_index, batch_logs)
            batch_cost = batch_logs['loss']
            if state['verbose']:
                iterator.set_postfix(cost=batch_cost)
            cost_sum += batch_cost

        generator.on_epoch_end()

        # average the metrics
        for lbl in self._model.metrics_names:
            epoch_logs[lbl] = epoch_logs.get(lbl) / len(iterator)

        cost = cost_sum / len(iterator)

        gradients = self._model.optimizer.get_gradient_values()
        x_grad = np.empty(x.shape)
        x_offset = 0
        for grad in gradients:
            w_size = grad.size
            x_grad[x_offset:x_offset + w_size] = grad.reshape(-1)
            x_offset += w_size
        assert x_offset == self._weights_size
        self._cost = cost
        self._gradients = x_grad
        return cost, x_grad

    def _validate(self, x, val_generator, state):
        # TODO: weight update should be optimized in the most common case
        self._update_weights(x)
        # test callback are in the newer version of the CallbackList API.
        # callbacks = state['callbacks']
        epoch_logs = state['epoch_logs']

        # callbacks.on_test_begin()

        val_outs = [0] * len(self._model.metrics_names)
        n_steps = len(val_generator)
        for batch_index in range(n_steps):
            inputs, outputs = val_generator[batch_index]
            batch_logs = {'batch': batch_index, 'size': inputs.shape[0]}

            # callbacks.on_test_batch_begin(batch_index, batch_logs)
            batch_outs = self._model.test_on_batch(inputs, outputs)
            if not isinstance(batch_outs, list):
                batch_outs = [batch_outs]
            for l, o in zip(self._model.metrics_names, batch_outs):
                batch_logs[l] = o
            # callbacks.on_test_batch_end(batch_index, batch_logs)
            for i, batch_out in enumerate(batch_outs):
                val_outs[i] += batch_out

        for l, o in zip(self._model.metrics_names, val_outs):
            o /= n_steps
            epoch_logs['val_' + l] = o

        # callbacks.on_test_end()

    def fit(self, inputs, outputs, **kwargs):
        return self.fit_generator(GeneratorWrapper(inputs, outputs), **kwargs)

    def fit_generator(self, generator, method='cg', epochs=1,
                      validation_data=None,
                      callbacks=None,
                      verbose=True):
        x0 = self._collect_weights()
        history = History()
        _callbacks = [BaseLogger(stateful_metrics=self._model.metrics_names)]
        _callbacks += (callbacks or []) + [history]
        callback_list = CallbackList(_callbacks)
        callback_list.set_model(self._model)
        callback_list.set_params({
            'epochs': epochs,
            'verbose': verbose,
            'metrics': list(self._model.metrics_names),
        })
        state = {
            'epoch': 0,
            'verbose': verbose,
            'callbacks': callback_list,
            'in_epoch': False,
            'epoch_logs': {},
        }
        min_options = {
            'maxiter': epochs,
        }

        val_generator = None
        if validation_data is not None:
            if isinstance(validation_data, keras.utils.Sequence):
                val_generator = validation_data
            elif isinstance(validation_data, tuple) and len(validation_data) == 2:
                val_generator = GeneratorWrapper(*validation_data)

        def on_iteration_end(xk):
            cb = state['callbacks']
            if val_generator is not None:
                self._validate(xk, val_generator, state)
            cb.on_epoch_end(state['epoch'], state['epoch_logs'])
            if state['verbose']:
                epoch_logs = state['epoch_logs']
                print('epoch: ', state['epoch'],
                      ', '.join(['{0}: {1}'.format(k, v) for k, v in epoch_logs.items()]))
            state['epoch'] += 1
            state['in_epoch'] = False
            state['epoch_logs'] = {}

        callback_list.on_train_begin()
        result = minimize(
            self._fun_generator, x0, method=method, jac=True, options=min_options,
            callback=on_iteration_end, args=(generator, state))
        self._update_weights(result['x'])
        callback_list.on_train_end()
        return result, history



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
def loss_logadd(y_true,y_pred):
    first_log = K.log(K.abs(y_pred))
    second_log = K.log(K.abs(y_true))
    return  K.mean(K.square(K.abs(y_pred - y_true) + first_log - second_log),axis=-1)

def loss_logadd2(y_true,y_pred):
    first_log = K.log(K.abs(y_pred))
    second_log = K.log(K.abs(y_true))
    #return  K.mean(K.square(K.abs(y_pred - y_true) + K.log(y_pred/y_true)),axis=-1)
    #return  K.mean(K.square(K.abs(y_pred - y_true) + tf.div(y_pred,y_true)),axis=-1)
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

    print(train_y_state)

    # define X-fold cross validation
    model = Sequential()

    #パラメータ設定
    model.add(Dense(30,input_dim=input_num))
    #model.add(BatchNormalization())
    model.add(Activation('relu'))
    #model.add(Dropout(0.1))
    #model.add(Dense(30))
    #model.add(BatchNormalization())
    #model.add(Activation('relu'))
    #model.add(Dropout(0.1))
    model.add(Dense(output_num))
    model.add(Activation('sigmoid'))
    #model.compile(optimizer='adam',loss='mean_squared_error',metrics=['mae'])
    #model.compile(optimizer='adam',loss=loss_logadd2,metrics=['mae'])
    model.compile(optimizer=GradientObserver(),loss=loss_logadd2)
    model.summary()

    opt = ScipyOptimizer(L-BFGS-B)

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
    epochs = 50000
    batch_size = 16
    verbose = 2
    #history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,
    #                    callbacks=cbks,validation_data=(test_x,test_y))
    #history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,callbacks=cbks)
    result, history = opt.fit_generator(train_x,train_y_state, epochs, verbose)

    ##学習実行
    #epochs = 1000
    #batch_size = 1
    #verbose = 2
    ##history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,
    ##                    callbacks=cbks,validation_data=(test_x,test_y))
    #history = model.fit(train_x,train_y_state,batch_size,epochs,verbose,callbacks=cbks)

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
#
#callbacks=cbks,validation_data=(X[test],Y[test])
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
