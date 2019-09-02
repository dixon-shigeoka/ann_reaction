import tensorflow as tf
import numpy as np
import matplotlib.pyplot as plt
import scipy.io
import time
import ctypes as c
import time
import os
import sys

np.random.seed(1234)
tf.set_random_seed(1234)
#numpy配列でprintした際に表示上の桁落ちを防ぐ
#np.set_printoptions(precision=15)

class PhysicsInformedNN:
    # Initialize the class
    def __init__(self, N, t_f, T0_f, p0_f, Yi0_f, layers):

        self.N = N
        self.layers = layers

        self.t0 = np.zeros([N,1])

        self.t_f_min = t_f[0]
        self.t_f_len = t_f[1] - t_f[0]
        self.T0_f_min = T0_f[0]
        self.T0_f_len = T0_f[1] - T0_f[0]
        self.p0_f_min = p0_f[0]
        self.p0_f_len = p0_f[1] - p0_f[0]

        self.t_f   = self.t_f_len * np.random.rand(self.N, 1) + self.t_f_min
        self.T0_f  = self.T0_f_len  * np.random.rand(self.N, 1)  + self.T0_f_min
        self.p0_f  = self.p0_f_len  * np.random.rand(self.N, 1)  + self.p0_f_min
        #self.Yi_f = self.Yi_f_len * nYi.random.rand(self.N, 1) + self.Yi_f_min
        #self.T0_f  = T0_f  * np.ones([N,1])
        #self.p0_f  = p0_f  * np.ones([N,1])
        self.Yi0_f = np.dot(np.ones([N,1]),Yi0_f)
        #print(self.t_f.dtype,self.t_f.dtype,self.p0_f.dtype)
        #print(self.t_f.shape,self.t_f.shape,self.p0_f.shape)
        #self.X = np.concatenate([self.t_f,self.T0_f,self.p0_f],axis=1)
        #self.X 1= np.hstack([self.t_f,self.T0_f,self.p0_f])

        X = np.concatenate([self.t_f,self.T0_f,self.p0_f],axis=1)
        self.lb = X.min(axis=0)
        self.ub = X.min(axis=0)


        self.omegai_f, self.totaldens = self.make_omegai(self.N, self.t_f, self.T0_f, self.p0_f, self.Yi0_f)
        self.rhoi0_f  = self.make_rhoi0(self.Yi0_f, self.totaldens)

        # Cut N2 mass fraction
        self.Yi0_f = np.delete(self.Yi0_f,-1,axis=1)
        self.rhoi0_f = np.delete(self.rhoi0_f,-1,axis=1)

        print(self.Yi0_f.shape, self.rhoi0_f.shape, self.p0_f.shape, self.T0_f.shape)

        # Initialize NNs
        self.layers = layers
        self.weights, self.biases = self.initialize_NN(layers)

        # tf placeholders and graph
        self.sess = tf.Session(config=tf.ConfigProto(allow_soft_placement=True,
                                                     log_device_placement=True))

        self.t0_tf = tf.placeholder(tf.float64, shape=[None, self.t0.shape[1]])
        self.T0_tf = tf.placeholder(tf.float64, shape=[None, self.T0_f.shape[1]])
        self.p0_tf = tf.placeholder(tf.float64, shape=[None, self.p0_f.shape[1]])
        self.Yi0_tf = tf.placeholder(tf.float64, shape=[None, self.Yi0_f.shape[1]])
        self.rhoi0_tf = tf.placeholder(tf.float64, shape=[None, self.rhoi0_f.shape[1]])

        self.t_f_tf = tf.placeholder(tf.float64, shape=[None, self.t_f.shape[1]])
        self.omegai_f_tf = tf.placeholder(tf.float64, shape=[None, self.omegai_f.shape[1]])

        self.rhoi0_pred = self.net_u(self.t0_tf, self.T0_tf, self.p0_tf, self.Yi0_tf)
        self.f_pred = self.net_f(self.t_f_tf, self.T0_tf, self.p0_tf, self.Yi0_tf, self.omegai_f_tf)

        self.loss = tf.reduce_mean(tf.square(self.rhoi0_tf - self.rhoi0_pred)) + \
                    tf.reduce_mean(tf.square(self.f_pred))


        self.optimizer = tf.contrib.opt.ScipyOptimizerInterface(self.loss,
                                                                method = 'L-BFGS-B',
                                                                options = {'maxiter': 50000,
                                                                           'maxfun': 50000,
                                                                           'maxcor': 50,
                                                                           'maxls': 50,
                                                                           'ftol' : 1.0 * np.finfo(float).eps})

        self.optimizer_Adam = tf.train.AdamOptimizer()
        self.train_op_Adam = self.optimizer_Adam.minimize(self.loss)

        init = tf.global_variables_initializer()
        self.sess.run(init)


    def initialize_NN(self, layers):
        weights = []
        biases = []
        num_layers = len(layers)
        for l in range(0,num_layers-1):
            W = self.xavier_init(size=[layers[l], layers[l+1]])
            b = tf.Variable(tf.zeros([1,layers[l+1]], dtype=tf.float64), dtype=tf.float64)
            weights.append(W)
            biases.append(b)
        return weights, biases

    def xavier_init(self, size):
        in_dim = size[0]
        out_dim = size[1]
        xavier_stddev = np.sqrt(2/(in_dim + out_dim))
        return tf.Variable(tf.truncated_normal([in_dim, out_dim], stddev=xavier_stddev), dtype=tf.float64)

    def make_omegai(self, N, t_f, T0_f, p0_f, Yi0_f) :
        omegai = np.zeros([self.N,self.Yi0_f.shape[1]-1])
        totaldens = np.zeros([self.N,1])

        N = c.c_int(self.N)
        delt = self.t_f[:]
        dtmp = self.T0_f[:]
        dprs = self.p0_f[:]
        aYi = self.Yi0_f[:,:]

        #byrefでポインタにして渡す，mtsの実行
        #多次元配列はFortranに渡した際に転置されるので注意
        fn1.imtss_omega_(c.byref(N),delt,dtmp,dprs,aYi,omegai,totaldens)
        return omegai, totaldens

    def make_rhoi0(self, Yi0_f, totaldens) :
        rhoi0 = self.Yi0_f * totaldens
        return rhoi0

    def neural_net(self, X, weights, biases):
        num_layers = len(weights) + 1
        #H[0:3] = 2.0*(X[0:3] - self.lb)/(self.ub - self.lb) - 1.0
        #H[3:11] = X[3:11]
        H = X
        for l in range(0,num_layers-2):
            W = weights[l]
            b = biases[l]
            H = tf.tanh(tf.add(tf.matmul(H, W), b))
        W = weights[-1]
        b = biases[-1]
        Y = tf.add(tf.matmul(H, W), b)
        return Y

    def net_u(self, t, T, p, Yi):
        rhoi = self.neural_net(tf.concat([t,T,p,Yi],axis=1), self.weights, self.biases)
        return rhoi

    def net_f(self, t, T, p, Yi, omegai):
        rhoi = self.neural_net(tf.concat([t,T,p,Yi],axis=1), self.weights, self.biases)
        #rhoi = self.net_u(t, T, p, Yi)
        rhoi_t = tf.gradients(rhoi, t)[0]
        f = rhoi_t - omegai
        return f

    def callback(self, loss):
        print('Loss:', loss)

    def train(self, nIter):
        tf_dict = {self.t0_tf: self.t0, self.T0_tf: self.T0_f, self.p0_tf: self.p0_f,
                   self.Yi0_tf: self.Yi0_f, self.rhoi0_tf: self.rhoi0_f,
                   self.t_f_tf: self.t_f, self.omegai_f_tf: self.omegai_f}

        start_time = time.time()

        # nIter = 0ならフルバッチL-BFGS-B, nIter /= 0ならAdam
        for it in range(nIter):
            self.sess.run(self.train_op_Adam, tf_dict)
            # Print
            if it % 10 == 0:
                elapsed = time.time() - start_time
                loss_value = self.sess.run(self.loss, tf_dict)
                print('It: %d, Loss: %.3e, Time: %.2f' %
                      (it, loss_value, elapsed))
                start_time = time.time()

        self.optimizer.minimize(self.sess,
                                feed_dict = tf_dict,
                                fetches = [self.loss],
                                loss_callback = self.callback)

    def predict(self, t, T0, p0, Yi0):
        t = np.expand_dims(t, axis=1)
        T0 = np.expand_dims(T0, axis=1)
        p0 = np.expand_dims(p0, axis=1)
        tf_dict = {self.t_f_tf:t, self.T0_tf:T0, self.p0_tf:p0,self.Yi0_tf:Yi0}
        rhoi_pred = self.sess.run(self.rhoi0_pred,tf_dict)
        Yi = rhoi/np.sum(rhoi,axis=1)
        return Yi




if __name__ == "__main__":

    #ctypesの引数設定
    fn1 = np.ctypeslib.load_library("mts_make_omega.so",".")
    fn1.imtss_omega_.argtypes = [
        c.POINTER(c.c_int),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
        np.ctypeslib.ndpointer(dtype=np.float64),
    ]
    fn1.imtss_omega_.restype = c.c_void_p


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

    # MTS training dataの読み込み
    #print('train_y',train_y)
    train_x = np.load(abspath_x)
    train_y = np.load(abspath_y)
    state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])
    #train_y = train_y[:,np.array([0,5])] #H mass fraction
    #train_x = np.delete(train_x,11,1) # delete delt
    #train_x = np.delete(train_x,10,1) # delete N2 mass fraction
    #train_y = np.delete(train_y,10,1) # delete N2 mass fraction
    train_y = np.delete(train_y,0,1) # delete Temperature
    train_y = np.delete(train_y,0,1) # delete Pressure


    #単精度と倍精度の変更
    #train_x = tf.cast(train_x,tf.float32)
    #train_y = tf.cast(train_y,tf.float32)


    #標準化または対数正規化
    #log conversion
    #train_x = np.log(train_x)
    #train_y = np.log(train_y)

    #standardization
    #mean_x = np.mean(train_x,axis=0)
    #mean_y = np.mean(train_y,axis=0)
    #std_x = np.std(train_x,axis=0)
    #std_y = np.std(train_y,axis=0)
    #train_x = (train_x - mean_x) / std_x
    #train_y = (train_y - mean_y) / std_y

    #min-max normalize
    #min_x = np.min(train_x,axis=0)
    #min_y = np.min(train_y,axis=0)
    #max_x = np.max(train_x,axis=0)
    #max_y = np.max(train_y,axis=0)
    #train_x = (train_x - min_x)/(max_x - min_x)
    #train_y = (train_y - min_y)/(max_y - min_y)


    #N_u = 1
    N = 30
    N_train = 0
    t_f = [0,3e-5]
    layers = [20] * 8
    layers = [11, *layers, 8]

    #T  = train_x[0,0:1 ]
    #p  = train_x[0,1:2 ]
    T  = [1190, 1210]
    p  = [1.01e5, 1.016e5 ]
    Yi = train_x[0,2:11]
    Yi = Yi.reshape(1,9)
    Exact = np.real(train_y)
    model = PhysicsInformedNN(N, t_f, T, p, Yi, layers)

    start_time = time.time()
    model.train(N_train)
    elapsed = time.time() - start_time
    print('Training time: %.4f' % (elapsed))

    t = np.array([3e-5])
    T = np.array([1200])
    p = np.array([1.01325e5])
    Yi = np.delete(Yi,-1,axis=1)
    Yi_pred = model.predict(t, T, p, Yi)

    error_u = np.linalg.norm(train_y-Yi_pred,2)/np.linalg.norm(train_y,2)
    print('Error u: %e' % (error_Yi))


##modelの保存
#abspath_weight = abspath + '/learned_model/stanford_weight.hdf5'
#open(abspath_model,"w").write(model.to_json())
#model.save_weights(abspath_weight)

