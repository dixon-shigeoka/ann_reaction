import numpy as np
import ctypes as c
import os
import time

#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_length = abspath + '/learning_data/data_length.npy'

#ctypesの引数設定
fn = np.ctypeslib.load_library("mtsdriver.so",".")
fn.imtss_.argtypes = [
    c.POINTER(c.c_double),
    c.POINTER(c.c_double),
    np.ctypeslib.ndpointer(dtype=np.float64),
    c.POINTER(c.c_double)
]
fn.imtss_.restype = c.c_void_p


# chemical_database

aemn = 1.e-10
delt_mts = 1.e-7 # base timeの時間刻み
list_num = 9    # 要素数
tcounter = 0    # 温度を変化させるループのカウンター
#total_num = 1  # 合計値
train_x = np.zeros([1,12],dtype=float)
train_y = np.zeros([1,11],dtype=float)
train_x_zeros = np.zeros([1,1])
train_y_zeros = np.zeros([1,1])
data_length = np.zeros([1,1])
dtmp = 0
counter = 0     # データ数カウンター

for i in range(1):  #質量分率，温度，密度を変化させるループ(base time)

    print(dtmp)
    data_length = np.append(data_length,counter)
    t1 = time.time()
    dtmp = 1200 + float(tcounter*50)
    dprs = 1.01325e5
    aMi = np.array([2.016,32.000,1.008,16.000,17.008,18.016,33.008,34.016,28.016])
    aMolarRatio = np.array([2,1,0,0,0,0,0,0,0])
    aMolarRatio = np.where(aMolarRatio == 0,aemn,aMolarRatio)
    aTotalMass = np.dot(aMi,aMolarRatio)
    aMixMolarRatio = np.multiply(aMi,aMolarRatio)
    aYi = aMixMolarRatio/aTotalMass
    equiv_error = 1.e-1
    tcounter = tcounter + 1
    print('tcounter is ', tcounter)

    while (equiv_error > 1.E-4 or dtmp < 2000) :     #時間方向にmts計算を行い訓練データを格納するループ(base timeからの変化)

        counter = counter + 1
        delt_mem = 0.e0

        train_x_append = np.append(train_x_zeros,dtmp)
        train_x_append = np.append(train_x_append,dprs)
        train_x_append = np.append(train_x_append,aYi)

        range_mts = 1 #次のfor loopを回す回数

        for ii in range(range_mts):    #deltをパラメータにする分だけmtsを回す

            #delt_rand = np.random.rand(1)*0.5*delt_mts
            delt_rand = delt_mts
            delt_mem = delt_mem + delt_rand #このループ内で回った時間

            #入力データを格納
            temp_before = dtmp
            train_x_append_in_mtsloop = np.append(train_x_append,delt_mem)
            train_x_moment = np.delete(train_x_append_in_mtsloop,0,0)
            train_x_moment = train_x_moment.reshape((1,12))
            train_x = np.concatenate([train_x,train_x_moment],axis=0)

            dtmp = c.c_double(dtmp)     #ctypes形式でラップ
            dprs = c.c_double(dprs)
            delt = c.c_double(delt_mem)

            #byrefでポインタにして渡す，mtsの実行
            fn.imtss_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(delt))
            dtmp = dtmp.value
            dprs = dprs.value
            delt = delt.value

            #MTS結果を教師データとして格納
            train_y_append = np.append(train_y_zeros,dtmp)
            train_y_append = np.append(train_y_append,dprs)
            train_y_append = np.append(train_y_append,aYi)
            train_y_moment = np.delete(train_y_append,0,0)
            train_y_moment = train_y_moment.reshape((1,11))
            train_y = np.concatenate([train_y,train_y_moment],axis=0)
            #print(train_x_append,dtmp)


        if (range_mts != 1) :
            dtmp = train_x_append[1]
            dprs = train_x_append[2]
            aYi  = train_x_append[3:12]
            aYi  = aYi.reshape((1,9))
            dtmp = c.c_double(dtmp)     #ctypes形式でラップ
            dprs = c.c_double(dprs)
            delt = c.c_double(delt_mts)
            #H2_before = aYi[0,0]

            #byrefでポインタにして渡す，mtsの実行
            fn.imtss_(c.byref(dtmp),c.byref(dprs),aYi,c.byref(delt))
            dtmp = dtmp.value
            dprs = dprs.value
            delt = delt.value

        #equiv_error = abs(aYi[0,0] - H2_before)
        equiv_error = abs(dtmp - temp_before)


t2 = time.time()
print('mts loop time is ', t2-t1)
counter = counter + 1
data_length = np.append(data_length,counter)


# OUTPUT
data_length = np.delete(data_length,0,0)
train_x = np.delete(train_x,0,0)
train_y = np.delete(train_y,0,0)
np.save(abspath_x,train_x)
np.save(abspath_y,train_y)
np.save(abspath_length,data_length)

