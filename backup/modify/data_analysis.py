import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os

#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/test_x.npy'
abspath_y = abspath + '/learning_data/test_y.npy'
abspath_model = abspath + '/learned_model/stanford_model.json'
abspath_weight = abspath + '/learned_model/stanford_weight.h5'
abspath_eval = abspath + '/output/eval.csv'
abspath_answer = abspath + '/output/answer.csv'


# MTS training dataの読み込み
#print('test_y',test_y)
test_x  = np.load(abspath_x)
test_y  = np.load(abspath_y)
print('x     ', test_x[0,:])
log_x   = np.log(test_x)
print('log_x ', log_x[0,:])
log_y   = np.log(test_y)
mean_x  = np.mean(log_x,axis=0)
mean_y  = np.mean(log_y,axis=0)
var_x   = np.var(log_x)
var_y   = np.var(log_y)
sqvar_x = np.sqrt(var_x)
sqvar_y = np.sqrt(var_y)
test_x = (log_x - mean_x) / sqvar_x
test_y = (log_y - mean_y) / sqvar_y
print('x mean', mean_x)

#プロット生成
temp = test_x[:,0]
aH2  = log_x[:,2]

plt.hist(aH2,bins=500)
plt.savefig("aH2.png")
plt.hist(temp,bins=100)
plt.savefig("temp.png")
