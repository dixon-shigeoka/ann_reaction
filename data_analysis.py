import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import os
from sklearn import preprocessing

#パスの指定
abspath = os.path.dirname(os.path.abspath(__file__))
abspath_x = abspath + '/learning_data/train_x.npy'
abspath_y = abspath + '/learning_data/train_y.npy'
abspath_length = abspath+ '/learning_data/data_length.npy'
abspath_model = abspath + '/learned_model/stanford_model.json'
abspath_eval = abspath + '/output/x.csv'
abspath_answer = abspath + '/output/y.csv'
abspath_randt = abspath + '/output/randt.csv'
print(abspath_x)

state_x = np.array(["temp","pres","H2","O2","H","O","OH","H2O","HO2","H2O2","N2","delt"])
select  = np.array([0,5])

# MTS training dataの読み込み
#print('train_y',train_y)
train_x = np.load(abspath_x)
train_y = np.load(abspath_y)
np.savetxt(abspath_eval,train_x,delimiter=',')
#np.savetxt(abspath_eval,train_data,delimiter=',')
np.savetxt(abspath_answer,train_y,delimiter=',')
#np.savetxt(abspath_randt,train_x,delimiter=',')
#train_x = train_x[0:100,select]
#train_y = train_y[0:100,select]
#train_x = train_x[:,select]
#train_y = train_y[:,select]

#print(train_x[:,0])
min_x = np.min(train_x,axis=0)
min_y = np.min(train_y,axis=0)
max_x = np.max(train_x,axis=0)
max_y = np.max(train_y,axis=0)
train_x = (train_x - min_x)/(max_x - min_x)
train_y = (train_y - min_y)/(max_y - min_y)

#プロット生成

plt.hist(train_x[:,1],bins=200)
plt.savefig("O.png")
plt.close()
plt.hist(train_x[:,0],bins=200)
plt.savefig("temp.png")
plt.close()

# Box-Cox conversion

#pt = preprocessing.PowerTransformer(method='box-cox')
#train_x = pt.fit_transform(train_x)

#プロット生成

#plt.hist(train_x[:,1],bins=200)
#plt.savefig("O_boxcox.png")
#plt.close()
#plt.hist(train_x[:,0],bins=200)
#plt.savefig("temp_boxcox.png")
#plt.close()
