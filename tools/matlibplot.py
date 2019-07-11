import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
from tensorflow.python.client import device_lib
import tensorflow as tf

def step_function(x):
  y = x > 0
  return y.astype(np.int)

x = np.arange(-5.0, 5.0, 0.1)
y = step_function(x)
plt.plot(x, y)
plt.ylim(-0.1, 1.1)
plt.savefig("step.png")

sess = tf.Session(config=tf.ConfigProto(log_device_placement=True))
