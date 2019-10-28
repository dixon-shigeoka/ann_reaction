node {
  name: "dense_1_input"
  op: "Placeholder"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: -1
        }
        dim {
          size: 10
        }
      }
    }
  }
}
node {
  name: "dense_1/random_uniform/shape"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 2
          }
        }
        tensor_content: "\n\000\000\000\036\000\000\000"
      }
    }
  }
}
node {
  name: "dense_1/random_uniform/min"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: -0.3872983346207417
      }
    }
  }
}
node {
  name: "dense_1/random_uniform/max"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.3872983346207417
      }
    }
  }
}
node {
  name: "dense_1/random_uniform/RandomUniform"
  op: "RandomUniform"
  input: "dense_1/random_uniform/shape"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "seed"
    value {
      i: 87654321
    }
  }
  attr {
    key: "seed2"
    value {
      i: 8974511
    }
  }
}
node {
  name: "dense_1/random_uniform/sub"
  op: "Sub"
  input: "dense_1/random_uniform/max"
  input: "dense_1/random_uniform/min"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_1/random_uniform/mul"
  op: "Mul"
  input: "dense_1/random_uniform/RandomUniform"
  input: "dense_1/random_uniform/sub"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_1/random_uniform"
  op: "Add"
  input: "dense_1/random_uniform/mul"
  input: "dense_1/random_uniform/min"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_1/kernel"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 10
        }
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "dense_1/kernel/Assign"
  op: "Assign"
  input: "dense_1/kernel"
  input: "dense_1/random_uniform"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "dense_1/kernel/read"
  op: "Identity"
  input: "dense_1/kernel"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/kernel"
      }
    }
  }
}
node {
  name: "dense_1/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "dense_1/bias"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "dense_1/bias/Assign"
  op: "Assign"
  input: "dense_1/bias"
  input: "dense_1/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "dense_1/bias/read"
  op: "Identity"
  input: "dense_1/bias"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/bias"
      }
    }
  }
}
node {
  name: "dense_1/MatMul"
  op: "MatMul"
  input: "dense_1_input"
  input: "dense_1/kernel/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: false
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: false
    }
  }
}
node {
  name: "dense_1/BiasAdd"
  op: "BiasAdd"
  input: "dense_1/MatMul"
  input: "dense_1/bias/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "data_format"
    value {
      s: "NHWC"
    }
  }
}
node {
  name: "activation_1/Relu"
  op: "Relu"
  input: "dense_1/BiasAdd"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_2/random_uniform/shape"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 2
          }
        }
        tensor_content: "\036\000\000\000\036\000\000\000"
      }
    }
  }
}
node {
  name: "dense_2/random_uniform/min"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: -0.31622776601683794
      }
    }
  }
}
node {
  name: "dense_2/random_uniform/max"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.31622776601683794
      }
    }
  }
}
node {
  name: "dense_2/random_uniform/RandomUniform"
  op: "RandomUniform"
  input: "dense_2/random_uniform/shape"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "seed"
    value {
      i: 87654321
    }
  }
  attr {
    key: "seed2"
    value {
      i: 3335364
    }
  }
}
node {
  name: "dense_2/random_uniform/sub"
  op: "Sub"
  input: "dense_2/random_uniform/max"
  input: "dense_2/random_uniform/min"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_2/random_uniform/mul"
  op: "Mul"
  input: "dense_2/random_uniform/RandomUniform"
  input: "dense_2/random_uniform/sub"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_2/random_uniform"
  op: "Add"
  input: "dense_2/random_uniform/mul"
  input: "dense_2/random_uniform/min"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_2/kernel"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "dense_2/kernel/Assign"
  op: "Assign"
  input: "dense_2/kernel"
  input: "dense_2/random_uniform"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "dense_2/kernel/read"
  op: "Identity"
  input: "dense_2/kernel"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/kernel"
      }
    }
  }
}
node {
  name: "dense_2/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "dense_2/bias"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "dense_2/bias/Assign"
  op: "Assign"
  input: "dense_2/bias"
  input: "dense_2/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "dense_2/bias/read"
  op: "Identity"
  input: "dense_2/bias"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/bias"
      }
    }
  }
}
node {
  name: "dense_2/MatMul"
  op: "MatMul"
  input: "activation_1/Relu"
  input: "dense_2/kernel/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: false
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: false
    }
  }
}
node {
  name: "dense_2/BiasAdd"
  op: "BiasAdd"
  input: "dense_2/MatMul"
  input: "dense_2/bias/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "data_format"
    value {
      s: "NHWC"
    }
  }
}
node {
  name: "activation_2/Relu"
  op: "Relu"
  input: "dense_2/BiasAdd"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_3/random_uniform/shape"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 2
          }
        }
        tensor_content: "\036\000\000\000\010\000\000\000"
      }
    }
  }
}
node {
  name: "dense_3/random_uniform/min"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: -0.39735970711951313
      }
    }
  }
}
node {
  name: "dense_3/random_uniform/max"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.39735970711951313
      }
    }
  }
}
node {
  name: "dense_3/random_uniform/RandomUniform"
  op: "RandomUniform"
  input: "dense_3/random_uniform/shape"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "seed"
    value {
      i: 87654321
    }
  }
  attr {
    key: "seed2"
    value {
      i: 3905091
    }
  }
}
node {
  name: "dense_3/random_uniform/sub"
  op: "Sub"
  input: "dense_3/random_uniform/max"
  input: "dense_3/random_uniform/min"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_3/random_uniform/mul"
  op: "Mul"
  input: "dense_3/random_uniform/RandomUniform"
  input: "dense_3/random_uniform/sub"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_3/random_uniform"
  op: "Add"
  input: "dense_3/random_uniform/mul"
  input: "dense_3/random_uniform/min"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "dense_3/kernel"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
        dim {
          size: 8
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "dense_3/kernel/Assign"
  op: "Assign"
  input: "dense_3/kernel"
  input: "dense_3/random_uniform"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "dense_3/kernel/read"
  op: "Identity"
  input: "dense_3/kernel"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/kernel"
      }
    }
  }
}
node {
  name: "dense_3/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 8
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "dense_3/bias"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 8
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "dense_3/bias/Assign"
  op: "Assign"
  input: "dense_3/bias"
  input: "dense_3/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "dense_3/bias/read"
  op: "Identity"
  input: "dense_3/bias"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/bias"
      }
    }
  }
}
node {
  name: "dense_3/MatMul"
  op: "MatMul"
  input: "activation_2/Relu"
  input: "dense_3/kernel/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: false
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: false
    }
  }
}
node {
  name: "dense_3/BiasAdd"
  op: "BiasAdd"
  input: "dense_3/MatMul"
  input: "dense_3/bias/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "data_format"
    value {
      s: "NHWC"
    }
  }
}
node {
  name: "activation_3/Sigmoid"
  op: "Sigmoid"
  input: "dense_3/BiasAdd"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "Adam/iterations/initial_value"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT64
        tensor_shape {
        }
        int64_val: 0
      }
    }
  }
}
node {
  name: "Adam/iterations"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "Adam/iterations/Assign"
  op: "Assign"
  input: "Adam/iterations"
  input: "Adam/iterations/initial_value"
  attr {
    key: "T"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/iterations"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "Adam/iterations/read"
  op: "Identity"
  input: "Adam/iterations"
  attr {
    key: "T"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/iterations"
      }
    }
  }
}
node {
  name: "Adam/lr/initial_value"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.001
      }
    }
  }
}
node {
  name: "Adam/lr"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "Adam/lr/Assign"
  op: "Assign"
  input: "Adam/lr"
  input: "Adam/lr/initial_value"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/lr"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "Adam/lr/read"
  op: "Identity"
  input: "Adam/lr"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/lr"
      }
    }
  }
}
node {
  name: "Adam/beta_1/initial_value"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.9
      }
    }
  }
}
node {
  name: "Adam/beta_1"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "Adam/beta_1/Assign"
  op: "Assign"
  input: "Adam/beta_1"
  input: "Adam/beta_1/initial_value"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_1"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "Adam/beta_1/read"
  op: "Identity"
  input: "Adam/beta_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_1"
      }
    }
  }
}
node {
  name: "Adam/beta_2/initial_value"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.999
      }
    }
  }
}
node {
  name: "Adam/beta_2"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "Adam/beta_2/Assign"
  op: "Assign"
  input: "Adam/beta_2"
  input: "Adam/beta_2/initial_value"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_2"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "Adam/beta_2/read"
  op: "Identity"
  input: "Adam/beta_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_2"
      }
    }
  }
}
node {
  name: "Adam/decay/initial_value"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "Adam/decay"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "Adam/decay/Assign"
  op: "Assign"
  input: "Adam/decay"
  input: "Adam/decay/initial_value"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/decay"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "Adam/decay/read"
  op: "Identity"
  input: "Adam/decay"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/decay"
      }
    }
  }
}
node {
  name: "activation_3_target"
  op: "Placeholder"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: -1
        }
        dim {
          size: -1
        }
      }
    }
  }
}
node {
  name: "activation_3_sample_weights"
  op: "Placeholder"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: -1
        }
      }
    }
  }
}
node {
  name: "loss/activation_3_loss/Abs"
  op: "Abs"
  input: "activation_3/Sigmoid"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Log"
  op: "Log"
  input: "loss/activation_3_loss/Abs"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Abs_1"
  op: "Abs"
  input: "activation_3_target"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Log_1"
  op: "Log"
  input: "loss/activation_3_loss/Abs_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/sub"
  op: "Sub"
  input: "loss/activation_3_loss/Log"
  input: "loss/activation_3_loss/Log_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Square"
  op: "Square"
  input: "loss/activation_3_loss/sub"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Mean/reduction_indices"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: -1
      }
    }
  }
}
node {
  name: "loss/activation_3_loss/Mean"
  op: "Mean"
  input: "loss/activation_3_loss/Square"
  input: "loss/activation_3_loss/Mean/reduction_indices"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "loss/activation_3_loss/Mean_1/reduction_indices"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "loss/activation_3_loss/Mean_1"
  op: "Mean"
  input: "loss/activation_3_loss/Mean"
  input: "loss/activation_3_loss/Mean_1/reduction_indices"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "loss/activation_3_loss/mul"
  op: "Mul"
  input: "loss/activation_3_loss/Mean_1"
  input: "activation_3_sample_weights"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/NotEqual/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "loss/activation_3_loss/NotEqual"
  op: "NotEqual"
  input: "activation_3_sample_weights"
  input: "loss/activation_3_loss/NotEqual/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Cast"
  op: "Cast"
  input: "loss/activation_3_loss/NotEqual"
  attr {
    key: "DstT"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "SrcT"
    value {
      type: DT_BOOL
    }
  }
  attr {
    key: "Truncate"
    value {
      b: false
    }
  }
}
node {
  name: "loss/activation_3_loss/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "loss/activation_3_loss/Mean_2"
  op: "Mean"
  input: "loss/activation_3_loss/Cast"
  input: "loss/activation_3_loss/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "loss/activation_3_loss/truediv"
  op: "RealDiv"
  input: "loss/activation_3_loss/mul"
  input: "loss/activation_3_loss/Mean_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "loss/activation_3_loss/Const_1"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "loss/activation_3_loss/Mean_3"
  op: "Mean"
  input: "loss/activation_3_loss/truediv"
  input: "loss/activation_3_loss/Const_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "loss/mul/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "loss/mul"
  op: "Mul"
  input: "loss/mul/x"
  input: "loss/activation_3_loss/Mean_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "metrics/mean_absolute_error/sub"
  op: "Sub"
  input: "activation_3/Sigmoid"
  input: "activation_3_target"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "metrics/mean_absolute_error/Abs"
  op: "Abs"
  input: "metrics/mean_absolute_error/sub"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "metrics/mean_absolute_error/Mean/reduction_indices"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: -1
      }
    }
  }
}
node {
  name: "metrics/mean_absolute_error/Mean"
  op: "Mean"
  input: "metrics/mean_absolute_error/Abs"
  input: "metrics/mean_absolute_error/Mean/reduction_indices"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "metrics/mean_absolute_error/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "metrics/mean_absolute_error/Mean_1"
  op: "Mean"
  input: "metrics/mean_absolute_error/Mean"
  input: "metrics/mean_absolute_error/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/Shape"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "training/Adam/gradients/grad_ys_0"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/Fill"
  op: "Fill"
  input: "training/Adam/gradients/Shape"
  input: "training/Adam/gradients/grad_ys_0"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Shape"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Shape_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/BroadcastGradientArgs"
  op: "BroadcastGradientArgs"
  input: "training/Adam/gradients/loss/mul_grad/Shape"
  input: "training/Adam/gradients/loss/mul_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Mul"
  op: "Mul"
  input: "training/Adam/gradients/Fill"
  input: "loss/activation_3_loss/Mean_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Sum"
  op: "Sum"
  input: "training/Adam/gradients/loss/mul_grad/Mul"
  input: "training/Adam/gradients/loss/mul_grad/BroadcastGradientArgs"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/mul_grad/Sum"
  input: "training/Adam/gradients/loss/mul_grad/Shape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Mul_1"
  op: "Mul"
  input: "loss/mul/x"
  input: "training/Adam/gradients/Fill"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Sum_1"
  op: "Sum"
  input: "training/Adam/gradients/loss/mul_grad/Mul_1"
  input: "training/Adam/gradients/loss/mul_grad/BroadcastGradientArgs:1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/mul_grad/Reshape_1"
  op: "Reshape"
  input: "training/Adam/gradients/loss/mul_grad/Sum_1"
  input: "training/Adam/gradients/loss/mul_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Reshape/shape"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/mul_grad/Reshape_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Reshape/shape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Shape"
  op: "Shape"
  input: "loss/activation_3_loss/truediv"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Tile"
  op: "Tile"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Shape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tmultiples"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Shape_1"
  op: "Shape"
  input: "loss/activation_3_loss/truediv"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Shape_2"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Const"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Prod"
  op: "Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Shape_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Const"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Const_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Prod_1"
  op: "Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Shape_2"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Const_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Maximum/y"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Maximum"
  op: "Maximum"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Prod_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Maximum/y"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/floordiv"
  op: "FloorDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Maximum"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Cast"
  op: "Cast"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/floordiv"
  attr {
    key: "DstT"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "SrcT"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Truncate"
    value {
      b: false
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/truediv"
  op: "RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Tile"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/Cast"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_3"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Shape"
  op: "Shape"
  input: "loss/activation_3_loss/mul"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Shape_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/BroadcastGradientArgs"
  op: "BroadcastGradientArgs"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/RealDiv"
  op: "RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/truediv"
  input: "loss/activation_3_loss/Mean_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Sum"
  op: "Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/BroadcastGradientArgs"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Shape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Neg"
  op: "Neg"
  input: "loss/activation_3_loss/mul"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/RealDiv_1"
  op: "RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Neg"
  input: "loss/activation_3_loss/Mean_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/RealDiv_2"
  op: "RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/RealDiv_1"
  input: "loss/activation_3_loss/Mean_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/mul"
  op: "Mul"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_3_grad/truediv"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/RealDiv_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Sum_1"
  op: "Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/mul"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/BroadcastGradientArgs:1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Reshape_1"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Sum_1"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/truediv"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Shape"
  op: "Shape"
  input: "loss/activation_3_loss/Mean_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Shape_1"
  op: "Shape"
  input: "activation_3_sample_weights"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/BroadcastGradientArgs"
  op: "BroadcastGradientArgs"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Mul"
  op: "Mul"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Reshape"
  input: "activation_3_sample_weights"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Sum"
  op: "Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Mul"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/BroadcastGradientArgs"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Shape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Mul_1"
  op: "Mul"
  input: "loss/activation_3_loss/Mean_1"
  input: "training/Adam/gradients/loss/activation_3_loss/truediv_grad/Reshape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Sum_1"
  op: "Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Mul_1"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/BroadcastGradientArgs:1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Reshape_1"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Sum_1"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/mul"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape"
  op: "Shape"
  input: "loss/activation_3_loss/Mean"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Size"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/add"
  op: "Add"
  input: "loss/activation_3_loss/Mean_1/reduction_indices"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Size"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/mod"
  op: "FloorMod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/add"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Size"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/range/start"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/range/delta"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/range"
  op: "Range"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/range/start"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Size"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/range/delta"
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Fill/value"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Fill"
  op: "Fill"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Fill/value"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/DynamicStitch"
  op: "DynamicStitch"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/range"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/mod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Fill"
  attr {
    key: "N"
    value {
      i: 2
    }
  }
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum/y"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum"
  op: "Maximum"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/DynamicStitch"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum/y"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/floordiv"
  op: "FloorDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/mul_grad/Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/DynamicStitch"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Tile"
  op: "Tile"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/floordiv"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tmultiples"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape_2"
  op: "Shape"
  input: "loss/activation_3_loss/Mean"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape_3"
  op: "Shape"
  input: "loss/activation_3_loss/Mean_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Const"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Prod"
  op: "Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape_2"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Const"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Const_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Prod_1"
  op: "Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Shape_3"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Const_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum_1/y"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum_1"
  op: "Maximum"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Prod_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum_1/y"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/floordiv_1"
  op: "FloorDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Maximum_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Cast"
  op: "Cast"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/floordiv_1"
  attr {
    key: "DstT"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "SrcT"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Truncate"
    value {
      b: false
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/truediv"
  op: "RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Tile"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/Cast"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean_1"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape"
  op: "Shape"
  input: "loss/activation_3_loss/Square"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Size"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 2
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/add"
  op: "Add"
  input: "loss/activation_3_loss/Mean/reduction_indices"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Size"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/mod"
  op: "FloorMod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/add"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Size"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
          }
        }
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/range/start"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/range/delta"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/range"
  op: "Range"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/range/start"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Size"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/range/delta"
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Fill/value"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Fill"
  op: "Fill"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Fill/value"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/DynamicStitch"
  op: "DynamicStitch"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/range"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/mod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Fill"
  attr {
    key: "N"
    value {
      i: 2
    }
  }
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum/y"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum"
  op: "Maximum"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/DynamicStitch"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum/y"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/floordiv"
  op: "FloorDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_1_grad/truediv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/DynamicStitch"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Tile"
  op: "Tile"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/floordiv"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tmultiples"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape_2"
  op: "Shape"
  input: "loss/activation_3_loss/Square"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape_3"
  op: "Shape"
  input: "loss/activation_3_loss/Mean"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Const"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Prod"
  op: "Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape_2"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Const"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Const_1"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Prod_1"
  op: "Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Shape_3"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Const_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum_1/y"
  op: "Const"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum_1"
  op: "Maximum"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Prod_1"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum_1/y"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/floordiv_1"
  op: "FloorDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Prod"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Maximum_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Cast"
  op: "Cast"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/floordiv_1"
  attr {
    key: "DstT"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "SrcT"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "Truncate"
    value {
      b: false
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/truediv"
  op: "RealDiv"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Tile"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/Cast"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Mean"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Const"
  op: "Const"
  input: "^training/Adam/gradients/loss/activation_3_loss/Mean_grad/truediv"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Square"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 2.0
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Mul"
  op: "Mul"
  input: "loss/activation_3_loss/sub"
  input: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Square"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Mul_1"
  op: "Mul"
  input: "training/Adam/gradients/loss/activation_3_loss/Mean_grad/truediv"
  input: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Mul"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Square"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Shape"
  op: "Shape"
  input: "loss/activation_3_loss/Log"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Shape_1"
  op: "Shape"
  input: "loss/activation_3_loss/Log_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
  attr {
    key: "out_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/BroadcastGradientArgs"
  op: "BroadcastGradientArgs"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Shape"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Sum"
  op: "Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Mul_1"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/BroadcastGradientArgs"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Reshape"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Shape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Sum_1"
  op: "Sum"
  input: "training/Adam/gradients/loss/activation_3_loss/Square_grad/Mul_1"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/BroadcastGradientArgs:1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tidx"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
  attr {
    key: "keep_dims"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Neg"
  op: "Neg"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Sum_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Reshape_1"
  op: "Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Neg"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Shape_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "Tshape"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/sub"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Log_grad/Reciprocal"
  op: "Reciprocal"
  input: "loss/activation_3_loss/Abs"
  input: "^training/Adam/gradients/loss/activation_3_loss/sub_grad/Reshape"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Log"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Log_grad/mul"
  op: "Mul"
  input: "training/Adam/gradients/loss/activation_3_loss/sub_grad/Reshape"
  input: "training/Adam/gradients/loss/activation_3_loss/Log_grad/Reciprocal"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Log"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Abs_grad/Sign"
  op: "Sign"
  input: "activation_3/Sigmoid"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Abs"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/loss/activation_3_loss/Abs_grad/mul"
  op: "Mul"
  input: "training/Adam/gradients/loss/activation_3_loss/Log_grad/mul"
  input: "training/Adam/gradients/loss/activation_3_loss/Abs_grad/Sign"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@loss/activation_3_loss/Abs"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/activation_3/Sigmoid_grad/SigmoidGrad"
  op: "SigmoidGrad"
  input: "activation_3/Sigmoid"
  input: "training/Adam/gradients/loss/activation_3_loss/Abs_grad/mul"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@activation_3/Sigmoid"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/dense_3/BiasAdd_grad/BiasAddGrad"
  op: "BiasAddGrad"
  input: "training/Adam/gradients/activation_3/Sigmoid_grad/SigmoidGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/BiasAdd"
      }
    }
  }
  attr {
    key: "data_format"
    value {
      s: "NHWC"
    }
  }
}
node {
  name: "training/Adam/gradients/dense_3/MatMul_grad/MatMul"
  op: "MatMul"
  input: "training/Adam/gradients/activation_3/Sigmoid_grad/SigmoidGrad"
  input: "dense_3/kernel/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/MatMul"
      }
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: false
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/gradients/dense_3/MatMul_grad/MatMul_1"
  op: "MatMul"
  input: "activation_2/Relu"
  input: "training/Adam/gradients/activation_3/Sigmoid_grad/SigmoidGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/MatMul"
      }
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: true
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/activation_2/Relu_grad/ReluGrad"
  op: "ReluGrad"
  input: "training/Adam/gradients/dense_3/MatMul_grad/MatMul"
  input: "activation_2/Relu"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@activation_2/Relu"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/dense_2/BiasAdd_grad/BiasAddGrad"
  op: "BiasAddGrad"
  input: "training/Adam/gradients/activation_2/Relu_grad/ReluGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/BiasAdd"
      }
    }
  }
  attr {
    key: "data_format"
    value {
      s: "NHWC"
    }
  }
}
node {
  name: "training/Adam/gradients/dense_2/MatMul_grad/MatMul"
  op: "MatMul"
  input: "training/Adam/gradients/activation_2/Relu_grad/ReluGrad"
  input: "dense_2/kernel/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/MatMul"
      }
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: false
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/gradients/dense_2/MatMul_grad/MatMul_1"
  op: "MatMul"
  input: "activation_1/Relu"
  input: "training/Adam/gradients/activation_2/Relu_grad/ReluGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/MatMul"
      }
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: true
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/gradients/activation_1/Relu_grad/ReluGrad"
  op: "ReluGrad"
  input: "training/Adam/gradients/dense_2/MatMul_grad/MatMul"
  input: "activation_1/Relu"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@activation_1/Relu"
      }
    }
  }
}
node {
  name: "training/Adam/gradients/dense_1/BiasAdd_grad/BiasAddGrad"
  op: "BiasAddGrad"
  input: "training/Adam/gradients/activation_1/Relu_grad/ReluGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/BiasAdd"
      }
    }
  }
  attr {
    key: "data_format"
    value {
      s: "NHWC"
    }
  }
}
node {
  name: "training/Adam/gradients/dense_1/MatMul_grad/MatMul"
  op: "MatMul"
  input: "training/Adam/gradients/activation_1/Relu_grad/ReluGrad"
  input: "dense_1/kernel/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/MatMul"
      }
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: false
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/gradients/dense_1/MatMul_grad/MatMul_1"
  op: "MatMul"
  input: "dense_1_input"
  input: "training/Adam/gradients/activation_1/Relu_grad/ReluGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/MatMul"
      }
    }
  }
  attr {
    key: "transpose_a"
    value {
      b: true
    }
  }
  attr {
    key: "transpose_b"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/AssignAdd/value"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT64
        tensor_shape {
        }
        int64_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/AssignAdd"
  op: "AssignAdd"
  input: "Adam/iterations"
  input: "training/Adam/AssignAdd/value"
  attr {
    key: "T"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/iterations"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/Cast"
  op: "Cast"
  input: "Adam/iterations/read"
  attr {
    key: "DstT"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "SrcT"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "Truncate"
    value {
      b: false
    }
  }
}
node {
  name: "training/Adam/add/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/add"
  op: "Add"
  input: "training/Adam/Cast"
  input: "training/Adam/add/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Pow"
  op: "Pow"
  input: "Adam/beta_2/read"
  input: "training/Adam/add"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub"
  op: "Sub"
  input: "training/Adam/sub/x"
  input: "training/Adam/Pow"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_1"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value/Minimum"
  op: "Minimum"
  input: "training/Adam/sub"
  input: "training/Adam/Const_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value"
  op: "Maximum"
  input: "training/Adam/clip_by_value/Minimum"
  input: "training/Adam/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt"
  op: "Sqrt"
  input: "training/Adam/clip_by_value"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Pow_1"
  op: "Pow"
  input: "Adam/beta_1/read"
  input: "training/Adam/add"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_1/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_1"
  op: "Sub"
  input: "training/Adam/sub_1/x"
  input: "training/Adam/Pow_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv"
  op: "RealDiv"
  input: "training/Adam/Sqrt"
  input: "training/Adam/sub_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul"
  op: "Mul"
  input: "Adam/lr/read"
  input: "training/Adam/truediv"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/zeros"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 10
          }
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 10
        }
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable/Assign"
  op: "Assign"
  input: "training/Adam/Variable"
  input: "training/Adam/zeros"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable/read"
  op: "Identity"
  input: "training/Adam/Variable"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_1"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_1"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_1/Assign"
  op: "Assign"
  input: "training/Adam/Variable_1"
  input: "training/Adam/zeros_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_1"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_1/read"
  op: "Identity"
  input: "training/Adam/Variable_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_1"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_2"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_2"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_2/Assign"
  op: "Assign"
  input: "training/Adam/Variable_2"
  input: "training/Adam/zeros_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_2"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_2/read"
  op: "Identity"
  input: "training/Adam/Variable_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_2"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_3"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_3"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_3/Assign"
  op: "Assign"
  input: "training/Adam/Variable_3"
  input: "training/Adam/zeros_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_3"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_3/read"
  op: "Identity"
  input: "training/Adam/Variable_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_3"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_4"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
          dim {
            size: 8
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_4"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
        dim {
          size: 8
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_4/Assign"
  op: "Assign"
  input: "training/Adam/Variable_4"
  input: "training/Adam/zeros_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_4"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_4/read"
  op: "Identity"
  input: "training/Adam/Variable_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_4"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_5"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 8
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_5"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 8
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_5/Assign"
  op: "Assign"
  input: "training/Adam/Variable_5"
  input: "training/Adam/zeros_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_5"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_5/read"
  op: "Identity"
  input: "training/Adam/Variable_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_5"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_6"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 10
          }
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_6"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 10
        }
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_6/Assign"
  op: "Assign"
  input: "training/Adam/Variable_6"
  input: "training/Adam/zeros_6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_6"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_6/read"
  op: "Identity"
  input: "training/Adam/Variable_6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_6"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_7"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_7"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_7/Assign"
  op: "Assign"
  input: "training/Adam/Variable_7"
  input: "training/Adam/zeros_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_7"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_7/read"
  op: "Identity"
  input: "training/Adam/Variable_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_7"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_8"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_8"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_8/Assign"
  op: "Assign"
  input: "training/Adam/Variable_8"
  input: "training/Adam/zeros_8"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_8"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_8/read"
  op: "Identity"
  input: "training/Adam/Variable_8"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_8"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_9"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_9"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_9/Assign"
  op: "Assign"
  input: "training/Adam/Variable_9"
  input: "training/Adam/zeros_9"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_9"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_9/read"
  op: "Identity"
  input: "training/Adam/Variable_9"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_9"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_10"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 30
          }
          dim {
            size: 8
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_10"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 30
        }
        dim {
          size: 8
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_10/Assign"
  op: "Assign"
  input: "training/Adam/Variable_10"
  input: "training/Adam/zeros_10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_10"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_10/read"
  op: "Identity"
  input: "training/Adam/Variable_10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_10"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_11"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
          dim {
            size: 8
          }
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Variable_11"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 8
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_11/Assign"
  op: "Assign"
  input: "training/Adam/Variable_11"
  input: "training/Adam/zeros_11"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_11"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_11/read"
  op: "Identity"
  input: "training/Adam/Variable_11"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_11"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_12/shape_as_tensor"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/zeros_12/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/zeros_12"
  op: "Fill"
  input: "training/Adam/zeros_12/shape_as_tensor"
  input: "training/Adam/zeros_12/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/Variable_12"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 1
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_12/Assign"
  op: "Assign"
  input: "training/Adam/Variable_12"
  input: "training/Adam/zeros_12"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_12"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_12/read"
  op: "Identity"
  input: "training/Adam/Variable_12"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_12"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_13/shape_as_tensor"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/zeros_13/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/zeros_13"
  op: "Fill"
  input: "training/Adam/zeros_13/shape_as_tensor"
  input: "training/Adam/zeros_13/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/Variable_13"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 1
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_13/Assign"
  op: "Assign"
  input: "training/Adam/Variable_13"
  input: "training/Adam/zeros_13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_13"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_13/read"
  op: "Identity"
  input: "training/Adam/Variable_13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_13"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_14/shape_as_tensor"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/zeros_14/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/zeros_14"
  op: "Fill"
  input: "training/Adam/zeros_14/shape_as_tensor"
  input: "training/Adam/zeros_14/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/Variable_14"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 1
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_14/Assign"
  op: "Assign"
  input: "training/Adam/Variable_14"
  input: "training/Adam/zeros_14"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_14"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_14/read"
  op: "Identity"
  input: "training/Adam/Variable_14"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_14"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_15/shape_as_tensor"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/zeros_15/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/zeros_15"
  op: "Fill"
  input: "training/Adam/zeros_15/shape_as_tensor"
  input: "training/Adam/zeros_15/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/Variable_15"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 1
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_15/Assign"
  op: "Assign"
  input: "training/Adam/Variable_15"
  input: "training/Adam/zeros_15"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_15"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_15/read"
  op: "Identity"
  input: "training/Adam/Variable_15"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_15"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_16/shape_as_tensor"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/zeros_16/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/zeros_16"
  op: "Fill"
  input: "training/Adam/zeros_16/shape_as_tensor"
  input: "training/Adam/zeros_16/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/Variable_16"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 1
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_16/Assign"
  op: "Assign"
  input: "training/Adam/Variable_16"
  input: "training/Adam/zeros_16"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_16"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_16/read"
  op: "Identity"
  input: "training/Adam/Variable_16"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_16"
      }
    }
  }
}
node {
  name: "training/Adam/zeros_17/shape_as_tensor"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_INT32
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_INT32
        tensor_shape {
          dim {
            size: 1
          }
        }
        int_val: 1
      }
    }
  }
}
node {
  name: "training/Adam/zeros_17/Const"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/zeros_17"
  op: "Fill"
  input: "training/Adam/zeros_17/shape_as_tensor"
  input: "training/Adam/zeros_17/Const"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "index_type"
    value {
      type: DT_INT32
    }
  }
}
node {
  name: "training/Adam/Variable_17"
  op: "VariableV2"
  attr {
    key: "container"
    value {
      s: ""
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "shape"
    value {
      shape {
        dim {
          size: 1
        }
      }
    }
  }
  attr {
    key: "shared_name"
    value {
      s: ""
    }
  }
}
node {
  name: "training/Adam/Variable_17/Assign"
  op: "Assign"
  input: "training/Adam/Variable_17"
  input: "training/Adam/zeros_17"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_17"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Variable_17/read"
  op: "Identity"
  input: "training/Adam/Variable_17"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_17"
      }
    }
  }
}
node {
  name: "training/Adam/mul_1"
  op: "Mul"
  input: "Adam/beta_1/read"
  input: "training/Adam/Variable/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_2/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_2"
  op: "Sub"
  input: "training/Adam/sub_2/x"
  input: "Adam/beta_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_2"
  op: "Mul"
  input: "training/Adam/sub_2"
  input: "training/Adam/gradients/dense_1/MatMul_grad/MatMul_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_1"
  op: "Add"
  input: "training/Adam/mul_1"
  input: "training/Adam/mul_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_3"
  op: "Mul"
  input: "Adam/beta_2/read"
  input: "training/Adam/Variable_6/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_3/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_3"
  op: "Sub"
  input: "training/Adam/sub_3/x"
  input: "Adam/beta_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Square"
  op: "Square"
  input: "training/Adam/gradients/dense_1/MatMul_grad/MatMul_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_4"
  op: "Mul"
  input: "training/Adam/sub_3"
  input: "training/Adam/Square"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_2"
  op: "Add"
  input: "training/Adam/mul_3"
  input: "training/Adam/mul_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_5"
  op: "Mul"
  input: "training/Adam/mul"
  input: "training/Adam/add_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const_2"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_3"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value_1/Minimum"
  op: "Minimum"
  input: "training/Adam/add_2"
  input: "training/Adam/Const_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value_1"
  op: "Maximum"
  input: "training/Adam/clip_by_value_1/Minimum"
  input: "training/Adam/Const_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt_1"
  op: "Sqrt"
  input: "training/Adam/clip_by_value_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_3/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1e-07
      }
    }
  }
}
node {
  name: "training/Adam/add_3"
  op: "Add"
  input: "training/Adam/Sqrt_1"
  input: "training/Adam/add_3/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv_1"
  op: "RealDiv"
  input: "training/Adam/mul_5"
  input: "training/Adam/add_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_4"
  op: "Sub"
  input: "dense_1/kernel/read"
  input: "training/Adam/truediv_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Assign"
  op: "Assign"
  input: "training/Adam/Variable"
  input: "training/Adam/add_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_1"
  op: "Assign"
  input: "training/Adam/Variable_6"
  input: "training/Adam/add_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_6"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_2"
  op: "Assign"
  input: "dense_1/kernel"
  input: "training/Adam/sub_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/mul_6"
  op: "Mul"
  input: "Adam/beta_1/read"
  input: "training/Adam/Variable_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_5/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_5"
  op: "Sub"
  input: "training/Adam/sub_5/x"
  input: "Adam/beta_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_7"
  op: "Mul"
  input: "training/Adam/sub_5"
  input: "training/Adam/gradients/dense_1/BiasAdd_grad/BiasAddGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_4"
  op: "Add"
  input: "training/Adam/mul_6"
  input: "training/Adam/mul_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_8"
  op: "Mul"
  input: "Adam/beta_2/read"
  input: "training/Adam/Variable_7/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_6/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_6"
  op: "Sub"
  input: "training/Adam/sub_6/x"
  input: "Adam/beta_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Square_1"
  op: "Square"
  input: "training/Adam/gradients/dense_1/BiasAdd_grad/BiasAddGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_9"
  op: "Mul"
  input: "training/Adam/sub_6"
  input: "training/Adam/Square_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_5"
  op: "Add"
  input: "training/Adam/mul_8"
  input: "training/Adam/mul_9"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_10"
  op: "Mul"
  input: "training/Adam/mul"
  input: "training/Adam/add_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const_4"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_5"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value_2/Minimum"
  op: "Minimum"
  input: "training/Adam/add_5"
  input: "training/Adam/Const_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value_2"
  op: "Maximum"
  input: "training/Adam/clip_by_value_2/Minimum"
  input: "training/Adam/Const_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt_2"
  op: "Sqrt"
  input: "training/Adam/clip_by_value_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_6/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1e-07
      }
    }
  }
}
node {
  name: "training/Adam/add_6"
  op: "Add"
  input: "training/Adam/Sqrt_2"
  input: "training/Adam/add_6/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv_2"
  op: "RealDiv"
  input: "training/Adam/mul_10"
  input: "training/Adam/add_6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_7"
  op: "Sub"
  input: "dense_1/bias/read"
  input: "training/Adam/truediv_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Assign_3"
  op: "Assign"
  input: "training/Adam/Variable_1"
  input: "training/Adam/add_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_1"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_4"
  op: "Assign"
  input: "training/Adam/Variable_7"
  input: "training/Adam/add_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_7"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_5"
  op: "Assign"
  input: "dense_1/bias"
  input: "training/Adam/sub_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/mul_11"
  op: "Mul"
  input: "Adam/beta_1/read"
  input: "training/Adam/Variable_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_8/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_8"
  op: "Sub"
  input: "training/Adam/sub_8/x"
  input: "Adam/beta_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_12"
  op: "Mul"
  input: "training/Adam/sub_8"
  input: "training/Adam/gradients/dense_2/MatMul_grad/MatMul_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_7"
  op: "Add"
  input: "training/Adam/mul_11"
  input: "training/Adam/mul_12"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_13"
  op: "Mul"
  input: "Adam/beta_2/read"
  input: "training/Adam/Variable_8/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_9/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_9"
  op: "Sub"
  input: "training/Adam/sub_9/x"
  input: "Adam/beta_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Square_2"
  op: "Square"
  input: "training/Adam/gradients/dense_2/MatMul_grad/MatMul_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_14"
  op: "Mul"
  input: "training/Adam/sub_9"
  input: "training/Adam/Square_2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_8"
  op: "Add"
  input: "training/Adam/mul_13"
  input: "training/Adam/mul_14"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_15"
  op: "Mul"
  input: "training/Adam/mul"
  input: "training/Adam/add_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const_6"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_7"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value_3/Minimum"
  op: "Minimum"
  input: "training/Adam/add_8"
  input: "training/Adam/Const_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value_3"
  op: "Maximum"
  input: "training/Adam/clip_by_value_3/Minimum"
  input: "training/Adam/Const_6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt_3"
  op: "Sqrt"
  input: "training/Adam/clip_by_value_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_9/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1e-07
      }
    }
  }
}
node {
  name: "training/Adam/add_9"
  op: "Add"
  input: "training/Adam/Sqrt_3"
  input: "training/Adam/add_9/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv_3"
  op: "RealDiv"
  input: "training/Adam/mul_15"
  input: "training/Adam/add_9"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_10"
  op: "Sub"
  input: "dense_2/kernel/read"
  input: "training/Adam/truediv_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Assign_6"
  op: "Assign"
  input: "training/Adam/Variable_2"
  input: "training/Adam/add_7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_2"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_7"
  op: "Assign"
  input: "training/Adam/Variable_8"
  input: "training/Adam/add_8"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_8"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_8"
  op: "Assign"
  input: "dense_2/kernel"
  input: "training/Adam/sub_10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/mul_16"
  op: "Mul"
  input: "Adam/beta_1/read"
  input: "training/Adam/Variable_3/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_11/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_11"
  op: "Sub"
  input: "training/Adam/sub_11/x"
  input: "Adam/beta_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_17"
  op: "Mul"
  input: "training/Adam/sub_11"
  input: "training/Adam/gradients/dense_2/BiasAdd_grad/BiasAddGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_10"
  op: "Add"
  input: "training/Adam/mul_16"
  input: "training/Adam/mul_17"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_18"
  op: "Mul"
  input: "Adam/beta_2/read"
  input: "training/Adam/Variable_9/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_12/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_12"
  op: "Sub"
  input: "training/Adam/sub_12/x"
  input: "Adam/beta_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Square_3"
  op: "Square"
  input: "training/Adam/gradients/dense_2/BiasAdd_grad/BiasAddGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_19"
  op: "Mul"
  input: "training/Adam/sub_12"
  input: "training/Adam/Square_3"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_11"
  op: "Add"
  input: "training/Adam/mul_18"
  input: "training/Adam/mul_19"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_20"
  op: "Mul"
  input: "training/Adam/mul"
  input: "training/Adam/add_10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const_8"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_9"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value_4/Minimum"
  op: "Minimum"
  input: "training/Adam/add_11"
  input: "training/Adam/Const_9"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value_4"
  op: "Maximum"
  input: "training/Adam/clip_by_value_4/Minimum"
  input: "training/Adam/Const_8"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt_4"
  op: "Sqrt"
  input: "training/Adam/clip_by_value_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_12/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1e-07
      }
    }
  }
}
node {
  name: "training/Adam/add_12"
  op: "Add"
  input: "training/Adam/Sqrt_4"
  input: "training/Adam/add_12/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv_4"
  op: "RealDiv"
  input: "training/Adam/mul_20"
  input: "training/Adam/add_12"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_13"
  op: "Sub"
  input: "dense_2/bias/read"
  input: "training/Adam/truediv_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Assign_9"
  op: "Assign"
  input: "training/Adam/Variable_3"
  input: "training/Adam/add_10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_3"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_10"
  op: "Assign"
  input: "training/Adam/Variable_9"
  input: "training/Adam/add_11"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_9"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_11"
  op: "Assign"
  input: "dense_2/bias"
  input: "training/Adam/sub_13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/mul_21"
  op: "Mul"
  input: "Adam/beta_1/read"
  input: "training/Adam/Variable_4/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_14/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_14"
  op: "Sub"
  input: "training/Adam/sub_14/x"
  input: "Adam/beta_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_22"
  op: "Mul"
  input: "training/Adam/sub_14"
  input: "training/Adam/gradients/dense_3/MatMul_grad/MatMul_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_13"
  op: "Add"
  input: "training/Adam/mul_21"
  input: "training/Adam/mul_22"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_23"
  op: "Mul"
  input: "Adam/beta_2/read"
  input: "training/Adam/Variable_10/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_15/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_15"
  op: "Sub"
  input: "training/Adam/sub_15/x"
  input: "Adam/beta_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Square_4"
  op: "Square"
  input: "training/Adam/gradients/dense_3/MatMul_grad/MatMul_1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_24"
  op: "Mul"
  input: "training/Adam/sub_15"
  input: "training/Adam/Square_4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_14"
  op: "Add"
  input: "training/Adam/mul_23"
  input: "training/Adam/mul_24"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_25"
  op: "Mul"
  input: "training/Adam/mul"
  input: "training/Adam/add_13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const_10"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_11"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value_5/Minimum"
  op: "Minimum"
  input: "training/Adam/add_14"
  input: "training/Adam/Const_11"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value_5"
  op: "Maximum"
  input: "training/Adam/clip_by_value_5/Minimum"
  input: "training/Adam/Const_10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt_5"
  op: "Sqrt"
  input: "training/Adam/clip_by_value_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_15/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1e-07
      }
    }
  }
}
node {
  name: "training/Adam/add_15"
  op: "Add"
  input: "training/Adam/Sqrt_5"
  input: "training/Adam/add_15/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv_5"
  op: "RealDiv"
  input: "training/Adam/mul_25"
  input: "training/Adam/add_15"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_16"
  op: "Sub"
  input: "dense_3/kernel/read"
  input: "training/Adam/truediv_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Assign_12"
  op: "Assign"
  input: "training/Adam/Variable_4"
  input: "training/Adam/add_13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_4"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_13"
  op: "Assign"
  input: "training/Adam/Variable_10"
  input: "training/Adam/add_14"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_10"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_14"
  op: "Assign"
  input: "dense_3/kernel"
  input: "training/Adam/sub_16"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/mul_26"
  op: "Mul"
  input: "Adam/beta_1/read"
  input: "training/Adam/Variable_5/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_17/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_17"
  op: "Sub"
  input: "training/Adam/sub_17/x"
  input: "Adam/beta_1/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_27"
  op: "Mul"
  input: "training/Adam/sub_17"
  input: "training/Adam/gradients/dense_3/BiasAdd_grad/BiasAddGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_16"
  op: "Add"
  input: "training/Adam/mul_26"
  input: "training/Adam/mul_27"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_28"
  op: "Mul"
  input: "Adam/beta_2/read"
  input: "training/Adam/Variable_11/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_18/x"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1.0
      }
    }
  }
}
node {
  name: "training/Adam/sub_18"
  op: "Sub"
  input: "training/Adam/sub_18/x"
  input: "Adam/beta_2/read"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Square_5"
  op: "Square"
  input: "training/Adam/gradients/dense_3/BiasAdd_grad/BiasAddGrad"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_29"
  op: "Mul"
  input: "training/Adam/sub_18"
  input: "training/Adam/Square_5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_17"
  op: "Add"
  input: "training/Adam/mul_28"
  input: "training/Adam/mul_29"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/mul_30"
  op: "Mul"
  input: "training/Adam/mul"
  input: "training/Adam/add_16"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Const_12"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 0.0
      }
    }
  }
}
node {
  name: "training/Adam/Const_13"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: inf
      }
    }
  }
}
node {
  name: "training/Adam/clip_by_value_6/Minimum"
  op: "Minimum"
  input: "training/Adam/add_17"
  input: "training/Adam/Const_13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/clip_by_value_6"
  op: "Maximum"
  input: "training/Adam/clip_by_value_6/Minimum"
  input: "training/Adam/Const_12"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Sqrt_6"
  op: "Sqrt"
  input: "training/Adam/clip_by_value_6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/add_18/y"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_DOUBLE
        tensor_shape {
        }
        double_val: 1e-07
      }
    }
  }
}
node {
  name: "training/Adam/add_18"
  op: "Add"
  input: "training/Adam/Sqrt_6"
  input: "training/Adam/add_18/y"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/truediv_6"
  op: "RealDiv"
  input: "training/Adam/mul_30"
  input: "training/Adam/add_18"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/sub_19"
  op: "Sub"
  input: "dense_3/bias/read"
  input: "training/Adam/truediv_6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "training/Adam/Assign_15"
  op: "Assign"
  input: "training/Adam/Variable_5"
  input: "training/Adam/add_16"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_5"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_16"
  op: "Assign"
  input: "training/Adam/Variable_11"
  input: "training/Adam/add_17"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_11"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/Adam/Assign_17"
  op: "Assign"
  input: "dense_3/bias"
  input: "training/Adam/sub_19"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "training/group_deps"
  op: "NoOp"
  input: "^loss/mul"
  input: "^metrics/mean_absolute_error/Mean_1"
  input: "^training/Adam/Assign"
  input: "^training/Adam/AssignAdd"
  input: "^training/Adam/Assign_1"
  input: "^training/Adam/Assign_10"
  input: "^training/Adam/Assign_11"
  input: "^training/Adam/Assign_12"
  input: "^training/Adam/Assign_13"
  input: "^training/Adam/Assign_14"
  input: "^training/Adam/Assign_15"
  input: "^training/Adam/Assign_16"
  input: "^training/Adam/Assign_17"
  input: "^training/Adam/Assign_2"
  input: "^training/Adam/Assign_3"
  input: "^training/Adam/Assign_4"
  input: "^training/Adam/Assign_5"
  input: "^training/Adam/Assign_6"
  input: "^training/Adam/Assign_7"
  input: "^training/Adam/Assign_8"
  input: "^training/Adam/Assign_9"
}
node {
  name: "group_deps"
  op: "NoOp"
  input: "^loss/mul"
  input: "^metrics/mean_absolute_error/Mean_1"
}
node {
  name: "IsVariableInitialized"
  op: "IsVariableInitialized"
  input: "dense_1/kernel"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/kernel"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_1"
  op: "IsVariableInitialized"
  input: "dense_1/bias"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/bias"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_2"
  op: "IsVariableInitialized"
  input: "dense_2/kernel"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/kernel"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_3"
  op: "IsVariableInitialized"
  input: "dense_2/bias"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/bias"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_4"
  op: "IsVariableInitialized"
  input: "dense_3/kernel"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/kernel"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_5"
  op: "IsVariableInitialized"
  input: "dense_3/bias"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/bias"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_6"
  op: "IsVariableInitialized"
  input: "Adam/iterations"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/iterations"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_INT64
    }
  }
}
node {
  name: "IsVariableInitialized_7"
  op: "IsVariableInitialized"
  input: "Adam/lr"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/lr"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_8"
  op: "IsVariableInitialized"
  input: "Adam/beta_1"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_9"
  op: "IsVariableInitialized"
  input: "Adam/beta_2"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_2"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_10"
  op: "IsVariableInitialized"
  input: "Adam/decay"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/decay"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_11"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_12"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_1"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_1"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_13"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_2"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_2"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_14"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_3"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_3"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_15"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_4"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_4"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_16"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_5"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_5"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_17"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_6"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_6"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_18"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_7"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_7"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_19"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_8"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_8"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_20"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_9"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_9"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_21"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_10"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_10"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_22"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_11"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_11"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_23"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_12"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_12"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_24"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_13"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_13"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_25"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_14"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_14"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_26"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_15"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_15"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_27"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_16"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_16"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "IsVariableInitialized_28"
  op: "IsVariableInitialized"
  input: "training/Adam/Variable_17"
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_17"
      }
    }
  }
  attr {
    key: "dtype"
    value {
      type: DT_DOUBLE
    }
  }
}
node {
  name: "init"
  op: "NoOp"
  input: "^Adam/beta_1/Assign"
  input: "^Adam/beta_2/Assign"
  input: "^Adam/decay/Assign"
  input: "^Adam/iterations/Assign"
  input: "^Adam/lr/Assign"
  input: "^dense_1/bias/Assign"
  input: "^dense_1/kernel/Assign"
  input: "^dense_2/bias/Assign"
  input: "^dense_2/kernel/Assign"
  input: "^dense_3/bias/Assign"
  input: "^dense_3/kernel/Assign"
  input: "^training/Adam/Variable/Assign"
  input: "^training/Adam/Variable_1/Assign"
  input: "^training/Adam/Variable_10/Assign"
  input: "^training/Adam/Variable_11/Assign"
  input: "^training/Adam/Variable_12/Assign"
  input: "^training/Adam/Variable_13/Assign"
  input: "^training/Adam/Variable_14/Assign"
  input: "^training/Adam/Variable_15/Assign"
  input: "^training/Adam/Variable_16/Assign"
  input: "^training/Adam/Variable_17/Assign"
  input: "^training/Adam/Variable_2/Assign"
  input: "^training/Adam/Variable_3/Assign"
  input: "^training/Adam/Variable_4/Assign"
  input: "^training/Adam/Variable_5/Assign"
  input: "^training/Adam/Variable_6/Assign"
  input: "^training/Adam/Variable_7/Assign"
  input: "^training/Adam/Variable_8/Assign"
  input: "^training/Adam/Variable_9/Assign"
}
node {
  name: "save/filename/input"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_STRING
        tensor_shape {
        }
        string_val: "model"
      }
    }
  }
}
node {
  name: "save/filename"
  op: "PlaceholderWithDefault"
  input: "save/filename/input"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
}
node {
  name: "save/Const"
  op: "PlaceholderWithDefault"
  input: "save/filename"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "shape"
    value {
      shape {
      }
    }
  }
}
node {
  name: "save/SaveV2/tensor_names"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_STRING
        tensor_shape {
          dim {
            size: 29
          }
        }
        string_val: "Adam/beta_1"
        string_val: "Adam/beta_2"
        string_val: "Adam/decay"
        string_val: "Adam/iterations"
        string_val: "Adam/lr"
        string_val: "dense_1/bias"
        string_val: "dense_1/kernel"
        string_val: "dense_2/bias"
        string_val: "dense_2/kernel"
        string_val: "dense_3/bias"
        string_val: "dense_3/kernel"
        string_val: "training/Adam/Variable"
        string_val: "training/Adam/Variable_1"
        string_val: "training/Adam/Variable_10"
        string_val: "training/Adam/Variable_11"
        string_val: "training/Adam/Variable_12"
        string_val: "training/Adam/Variable_13"
        string_val: "training/Adam/Variable_14"
        string_val: "training/Adam/Variable_15"
        string_val: "training/Adam/Variable_16"
        string_val: "training/Adam/Variable_17"
        string_val: "training/Adam/Variable_2"
        string_val: "training/Adam/Variable_3"
        string_val: "training/Adam/Variable_4"
        string_val: "training/Adam/Variable_5"
        string_val: "training/Adam/Variable_6"
        string_val: "training/Adam/Variable_7"
        string_val: "training/Adam/Variable_8"
        string_val: "training/Adam/Variable_9"
      }
    }
  }
}
node {
  name: "save/SaveV2/shape_and_slices"
  op: "Const"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_STRING
        tensor_shape {
          dim {
            size: 29
          }
        }
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
      }
    }
  }
}
node {
  name: "save/SaveV2"
  op: "SaveV2"
  input: "save/Const"
  input: "save/SaveV2/tensor_names"
  input: "save/SaveV2/shape_and_slices"
  input: "Adam/beta_1"
  input: "Adam/beta_2"
  input: "Adam/decay"
  input: "Adam/iterations"
  input: "Adam/lr"
  input: "dense_1/bias"
  input: "dense_1/kernel"
  input: "dense_2/bias"
  input: "dense_2/kernel"
  input: "dense_3/bias"
  input: "dense_3/kernel"
  input: "training/Adam/Variable"
  input: "training/Adam/Variable_1"
  input: "training/Adam/Variable_10"
  input: "training/Adam/Variable_11"
  input: "training/Adam/Variable_12"
  input: "training/Adam/Variable_13"
  input: "training/Adam/Variable_14"
  input: "training/Adam/Variable_15"
  input: "training/Adam/Variable_16"
  input: "training/Adam/Variable_17"
  input: "training/Adam/Variable_2"
  input: "training/Adam/Variable_3"
  input: "training/Adam/Variable_4"
  input: "training/Adam/Variable_5"
  input: "training/Adam/Variable_6"
  input: "training/Adam/Variable_7"
  input: "training/Adam/Variable_8"
  input: "training/Adam/Variable_9"
  attr {
    key: "dtypes"
    value {
      list {
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_INT64
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
      }
    }
  }
}
node {
  name: "save/control_dependency"
  op: "Identity"
  input: "save/Const"
  input: "^save/SaveV2"
  attr {
    key: "T"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@save/Const"
      }
    }
  }
}
node {
  name: "save/RestoreV2/tensor_names"
  op: "Const"
  device: "/device:CPU:0"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_STRING
        tensor_shape {
          dim {
            size: 29
          }
        }
        string_val: "Adam/beta_1"
        string_val: "Adam/beta_2"
        string_val: "Adam/decay"
        string_val: "Adam/iterations"
        string_val: "Adam/lr"
        string_val: "dense_1/bias"
        string_val: "dense_1/kernel"
        string_val: "dense_2/bias"
        string_val: "dense_2/kernel"
        string_val: "dense_3/bias"
        string_val: "dense_3/kernel"
        string_val: "training/Adam/Variable"
        string_val: "training/Adam/Variable_1"
        string_val: "training/Adam/Variable_10"
        string_val: "training/Adam/Variable_11"
        string_val: "training/Adam/Variable_12"
        string_val: "training/Adam/Variable_13"
        string_val: "training/Adam/Variable_14"
        string_val: "training/Adam/Variable_15"
        string_val: "training/Adam/Variable_16"
        string_val: "training/Adam/Variable_17"
        string_val: "training/Adam/Variable_2"
        string_val: "training/Adam/Variable_3"
        string_val: "training/Adam/Variable_4"
        string_val: "training/Adam/Variable_5"
        string_val: "training/Adam/Variable_6"
        string_val: "training/Adam/Variable_7"
        string_val: "training/Adam/Variable_8"
        string_val: "training/Adam/Variable_9"
      }
    }
  }
}
node {
  name: "save/RestoreV2/shape_and_slices"
  op: "Const"
  device: "/device:CPU:0"
  attr {
    key: "dtype"
    value {
      type: DT_STRING
    }
  }
  attr {
    key: "value"
    value {
      tensor {
        dtype: DT_STRING
        tensor_shape {
          dim {
            size: 29
          }
        }
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
        string_val: ""
      }
    }
  }
}
node {
  name: "save/RestoreV2"
  op: "RestoreV2"
  input: "save/Const"
  input: "save/RestoreV2/tensor_names"
  input: "save/RestoreV2/shape_and_slices"
  device: "/device:CPU:0"
  attr {
    key: "dtypes"
    value {
      list {
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_INT64
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
        type: DT_DOUBLE
      }
    }
  }
}
node {
  name: "save/Assign"
  op: "Assign"
  input: "Adam/beta_1"
  input: "save/RestoreV2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_1"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_1"
  op: "Assign"
  input: "Adam/beta_2"
  input: "save/RestoreV2:1"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/beta_2"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_2"
  op: "Assign"
  input: "Adam/decay"
  input: "save/RestoreV2:2"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/decay"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_3"
  op: "Assign"
  input: "Adam/iterations"
  input: "save/RestoreV2:3"
  attr {
    key: "T"
    value {
      type: DT_INT64
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/iterations"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_4"
  op: "Assign"
  input: "Adam/lr"
  input: "save/RestoreV2:4"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@Adam/lr"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_5"
  op: "Assign"
  input: "dense_1/bias"
  input: "save/RestoreV2:5"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_6"
  op: "Assign"
  input: "dense_1/kernel"
  input: "save/RestoreV2:6"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_1/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_7"
  op: "Assign"
  input: "dense_2/bias"
  input: "save/RestoreV2:7"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_8"
  op: "Assign"
  input: "dense_2/kernel"
  input: "save/RestoreV2:8"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_2/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_9"
  op: "Assign"
  input: "dense_3/bias"
  input: "save/RestoreV2:9"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/bias"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_10"
  op: "Assign"
  input: "dense_3/kernel"
  input: "save/RestoreV2:10"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@dense_3/kernel"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_11"
  op: "Assign"
  input: "training/Adam/Variable"
  input: "save/RestoreV2:11"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_12"
  op: "Assign"
  input: "training/Adam/Variable_1"
  input: "save/RestoreV2:12"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_1"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_13"
  op: "Assign"
  input: "training/Adam/Variable_10"
  input: "save/RestoreV2:13"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_10"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_14"
  op: "Assign"
  input: "training/Adam/Variable_11"
  input: "save/RestoreV2:14"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_11"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_15"
  op: "Assign"
  input: "training/Adam/Variable_12"
  input: "save/RestoreV2:15"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_12"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_16"
  op: "Assign"
  input: "training/Adam/Variable_13"
  input: "save/RestoreV2:16"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_13"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_17"
  op: "Assign"
  input: "training/Adam/Variable_14"
  input: "save/RestoreV2:17"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_14"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_18"
  op: "Assign"
  input: "training/Adam/Variable_15"
  input: "save/RestoreV2:18"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_15"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_19"
  op: "Assign"
  input: "training/Adam/Variable_16"
  input: "save/RestoreV2:19"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_16"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_20"
  op: "Assign"
  input: "training/Adam/Variable_17"
  input: "save/RestoreV2:20"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_17"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_21"
  op: "Assign"
  input: "training/Adam/Variable_2"
  input: "save/RestoreV2:21"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_2"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_22"
  op: "Assign"
  input: "training/Adam/Variable_3"
  input: "save/RestoreV2:22"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_3"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_23"
  op: "Assign"
  input: "training/Adam/Variable_4"
  input: "save/RestoreV2:23"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_4"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_24"
  op: "Assign"
  input: "training/Adam/Variable_5"
  input: "save/RestoreV2:24"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_5"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_25"
  op: "Assign"
  input: "training/Adam/Variable_6"
  input: "save/RestoreV2:25"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_6"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_26"
  op: "Assign"
  input: "training/Adam/Variable_7"
  input: "save/RestoreV2:26"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_7"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_27"
  op: "Assign"
  input: "training/Adam/Variable_8"
  input: "save/RestoreV2:27"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_8"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/Assign_28"
  op: "Assign"
  input: "training/Adam/Variable_9"
  input: "save/RestoreV2:28"
  attr {
    key: "T"
    value {
      type: DT_DOUBLE
    }
  }
  attr {
    key: "_class"
    value {
      list {
        s: "loc:@training/Adam/Variable_9"
      }
    }
  }
  attr {
    key: "use_locking"
    value {
      b: true
    }
  }
  attr {
    key: "validate_shape"
    value {
      b: true
    }
  }
}
node {
  name: "save/restore_all"
  op: "NoOp"
  input: "^save/Assign"
  input: "^save/Assign_1"
  input: "^save/Assign_10"
  input: "^save/Assign_11"
  input: "^save/Assign_12"
  input: "^save/Assign_13"
  input: "^save/Assign_14"
  input: "^save/Assign_15"
  input: "^save/Assign_16"
  input: "^save/Assign_17"
  input: "^save/Assign_18"
  input: "^save/Assign_19"
  input: "^save/Assign_2"
  input: "^save/Assign_20"
  input: "^save/Assign_21"
  input: "^save/Assign_22"
  input: "^save/Assign_23"
  input: "^save/Assign_24"
  input: "^save/Assign_25"
  input: "^save/Assign_26"
  input: "^save/Assign_27"
  input: "^save/Assign_28"
  input: "^save/Assign_3"
  input: "^save/Assign_4"
  input: "^save/Assign_5"
  input: "^save/Assign_6"
  input: "^save/Assign_7"
  input: "^save/Assign_8"
  input: "^save/Assign_9"
}
versions {
  producer: 38
}