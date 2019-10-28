import tensorflow as tf
import os

abspath = os.path.dirname(os.path.abspath(__file__))
model_dir = abspath + '/learned_model/model_for_c++/'
#output_node_names = ['sigmoid/Sigmoid']
output_node_names = "activation_3/Sigmoid"

def main():
    with tf.Session() as sess:
        # Restore the graph
        saver = tf.train.import_meta_graph(model_dir + "model.meta")

        # Load weights
        saver.restore(sess, model_dir + 'model.ckpt')

        # Freeze the graph
        frozen_graph_def = tf.graph_util.convert_variables_to_constants(
            sess,
            sess.graph_def,
            output_node_names.split(",")
            )

        # Save the frozen graph
        with open(model_dir + 'frozen_graph.pb', 'wb') as f:
            f.write(frozen_graph_def.SerializeToString())


if __name__ == '__main__':
    main()
