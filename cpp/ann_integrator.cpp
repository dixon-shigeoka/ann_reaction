#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
// #include <opencv2/opencv.hpp> //
#include <tensorflow/c/c_api.h>
#include "tf_utils.hpp"

//#define MODEL_FILENAME RESOURCE_DIR"frozen_graph.pb"
#define MODEL_FILENAME RESOURCE_DIR"stanford_model.pb"


static int displayGraphInfo()
{
    TF_Graph *graph = tf_utils::LoadGraphDef(MODEL_FILENAME);
    if (graph == nullptr) {
        std::cout << "Can't load graph" << std::endl;
        return 1;
    }

    size_t pos = 0;
    TF_Operation* oper;
    printf("--- graph info ---\n");
    while ((oper = TF_GraphNextOperation(graph, &pos)) != nullptr) {
        printf("%s\n", TF_OperationName(oper));
    }
    printf("--- graph info ---\n");

    TF_DeleteGraph(graph);
    return 0;
}



void ann_integrator(double atmp, double aprs, double aYi[])
{
    // printf("Hello from TensorFlow C library version %s\n", TF_Version());


    // read the statistic data from csv files //
    std::ifstream ifs("./resource/statistics.csv");
    std::string line;
    const std::string delim = ",";

    int row = 0;
    int col = 0;
    std::string field;
    double statistics[2][10];

    while ( std::getline(ifs, line) ) {
      std::istringstream stream(line);
      std::vector<std::string> result;
      std::vector<std::string> strvec;
      while (getline(stream,field,',')){
        result.push_back(field);
        strvec = result;
      }
      for (int i=0; i<strvec.size();i++){
        statistics[col][i] = std::stod(strvec.at(i));
      }
      col++;
    }

    //std::vector<double> test_x = {1360.0, 151987.5, 1.11900533e-01, 8.88099467e-01,
    //                              2.79751332e-14, 4.44049734e-13, 4.72024867e-13,
    //                              5.00000000e-13, 9.16074600e-13, 9.44049734e-13};

    std::vector<double> test_x = {atmp, aprs, aYi}

    // preprocessing for input data //
    for (int i=0; i<test_x.size(); i++){
      test_x[i] = log(test_x[i]);
      test_x[i] = (test_x[i] - statistics[0][i]) / statistics[1][i];
    }

    // get graph info //
    displayGraphInfo();

    TF_Graph *graph = tf_utils::LoadGraphDef(MODEL_FILENAME);
    if (graph == nullptr) {
        std::cout << "Can't load graph" << std::endl;
        return 1;
    }

    // prepare input tensor //
    TF_Output input_op = { TF_GraphOperationByName(graph, "dense_1_input"), 0 };
    if (input_op.oper == nullptr) {
        std::cout << "Can't init input_op" << std::endl;
        return 2;
    }


    const std::vector<int64_t> input_dims = { 1, 10 };
    std::vector<double> input_vals;
    input_vals = test_x;
    //image.reshape(0, 1).copyTo(input_vals); // Mat to vector

    TF_Tensor* input_tensor = tf_utils::CreateTensor(TF_DOUBLE,
        input_dims.data(), input_dims.size(),
        input_vals.data(), input_vals.size() * sizeof(double));

    // prepare output tensor //
    TF_Output out_op = { TF_GraphOperationByName(graph, "activation_3/Sigmoid"), 0 };
    if (out_op.oper == nullptr) {
        std::cout << "Can't init out_op" << std::endl;
        return 3;
    }

    TF_Tensor* output_tensor = nullptr;

    // prepare session //
    TF_Status* status = TF_NewStatus();
    TF_SessionOptions* options = TF_NewSessionOptions();
    TF_Session* sess = TF_NewSession(graph, options, status);
    TF_DeleteSessionOptions(options);

    if (TF_GetCode(status) != TF_OK) {
        TF_DeleteStatus(status);
        return 4;
    }

    //run session //
    TF_SessionRun(sess,
        nullptr, // Run options.
        &input_op, &input_tensor, 1, // Input tensors, input tensor values, number of inputs.
        &out_op, &output_tensor, 1, // Output tensors, output tensor values, number of outputs.
        nullptr, 0, // Target operations, number of targets.
        nullptr, // Run metadata.
        status // Output status.
    );

    if (TF_GetCode(status) != TF_OK) {
      std::cout << TF_GetCode(status) << std::endl;
      std::cout << "Error run session";
      TF_DeleteStatus(status);
      return 5;
    }

    TF_CloseSession(sess, status);
    if (TF_GetCode(status) != TF_OK) {
        std::cout << "Error close session";
        TF_DeleteStatus(status);
        return 6;
    }

    TF_DeleteSession(sess, status);
    if (TF_GetCode(status) != TF_OK) {
        std::cout << "Error delete session";
        TF_DeleteStatus(status);
        return 7;
    }

    const auto probs = static_cast<double*>(TF_TensorData(output_tensor));
    std::cout << TF_TensorData(output_tensor) << std::endl;

    //for (int i = 0; i < 8; i++) {
    //    printf("prob of %d: %.19lf\n", i, probs[i]);
    //}

    TF_DeleteTensor(input_tensor);
    TF_DeleteTensor(output_tensor);
    TF_DeleteGraph(graph);
    TF_DeleteStatus(status);

    //return probs[];
}
