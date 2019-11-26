cmake_minimum_required(VERSION 2.8)
project(Predictor)

# Create executable file
#add_executable(Predictor
#    predict.cpp
#    tf_utils.cpp
#    tf_utils.hpp
#)

# Create object file without linking
add_library(Predictor STATIC
    predict.cpp
    tf_utils.cpp
    tf_utils.hpp
)

# For OpenCV
#find_package(OpenCV REQUIRED)
#if(OpenCV_FOUND)
#    target_include_directories(Predictor PUBLIC ${OpenCV_INCLUDE_DIRS})
#    target_link_libraries(Predictor ${OpenCV_LIBS})
#endif()

# For Tensorflow
target_link_libraries(Predictor ${PROJECT_SOURCE_DIR}/libtensorflow/libtensorflow-cpu-linux-x86_64-1.14.0/lib/libtensorflow.so.1)
target_include_directories(Predictor PUBLIC ${PROJECT_SOURCE_DIR}/libtensorflow/libtensorflow-cpu-linux-x86_64-1.14.0/include)
file(COPY ${PROJECT_SOURCE_DIR}/libtensorflow/libtensorflow-cpu-linux-x86_64-1.14.0/lib/libtensorflow.so.1 DESTINATION ${PROJECT_BINARY_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lstdc++ -c")

#Copy resouce
file(COPY ${CMAKE_SOURCE_DIR}/resource/ DESTINATION ${PROJECT_BINARY_DIR}/resource/)
add_definitions(-DRESOURCE_DIR="${PROJECT_BINARY_DIR}/resource/")
