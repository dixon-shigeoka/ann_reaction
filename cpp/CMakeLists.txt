cmake_minimum_required(VERSION 2.8)
project(1dsteady)
#enable_language(Fortran)

# Create executable file
#add_executable(1dsteady
#    predict.cpp
#    tf_utils.cpp
#    tf_utils.hpp
#    1dsteady.f90
#)

# Create object file without linking
add_library(Predictor
  ann_integrator.cpp  #for fortran main
  ann_integrator.hpp
  #main.cpp               #for c++ main
  #main.hpp
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
#target_link_libraries(Predictor ${PROJECT_SOURCE_DIR}/libtensorflow/libtensorflow-cpu-linux-x100_64-1.14.0/lib/libtensorflow.so.1)
target_include_directories(Predictor PUBLIC ${PROJECT_SOURCE_DIR}/libtensorflow/libtensorflow-cpu-linux-x86_64-1.14.0/include)
file(COPY ${PROJECT_SOURCE_DIR}/libtensorflow/libtensorflow-cpu-linux-x86_64-1.14.0/lib/libtensorflow.so.1 DESTINATION ${PROJECT_BINARY_DIR})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -c")

#Copy resouce
file(COPY ${CMAKE_SOURCE_DIR}/resource/ DESTINATION ${PROJECT_BINARY_DIR}/resource/)
add_definitions(-DRESOURCE_DIR="${PROJECT_BINARY_DIR}/resource/")
