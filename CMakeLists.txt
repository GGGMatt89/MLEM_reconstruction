PROJECT(reconstruction)
cmake_minimum_required(VERSION 2.8)
# add the folder containing the dlib folder to the include path
#INCLUDE_DIRECTORIES(./lib)
#include_directories(./src)

ADD_EXECUTABLE(./reconstruction main.cpp bgeot_ftool.cpp aqphy.cpp image.cpp mlem.cpp)

#ADD_EXECUTABLE(./reconstruction main.cpp bgeot_ftool.cpp aqphy.cpp traitement.cpp image.cpp Point.cpp FSV/cone.cpp FSV/ellipse.cpp FSV/intell.cpp FSV/plan.cpp FSV/RayCasting.cpp FSV/utils.cpp FSV/volume.cpp FSV/line_cone.cpp)
#${CMAKEDEMO_SOURCE_DIR}
