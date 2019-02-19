TEMPLATE = app
TARGET = 
DEPENDPATH += . ../src2 
INCLUDEPATH += . ../src2 e:/work/qwt6/src e:/work/boost
QMAKE_CXXFLAGS += -O3
QT += opengl
LIBS += -Le:/work/qwt6/lib -Le:/work/boost/lib -L/usr/lib -lgsl -lfftw3
OBJECTS_DIR = ../tmp

unix: LIBS += -lqwt6
win32: LIBS += -lqwt -Le://work//GnuWin32//lib

win32: INCLUDEPATH += e://work//GnuWin32//include

# Input
SOURCES += main.cpp\
           poly.cpp\
           util.cpp\
           vec3.cpp\
           povray_export.cpp\
           mat3.cpp\
           piezo_tensor.cpp\

HEADERS += errors.h\
           poly.h\
           vec3.h\
           povray_export.h\
           util.h\
           plane_wave.h\
           piezo_tensor.h\
           matrix_fftw.h
