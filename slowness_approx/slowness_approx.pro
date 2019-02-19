######################################################################
# Automatically generated by qmake (2.01a) Fri Jan 28 19:25:48 2011
######################################################################

TEMPLATE = app
TARGET = slowness_approx
DEPENDPATH += . ../src2 ../src
INCLUDEPATH += . ../src2 ../src e:/work/boost e:/work/gnuwin32/include/ /usr/include/eigen3 d:/work/boost d:/work/gnuwin32/include/ d:/work/eigen3 d:/work/fftw/include
QT += opengl
QMAKE_CXXFLAGS += -O3 -std=c++0x
OBJECTS_DIR = ../tmp
LIBS += -Le:/work/gnuwin32/bin -Ld:/work/GnuWin32/bin -Ld:/work/fftw/lib -lgsl -lfftw3 -lgslcblas
DEFINES += SLOWNESS_WORK
CONFIG += depend_includepath

# Input
HEADERS += ../src2/composite_wave.h\
           data_sphere.h \
           gl_view.h \
           main_dialog.h \
           triang_ind.h \
           gl_utils.h \
           world.h \
           ../src2/util.h \
           ../src2/plan_fftw.h\
           ../src2/plane_wave.h\
           ../src2/plane_wave_c.h\
           ../src2/poly.h \
           ../src2/vec3.h \
           height_calculator.h\
           ../src2/mat3.h\
           GraphWidget.h\
           ../src2/piezo_tensor.h\
           ../src2/matrix_fftw.h\
           ../src2/wave_matrix.h

FORMS += main_dialog.ui

SOURCES += GraphWidget.cpp\
           data_sphere.cpp \
           gl_utils.cpp \
           gl_view.cpp \
           height_calculator.cpp \
           main.cpp \
           main_dialog.cpp \
           ../src2/mat3.cpp\
           ../src2/matrix_fftw.cpp\
           ../src2/piezo_tensor.cpp\
           ../src2/plan_fftw.cpp \
           ../src2/plane_wave.cpp \
           ../src2/plane_wave_c.cpp\
           ../src2/poly.cpp \
           ../src2/povray_export.cpp \
           ../src2/storage.cpp \
           triang_ind.cpp \
           ../src2/util.cpp \
           ../src2/vec3.cpp \
           ../src2/vec3c.cpp\
           ../src2/wave_matrix.cpp \
           world.cpp \
           ../src2/composite_wave.cpp