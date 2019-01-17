#/***************************************************************************
#**                                                                        **
#**  Maurice - a 2D deformable mapping/animation program                   **
#**  Copyright (C) 2014-2016 Graphics Research Group, IIIT Delhi           **
#**                                                                        **
#**  This program is free software: you can redistribute it and/or modify  **
#**  it under the terms of the GNU General Public License as published by  **
#**  the Free Software Foundation, either version 3 of the License, or     **
#**  (at your option) any later version.                                   **
#**                                                                        **
#**  This program is distributed in the hope that it will be useful,       **
#**  but WITHOUT ANY WARRANTY; without even the implied warranty of        **
#**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         **
#**  GNU General Public License for more details.                          **
#**                                                                        **
#**  You should have received a copy of the GNU General Public License     **
#**  along with this program.  If not, see http://www.gnu.org/licenses/.   **
#**                                                                        **
#****************************************************************************
#**           Author: Ojaswa Sharma                                        **
#**           E-mail: ojaswa@iiitd.ac.in                                   **
#**           Date  : 16.07.2015                                           **
#****************************************************************************/

#QMAKE_CC = /usr/local/bin/clang-omp
#QMAKE_CXX = /usr/local/bin/clang-omp++
#QMAKE_LINK = /usr/local/bin/clang-omp
#QMAKE_LINK_SHLIB = /usr/local/bin/clang-omp

#QMAKE_CFLAGS += -fopenmp -Ofast
#QMAKE_CXXFLAGS += -fopenmp -Ofast
#QMAKE_LFLAGS += -fopenmp
QMAKE_CFLAGS += -Ofast
QMAKE_CXXFLAGS += -Ofast
QMAKE_LFLAGS +=

QT += core gui opengl
TEMPLATE = app
TARGET = Maurice

#CONFIG += debug
#DEFINES += NDEBUG

#INCLUDEPATH += /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.11.sdk/usr/include
INCLUDEPATH += . src/ui src/utils src/algorithm
INCLUDEPATH += /usr/local/include/eigen3 /usr/local/include/opencv /usr/local/include 
#INCLUDEPATH += /usr/local/include/libiomp/
INCLUDEPATH += extern/poly2tri

LIBS += `/usr/local/bin/pkg-config --libs opencv`

OBJECTS_DIR=build
MOC_DIR=build
UI_DIR=build
RCC_DIR=build

# Input
HEADERS +=  src/ui/mainwindow.h \
            src/ui/openglwidget.h \
            src/algorithm/shape.h \
            extern/poly2tri/poly2tri.h \
            extern/poly2tri/common/shapes.h \
            extern/poly2tri/common/utils.h \
            extern/poly2tri/sweep/advancing_front.h \
            extern/poly2tri/sweep/cdt.h \
            extern/poly2tri/sweep/sweep_context.h \
            extern/poly2tri/sweep/sweep.h \
	    src/utils/gl_utils.h \
	    src/utils/shader_utils.h \
    src/ui/cagedialog.h

FORMS +=    src/ui/mainwindow.ui \
    src/ui/cagedialog.ui

SOURCES +=  src/main.cpp \
            src/ui/mainwindow.cpp \
            src/ui/openglwidget.cpp \
            src/algorithm/shape.cpp \
            extern/poly2tri/common/shapes.cc \
            extern/poly2tri/sweep/advancing_front.cc \
            extern/poly2tri/sweep/cdt.cc \
            extern/poly2tri/sweep/sweep_context.cc \
            extern/poly2tri/sweep/sweep.cc \
	    src/utils/gl_utils.cpp \
	    src/utils/shader_utils.cpp \
    src/ui/cagedialog.cpp

RESOURCES += res/icons.qrc res/shaders.qrc

macx{
    CONFIG-=app_bundle
}
