cmake_minimum_required(VERSION 3.5.0)

project(Maurice)

set(TARGET ${CMAKE_PROJECT_NAME})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR})

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)


#Dependencies
find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Core)
if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
	find_package(Qt6 COMPONENTS Core Widgets OpenGLWidgets Gui)
else()
	find_package(Qt5 COMPONENTS Core Widgets OpenGL Gui REQUIRED)
endif()



find_package(OpenGL REQUIRED)
find_package(GLM REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP)


set(SOURCES
	src/main.cpp
	src/ui/mainwindow.cpp
	src/ui/openglwidget.cpp
	src/algorithm/shape.cpp
	extern/poly2tri/common/shapes.cc
	extern/poly2tri/sweep/advancing_front.cc
	extern/poly2tri/sweep/cdt.cc
	extern/poly2tri/sweep/sweep_context.cc
	extern/poly2tri/sweep/sweep.cc
	src/ui/cagedialog.cpp
	)

set(HEADERS
	src/ui/mainwindow.h
	src/ui/openglwidget.h
	src/algorithm/shape.h
	extern/poly2tri/poly2tri.h
	extern/poly2tri/common/shapes.h
	extern/poly2tri/common/utils.h
	extern/poly2tri/sweep/advancing_front.h
	extern/poly2tri/sweep/cdt.h
	extern/poly2tri/sweep/sweep_context.h
	extern/poly2tri/sweep/sweep.h
	src/ui/cagedialog.h
	)

set(UI_SOURCES
	src/ui/mainwindow.ui
	src/ui/cagedialog.ui
	)

set(RESOURCES
	res/icons.qrc
	res/shaders.qrc
	)

add_executable(${TARGET} ${SOURCES} ${HEADERS} ${UI_SOURCES} ${RESOURCES})

target_include_directories(${TARGET} PRIVATE
	${CMAKE_SOURCE_DIR}/src/ui
	${CMAKE_SOURCE_DIR}/src/algorithm
	${CMAKE_SOURCE_DIR}/extern/poly2tri
	${GLM_INCLUDE_DIRS}/../include
	)

target_link_libraries(${TARGET} LINK_PUBLIC ${OpenCV_LIBS} Eigen3::Eigen glm::glm)
if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
	target_link_libraries(${TARGET} LINK_PUBLIC Qt6::Core Qt6::Widgets Qt6::Gui Qt6::OpenGLWidgets)
else()
	target_link_libraries(${TARGET} LINK_PUBLIC Qt5::Core Qt5::Widgets Qt5::Gui Qt5::OpenGL)
endif()
if(OpenMP_CXX_FOUND)
	target_link_libraries(${TARGET} PRIVATE $<$<CONFIG:Release>:OpenMP::OpenMP_CXX>)
endif()

