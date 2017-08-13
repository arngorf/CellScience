PROJECT( mt )

##############################################################################
# POLICIES ###################################################################
##############################################################################

cmake_minimum_required(VERSION 2.8.8 FATAL_ERROR)

if (POLICY CMP0053)
  cmake_policy(SET CMP0053 NEW)
endif()

if (POLICY CMP0057)
  cmake_policy(SET CMP0057 NEW)
endif()

##############################################################################
# OPTIONS ####################################################################
##############################################################################

#OPTION(useFFTW "Turn FTTW on and off" ON)
#OPTION(usePython2 "Turn Python2.7 on and off" ON)
#MESSAGE("useFFTW = ${useFFTW}")
#MESSAGE("usePython2 = ${usePython2}")

##############################################################################
# SEARCH FOR SOURCE FILES ####################################################
##############################################################################

file(GLOB_RECURSE SRC_PATH "src/*.cpp")
MESSAGE("Sources: \n ${SRC_PATH}")
file(GLOB_RECURSE QT_RESOURCES "src/*.qrc")
MESSAGE("Sources: \n ${QT_RESOURCES}")

# sources files
Set (PROJECT_SRCS
    ${SRC_PATH}
)

##############################################################################
# FIND LIBARIES ##############################################################
##############################################################################

SET(EIGEN3_INCLUDE_DIR SYSTEM /usr/include/eigen3)

find_package(VTK REQUIRED)
if (VTK_FOUND)
        include(${VTK_USE_FILE})
        message("VTK found, VTK_USE_FILE: ${VTK_USE_FILE}")
else()
        message(ERROR "VTK not found")
endif()

#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)

#find_package(Boost 1.45.0 COMPONENTS system filesystem date_time program_options serialization thread unit_test_framework)

find_package(Qt5Widgets REQUIRED)

set(CMAKE_AUTOMOC ON)

qt5_add_resources(UI_RESOURCES src/application.qrc)

ADD_EXECUTABLE( mt ${PROJECT_SRCS} ${UI_RESOURCES} )

target_compile_features(mt PRIVATE cxx_range_for)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

##############################################################################
# SET FLAGS DIRECTORIES AND DEFINITIONS ######################################
##############################################################################

set(CMAKE_CXX_FLAGS " -O3 -Wall -Wextra -fopenmp -Wno-sign-compare")
SET(CMAKE_EXE_LINKER_FLAGS "-lm -ltiff")

MESSAGE("Linker flags: ${CMAKE_EXE_LINKER_FLAGS}")

INCLUDE_DIRECTORIES (
    /usr/include
    /usr/local/include
    ${EIGEN3_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

##############################################################################
# TARGET LINK LIBARIES #######################################################
##############################################################################

TARGET_LINK_LIBRARIES( mt ${LINK_LIBRARIES} ${VTK_LIBRARIES} Qt5::Widgets ) #${Boost_LIBRARIES}
