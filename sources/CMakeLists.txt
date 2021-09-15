cmake_minimum_required (VERSION 3.8)

project ("Focusinterpolation")


cmake_minimum_required (VERSION 3.8)


###################################################################

# Change this path to your Open-CV folder in which the cmake file OpenCVConfig.cmake can be found.
# Note that even if you use windows, cmake wants / slashes, like on Unixoid systems.
# A UNIX example for a path would be
set(OpenCV_DIR "/opt/local/libexec/opencv4/cmake")

# A Windows example would be
#set(OpenCV_DIR "C:/opencv/build")

# Change this path to the folder where the fitsio.h file from cfitsio can be found
# Note that even if you use windows, cmake wants / slashes, like on Unixoid systems.
# A UNIX example for a path would be
set(Cfitsio_INCLUDE_DIR "/opt/local/include")
# A Windows example would be

#set(Cfitsio_INCLUDE_DIR "C:/Users/benja/source/repos/cfitsio/cfitsio-3.49")


# Change this path to the folder where the fitsio.lib file can be found
# Note that even if you use windows, cmake wants / slashes, like on Unixoid systems.
# A UNIX example for a path would be
set(Cfitsio_LIB_DIR "/opt/local/lib")

# A Windows example would be
#set(Cfitsio_LIB_DIR "C:/Users/benja/source/repos/cfitsio/cfitsio-3.49/cfitsio.build/Release")



# No further changes to the cmake file should be necessary.
# However, do not forget to copy the libraries Open-CV and fitsio into the folder of the executables.
# the focusinterpolation library must also be copied into the same
# directory where the executables are.


###################################################################
if (${UNIX})
    set(UNIXOID TRUE)
    add_definitions(-DUNIXOID)
endif()

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    add_definitions(-DMACOSX)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


find_path(CFITSIO_INCLUDES fitsio.h ${Cfitsio_INCLUDE_DIR})

find_package(OpenCV REQUIRED )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${Cfitsio_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS})

find_library (CFITSIO_LIBRARY NAMES cfitsio PATHS ${Cfitsio_LIB_DIR})

include(GenerateExportHeader)

ADD_LIBRARY(focusinterpolation  SHARED focusinterpolation.cpp focusinterpolation.h focusinterpolation_exports.h)
TARGET_LINK_LIBRARIES(focusinterpolation ${CFITSIO_LIBRARY} ${OpenCV_LIBRARIES})


GENERATE_EXPORT_HEADER(focusinterpolation           
    BASE_NAME focusinterpolation
    EXPORT_MACRO_NAME FOCUSINTERPOLATION_API
    EXPORT_FILE_NAME focusinterpolation_exports.h
    STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)

ADD_EXECUTABLE(FM FM.cpp)
TARGET_LINK_LIBRARIES(FM focusinterpolation ${CFITSIO_LIBRARY})

ADD_EXECUTABLE(HFDM HFDM.cpp)
TARGET_LINK_LIBRARIES(HFDM focusinterpolation ${CFITSIO_LIBRARY})

