cmake_minimum_required (VERSION 3.12)

project ("Focusinterpolation")



###################################################################

## comment out if you do not want to use the beta variable in the Levenberg-Marquardt fit

# Change this path to your Open-CV folder in which the cmake file OpenCVConfig.cmake can be found.
# Note that even if you use windows, cmake wants / slashes, like on Unixoid systems.
# A UNIX example for a path would be
#set(OpenCV_DIR "/opt/local/libexec/opencv4/cmake")

# A Windows example would be
set(OpenCV_DIR "C:/Users/benja/source/repos/OpenCV/build/")

# Change this path to the folder where the fitsio.h file from cfitsio can be found
# Note that even if you use windows, cmake wants / slashes, like on Unixoid systems.
# A UNIX example for a path would be
# set(Cfitsio_INCLUDE_DIR "/opt/local/include")
# A Windows example would be

set(CFITSIO_INCLUDE_DIR "C:/Users/benja/source/repos/cfitsio/cfitsio-4.0.0/")


# Change this path to the folder where the fitsio.lib file can be found
# Note that even if you use windows, cmake wants / slashes, like on Unixoid systems.
# A UNIX example for a path would be
#set(Cfitsio_LIB_DIR "/opt/local/lib")

# A Windows example would be
set(CFITSIO_LIB_DIR "C:/Users/benja/source/repos/cfitsio/cfitsio-4.0.0/cfitsio.build/Release/")


# No further changes to the cmake file should be necessary.
# However, do not forget to copy the libraries Open-CV and fitsio into the folder of the executables.
# the focusinterpolation library must also be copied into the same
# directory where the executables are.


###################################################################





find_path(CFITSIO_INCLUDES fitsio.h ${CFITSIO_INCLUDE_DIR})

find_package(OpenCV REQUIRED )

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR} ${CFITSIO_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS})

find_library(CFITSIO_LIBRARY NAMES cfitsio PATHS ${CFITSIO_LIB_DIR})

include(GenerateExportHeader)


ADD_LIBRARY(focusinterpolation  "focusinterpolation.cpp" focusinterpolation.h focusinterpolation_exports.h)



if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    add_definitions(-DMACOSX)
    set_property(TARGET focusinterpolation PROPERTY CXX_STANDARD 14)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
   set(CLANGCOMPILER TRUE)
   add_definitions(-DCLANGCOMPILER)
   set_property(TARGET focusinterpolation PROPERTY CXX_STANDARD 14)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
   set(GNUCOMPILER TRUE)
   add_definitions(-DGNUCOMPILER)
   set_property(TARGET focusinterpolation PROPERTY CXX_STANDARD 14)
endif()

if(MSVC)
    set_property(TARGET focusinterpolation PROPERTY CXX_STANDARD 17)
    set_property(TARGET focusinterpolation PROPERTY CXX_STANDARD_REQUIRED ON)
endif()

if(MSVC)
    target_compile_options(focusinterpolation PUBLIC "/Zc:__cplusplus")
endif()


TARGET_LINK_LIBRARIES(focusinterpolation ${CFITSIO_LIBRARY} ${OpenCV_LIBRARIES})


GENERATE_EXPORT_HEADER(focusinterpolation           
    BASE_NAME focusinterpolation
    EXPORT_MACRO_NAME FOCUSINTERPOLATION_API
    EXPORT_FILE_NAME focusinterpolation_exports.h
    STATIC_DEFINE SHARED_EXPORTS_BUILT_AS_STATIC)

ADD_EXECUTABLE(FM FM.cpp)

if(MSVC)
    target_compile_options(FM PUBLIC "/Zc:__cplusplus")
endif()

set_property(TARGET FM PROPERTY CXX_STANDARD 17)
set_property(TARGET FM PROPERTY CXX_STANDARD_REQUIRED ON)

TARGET_LINK_LIBRARIES(FM focusinterpolation ${CFITSIO_LIBRARY})

ADD_EXECUTABLE(HFDM HFDM.cpp)

if(MSVC)
    target_compile_options(HFDM PUBLIC "/Zc:__cplusplus")
endif()

TARGET_LINK_LIBRARIES(HFDM focusinterpolation ${CFITSIO_LIBRARY})


