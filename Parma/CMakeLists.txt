cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

project(Parma4 C CXX Fortran)

include(FortranCInterface)
FortranCInterface_VERIFY(CXX)

set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS .f90)
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-range-check")

set(CMAKE_CXX_STANDARD 17)

set(FORTRAN_GENERATOR "src/main-generator_forC.f90 src/subroutines.f90")
set(CPP_LIBFILE " src/Parma.cc ")

set(FSOURCES ${FORTRAN_GENERATOR} ${CPP_LIBFILE})

add_library(Parma4 STATIC)

target_sources(Parma4
    PRIVATE
        src/main-generator_forC.f90
        src/subroutines.f90
        src/Parma.cc)

target_include_directories(Parma4 PUBLIC "include")
