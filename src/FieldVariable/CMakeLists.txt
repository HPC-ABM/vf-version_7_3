set(FV_SRC FieldVariable.cpp)
file( GLOB USR_FV_SRC Usr_FieldVariables/*.cpp)

set(fv_lib_type SHARED)
message("glob of FV is ${USR_FV_SRC}")
add_library(FieldVariable ${fv_lib_type} ${FV_SRC} ${USR_FV_SRC})

# Added by Nuttiiya
find_package(OpenMP)
if (OPENMP_FOUND)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  add_definitions(-D_OMP)
endif()

## Added by Nuttiiya for atomic flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")