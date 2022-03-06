project("molecular_dynamics_openmp")

# file (GLOB SOURCES "src/*.cpp")
# find_package(UNIVERSAL CONFIG REQUIRED)
# Include universal (posits) ##############################################
include_directories("/homelocal/raul_local/Documents/libs/universal/include")
include_directories("~/universal-3.48/include")
find_package(OpenMP)


add_executable(${PROJECT_NAME} src/molecular_dynamics.cpp)
# target_link_libraries(${PROJECT_NAME} UNIVERSAL::UNIVERSAL)
target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
set_property(TARGET  ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

