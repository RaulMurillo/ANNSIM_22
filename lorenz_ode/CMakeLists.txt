project("lorenz_ode")

# file (GLOB SOURCES "src/*.cpp")
# find_package(UNIVERSAL CONFIG REQUIRED)
# Include universal (posits) ##############################################
include_directories("/homelocal/raul_local/Documents/libs/universal/include")

add_executable(${PROJECT_NAME} src/lorenz_ode.cpp)
# target_link_libraries(${PROJECT_NAME} UNIVERSAL::UNIVERSAL)
set_property(TARGET  ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

