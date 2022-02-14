project("molecular_dynamics")

# file (GLOB SOURCES "src/*.cpp")
# find_package(UNIVERSAL CONFIG REQUIRED)
# Include universal (posits) ##############################################
include_directories("~/universal-3.48/include")


add_executable(${PROJECT_NAME} src/molecular_dynamics.cpp)
# target_link_libraries(${PROJECT_NAME} UNIVERSAL::UNIVERSAL)
set_property(TARGET  ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

