project("brownian_motion")

# file (GLOB SOURCES "src/*.cpp")
# find_package(UNIVERSAL CONFIG REQUIRED)
# Include universal (posits) ##############################################
include_directories("~/universal-3.52/include")


#add_executable(${PROJECT_NAME} src/brownian_motion_test.cpp src/brownian_motion.cpp)
add_executable(${PROJECT_NAME} src/brownian_motion_test.cpp)
# target_link_libraries(${PROJECT_NAME} UNIVERSAL::UNIVERSAL)
set_property(TARGET  ${PROJECT_NAME} PROPERTY CXX_STANDARD 17)

