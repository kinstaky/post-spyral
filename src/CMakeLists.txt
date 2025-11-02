add_subdirectory(event)

add_library(path_manager path_manager.cpp)
target_include_directories(path_manager PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(path_manager PUBLIC toml11::toml11)

add_library(nuclear_data nuclear_data.cpp)
target_include_directories(nuclear_data PUBLIC ${PROJECT_SOURCE_DIR})
target_link_libraries(nuclear_data PUBLIC path_manager)

add_library(particle particle.cpp)
target_include_directories(particle PUBLIC ${ROOT_INCLUDE_DIRS})
target_link_libraries(particle PUBLIC nuclear_data catima ROOT::Core)
