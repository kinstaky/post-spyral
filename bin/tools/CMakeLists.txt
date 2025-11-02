add_executable(convert_cut convert_cut.cpp)
target_include_directories(convert_cut PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(convert_cut PRIVATE ROOT::Core ROOT::Graf ROOT::RIO)
