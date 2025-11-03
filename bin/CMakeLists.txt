add_subdirectory(tools)

add_executable(ic_fit ic_fit.cpp)
target_include_directories(ic_fit PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(ic_fit PRIVATE ROOT::Core ROOT::Graf)


add_executable(generate_hardcode_mass generate_hardcode_mass.cpp)
target_link_libraries(generate_hardcode_mass PRIVATE nuclear_data path_manager)

add_executable(missing_mass missing_mass.cpp)
target_include_directories(missing_mass PRIVATE ${ROOT_INCLUDE_DIRS})
target_link_libraries(missing_mass PRIVATE
	solved_event particle
	ROOT::Core ROOT::RIO ROOT::Tree
)

add_executable(search_beam_energy search_beam_energy.cpp)
target_link_libraries(search_beam_energy PRIVATE
	particle solved_event cxxopts
	ROOT::RIO ROOT::Graf ROOT::Hist
)

add_executable(backtrack_beam_energy backtrack_beam_energy.cpp)
target_link_libraries(backtrack_beam_energy PRIVATE
	particle cxxopts
)
