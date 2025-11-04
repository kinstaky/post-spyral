add_executable(test_energy_loss energy_loss.cpp)
target_link_libraries(test_energy_loss PRIVATE particle)
add_test(NAME energy_loss_test COMMAND test_energy_loss)

add_executable(test_scatter scatter.cpp)
target_link_libraries(test_scatter PRIVATE particle)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_link_libraries(test_scatter PRIVATE ROOT::RIO ROOT::Graf)
endif()
add_test(NAME scatter_test COMMAND test_scatter)

add_executable(test_breakup breakup.cpp)
target_link_libraries(test_breakup PRIVATE particle)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_link_libraries(test_breakup PRIVATE ROOT::RIO ROOT::Graf)
endif()
add_test(NAME breakup_test COMMAND test_breakup)


file(
	COPY ${PROJECT_SOURCE_DIR}/test/data/
	DESTINATION ${PROJECT_BINARY_DIR}/test/data/
)