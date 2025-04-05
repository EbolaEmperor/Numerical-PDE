add_executable(solve1D solve1D.cpp)
target_link_directories(solve1D PRIVATE ${JSONCPP_LIBRARY_DIRS})
target_link_libraries(solve1D ${MAIN_LIB} ${JSONCPP_LIBRARIES})

add_executable(solve2D solve2D.cpp)
target_link_directories(solve2D PRIVATE ${JSONCPP_LIBRARY_DIRS})
target_link_libraries(solve2D ${MAIN_LIB} ${JSONCPP_LIBRARIES})
