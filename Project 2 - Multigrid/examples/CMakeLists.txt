add_executable(solve1D solve1D.cpp)
target_link_libraries(solve1D ${MAIN_LIB} jsoncpp)

add_executable(solve2D solve2D.cpp)
target_link_libraries(solve2D ${MAIN_LIB} jsoncpp)
