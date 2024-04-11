add_subdirectory(Core)
add_subdirectory(IntergridOp)
add_subdirectory(Solver)
add_library(${MAIN_LIB}
  ${Core_Files}
  ${IntergridOp_Files}
  ${Solver_Files}
)

target_include_directories(${MAIN_LIB} PUBLIC ./)