cmake_minimum_required(VERSION 3.20)
project(parser VERSION 1.0)
set(CMAKE_OSX_ARCHITECTURES "arm64")

# Percorso di CPLEX
set(CPLEXDIR "/Applications/CPLEX_Studio2211/cplex")

# Include directory (per permettere a CMake di trovare cplex.h)
include_directories("${CPLEXDIR}/include" "${CPLEXDIR}/include/ilcplex")

# Flags di compilazione e linking
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -I${CPLEXDIR}/include -I${CPLEXDIR}/include/ilcplex")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L${CPLEXDIR}/lib/arm64_osx/static_pic -lcplex -lpthread -ldl")

# Crea l'eseguibile
add_executable(parser src/main.c src/tsp_utils.c src/chrono.c src/heuristics.c src/solver.c src/callback.c src/cpx_utils.c)

# Link a CPLEX
target_link_libraries(parser 
    cplex
    libcplex.a  
    pthread 
    dl 
    m
)


include(CPack)