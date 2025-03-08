cmake_minimum_required(VERSION 3.20)
project(parser VERSION 1.0)

# Impostazione del percorso di CPLEX
set(CPLEXDIR "/opt/ibm/ILOG/CPLEX_Studio2212/cplex")

# Verifica che CPLEX esista
if (NOT EXISTS "${CPLEXDIR}")
    message(FATAL_ERROR "Il percorso a CPLEX non è valido: ${CPLEXDIR}")
endif()

# Include directory
include_directories("${CPLEXDIR}/include")

# Aggiungi la directory delle librerie (statiche o condivise)
link_directories("${CPLEXDIR}/lib/x86-64_linux/static_pic")  # Usa "shared" se preferisci librerie dinamiche

# Crea l'eseguibile
add_executable(parser src/main.c src/tsp_functions.c src/chrono.c)

# Aggiungi il link alle librerie di CPLEX (statiche o dinamiche)
target_link_libraries(parser 
    cplex
    libcplex.a  
    pthread 
    dl 
    m
)

# Se usi librerie dinamiche (shared), puoi aggiungere:
# target_link_libraries(parser PRIVATE "${CPLEXDIR}/lib/x86-64_linux/shared/libcplex.so")

include(CPack)