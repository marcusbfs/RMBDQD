
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# Optimizações para UNIX OS e compilador GNU GFortran

# Cuidado ao atualizar gfortran!!!
# Pode que quebrar alguma flag de otimização

if (UNIX)
    if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
        message(STATUS "Unix OS and GNU GFortran compiler")
        set(common "-Wall -Wextra -Wimplicit-interface -fmax-errors=1")
        set(par "-cpp -floop-nest-optimize")
        set (CMAKE_Fortran_FLAGS "-O3 ${par} -fno-f2c -march=native -ffast-math -funroll-loops")
        set (CMAKE_Fortran_FLAGS_RELEASE "-fno-f2c -O3")
        set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g")
        message(STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
    endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
endif (UNIX)
message(STATUS "The host processor is ${CMAKE_HOST_SYSTEM_PROCESSOR}")


# Diretório com as subrotinas
set(libdir "mod")

ADD_LIBRARY(lib ${libdir}/decomp.f90     ${libdir}/jacobi.f90    ${libdir}/mult.f90
                ${libdir}/redbak.f90     ${libdir}/scheck.f90    ${libdir}/commonutils.f90
                ${libdir}/subrotinas.f90 ${libdir}/subspace.f90  ${libdir}/help.f90
                ${libdir}/f90getopt.f90 )

ADD_LIBRARY(lcpp ${libdir}/get_os.cpp)
ADD_LIBRARY(libf ${libdir}/frequencia.f90)
ADD_EXECUTABLE ( ${PROJECT_NAME} main.f90 )

TARGET_LINK_LIBRARIES( libf lib lcpp)
TARGET_LINK_LIBRARIES( ${PROJECT_NAME} lib libf lcpp )

install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin)
