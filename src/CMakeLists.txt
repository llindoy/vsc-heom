cmake_minimum_required (VERSION 3.0)
project(vsc-heom LANGUAGES C CXX)

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake_modules/ ${CMAKE_MODULE_PATH})

option(USE_OPENMP "Compile with openmp" ON)
option(BUILD_THERMALISE "Compile the executables that thermalise the ADO tensor before evaluating the correlation functions used in the evaluation of the rate constant." OFF)
option(BUILD_IR "Compile the executables for evaluating the dipole-dipole correlation function" OFF)

if(USE_OPENMP)
    add_compile_definitions(USE_OPENMP)
endif()


if(COMPILER_SET)
else()
    if(MSVC)
        message(FATAL_ERROR "Request MSVC compiler.  I don't know how this works so exiting.")
    else()
        set(CMAKE_CXX_FLAGS "-O3 -Wall -g -Wextra -Wshadow -pedantic -Wold-style-cast ${CMAKE_CXX_FLAGS}")
        option(WITH_WERROR "Compile with '-Werror' C++ compiler flag" ON)
        if(WITH_WERROR)
            #set(CMAKE_CXX_FLAGS "-Werror ${CMAKE_CXX_FLAGS}")
        endif(WITH_WERROR)


        #additional compilation flags for GNU compilers
        if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto -fwhole-program -fcx-fortran-rules")
            if(USE_OPENMP)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
            endif()
        #additional compilation flags for inteal compilers
        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report=3 -ipo -O3 -xHost -restrict")
            message(STATUS "using intel compiler.")
            if(USE_OPENMP)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
            endif()

        elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "IntelLLVM")
            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -qopt-report=3 ")
            if(USE_OPENMP)
                message(STATUS "Using Openmp")
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -openmp")
            endif()

        elseif ( ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang") )
            find_package(LLVM REQUIRED)
            find_package(Clang REQUIRED)

            set(LINALG_LINALG ${LIBS_LIBS} ${LLVM_LIBRARIES} ${CLANG_LIBS})
            

            #if we are using clang then find threads
            set(CMAKE_THREAD_LIBS_INIT "-lpthread")
            set(CMAKE_HAVE_THREADS_LIBRARY 1)
            set(CMAKE_USE_WIN32_THREADS_INIT 0)
            set(CMAKE_USE_PTHREADS_INIT 1)
            set(THREADS_PREFER_PTHREAD_FLAG ON)

            set(CMAKE_CXX_FLAGS "-stdlib=libc++ ${CMAKE_CXX_FLAGS}")
            set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -stdlib=libc++ ")


            if(USE_OPENMP)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp=libomp")
            endif()
        endif()
        

        if(CMAKE_VERSION VERSION_LESS 3.1)
            set(CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS}")
        else()
          if(NOT DEFINED CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD STREQUAL "98")
            set(CMAKE_CXX_STANDARD 17)
          endif()

          if(CMAKE_CXX_STANDARD GREATER 17)
            cmake_minimum_required(VERSION 3.8)
          endif()

          set(CMAKE_CXX_STANDARD_REQUIRED ON)
        endif()
    endif()
    set(COMPILER_SET TRUE INTERNAL BOOL "")
endif()


add_executable(heom_mol.x rate_theory/main.cpp)
add_executable(heom_cavity.x rate_theory/main_cavity.cpp)
add_executable(heom_cavity_loss.x rate_theory/main_cavity_loss.cpp)
add_executable(heom_mol_2d.x rate_theory/main_2d.cpp)
add_executable(heom_cavity_loss_2d.x rate_theory/main_cavity_loss_2d.cpp)

if(BUILD_THERMALISE)
add_executable(heom_cavity_loss_thermalise.x rate_theory/main_cavity_loss_thermalise.cpp)
add_executable(heom_mol_thermalise.x rate_theory/main_thermalise.cpp)
endif()

if(BUILD_IR)
add_executable(heom_mol_ir.x ir_spectrum/main_ir.cpp)
add_executable(heom_cavity_loss_ir.x ir_spectrum/main_cavity_loss_ir.cpp)
add_executable(heom_mol_2d_ir.x ir_spectrum/main_2d_ir.cpp)
add_executable(heom_cavity_loss_2d_ir.x ir_spectrum/main_cavity_loss_ir_2d.cpp)

if(BUILD_THERMALISE)
add_executable(heom_mol_ir_thermalise.x ir_spectrum/main_ir_therm.cpp)
add_executable(heom_cavity_loss_ir_thermalise.x ir_spectrum/main_cavity_loss_ir_therm.cpp)
endif()
endif()

include(FetchContent)
#now make sure we have the rapidjson library 
#now we can also set up rapidjson
set(RAPIDJSON_ROOT_DIR ${CMAKE_SOURCE_DIR}/external/rapidjson)
if(EXISTS ${RAPIDJSON_ROOT_DIR})
    SET(rapidjson_repo "file://${RAPIDJSON_ROOT_DIR}/")
else()
    SET(rapidjson_repo "https://github.com/Tencent/rapidjson")
endif()

FetchContent_Declare(rapidjson
    GIT_REPOSITORY ${rapidjson_repo}
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
    EXCLUDE_FROM_ALL
)

FetchContent_Populate(rapidjson)

add_definitions(-DRAPIDJSON_LIBRARY_FOUND)
set(EOS_INCLUDE_DIRECTORIES "${rapidjson_SOURCE_DIR}/include/;${EOS_INCLUDE_DIRECTORIES}" CACHE INTERNAL "")


if(DEFINED BLA_VENDOR)
    if(BLA_VENDOR MATCHES "Intel")
        find_package(LAPACK)
        find_package(BLAS)
        if(BLAS_FOUND)
            set(MKL_INCLUDE_DIRS ${MKLROOT}/include)
            #FIND_PATH(MKL_INCLUDE_DIR "mkl_cblas.h" PATHS ${MKLROOT}/include)
            message(STATUS ${MKL_INCLUDE_DIRS})
            add_compile_definitions(USE_MKL)
            add_compile_definitions(BLAS_HEADER_INCLUDED)
            add_compile_definitions(NOT_ALLOWS_RETURN_TYPE )
            add_compile_definitions(BLAS_NO_TRAILING_UNDERSCORE )
            if(DLA_VENDOR MATCHES "64ilp")
                add_compile_definitons(BLAS_64_BIT)
            endif()
            set(LINALG_LIBS ${LINALG_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
            include_directories(SYSTEM ${MKL_INCLUDE_DIRS})
            set(USE_SPARSE_BLAS OFF)
        else()
            set(USE_DEFAULT_FIND_BLAS ON)
        endif()
    elseif(BLA_VENDOR MATCHES "OpenBLAS")   
        add_compile_definitions(USE_OPENBLAS)
        set(USE_DEFAULT_FIND_BLAS ON)
    elseif(BLA_VENDOR MATCHES "FLAME")
        #here we allow for us to handle the multithreaded case.  This is not dealt with by the standard FindBLAS and FindLapack files.  
        #this has been setup to work with the aocl versions of blis and flame
        if(USE_OPENMP)
            #first we attempt to find the standard FLAME LIBRARY
            set(BLA_VENDOR "FLAME")
            include(FindBLAS)
            if(NOT BLAS_FOUND)
                message(FATAL_ERROR "BLAS NOT FOUND")
            endif()
            check_blas_libraries(
              BLAS_LIBRARIES
              BLAS
              sgemm
              ""
              "blis-mt"
              ""
              ""
              ""
              )
            message(STATUS ${BLAS_LIBRARIES})
            set(BLA_VENDOR "FLAME")
            find_package(LAPACK REQUIRED)
            set(LINALG_LIBS ${LINALG_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
            set(USE_DEFAULT_FIND_BLAS OFF)
        else()
            add_compile_definitions(USE_FLAME)
            find_package(BLAS REQUIRED)
            find_package(LAPACK REQUIRED)
            message(STATUS ${BLAS_LIBRARIES})
            set(LINALG_LIBS ${LINALG_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
            set(USE_DEFAULT_FIND_BLAS OFF)
        endif()
    else()
        find_package(BLAS REQUIRED)
        find_package(LAPACK REQUIRED)
        message(STATUS ${BLAS_LIBRARIES})
        set(LINALG_LIBS ${LINALG_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
        set(USE_DEFAULT_FIND_BLAS OFF)
    endif()
else()
    set(USE_DEFAULT_FIND_BLAS ON)
endif()

#if we aren't dealing with one of the vendors we treat specially then we will just go ahead and use the default
#findBLAS and findLAPACK calls
if(USE_DEFAULT_FIND_BLAS)
    find_package(BLAS)
    find_package(LAPACK)

    if(NOT BLAS_FOUND)
        set(BLAS_LIBRARIES "-lblas")
    endif()



    if(NOT LAPACK_FOUND)
        set(LAPACK_LIBRARIES "-llapack -lm")
    endif()
    set(LINALG_LIBS ${LINALG_LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()



target_include_directories(heom_mol_2d.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_mol_2d.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_mol_2d.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_mol_2d.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_cavity_loss_2d.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity_loss_2d.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity_loss_2d.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity_loss_2d.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_mol.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_mol.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_mol.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_mol.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_cavity.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)
 
target_include_directories(heom_cavity_loss.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity_loss.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity_loss.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity_loss.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)



if(BUILD_THERMALISE)
target_include_directories(heom_mol_thermalise.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_mol_thermalise.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_mol_thermalise.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_mol_thermalise.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_cavity_loss_thermalise.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity_loss_thermalise.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity_loss_thermalise.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity_loss_thermalise.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)
endif()

if(BUILD_IR)

target_include_directories(heom_mol_ir.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_mol_ir.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_mol_ir.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_mol_ir.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_cavity_loss_2d_ir.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity_loss_2d_ir.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity_loss_2d_ir.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity_loss_2d_ir.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_mol_2d_ir.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_mol_2d_ir.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_mol_2d_ir.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_mol_2d_ir.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_cavity_loss_ir.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity_loss_ir.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity_loss_ir.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity_loss_ir.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

if(BUILD_THERMALISE)
target_include_directories(heom_mol_ir_thermalise.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_mol_ir_thermalise.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_mol_ir_thermalise.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_mol_ir_thermalise.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)

target_include_directories(heom_cavity_loss_ir_thermalise.x PRIVATE ${EOS_INCLUDE_DIRECTORIES})
target_include_directories(heom_cavity_loss_ir_thermalise.x PRIVATE ${LINALG_INCLUDE_DIRECTORIES})
target_link_libraries(heom_cavity_loss_ir_thermalise.x linalg  ${LINALG_LIBS} ${LIBS}) 
install(TARGETS heom_cavity_loss_ir_thermalise.x 
        DESTINATION ${CMAKE_SOURCE_DIR}/bin)
endif()
endif()


