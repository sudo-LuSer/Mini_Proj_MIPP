cmake_minimum_required(VERSION 3.10)
project(MINI_Proj_Matrix_Operation)

set(CMAKE_CXX_STANDARD 17)

add_executable(Mat_op 
    main.cpp 
    MATH__.cpp
)

target_include_directories(Mat_op PRIVATE 
    /mnt/c/Users/Nasr-allah/Desktop/Mini_Proj_MIPP/include/MIPP/include
)