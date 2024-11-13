cmake_minimum_required(VERSION 3.10)
project(MyFortranProject Fortran)
add_executable(myprogram global_parameters.f90 utility_routines.f90 main.f90)
