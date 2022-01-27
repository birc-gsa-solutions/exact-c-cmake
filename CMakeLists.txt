cmake_minimum_required (VERSION 3.21)
project (PROJECT)
enable_testing()

add_executable(naive src/naive.c src/parsers.c src/sam.c)
add_executable(lin src/lin.c src/parsers.c src/sam.c)
