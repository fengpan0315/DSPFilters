cmake_minimum_required(VERSION 3.0)

project (DSPFilters CXX)

set (CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_BUILD_TYPE Release)

if((${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC))
  set(MYFLAGS "/O2 /WX- /MT")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MYFLAGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${MYFLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${MYFLAGS}")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${MYFLAGS}")
  include(CheckSymbolExists)
  check_symbol_exists(snprintf "stdio.h" HAVE_SNPRINTF)
  if(NOT HAVE_SNPRINTF)
  add_definitions(-Dsnprintf=_snprintf)
  endif()
endif()

aux_source_directory(${CMAKE_CURRENT_SOURCE_DIR}/source SOURCE_LIB)

add_library(${PROJECT_NAME} STATIC ${SOURCE_LIB})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set_property(TARGET ${PROJECT_NAME} PROPERTY POSITION_INDEPENDENT_CODE ON)
