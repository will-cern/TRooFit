
if( XAOD_ANALYSIS )

atlas_subdir( TRooFit )

find_package( ROOT COMPONENTS MathCore RIO Core Tree Hist TreePlayer RooFit RooFitCore Matrix Gpad Graf RooStats )

atlas_add_root_dictionary( TRooFitLib TRooFitDictSource 
                           ROOT_HEADERS TRooFit/TRoo*.h TRooFit/apps/*.h TRooFit/Utils.h TRooFit/LinkDef.h )

atlas_add_library( TRooFitLib src/*.cxx src/apps/*.cxx ${TRooFitDictSource}
                    PUBLIC_HEADERS TRooFit
                    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                    LINK_LIBRARIES ${ROOT_LIBRARIES}
)

else()

#standalone project compilation

cmake_minimum_required( VERSION 3.2 FATAL_ERROR )
set(CMAKE_CXX_STANDARD 14)
project(TRooFit)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package( ROOT COMPONENTS MathCore RIO Core Tree Hist TreePlayer RooFit RooFitCore Matrix Gpad Graf RooStats )
include(${ROOT_USE_FILE})

include_directories(${ROOT_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR})
add_definitions(${ROOT_CXX_FLAGS})

file(GLOB HEADERS TRooFit/TRoo*.h TRooFit/apps/*.h TRooFit/Utils.h)

ROOT_GENERATE_DICTIONARY(G__TRooFit ${HEADERS} LINKDEF TRooFit/LinkDef.h)

file(GLOB SOURCES src/*.cxx src/apps/*.cxx)

add_library(TRooFit SHARED ${SOURCES} G__TRooFit.cxx)
target_link_libraries(TRooFit ${ROOT_LIBRARIES})

install( DIRECTORY TRooFit DESTINATION include/TRooFit FILES_MATCHING PATTERN "*.h" )
install( DIRECTORY ${CMAKE_BINARY_DIR}/ DESTINATION lib FILES_MATCHING PATTERN "libTRooFit*" )






endif()