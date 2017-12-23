
atlas_subdir( TRooFit )

find_package( ROOT COMPONENTS MathCore RIO Core Tree Hist TreePlayer RooFit RooFitCore Matrix Gpad Graf RooStats )

atlas_add_root_dictionary( TRooFitLib TRooFitDictSource 
                           ROOT_HEADERS TRooFit/TRoo*.h TRooFit/apps/*.h TRooFit/Utils.h TRooFit/LinkDef.h )

atlas_add_library( TRooFitLib src/*.cxx src/apps/*.cxx ${TRooFitDictSource}
                    PUBLIC_HEADERS TRooFit
                    PRIVATE_INCLUDE_DIRS ${ROOT_INCLUDE_DIRS}
                    LINK_LIBRARIES ${ROOT_LIBRARIES}
)
