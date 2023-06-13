PROJECT(VESSELFILTER)

SET( CONFIG_ITK_FILE /home/yliu/Imaging/Libraries/ITK/itk-3.20.1/install/lib/InsightToolkit/ITKConfig.cmake)
INCLUDE( ${CONFIG_ITK_FILE} )

INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIRS})
LINK_DIRECTORIES(${ITK_LIBRARY_DIRS})

SET(TARGET_LIBRARIES ITKIO ITKAlgorithms ITKBasicFilters ITKCommon)


ADD_EXECUTABLE(segment segment.cpp  )
ADD_EXECUTABLE(segmentVesselness segmentVesselness.cpp  )
ADD_EXECUTABLE(segmentFastMarching segmentFastMarching.cpp  )
ADD_EXECUTABLE(segmentFMWithInitialGraph segmentFMWithInitialGraph.cpp  )
ADD_EXECUTABLE(phantomGenerate phantomGenerate.cpp  )
ADD_EXECUTABLE(phantomPerturb phantomPerturb.cpp  )
ADD_EXECUTABLE(phantomCorrupt phantomCorrupt.cpp  )
ADD_EXECUTABLE(explant explant.cpp  )
ADD_EXECUTABLE(convertVesselToGML convertVesselToGML.cpp  )
ADD_EXECUTABLE(skeleton skeleton.cpp  )
ADD_EXECUTABLE(vesselTree vesselTree.cpp  )
ADD_EXECUTABLE(labelToGraph LabelToGraph.cpp )


TARGET_LINK_LIBRARIES (segment ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (segmentVesselness ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (segmentFastMarching ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (segmentFMWithInitialGraph ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (phantomGenerate ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (phantomPerturb ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (phantomCorrupt ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (explant ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (convertVesselToGML ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (skeleton ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (vesselTree ${TARGET_LIBRARIES})
TARGET_LINK_LIBRARIES (labelToGraph ${TARGET_LIBRARIES})
