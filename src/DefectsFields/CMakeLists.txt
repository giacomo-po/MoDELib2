include_directories(${MoDELibRootDir}/include/Utilities/)
include_directories(${MoDELibRootDir}/include/DislocationMobilities)
include_directories(${MoDELibRootDir}/include/DislocationDynamicsIO)
include_directories(${MoDELibRootDir}/include/DislocationMicrostructure/)
include_directories(${MoDELibRootDir}/include/DefectsFields/)

target_sources(${PROJECT_NAME} PRIVATE 
	${CMAKE_CURRENT_SOURCE_DIR}/DefectsFieldsExtractor.cpp
			)

