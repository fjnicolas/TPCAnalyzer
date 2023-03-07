art_make(  
 MODULE_LIBRARIES    
 ${ART_FRAMEWORK_CORE}    
 ${ART_FRAMEWORK_PRINCIPAL}    
 ${ART_ROOT_IO_TFILE_SUPPORT}      
${ART_FRAMEWORK_SERVICES_REGISTRY}    
 ${ART_ROOT_IO_TFILESERVICE_SERVICE}    
 ${FHICLCPP}    
 ${ROOT_BASIC_LIB_LIST}    
 canvas    
 art_Persistency_Common    
 art_Persistency_Provenance     
 art_Utilities    
 nusimdata_SimulationBase
 lardataobj_RecoBase
 larsim_Utils
 larsim_MCCheater_ParticleInventoryService_service
 larsim_MCCheater_BackTrackerService_service
 larcorealg_Geometry
 larcore_Geometry_Geometry_service
 larsim_Simulation lardataobj_Simulation
 larcore_Geometry_Geometry_service
 larcorealg_Geometry
 sbndcode_RecoUtils
 larreco_SpacePointSolver

 )
 install_source()
 install_fhicl()

add_subdirectory(job)
