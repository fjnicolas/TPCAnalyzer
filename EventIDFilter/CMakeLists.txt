art_make(
MODULE_LIBRARIES
  ${ART_FRAMEWORK_CORE}
  ${ART_FRAMEWORK_PRINCIPAL}
  ${ART_FRAMEWORK_SERVICES_REGISTRY}
  ${ART_FRAMEWORK_SERVICES_OPTIONAL}
  ${ART_ROOT_IO_TFILE_SUPPORT}
  ${ART_ROOT_IO_TFILESERVICE_SERVICE}
  ${MF_MESSAGELOGGER}
  ${FHICLCPP}
  ${ROOT_BASIC_LIB_LIST}
  nusimdata::SimulationBase

)
install_fhicl()
install_source()
FILE(GLOB fcl_files *.fcl)

add_subdirectory(job)
