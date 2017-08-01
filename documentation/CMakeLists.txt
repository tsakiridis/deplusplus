# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set root documentation folder and doxygen file
    set(DOXYGEN_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
    set(DOXYGEN_IN ${DOXYGEN_ROOT}/dox)

    add_custom_target( documentation
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_IN}
        WORKING_DIRECTORY ${DOXYGEN_ROOT}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)
