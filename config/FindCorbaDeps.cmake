#
# CORBA
#

OPTION( ENABLE_CORBA "Enable CORBA" OFF)
#DEPENDENT_OPTION( ENABLE_CORBA "Enable CORBA (using TAO)" ON "ACE_CONFIG AND TAO_ORB AND TAO_ORBSVCS AND CORBA_ENABLED" OFF)

# default to ACE/TAO if a corba implementation is not specified
IF (NOT CORBA_IMPLEMENTATION)
    SET( CORBA_IMPLEMENTATION "TAO" CACHE STRING "The implementation of CORBA to use (allowed values: TAO or OMNIORB )" )
ENDIF (NOT CORBA_IMPLEMENTATION)

if (ENABLE_CORBA)
    IF(${CORBA_IMPLEMENTATION} STREQUAL "TAO")
        # Look for TAO and ACE
        INCLUDE(${PROJ_SOURCE_DIR}/config/FindTAO.cmake)
        IF(NOT FOUND_TAO)
            MESSAGE(FATAL_ERROR "cannot find TAO")
        ELSE(NOT FOUND_TAO)
            MESSAGE(STATUS "CORBA enabled: TAO")
        ENDIF(NOT FOUND_TAO)
    ELSEIF(${CORBA_IMPLEMENTATION} STREQUAL "OMNIORB")
        INCLUDE(${PROJ_SOURCE_DIR}/config/FindOmniORB.cmake)
        IF(NOT OMNIORB4_FOUND)
            MESSAGE(FATAL_ERROR "cannot find OmniORB4")
        ELSE(NOT OMNIORB4_FOUND)
            MESSAGE(STATUS "CORBA enabled: OMNIORB")
        ENDIF(NOT OMNIORB4_FOUND)
    ELSE(${CORBA_IMPLEMENTATION} STREQUAL "TAO")
        MESSAGE(FATAL_ERROR "Unknown CORBA implementation '${CORBA_IMPLEMENTATION}': must be TAO or OMNIORB.")
    ENDIF(${CORBA_IMPLEMENTATION} STREQUAL "TAO")
endif (ENABLE_CORBA)
