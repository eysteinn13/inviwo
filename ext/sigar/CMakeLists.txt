#--------------------------------------------------------------------
# Sigar Library
ivw_project(sigar)

#--------------------------------------------------------------------
# Define defintions
if(WIN32) 
    ## make sure we only use the smallest set of 
    ## headers on win32. Otherwise we get clashes 
    ## between winsock2.h and winsock.h
    add_definitions(-DWIN32_LEAN_AND_MEAN)
endif()

add_definitions(-D_CRT_SECURE_NO_WARNINGS)
ivw_add_definition(IVW_SIGAR)
    
source_group("CMake Files" FILES ${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt)

set(HEADER_FILES include/sigar.h 
    include/sigar_fileinfo.h 
    include/sigar_format.h 
    include/sigar_getline.h 
    include/sigar_log.h 
    include/sigar_private.h 
    include/sigar_ptql.h 
    include/sigar_util.h
)

include_directories(include/)

## linux
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    set(SIGAR_SRC src/os/linux/linux_sigar.c)
    include_directories(src/os/linux/)
endif(CMAKE_SYSTEM_NAME STREQUAL "Linux")

## macosx, freebsd
if(CMAKE_SYSTEM_NAME MATCHES "(Darwin|FreeBSD)")
    set(SIGAR_SRC src/os/darwin/darwin_sigar.c)
    include_directories(src/os/darwin/)
    if(CMAKE_SYSTEM_NAME MATCHES "(Darwin)")
        ivw_include_directories(/Developer/Headers/FlatCarbon/)
        add_definitions(-DDARWIN)
        set(SIGAR_LINK_FLAGS "-framework CoreServices -framework IOKit")
    else()
        ## freebsd needs libkvm
        set(SIGAR_LINK_FLAGS "-lkvm")
    endif()
endif()

## solaris
if(CMAKE_SYSTEM_NAME MATCHES "(Solaris|SunOS)" )
    set(SIGAR_SRC
        src/os/solaris/solaris_sigar.c
        src/os/solaris/get_mib2.c
        src/os/solaris/kstats.c
        src/os/solaris/procfs.c
    )

    include_directories(src/os/solaris/)
    add_definitions(-DSOLARIS)
    set(SIGAR_LINK_FLAGS -lkstat -ldl -lnsl -lsocket -lresolv)
endif()

## solaris
if(CMAKE_SYSTEM_NAME MATCHES "(hpux)" )
    set(SIGAR_SRC src/os/hpux/hpux_sigar.c)
    include_directories(src/os/hpux/)
    ivw_add_definition(-DSIGAR_HPUX)
    set(SIGAR_LINK_FLAGS -lnm)
endif()

## aix
if(CMAKE_SYSTEM_NAME MATCHES "(AIX)" )
    set(SIGAR_SRC src/os/aix/aix_sigar.c)

    include_directories(src/os/aix/)
    set(SIGAR_LINK_FLAGS -lodm -lcfg)
endif()

if(WIN32)
    if(BUILD_SHARED_LIBS)
        ivw_add_definition(SIGAR_SHARED)
        add_definitions(-DSIGAR_EXPORTS)
    endif()
    set(SIGAR_SRC src/os/win32/peb.c src/os/win32/win32_sigar.c src/os/win32/wmi.cpp)
    include_directories(src/os/win32)
else()
    add_definitions(-DHAVE_CONFIG_H)
endif()

set(SIGAR_SRC ${SIGAR_SRC}
    src/sigar.c
    src/sigar_cache.c
    src/sigar_fileinfo.c
    src/sigar_format.c
    src/sigar_getline.c
    src/sigar_ptql.c
    src/sigar_signal.c
    src/sigar_util.c
)

add_library(sigar ${SIGAR_SRC} ${HEADER_FILES})

if(WIN32)
    target_link_libraries(sigar PRIVATE ws2_32 netapi32 version)
endif(WIN32)

if(SIGAR_LINK_FLAGS)
    ivw_add_link_flags(${SIGAR_LINK_FLAGS})
endif()

#--------------------------------------------------------------------
# Creates VS folder structure
ivw_folder(sigar ext)

#--------------------------------------------------------------------
# Make package (for other projects to find)
ivw_default_install_comp_targets(core sigar)
ivw_make_package(Sigar sigar)

ivw_suppress_compiler_warnings(sigar)

