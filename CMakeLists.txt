cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(SlicerOMEZarr)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/fideus-labs/SlicerOMEZarr")
set(EXTENSION_CATEGORY "Informatics")
set(EXTENSION_CONTRIBUTORS "Valentin Boussot (Fideus Labs), Matt McCormick (Fideus Labs)")
set(EXTENSION_DESCRIPTION "Read OME-Zarr (OME-NGFF) images into 3D Slicer: multiscale level selection, spatial metadata and RFC-4 anatomical orientation mapped to RAS, region-of-interest loading from local or remote stores.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/fideus-labs/SlicerOMEZarr/main/SlicerOMEZarr.png")
set(EXTENSION_SCREENSHOTURLS "")
set(EXTENSION_DEPENDS "NA")

#-----------------------------------------------------------------------------
# Extension dependencies
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
# Extension modules
add_subdirectory(OMEZarr)
## NEXT_MODULE

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_GENERATE_CONFIG})
include(${Slicer_EXTENSION_CPACK})
