cmake_minimum_required(VERSION 3.16.3...3.19.7 FATAL_ERROR)

project(OMEZarr)

#-----------------------------------------------------------------------------
# Extension meta-information
set(EXTENSION_HOMEPAGE "https://github.com/fideus-labs/SlicerOMEZarr")
set(EXTENSION_CATEGORY "Informatics")
set(EXTENSION_CONTRIBUTORS "Valentin Boussot (Fideus Labs), Matt McCormick (Fideus Labs)")
set(EXTENSION_DESCRIPTION "Open and save OME-Zarr (OME-NGFF) images, the chunked multi-resolution format used for large microscopy and medical volumes, from a local folder or a web address. Large images open at a resolution that fits in memory, and slice views can then be refined to full resolution.")
set(EXTENSION_ICONURL "https://raw.githubusercontent.com/fideus-labs/SlicerOMEZarr/main/SlicerOMEZarr.png")
set(EXTENSION_SCREENSHOTURLS "https://raw.githubusercontent.com/fideus-labs/SlicerOMEZarr/main/Screenshots/main.png")
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
