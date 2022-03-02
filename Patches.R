

library(Morpho)
library(Rvcg)

#---------------#
# APPLY PATCHES #
#---------------#
# This is the code for the semi-automated process used to apply patch landmarks to
# all specimens, based on a template specimen.
# It requires a file with the patch landmarks only, the corresponding mesh for the
# specimen on which the patches are placed and a file containing meshes for all the
# specimens in the analysis

# 1. Import patch template points

# lower jaw
patch.coords <- as.matrix( read.table( paste( "patch/mandible_patchpoints.pts",
                                              sep = "") )[,-1] )

# 2. Extract initial specimen (template) to be patched from your list of specimens
template.specimen <- arr.raw[,,89]

# 3. Import mesh file for template specimen (should match template specimen)

template.mesh <- vcgImport( paste( "ply ASCII/Thunnus_tonggol_lower_jaw.ply", 
                                   sep = "" ) )

# list all anatomical (fixed) and curves landmarks
points <-(1:8)
curves <-c(9:13,14:43,44:68,69:98,99:128,129:153,154:168,169:173,174:183,184:223)

# 4. Plot template mesh
shade3d( template.mesh, col = "white" )

# 5. Add initial landmarks to template mesh (check landmark placement)
spheres3d( template.specimen[points,], radius = 0.5, col = "red" )
spheres3d( template.specimen[curves,], radius = .2, col = "gold" )

# 6. Specify semilandmark curves
curvein <- list( curves )

# 7. Create atlas to apply patch to all specimens
atlas <- createAtlas( mesh = template.mesh, landmarks = template.specimen,
                      patch = patch.coords, corrCurves = curvein,
                      patchCurves = NULL, keep.fix = c(1:8) )

# 8. Create patch image in 3D to check patch landmarks
# mesh
shade3d( atlas$mesh, col = bone1 )
# landmarks
spheres3d( atlas$landmarks, radius = .1, color = "red" )
# patch points
spheres3d( atlas$patch, radius = .1, color = "blue" )

# 9. Apply patch to all specimens, using mesh files for each
# link to folder containing meshes
patched.specimens <- placePatch( atlas = atlas, dat.array = arr.raw, inflate = 1,
                                 path = paste( "ply ASCII",
                                               sep = "" ),
                                 fileext = ".ply" )
