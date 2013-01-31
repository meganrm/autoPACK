#include as follow : execfile('pathto/POP56.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP56= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP56.sph',
radii = [[3.79, 5.2000000000000002, 5.2400000000000002, 4.54]],
cutoff_boundary = 0,
Type = 'MultiSphere',
cutoff_surface = 0,
gradient = '',
jitterMax = [0.5, 0.5, 0.10000000000000001],
packingPriority = 0,
rotAxis = [0.0, 2.0, 1.0],
nbJitter = 5,
molarity = 1.0,
rotRange = 6.2831,
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP56.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP56',
positions = [[(-0.080000000000000002, 1.0700000000000001, -8.4900000000000002), (2.3599999999999999, 1.1000000000000001, -21.77), (-0.78000000000000003, -2.4900000000000002, -16.390000000000001), (7.6699999999999999, -1.6299999999999999, -11.81)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP56)
