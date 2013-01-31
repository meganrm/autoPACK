#include as follow : execfile('pathto/POP12.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP12= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP12.sph',
radii = [[7.8200000000000003, 4.6600000000000001, 4.4800000000000004, 4.7999999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP12.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP12',
positions = [[(1.48, -0.33000000000000002, -4.6900000000000004), (-1.6100000000000001, 0.94999999999999996, -20.68), (-3.0899999999999999, -2.7999999999999998, -11.31), (4.8499999999999996, 1.8600000000000001, -15.15)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP12)
