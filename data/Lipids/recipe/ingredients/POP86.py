#include as follow : execfile('pathto/POP86.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP86= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP86.sph',
radii = [[4.0599999999999996, 7.4699999999999998, 5.0300000000000002, 6.0499999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP86.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP86',
positions = [[(-0.40999999999999998, -3.2599999999999998, 8.0199999999999996), (-6.9400000000000004, 4.79, 9.3900000000000006), (1.8600000000000001, -1.2, 21.420000000000002), (2.1699999999999999, 0.11, 15.32)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP86)
