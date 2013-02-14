#include as follow : execfile('pathto/POP20.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP20= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP20.sph',
radii = [[1.0, 4.8200000000000003, 4.0899999999999999, 3.8100000000000001, 2.0600000000000001, 0.76000000000000001, 3.52, 2.5499999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP20.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP20',
positions = [[(1.02, -0.80000000000000004, 25.370000000000001), (5.2199999999999998, 2.4900000000000002, 6.1500000000000004), (1.54, -0.35999999999999999, 12.65), (-2.1200000000000001, -1.0800000000000001, 19.0), (-1.75, 0.38, 23.690000000000001), (-0.60999999999999999, -1.4399999999999999, 25.870000000000001), (-4.8600000000000003, -0.82999999999999996, 14.43), (0.40000000000000002, 0.25, 6.5)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP20)