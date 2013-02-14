#include as follow : execfile('pathto/POP90.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP90= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP90.sph',
radii = [[3.5, 2.6499999999999999, 2.5299999999999998, 4.0, 3.6899999999999999, 2.2599999999999998, 1.45, 1.6000000000000001]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP90.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP90',
positions = [[(0.54000000000000004, 4.6600000000000001, 14.08), (-1.5600000000000001, -5.5700000000000003, 9.75), (-4.9699999999999998, -4.6100000000000003, 5.0300000000000002), (-0.72999999999999998, 2.9700000000000002, 5.5300000000000002), (1.3200000000000001, 0.66000000000000003, 20.260000000000002), (0.55000000000000004, -3.71, 13.73), (-0.47999999999999998, -2.6400000000000001, 18.16), (3.1800000000000002, 3.25, 22.039999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP90)