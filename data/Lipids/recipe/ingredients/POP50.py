#include as follow : execfile('pathto/POP50.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP50= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP50.sph',
radii = [[3.29, 1.98, 1.6799999999999999, 3.0899999999999999, 2.98, 3.1499999999999999, 3.4700000000000002, 2.9399999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP50.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP50',
positions = [[(3.4900000000000002, 1.3700000000000001, -16.109999999999999), (2.7200000000000002, -0.66000000000000003, -22.449999999999999), (2.8100000000000001, -4.0300000000000002, -21.82), (-1.24, -0.17999999999999999, -18.93), (-6.8099999999999996, -1.8600000000000001, -5.6100000000000003), (1.0, 2.1800000000000002, -2.6299999999999999), (-4.7000000000000002, -0.53000000000000003, -12.68), (2.9399999999999999, 3.79, -9.0399999999999991)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP50)