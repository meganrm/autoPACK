#include as follow : execfile('pathto/POP62.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP62= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP62.sph',
radii = [[2.8599999999999999, 2.5699999999999998, 1.5800000000000001, 3.73, 1.95, 3.9700000000000002, 3.5299999999999998, 1.9099999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP62.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP62',
positions = [[(0.10000000000000001, -2.6499999999999999, 21.539999999999999), (-2.77, -1.9299999999999999, 9.4800000000000004), (-7.3099999999999996, -2.5899999999999999, 1.8799999999999999), (4.8300000000000001, 2.79, 13.65), (-2.2000000000000002, -1.27, 15.15), (3.9900000000000002, 5.29, 5.3099999999999996), (0.70999999999999996, 0.85999999999999999, 18.760000000000002), (-5.3399999999999999, -3.8399999999999999, 5.4199999999999999)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP62)