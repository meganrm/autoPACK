#include as follow : execfile('pathto/POP55.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP55= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP55.sph',
radii = [[4.0599999999999996, 3.5, 5.5999999999999996, 5.3799999999999999]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP55.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP55',
positions = [[(0.12, -2.1299999999999999, -17.210000000000001), (2.1299999999999999, 1.0800000000000001, -21.079999999999998), (-0.76000000000000001, 5.8099999999999996, -7.3700000000000001), (2.0299999999999998, -1.8999999999999999, -9.5700000000000003)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP55)
