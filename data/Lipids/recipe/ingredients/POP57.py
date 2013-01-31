#include as follow : execfile('pathto/POP57.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP57= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP57.sph',
radii = [[6.0899999999999999, 3.2599999999999998, 6.0700000000000003, 6.5099999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP57.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'POP57',
positions = [[(4.6500000000000004, 2.02, 6.5099999999999998), (-3.8799999999999999, -0.42999999999999999, 21.190000000000001), (-6.8700000000000001, -2.2200000000000002, 7.1500000000000004), (0.059999999999999998, 0.089999999999999997, 17.34)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP57)
