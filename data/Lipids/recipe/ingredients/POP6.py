#include as follow : execfile('pathto/POP6.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
POP6= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/POP6.sph',
radii = [[7.2599999999999998, 4.4800000000000004, 5.5, 6.5300000000000002]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/POP6.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, 1.0],
name = 'POP6',
positions = [[(0.89000000000000001, -2.7000000000000002, -2.9900000000000002), (-0.54000000000000004, -0.29999999999999999, -21.449999999999999), (4.4500000000000002, -0.78000000000000003, -12.68), (-4.2599999999999998, -1.6000000000000001, -12.06)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(POP6)
