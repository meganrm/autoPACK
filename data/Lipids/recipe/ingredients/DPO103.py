#include as follow : execfile('pathto/DPO103.py',globals(),{'recipe':recipe_variable_name})
from AutoFill.Ingredient import SingleSphereIngr, MultiSphereIngr
DPO103= MultiSphereIngr( 
packingMode = 'random',
color = [1, 0, 0],
sphereFile = 'http://autofill.googlecode.com/svn/data//Lipids/spheres/DPO103.sph',
radii = [[3.73, 1.24, 3.8700000000000001, 1.48, 1.1200000000000001, 3.0800000000000001, 3.48, 4.7599999999999998]],
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
meshFile = 'http://autofill.googlecode.com/svn/data//Lipids/geoms/0/DPO103.c4d',
perturbAxisAmplitude = 0.1,
principalVector = [0.0, 0.0, -1.0],
name = 'DPO103',
positions = [[(-2.1699999999999999, -2.2400000000000002, -13.630000000000001), (-5.9699999999999998, 1.74, -21.75), (4.4900000000000002, 1.72, -22.98), (-3.5499999999999998, 2.1499999999999999, -23.440000000000001), (-6.3899999999999997, 1.97, -18.850000000000001), (-0.96999999999999997, 0.44, -20.149999999999999), (7.3399999999999999, 0.81000000000000005, -15.74), (-1.1299999999999999, -3.7599999999999998, -4.9299999999999997)]],
placeType = 'jitter',
useRotAxis = 1,
nbMol = 0,
)
recipe.addIngredient(DPO103)