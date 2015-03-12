import sys
import numpy as np
from boundedboxes import *
import binarysearchmr2 as bs

global psdPttoFace

"""includes all the Data Structures Present and used in boundedboxes.py"""

""" return the average of two numbers """
def average(x, y):
	return (x + y) / 2.0

""" returns a vector which points in the direction V (given in list format) with a new magnitude of M. 
	Is not possible/doesn't work if V is the zero vector. """
def scaleVector(v, m):
	antiscale = 1.0 / np.linalg.norm(v)
	scale = m * antiscale
	scaleVector = np.multiply(scale, v)
	return list(scaleVector)

""" Returns the bounded box that PT is inside as a tuple of variable length.
	ASSUMES runIteration() WAS CALLED BEFORE THIS METHOD IS INVOKED.
	If a point is along the boundary of two bounded boxes, getBox(pt, startBox) will return the smaller of the two,
	or both if the boxes are the same size. """
def getBox(pt, startBox):
	if (not startBox.isPtinBox(pt) or startBox == None):
		print "There is a problem"
		print "Pt must be in it's starting box"
		#throw an exception
		sys.exit()
		return None

	leftBox = startBox.child1
	rightBox = startBox.child2

	if (leftBox == None): #At the bottom level of a bounding box tree
		if (not (rightBox == None)):
			print "A box cannot have only one child"
			return
		return startBox
	elif (leftBox.isPtinBox(pt) and rightBox.isPtinBox(pt)): #Pt. is in both children
		goLeft = getBox(pt, leftBox)
		goRight = getBox(pt, rightBox)

		leftArea = goLeft.getBoxVol()
		rightArea = goRight.getBoxVol()

		if (leftArea > rightArea):
			return goLeft
		elif (leftArea == rightArea):
			return goLeft #, goRight #both of them
		else:
			return goRight
	if (leftBox.isPtinBox(pt)):
		return getBox(pt, leftBox)
	if (rightBox.isPtinBox(pt)):
		return getBox(pt, rightBox)

""" A Two Dimensional test which returns "inside" if inputPt is inside of the geometry, and else "outside".
	Helper function used in get2DGeometry method. If CURPSDPT belonds to more than one polygon, then
	no calcuations are needed to determine that INPUTPT is inside (since CURPSDPT is already determined
	to belong to the edge closest to INPUTPT. """
def insideOutsideTest2D(inputPt, curPsdPt, psdptDict):
	insideOrOutside = ""
	curPsdTuple = curPsdPt.ptToTuple()
	betweenPoints = curPsdPt.psdptBetween() #what points is my current psuedoPoint between ?
	curPolys = psdptDict[curPsdPt.ptToTuple()] #get the polygons associated with the current PsuedoPoint
	howManyPolys = len(curPolys)

	if (howManyPolys == 1):
		#check the Face associated with betweenPoints and curPolys.pop() for inside/outside
		myShape = curPolys.pop()
		curPolys.add(myShape)
		thisFace = None
		for face in myShape.boundaries:
			if (face.inPlanePoints() == betweenPoints):
				thisFace = face

		if (thisFace == None):
			print "no such face exists with this polygon...."
			return

		psdptVector = curPsdPt.getPointArray()
		normal = thisFace.getNormal()

		dot1 = normal
		myPoint = inputPt.getPointArray()

		dot2 = np.subtract(myPoint, psdptVector)

		dotProduct = np.dot(dot1, dot2)

		if (dotProduct > 0):
			insideOrOutside = "outside"
		else:
			insideOrOutside = "inside"
	elif (howManyPolys == 2):
		insideOrOutside = "inside"
	else:
		print "Impossible to an edge shared by " + str(howManyPolys) + " polygons."
	return insideOrOutside

""" Return -1 if PT is inside the geometry, and (+)1 if PT is outside. Either PT is closest to the face of a polygon,
		(MINPOLY), or an edge (shared by two Polygons) given in a tuple MINEDGE which contains two Polygon instances.
		Use MINPOLY for the test iff MINEDGE is None, else use MINEDGE. 
		The case in which the dot product is zero will not be an issue, unless the given geometry is not a closed space. """
def insideOutsideTest3D(pt, minPoly, minEdge):
		curPoint = pt.getPointArray()
		inPlanePt = minPoly.p1.getPointArray()
		if (minEdge == None):
			normal = minPoly.getNormal()
			planeToPt = np.subtract(curPoint, inPlanePt)
			if (np.dot(normal, planeToPt) >= 0):
				return 1
			else:
				return -1
		#so if you use minEdge
		#minEdge[0].printPolygon()
		#minEdge[1].printPolygon()
		edge1 = insideOutsideTest3D(pt, minEdge[0], None)
		edge2 = insideOutsideTest3D(pt, minEdge[1], None)
		# print edge1
		# print edge2
		if ( (edge1 == -1 and edge2 == 1) or (edge1 == 1 and edge2 == -1) ):
			print "There is a major problem, insideOutSide3D test to closest edges indicates both inside and outside. Check if polygon set up properly."
			return 0
		if (edge1 + edge2 > 0):
			return 1
		return -1

""" Return the norm of the vector perpendicular to the orthogonal projection (as a tuple) of Q onto the line between P1 and P2.
	(Hence, a postive number is always returned)
	Although, if the norm of the vector projection is greater than the length of the line segment, then just return the norm of
	the distance vector between newQ and (P2 - P1)
	(Also returning the perpendicular distance from Q to line from P1 to P2).
	P1, P2, and Q are lists of coordinates (and therefore are vectors all in Rn.
	MAKE SURE ALL INTS COERCED INTO FLOATS. """
def getDistofProjection(P1, P2, Q):
	basis = np.subtract(P2, P1)
	newQ = np.subtract(Q, P1)
	
	projQ = getnewQprojBasis(newQ, basis)

	if ((np.linalg.norm(projQ) > np.linalg.norm(basis)) or (np.linalg.norm(np.subtract(projQ, basis)) > np.linalg.norm(basis))):
		notOrthogonalDist1 = np.linalg.norm(np.subtract(Q, P1))
		notOrthogonalDist2 = np.linalg.norm(np.subtract(Q, P2))
		return min(notOrthogonalDist1, notOrthogonalDist2)

	z = np.subtract(newQ, projQ)

	return np.linalg.norm(z)

""" Given a PLANE and a PT, return the signed distance of the pt to that plane.
	Sign is positive if PT is on the same side as PLANE's normal vector.
	if PT is on the opposite side, sign is negative. """
def getSigned3DProj(plane, pt):
	a = plane.a
	b = plane.b
	c = plane.c
	d = plane.d
	d = d * -1

	x = pt.x
	y = pt.y
	z = pt.z

	top = a*x + b*y + c*z + d
	bottom = (a**2 + b**2 + c**2)**0.5

	if bottom==0: #one of the shapes I tested threw an divide by zero error here. 
		return 0
	else:
		return top / bottom

def getnewQprojBasis(newQ, basis):
	top = (np.dot(basis, newQ))

	bottom = np.dot(basis, basis)

	scale = float(top) / float(bottom)

	proj = np.multiply(scale, basis)

	return proj

""" return whether int/float x is between a and b numerically, x or y could be bigger number"""
def between(a, x, y):
	bigger = 0
	smaller = 0
	if (x >= y):
		bigger, smaller = x, y
	else:
		bigger, smaller = y, x
	if (a <= bigger and a >= smaller):
		return True
	return False

def ptBetween(a, pt1, pt2):
	if (between(a.x, pt1.x, pt2.x) and between(a.y, pt1.y, pt2.y) and between(a.z, pt1.z, pt2.z)):
		return True
	return False

""" Returns the Point DIST away iff pt plus/minus DIST (a scalar) in DIRECTION can leave BOX, otherwise returns 0. 
	DIR is a string, possibilities: '-x', '+x', '-y', '+y', '-z', '+z'. """
def canLeaveBox(pt, dist, box, direction):

	if (direction == '-x'):
		movedPt1 = (pt.x - dist, pt.y, pt.z)
		if (box.isPtTupinBox(movedPt1)):
			return 0
		return movedPt1
	elif (direction == '+x'):
		movedPt2 = (pt.x + dist, pt.y, pt.z)
		if (box.isPtTupinBox(movedPt2)):
			return 0
		return movedPt2
	elif (direction == '-y'):
		movedPt3 = (pt.x, pt.y - dist, pt.z)
		if (box.isPtTupinBox(movedPt3)):
			return 0
		return movedPt3
	elif (direction == '+y'):
		movedPt4 = (pt.x, pt.y + dist, pt.z)
		if (box.isPtTupinBox(movedPt4)):
			return 0
		return movedPt4
	elif (direction == '-z'):
		movedPt5 = (pt.x, pt.y, pt.z - dist)
		if (box.isPtTupinBox(movedPt5)):
			return 0
		return movedPt5
	elif (direction == '+z'):
		movedPt6 = (pt.x, pt.y, pt.z + dist)
		if (box.isPtTupinBox(movedPt6)):
			return 0
		return movedPt6
	else:
		return "not a valid direction"

""" POINT DATA STRUCTURE, has a label#, and coordinates."""
class Point:
	""" constructor takes in a list of strings (unparsed) of form
		[label, x coord, y coord, z coord] """
	def __init__(self, index=0, x=0,y=0,z=0):
		self.number = float(index)
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.sisterPoints = set() #DOES NOT include self in list of sister points
	
	def parse(self, unparsed):
		self.number=unparsed[0]
		self.x = unparsed[1]
		self.y = unparsed[2]
		self.z = unparsed[3]

	""" prints the point in the form: Point #, (x, y, z) """
	def printPoint(self):
		print "Point " + str(self.number) + ": (" + str(self.x) + ", " + str(self.y) + ", " + str(self.z) + ")"
	
	""" turns the coordinates into array, returning as [x, y, z]"""
	def getPointArray(self):
		return [self.x, self.y, self.z]

	""" changes PsdPt number to n"""
	def changePtNumber(self, n):
		self.number = n

	""" Returns the point number of this point. """
	def getPtNumber(self):
		return self.number

	""" returns the coordinates as an immutable tuple in form (x, y, z). """
	def ptToTuple(self):
		return (self.x, self.y, self.z)

	""" will fill the set of sisterPoints using a list of Polygons. Finds the polygons that self is part of,
		and stores the points belonging to the same Polygon(s) as self does.
		NEW: INCLUDES: PSUEDOPOINTS as sister pts"""
	def getSisterPoints(self, polyList):

		""" Returns True iff this Point is one of the points defining the Polygon instance POLY. """
		def ptinPoly(poly):
			return poly.isPointinPoly(self)

		filteredPolys = filter(ptinPoly, polyList)
		for poly in filteredPolys:
			for pt in poly.getPtsList():
				if (pt.number != self.number):
					self.sisterPoints.add(pt)

""" PSEUDO POINT DATA STRUCTURE, exactly like points but exclusively along the edge and are NOT
	included in the geometry of the inputted points and polygons. All haracterized by self.number = -1.
	Attribute btween stores a sorted list of length 2 of the point numbers which THIS (self) PseudoPoint
	resides between. Used for accessing sister PsuedoPoints. """
class PsuedoPoint(Point):
	""" constructor takes in a list with coords [x, y, z], and two ints of the points it is between, and then calls super constructor
		with -1 added at the begining so self.number = 1. """
	def __init__(self, xyzList, pNum1, pNum2):
		Point.__init__(self, index="-1", x=xyzList[0],y=xyzList[1],z=xyzList[2])
		self.btween = [pNum1, pNum2]
		self.btween.sort()
		self.sisterpoints = set()

	""" Alternate printPoint() method. """
	def printPsdPt(self):
		Point.printPoint(self)
		print "between Point " + str(self.btween[0]) + " and Point " + str(self.btween[1])

	""" Returns True if values in self.btween are same as in other.btween. Indices should algin bc both lists are sorted. """
	def sameBetween(self, other):
		return (self.btween[0] == other.btween[0] and self.btween[1] == other.btween[1])

	""" Returns the point numbers in a list which this (self) PsuedoPoint is place between. """
	def psdptBetween(self):
		return self.btween

	""" Returns a list of sister psuedo points of THIS (self) PsuedoPoint found in PSDPTLIST. sister
		PsuedoPoints are defined as other psuedo points which are in the same polygon and
		have the same values in sispt.btween as THIS PsuedoPoint. PSDPTLIST is a list of
		PsuedoPoints. It is a very good idea to filter out the list before
		invoking this method, such as only constructing the list as to have PsuedoPoints
		contained in teh same polygon (although in this method there is no way to check this condition
		directly bc lack of information). """
	def getSisterPsdPts(self, psdptList):
		fillThis = set()
		for other in psdptList:
			if (self.sameBetween(other) and not (self.ptToTuple() == other.ptToTuple())):
				fillThis.add(other)
		return fillThis

class Plane:
	""" Constructs a plane, inputs are floats/ints and follow equation ax + by + cz = d, and where PT is
		a Point instance in this plane. """
	def __init__(self, a, b, c, d, pt):
		self.a = float(a)
		self.b = float(b)
		self.c = float(c)
		self.d = float(d)
		self.pt = pt
	
	""" prints the plane equation"""
	def printEqn(self):
		print str(self.a) + "x + " + str(self.b) + "y +" + str(self.c) + "z = " + str(self.d)

	""" Returns the normal vector for this (self) plane. """
	def getNormal(self):
		return [self.a, self.b, self.c]

	""" Returns the normal vector, scaled to size N of this (self) plane. """
	def getNormalSizeN(self, n):
		original = self.getNormal()
		magnitude = np.linalg.norm(original)
		antiScale = 1.0 / magnitude
		scale = n * antiScale
		return np.multiply(scale, [self.a, self.b, self.c])

	""" returns True iff PT, a Point Instance, is in THIS(self) plane.
		ASSUMES this(self) plane does NOT carry boundaries as the edges of a polygon indicate. """
	def isPtinPlane(self, pt):
		ax = self.a * pt.x
		by = self.b * pt.y
		cz = self.c * pt.z

		return abs(self.d - (ax + by + cz)) < 0.00001

	""" Returns a plane parallel to THIS (self) Plane. Distance between the two planes is given by N (scalar). """
	def createParallelPlane(self, n):
		n = float(n)
		extrudeFrom = self.p1.getPointArray()
		normal = self.getNormal()
		A = normal[0]
		B = normal[1]
		C = normal[2]
		addedVector = np.add(normal, extrudeFrom)
		addedNorm = np.linalg.norm(addVector)
		toNormalized = n / addedNorm

		scaledInPlane = np.multiply(toNormalized, addedVector)
		D = (A * scaledInPlane[0]) + (B * scaledInPlane[1]) + (C * scaledInPlane[2])
		ptInPlane = Point()
		ptInPlane.parse([-2] + scaledInPlane)
		return Plane(A, B, C, D, ptInPlane)

class Polygon:
	"""constructor, passes in a list like [polynum, p1num, p2num, p3num] which are all strings (hence unparsed).
		ptNumbersList is a list of points from which it picks the indexes p1, p2, p3 to be the points
		in the polygon (are in CCW order).
		self.psuedoPts is a list of PsuedoPoint instances along the edges of this polygon, initially empty (to be filled
			in method invocation to makeBarrierPseudoPts (in boundboxes.py) which is called by iterateBoxes which is called by
			runIteration which is run by the main hopefully).
		self.boundaries is a set of Face instances representing extruding walls from the edges of this polygon (to be filled
			by the addBoundary method of this class, which is called in makeBoundaries (in extrusion2D.py)).
		self.plane is the Plane instance in which this polygon lies. Initially set to None, but set to the Plane instance
			returned by the getPlaneEqn instance method. """
	def __init__(self, unparsed, ptNumbersList):
		self.number = int(unparsed[0])
		self.p1 = ptNumbersList[int(unparsed[1])]
		self.p2 = ptNumbersList[int(unparsed[2])]
		self.p3 = ptNumbersList[int(unparsed[3])]
		self.psuedoPts = []
		self.boundaries = set()
		self.plane = None
	"""prints the polygon so it looks nice."""
	def printPolygon(self):
		print "Polygon " + str(self.number) + ": (" + str(self.p1.number) + ", " + str(self.p2.number) + ", " + str(self.p3.number) + ")"

	"""Adds a boundary to the set of boundaries associated/extruding from this polygon. """
	def addBoundary(self, face):
		self.boundaries.add(face)
	
	""" returns the pts of each of its defining points (p1, p2, p3) in an array. """
	def getPtsList(self):
		return [self.p1, self.p2, self.p3]

	""" returns the pt numbers of each of its defining points (p1, p2, p3) in an array. """
	def getPtsNumList(self):
		return [self.p1.number, self.p2.number, self.p3.number]

	""" Return a list of the polygons edges. (Edges represent by a tuple containing two Point instances) """
	def getEdges(self):
		return [(self.p1, self.p2), (self.p2, self.p3), (self.p3, self.p1)]

	""" Returns True iff this Polygon shares an edge with POLY (another Polygon instance). """
	def sharesEdge(self, poly):
		pts1 = self.getPtsNumList()
		pts2 = poly.getPtsNumList()
		pts1.sort()
		pts2.sort()

		i = 0
		countSame = 0
		while (i < len(pts1)):
			if (pts1[i] == pts2[i]):
				countSame += 1
			i += 1
		return countSame > 1

	""" Return the list of PsuedoPoint instances lying on the edges of this Polygon. """
	def getPsdPts(self):
		return self.psuedoPts

	""" Returns the last Point in the polygon given the point numbers of the two other points. """
	def getOtherPt(self, p1num, p2num):
		if ((p1num == self.p1.number and p2num == self.p2.number) or (p1num == self.p2.number and p2num == self.p1.number)):
			return self.p3
		elif ((p1num == self.p1.number and p2num == self.p3.number) or (p1num == self.p3.number and p2num == self.p1.number)):
			return self.p2
		elif ((p1num == self.p3.number and p2num == self.p2.number) or (p1num == self.p2.number and p2num == self.p3.number)):
			return self.p1
		else:
			return None

	""" returns true if the given pt is one of the poitns that defines SELF polygon. """
	def isPointinPoly(self, pt):
		return (pt.number == self.p1.number or pt.number == self.p2.number or pt.number == self.p3.number)

	""" Returns True iff PT (a Point instance) is both within the same plan as this SELF polygon, and if PT is inside the boundaries
		defining this SELF polygon. """
	def isInside(self, pt):
		if (self.plane == None):
			self.plane = self.getPlaneEqn()
		if (not self.plane.isPtinPlane(pt)): #if the point is not in the same plane as THIS/self polygon, then return False
			return False

		""" Return True iff vectors (in list format) U and V point in the same direction.
		Assumes that U and V are cross products from related vectors, (so they
		are either parallel or perpendicular). (In reality, returns whether the
		angle between the two vectors is between zero and than pi/2. """
		def pointsInSameDir(u, v):
			return np.dot(u, v) >= 0

		a = self.p1.getPointArray()
		b = self.p2.getPointArray()
		c = self.p3.getPointArray()
		curPt = pt.getPointArray()

		#preserving CCW order
		aToB = np.subtract (b, a)
		bToC = np.subtract(c, b)
		cToA = np.subtract(a, c)

		aToPt = np.subtract(curPt, a)
		bToPt = np.subtract(curPt, b)
		cToPt = np.subtract(curPt, c)

		cross1 = np.cross(aToB, aToPt)
		cross2 = np.cross(bToC, bToPt)
		cross3 = np.cross(cToA, cToPt)

		referenceCross = np.cross(aToB, np.subtract(c, a))

		return pointsInSameDir(referenceCross, cross1) and pointsInSameDir(referenceCross, cross2) and pointsInSameDir(referenceCross, cross3)

	"""returns a Plane, for plane eqn ax + by + cz = d"""
	def getPlaneEqn(self):
		point1 = self.p1.getPointArray()
		point2 = self.p2.getPointArray()
		point3 = self.p3.getPointArray()
		
		vector1 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
		vector2 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
		normal = np.cross(vector1, vector2)
		normal = self.getNormal()
		A = normal[0]
		B = normal[1]
		C = normal[2]
		D = (A * point1[0]) + (B * point1[1]) + (C * point1[2])
		
		planeEqn = Plane(A, B, C, D, self.p1)
		self.plane = planeEqn
		return planeEqn

	""" returns the normal vector of the plane in which this Polygon lies. """
	def getNormal(self):
		point1 = self.p1.getPointArray()
		point2 = self.p2.getPointArray()
		point3 = self.p3.getPointArray()
		
		vector1 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
		vector2 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
		normal = np.cross(vector1, vector2)
		return list(normal)

""" A face (used for extrusion barriers) is a special Polygon made with 4 Points. If need be,
	it can be split along a diagonal into 2 Polygons each composed of 3 Points. """
class Face(Polygon):
	""" constructs a face using points p1, p2, p3, and p4 which are all Point Instances.
		Directionality follows: p1 ---> p2 ----> p3 ----> p4.... CCW. 
		POLYGON is the polygon from which this Face is extruding. Because of this, it is possible
		to have a boundary covered by 2 Face instances (bc face may be shared between two polygons),
		but SELF.POLYGON (as well as directionality) will distinguish this two faces from one another. """
	def __init__(self, p1, p2, p3, p4, polygon):
		self.p1 = p1
		self.p2 = p2
		self.p3 = p3
		self.p4 = p4
		self.polygon = polygon

	""" Prints this face out. """
	def printFace(self):
		print "Face with points " + str(self.p1.number) + " and " + str(self.p2.number) + " with polygon " + str(self.polygon.number)

	""" Returns the two points of this face associated with the plane in which all original points lie in a sorted list. """
	def inPlanePoints(self):
		make = [self.p1.number, self.p2.number]
		make.sort()
		return make

	""" Returns the normal of the face, calculated using the same CCW directionality as indicated
		in the constructor method. cross p2 - p1 with p3 - p1 to get normal.
		NORMAL always Faces Outwards. """
	def getNormal(self):
		point1 = self.p1.getPointArray()
		point2 = self.p2.getPointArray()
		point3 = self.p3.getPointArray()
		point4 = self.p4.getPointArray()
		
		vector1 = [point2[0] - point1[0], point2[1] - point1[1], point2[2] - point1[2]]
		vector2 = [point3[0] - point1[0], point3[1] - point1[1], point3[2] - point1[2]]
		normal = np.cross(vector1, vector2)

		return list(normal)

class BoundedBox:
	"""constructs a bounded box given dimensions x y z about Point instance centerpt"""
	def __init__(self, x, y, z, centerpt, parent = None):
		self.x = x
		self.y = y
		self.z = z
		self.centerpt = centerpt
		self.pointsInside = []
		self.parent = parent
		self.longEdge = self.getLongestEdge()
		self.psuedoPts = set()
		self.child1 = None
		self.child2 = None
		self.facesInside = set()
# i think boxtype should be an attribute of the box. 
		self.boxType=-1 

	def getTopRight(self):
		return (self.getMaxXcoord(), self.getMaxYcoord(), self.getMaxZcoord())

	def getBottomRight(self):
		return

	def getMaxXcoord(self):
		return self.centerpt.x + self.x / 2.0

	def getMinXccord(self):
		return self.centerpt.x - self.x / 2.0

	def getMaxYcoord(self):
		return self.centerpt.y + self.y / 2.0

	def getMinYccord(self):
		return self.centerpt.y - self.y / 2.0

	def getMaxZcoord(self):
		return self.centerpt.z + self.z / 2.0

	def getMinZccord(self):
		return self.centerpt.z - self.z / 2.0

	""" returns the volume of THIS (self) bounding box (according to the predefined dimensions).
		IF THE BOX IS TWO DIMENSIONAL.... will ignore the z-coordinate. """
	def getBoxVol(self):
		if (self.z == 0):
			return self.x * self.y
		return self.x * self.y * self.z

	def isPtinBox(self, pt):
		if (pt.x <= self.getMaxXcoord() and pt.x >= self.getMinXccord()):
			if (pt.y <= self.getMaxYcoord() and pt.y >=  self.getMinYccord()):
				if (pt.z <= self.getMaxZcoord() and pt.z >= self.getMinZccord()):
					return True
		return False

	def isPtTupinBox(self, ptTup):
		if (ptTup[0] <= self.getMaxXcoord() and ptTup[0] >= self.getMinXccord()):
			if (ptTup[1]<= self.getMaxYcoord() and ptTup[1] >=  self.getMinYccord()):
				if (ptTup[2] <= self.getMaxZcoord() and ptTup[2] >= self.getMinZccord()):
					return True
		return False

	""" will return the x coor of the centerprt if dir = "x", y coor if dir = "y", and z coor if dir = "z" """
	def getDirCenterPt(self, dir):
		if (dir == 'x'):
			return self.centerpt.x
		elif (dir == 'y'):
			return self.centerpt.y
		else:
			return self.centerpt.z

	def printBox(self):
		print "this box has dimensions x = " + str(self.x) + ", y = " + str(self.y) + ", z = " + str(self.z) + " with centerpt: (" + str(self.centerpt.x) + ", " + str(self.centerpt.y) + ", " + str(self.centerpt.z) + ")"

	""" will return the type of the BoundedBox as an integer, MUST getPointsinBox before asking for the type"""
	def getBoxType(self):
		#i forgot the types
		numPts = len(self.pointsInside)
		if (numPts == 0):
			#do not need to check for number of edges
			return 0
		elif (numPts == 1):
			return 1
		else:
			return -1

	def setBoxType(self):
		#i forgot the types
		numPts = len(self.pointsInside)
		if (numPts == 0):
			#do not need to check for number of edges
			self.boxType= 0
		elif (numPts == 1):
			self.boxType= 1


	""" takes in a list of Point instances, resets the attribute pointsInside and then
		will refill the list pointsInside with only the Point instances within the box.
		Does mutation ONLY on self.pointsInside and does not return any new list or value """
	def getPointsinBox(self, ptsList):
		self.pointsInside = []
		
		for pt in ptsList:
			if ( (between(pt.x, self.getMaxXcoord(), self.getMinXccord())) and 
				(between(pt.y, self.getMaxYcoord(), self.getMinYccord())) and 
				(between(pt.z, self.getMaxZcoord(), self.getMinZccord())) ):
				self.pointsInside.append(pt)
		return

	def sortPoints(self, axis):
		self.pointsInside.sort(key=lambda point: point.__dict__[axis], reverse=False)

	""" Prints the points inside of the box. RUN getPointsinBox BEFORE calling this """
	def printPoints(self):
		for pt in self.pointsInside:
			print(pt.number)

	""" REQUIRED THAT runIteration(ptsLst, plyList) IS CALLED PREVIOUSLY TO CALLING THIS FUNCTION.
		the idea is that this method will fill THIS bounded box (self) with the faces inside of it, coming from
		it's two child boxes (who should already have their facesInside filled in the making of psuedo points). 
		In this way, fills facesInside with the union of self.child1.facesInside and self.child2.facesInside"""
	def getFacesinBox(self, psdptDict):
		if (self.getBoxType != -1):
			#this condition mandates that in spiltting the boxes (at the final level, should plainly look at box.psdPts and use the
			#psdptDict dictionary to look up the faces inside the box
			for psdpt in self.psuedoPts:
				curPolySet = psdptDict[psdpt.ptToTuple()]
				self.facesInside = self.facesInside | curPolySet
			return
		else:
			self.facesInside = getFacesinBox(self.child1.facesInside) | getFacesinBox(self.child2.facesInside)

	"""" 2D TEST
		Given INPUTPT (Assumed to be a Point instance inside of THIS (self) box), but it does not have to be inside this(self) box.
		return the nearest distance from PT to a face inside
		of THIS (self) box. Also takes in a psdptDict to call getFacesinBox, if THIS (self) box does have
		children boxes. REQUIRED THAT YOU CALL runIteration(ptsLst, plyList) PREVIOUS TO INVOKING THIS METHOD. 
		STARTBOX needed so if you need to acquire neighbor boxes. CHECKOUTSIDE indicates whether you should check
		if the distance calculated within the box needes to be checked if it leave this (self) box.
		------------> to implement: return a SIGNED distance, negative if inside, positive if outside. """
	def get2DGeometry(self, inputPt, psdptDict, startBox, checkOutside):
		myfaces = set()
		if (self.child1 == None): #no children
			myfaces = self.facesInside
		else:
			myfaces = self.getFacesinBox(psdptDict)

		sisters = set()
		alreadyChecked = set()
		distances = set()
		runningMinPt = None
		directions = ['-x', '+x', '-y', '+y']# '-z', '+z']
		
		for pt in self.psuedoPts: #optimized so you don't double check edges.
			#Get the psuedo sister points of PT which are inside THIS (self) box
			if (not pt in alreadyChecked):
				alreadyChecked.add(pt)
				sisters = pt.getSisterPsdPts(self.psuedoPts)
				numSisters = len(sisters)
				if (numSisters == 1):
					#calculate edge with this sister point
					mySister = sisters.pop()
					alreadyChecked.add(mySister)

					dist = getDistofProjection(pt.getPointArray(), mySister.getPointArray(), inputPt.getPointArray())

					distances.add(dist)

					if (dist == min(distances)):
						runningMinPt = pt

				elif (numSisters == 0):
					#calculate edge with the real point in this box
					soloRealPoint = self.pointsInside[0].getPointArray()

					dist = getDistofProjection(pt.getPointArray(), soloRealPoint, inputPt.getPointArray())

					distances.add(dist)

					if (dist == min(distances)):
						runningMinPt = pt 

				else:
					print "there is a major problem, expected 1 or 0 sisterpoints and got " + str(numSisters)
					for i in sisters:
						i.printPoint()
					return
		#REORGANIZE. if min(dist = 0, you don't need to checkOutside or do an insideoutside test)
		if (len(distances) == 0): #no geometry within box
			if (checkOutside):
				return self.parent.get2DGeometry(inputPt, psdptDict, startBox, True)
			return sys.maxint
		#Distance intially positive. then, check if you can go outside of the box.
		curMin = abs(min(distances))

		if (curMin == 0):
			return curMin

		#Do insideOutsideTest, curMin > 0
		if (checkOutside):
			if (runningMinPt == None):
				print "there is no psuedoPoint with closest geometry"
				return

			""" TEMPORARY!!! CHANGE BACK TO # COMMENT IF THIS DOES NOT WORK."""
			soFar = False
			insideOrOutside = "outside"
			for poly in self.facesInside:
				soFar = poly.isInside(inputPt) or soFar
				if (soFar):
					insideOrOutside = "inside"
					break
			#insideOrOutside = insideOutsideTest2D(inputPt, runningMinPt, psdptDict)
			if (insideOrOutside == "inside"):
				sign = -1
			elif (insideOrOutside == "outside"):
				sign = 1
			else:
				print "problem with insideOutsideTest2D"
				return sys.maxint
		else:
			sign = 1

		outsideDist = set()
		
		if (checkOutside): #if you need to check outside, ONLY called on initial box containing inputPt call basically. FIX FIX FIX FIX
			for direct in directions:
				left = canLeaveBox(inputPt, curMin, self, direct)
				#if left is not zero, then get the box left is in, then test inputPt for geometry in the newly acquired Box.
				if (left != 0):
					if (startBox.isPtTupinBox(left)):
						leftPoint = Point(index= -1, x=left[0], y=left[1], z=left[2])
						neighborBox = getBox(leftPoint, startBox)
						outsideDist.add(abs(neighborBox.get2DGeometry(inputPt, psdptDict, startBox, False)))
						#add absolute value will ensure all distances are initially positive.
			if (len(outsideDist) != 0):
				return sign * min(curMin, min(outsideDist)) #both curMin and all entries in outsideDist are positive by abs() function
		return sign * curMin

	"""" 3D TEST.
		Given PT (Assumed to be inside of THIS (self) box), but it does not have to be inside this(self) box.
		return the nearest distance from PT to a face inside
		of THIS (self) box. Also takes in a psdptDict to call getFacesinBox, if THIS (self) box does have
		children boxes. REQUIRED THAT YOU CALL runIteration(ptsLst, plyList) PREVIOUS TO INVOKING THIS METHOD. 
		STARTBOX needed so if you need to acquire neighbor boxes. CHECKOUTSIDE indicates whether you should check
		if the distance calculated within the box needes to be checked if it leave this (self) box.
		------------> to implement: return a SIGNED distance, negative if inside, positive if outside. """
	def get3DGeometry(self, inputPt, psdptDict, startBox, checkOutside):
		myfaces = set()
		if (self.child1 != None): #no children
			self.getFacesinBox(psdptDict)

		distances = {}
		runningMin = sys.maxint #MUST be positive
		runningMinPoly = None
		runningMinPolyTuple = None
		curMinSign = 1
		edgeisMin = False
		directions = ['-x', '+x', '-y', '+y', '-z', '+z']
		
		for poly in self.facesInside: #check distance to each face

			useEdge = False

			if (runningMinPoly == None):
				runningMinPoly = poly

			curPlane = poly.getPlaneEqn()

			dist = getSigned3DProj(curPlane, inputPt)

			backToPlane = curPlane.getNormalSizeN(dist)

			backInPlane1 = np.subtract(inputPt.getPointArray(), backToPlane)
			backInPlane2 = np.add(inputPt.getPointArray(), backToPlane)
			testBack1 = list(backInPlane1)
			testBack2 = list(backInPlane2)
			testBackPt1 = Point()
			testBackPt1.parse([0.5]+testBack1)
			testBackPt2 = Point()
			testBackPt2.parse([0.5] + testBack2)
			if ( (not poly.isInside(testBackPt1)) and (not poly.isInside(testBackPt2)) ):
				#poly.printPolygon()
				#print "help caroline"
				useEdge = True
				#look at psuedoPoints in polygon, form
				edgeDistances = set()
				for edge in poly.getEdges():
					
					p1 = edge[0].getPointArray()
					p2 = edge[1].getPointArray()
					#check the edges of the current polygon.!!!
					edgeProj = getDistofProjection(p1, p2, inputPt.getPointArray())
					#print edgeProj
					#print "( " + str(edge[0].number) + ", " + str(edge[1].number) + ")"
					edgeDistances.add(edgeProj)

				dist = min(edgeDistances)

			if (dist == 0):
				return dist, runningMinPoly

			if (useEdge and dist == runningMin and runningMinPoly.sharesEdge(poly)):
				#the case in which inputPt is closest to an edge shared by two polygons
				runningMinPolyTuple = (runningMinPoly, poly)
				edgeisMin = True

			if (abs(dist) < runningMin):
				runningMin = abs(dist)
				runningMinPoly = poly
				runningMinPolyTuple = None
				if (useEdge):
					edgeisMin = True
				else:
					edgeisMin = False
				#curMinSign = dist / abs(dist)
			#only need to test a close face for inside Outside

		#if min(dist = 0, you don't need to checkOutside or do an insideoutside test)
		if (len(self.facesInside) == 0): #no geometry within box
			if (checkOutside and self.parent != None):
				return self.parent.get3DGeometry(inputPt, psdptDict, startBox, True)
			return sys.maxint, None
		#Distance intially positive. then, check if you can go outside of the box.

		if (runningMin == 0):
			return runningMin, runningMinPoly

		#find the sign
		if (runningMinPolyTuple == None and edgeisMin): #case where closest geometry is not closed
			if (runningMinPoly.isInside(inputPt)):
				curMinSign = -1
			else:
				curMinSign = 1
		else:
			curMinSign = insideOutsideTest3D(inputPt, runningMinPoly, runningMinPolyTuple)

		runningMinOutside = (sys.maxint, None)
		
		if (checkOutside): #if you need to check outside, ONLY called on initial box containing inputPt call basically. FIX FIX FIX FIX
			for direct in directions:
				left = canLeaveBox(inputPt, runningMin, self, direct)
				#if left is not zero, then get the box left is in, then test inputPt for geometry in the newly acquired Box.
				if (left != 0):
					if (startBox.isPtTupinBox(left)):
						leftPoint = Point(index=-1, x=left[0], y=left[1], z=left[2])
						neighborBox = getBox(leftPoint, startBox)
						curOutside = neighborBox.get3DGeometry(inputPt, psdptDict, startBox, False)
						if (runningMinOutside[0] > abs(curOutside[0])):
							runningMinOutside = (abs(curOutside[0]), curOutside[1])
						#add absolute value will ensure all distances are initially positive.
			if (runningMinOutside[1] != None and runningMinOutside[0] < runningMin):
				return curMinSign * runningMinOutside[0], runningMinOutside[1]
		return curMinSign * runningMin, runningMinPoly

	""" will NOT modify the original BoundedBox instance. Will return two new bounded boxes
		in a tuple, each of which is half by of the original BoundedBox instance cut along its longest edge, """
	def spliceBoxLong(self):
		longest = self.longEdge #you already call the function in the initation of the box
		if (longest == 'x'):
			return self.spliceBoxAlongX()
		elif (longest == 'y'):
			return self.spliceBoxAlongY()
		else:
			return self.spliceBoxAlongZ()

	def splitPoints(self, axis):
		self.sortPoints(axis)
		max_index= len(self.pointsInside)-1
		leftPoints=[]
		rightPoints=[]
		if len(self.pointsInside) <=1:
			raise Exception( "shouldn't be spliting this box")

		#if all the points are below the cutpoint
		elif self.pointsInside[max_index].__dict__[axis] < self.centerpt.__dict__[axis]:
			leftPoints=self.pointsInside
		#if all the points are above the cutpoint
		elif self.pointsInside[0].__dict__[axis] > self.centerpt.__dict__[axis]:
			rightPoints=self.pointsInside
		
		#if there are exactly 2 points, each on either side of the cutpoint
		elif len(self.pointsInside) ==2:
			leftPoints.append(self.pointsInside[0])
			rightPoints.append(self.pointsInside[1])

		else:
		#binary search for the index of the last point below the cutpoint
			cut_index=bs.binary_search(self.pointsInside, self.centerpt.__dict__[axis], axis)
		#if the point lies on the cutting line, add point to both sides
			if self.pointsInside[cut_index].__dict__[axis] == self.centerpt.__dict__[axis]:
				leftPoints=self.pointsInside[:cut_index+1]
				rightPoints=self.pointsInside[cut_index:]

			else:
				leftPoints=self.pointsInside[:cut_index+1]
				rightPoints=self.pointsInside[cut_index+1:]
				# if self.pointsInside[cut_index] in same poly as self.pointsInside[cut_index+1]:
				# 	make new PseudoPoint()
		
		return leftPoints, rightPoints


	def spliceBoxAlongX(self):

		centerLeftShift = self.centerpt.x - 0.25 * self.x
		centerRightShift = self.centerpt.x + 0.25 * self.x

		leftcenter = Point(index=-1, x=centerLeftShift, y=self.centerpt.y, z=self.centerpt.z)
		rightcenter = Point(index=-1, x=centerRightShift, y=self.centerpt.y, z=self.centerpt.z)

		leftBox = BoundedBox(self.x / 2.0, self.y, self.z, leftcenter, self)
		rightBox = BoundedBox(self.x / 2.0, self.y, self.z, rightcenter, self)

		self.child1 = leftBox
		self.child2 = rightBox

		# leftBox.getPointsinBox(self.pointsInside)
		# rightBox.getPointsinBox(self.pointsInside)

		leftBox.pointsInside, rightBox.pointsInside= self.splitPoints('x')

		
		#look at the x coordinates of the psuedo points
		for psdpt in self.psuedoPts:
			if (psdpt.x > self.centerpt.x): #if the x coordinate of the psuedo pt is greater than middle x
				#add the psuedo pt to the right box
				rightBox.psuedoPts.add(psdpt)
			elif (psdpt.x < self.centerpt.x): #if x coor is less than x coor of center pt, add psd pt to leftBox
				leftBox.psuedoPts.add(psdpt)
			else: #otherwise, the x coor of the psd pt is equal to the x coor of the center pt
				rightBox.psuedoPts.add(psdpt)
				leftBox.psuedoPts.add(psdpt)

		return leftBox, rightBox

	def spliceBoxAlongY(self):
		centerLeftShift = self.centerpt.y - 0.25 * self.y
		centerRightShift = self.centerpt.y + 0.25 * self.y

		leftcenter = Point(-1, self.centerpt.x, centerLeftShift, self.centerpt.z)
		rightcenter = Point(-1, self.centerpt.x, centerRightShift, self.centerpt.z)

		leftBox = BoundedBox(self.x, self.y / 2.0, self.z, leftcenter, self)
		rightBox = BoundedBox(self.x, self.y / 2.0, self.z, rightcenter, self)

		self.child1 = leftBox
		self.child2 = rightBox

		# leftBox.getPointsinBox(self.pointsInside)
		# rightBox.getPointsinBox(self.pointsInside)

		leftBox.pointsInside, rightBox.pointsInside= self.splitPoints('y')


		for psdpt in self.psuedoPts:
			if (psdpt.y > self.centerpt.y):
				rightBox.psuedoPts.add(psdpt)
			elif (psdpt.y < self.centerpt.y):
				leftBox.psuedoPts.add(psdpt)
			else:
				rightBox.psuedoPts.add(psdpt)
				leftBox.psuedoPts.add(psdpt)

		return leftBox, rightBox

	def spliceBoxAlongZ(self):
		centerLeftShift = self.centerpt.z - 0.25 * self.z
		centerRightShift = self.centerpt.z + 0.25 * self.z

		leftcenter = Point(-1, self.centerpt.x, self.centerpt.y, centerLeftShift)
		rightcenter = Point(-1, self.centerpt.x, self.centerpt.y, centerRightShift)

		leftBox = BoundedBox(self.x, self.y, self.z / 2.0, leftcenter, self)
		rightBox = BoundedBox(self.x, self.y, self.z / 2.0, rightcenter, self)

		self.child1 = leftBox
		self.child2 = rightBox

		# leftBox.getPointsinBox(self.pointsInside)
		# rightBox.getPointsinBox(self.pointsInside)
		
		leftBox.pointsInside, rightBox.pointsInside= self.splitPoints('z')


		for psdpt in self.psuedoPts:
			if (psdpt.z > self.centerpt.z):
				rightBox.psuedoPts.add(psdpt)
			elif (psdpt.z < self.centerpt.z):
				leftBox.psuedoPts.add(psdpt)
			else:
				rightBox.psuedoPts.add(psdpt)
				leftBox.psuedoPts.add(psdpt)

		return leftBox, rightBox
	
	"""returns the longest edge of the box as a string 'x', 'y', or 'z' ."""
	def getLongestEdge(self):
		longest = max(self.x, self.y, self.z)
		if (longest == 0):
			raise Exception("you can't have a box with zero dimensions....")
		elif (longest == self.x):
			return 'x'
		elif (longest == self.y):
			return 'y'
		elif (longest == self.z):
			return 'z'


