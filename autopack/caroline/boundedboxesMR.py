import sys
import numpy as np
from DataStructuresMR import *

"""Returns a functional list of Point instances from a list VERTICES which contains coordinates in numpy format"""
def parseVerticies(vertices):
	fillThis = []
	count = 0
	for v in vertices:
		fillThis.append(Point(index=count,x=v[0], y=v[1], z=v[2]))
		count += 1
	return fillThis

"""Returns a functional list of Polygon instances from lists FACES which contains coordinates in numpy format
	and PTSLIST which is a list of point instances."""
def parseFaces(faces, ptsList):
	fillThis = []
	count = 0
	for f in faces:
		#print f
		fillThis.append(Polygon([count] + f, ptsList))
		count +=1
	return fillThis

""" takes in list, returns list with 'cm' removed "" """
def removeCM(list):
	while (list.count('cm') > 0):
		list.remove('cm')
	return list

""" Returns true if there are duplicate points in PTSLIST. """
def duplicates(ptsList):
	#i want a hash table
	return False

""" Returns true if there is coordinately identical point to PT in PTSLIST """
def alreadyExists(pt, ptsList):
	for curPt in ptsList:
		if (curPt.x == pt.x and curPt.y == pt.y and curPt.z == pt.z):
			return True
	return False

"""READING/CREATING POINTS WITH COORDINATES, will return the full list of Point instances corresponding
	to the points given"""
def readPoints(pointsFile):
	pointsTXT = open(str(pointsFile), "r")

	ptdata = pointsTXT.readlines()

	pointsTXT.close()

	if (len(ptdata) == 1):
		for line in ptdata:
			coordinates = line.split('\r')
	else:
		coordinates = ptdata

	#the index of each Point in the list corresponds to its Points.number value
	pointsList = []

	coordinates.pop(0)

	for line in coordinates:
		pointData = line.split()
		pointData = removeCM(pointData)
		if (len(pointData) == 0):
			continue
		if (not len(pointData) == 4):
			print "Error reading points file."
			return
		makePoint = Point()
		makePoint.parse(pointData)
		pointsList.append(makePoint)

	return pointsList

"""READING/CREATING POLYGONS, will return a list of Polygon Instances, needs a list of points to
	create a polygon instance"""
def readPolygons(polygonFile, ptsList):
	polygonsTXT = open(str(polygonFile), "r")

	plydata = polygonsTXT.readlines()

	polygonsTXT.close()

	if (len(plydata) == 1):
		for line in plydata:
			polygons = line.split('\r')
	else:
		polygons = plydata

	polygonList = []

	polygons.pop(0)

	for line in polygons:
		polyData = line.split()
		if (len(polyData) == 0):
			continue
		if (not len(polyData) == 4):
			print "Error in reading Polygon file."
			return
		makePolygon = Polygon(polyData, ptsList)
		polygonList.append(makePolygon)

	return polygonList

""" Make Initial Bounded Box given a list PTSLIST, with a buffer of B, default is  0.1."""
def makeInitialBox(vertList, gridPoints, b=0.1):
	ptsList=vertList+gridPoints
	maxXcoord = ptsList[0].x
	minXcoord = ptsList[0].x

	maxYcoord = ptsList[0].y
	minYcoord = ptsList[0].y

	maxZcoord = ptsList[0].z
	minZcoord = ptsList[0].z

	for pt in ptsList:
		maxXcoord = max(maxXcoord, pt.x)
		minXcoord = min(minXcoord, pt.x)

		maxYcoord = max(maxYcoord, pt.y)
		minYcoord = min(minYcoord, pt.y)

		maxZcoord = max(maxZcoord, pt.z)
		minZcoord = min(minZcoord, pt.z)

	xDim = maxXcoord - minXcoord
	yDim = maxYcoord - minYcoord
	zDim = maxZcoord - minZcoord

	xCenter = average(maxXcoord, minXcoord)
	yCenter = average(maxYcoord, minYcoord)
	zCenter = average(maxZcoord, minZcoord)

	#print str(xCenter) + " " + str(yCenter) + " " + str(zCenter)

	centerpt = Point(index=-1, x=xCenter, y=yCenter, z=zCenter)
	bb= BoundedBox(xDim + 2*b, yDim + 2*b, zDim + 2*b, centerpt)
	bb.pointsInside=vertList
	bb.gridPointsInside=gridPoints
	return bb

""" given two boxes of the same parent, will create Psuedo Points on their newly created 
	and shared Edge, according to the edges as defined by the REAL points in the input PTSLIST. 
	Ideally, PTSLIST should be the box1.parent.pointsInside. the function will modify
	plyList whose respective Polygons will gain a created PsuedoPoint if needed.
	Returns a dictionary of psuedo Points (coordinates represented by tuples) mapped to a set of their
	associated polygons"""
def makeBarrierPseudoPts(box1, box2, ptsList, plyList):

	psdptDict = {}

	keepEdge = box1.parent.longEdge
	if (keepEdge != box2.parent.longEdge):
			print "the two boxes must have the same parent"

	fixedCoor = box1.parent.getDirCenterPt(keepEdge)

	for pt in ptsList:
		pt.getSisterPoints(plyList)
		for sispt in pt.sisterPoints:
			#MAKE A PSUEDO POINT WITH LONGEST EDGE IN TACT BETWEEN pt AND sispt
			newpt = makePsuedoPointBtweenEdge(pt, sispt, keepEdge, fixedCoor)

			#ADD THIS PSUEDO POINT TO BOTH BOXES LISTS OF PSUEDO POINTS
			""" satifies conditions:
				1) newpt must be a new psuedo point in box1
				2) newpt must be a new psuedo point in box2 
				3) must not already be a REAl Point
				4) must be in box1 (actually checking if the pt is along the border between box1 and box2
				5) newpt must be between the pt and its corresponding sister point
			"""
			if ((not containsPt(box1.psuedoPts, newpt)) and (not containsPt(box2.psuedoPts, newpt)) and not(isAlreadypt(newpt, ptsList)) and box1.isPtinBox(newpt) and ptBetween(newpt, pt, sispt)):
				box1.psuedoPts.add(newpt)
				box2.psuedoPts.add(newpt)

				#ADD THIS PSUEDO POINT TO THE POLYGON(s) ASSOCIATED WITH THE EDGES OF sispt AND pt
				#ADD THIS PSUEDO POINT TO THE GLOBAL PSD PT DICTIONARY WITH THE ASSOCIATED POLYGON ADDED TO VALUE
				#ADD THE POLYGON(s) ASSOCIATED WITH THIS PSUEDO POINT TO BOX (child) BOXES
				#SINCE THE PSUEDO POINT IS AN INTERSECTION OF AN EDGE/FACE WITH THE BOUNDARY BETWEEN box1 and box2,
				#THE POLYGON(s) ASSOCIATED WITH THIS NEW PSUEDO POINT BELONG TO BOTH box1 AND box2
				for polygon in plyList:
					plyPts = polygon.getPtsList()
					if (containsPt(plyPts, pt) and containsPt(plyPts, sispt)):
						curTuple = newpt.ptToTuple()
						if (curTuple in psdptDict and not (polygon in psdptDict[curTuple])):
							#if an entry in the dictionary already corresponds to the new psuedo point (but the polygon is new)
							#this also means that this Psuedo Point lies in the same place as real point.
							psdptDict[curTuple].add(polygon)
						else: #need to create a new entry in the dictionary
							psdptDict[curTuple] = set()
							psdptDict[curTuple].add(polygon)
						#add the psuedoPoint the polygons it belongs to.
						polygon.psuedoPts.append(newpt)
	return psdptDict

""" Returns True if L has the desired point, otherwise returns false. All items of L considered to be Points,
	or Psuedo Points."""
def containsPt(L, pt):
	for item in L:
		if (item.number == -1):
			return isAlreadypt(pt, L)
		else:
			if (item.number == pt.number and item.x == pt.x and item.y == pt.y and item.z == pt.z):
				return True
	return False

def isAlreadypt(pt, ptsList):
	for item in ptsList:
		if (abs(item.x - pt.x) < 0.00001 and abs(item.y - pt.y) < 0.00001 and abs(item.z - pt.z) < 0.000001):
			return True
	return False


""" creates a psuedo point between pt1 and pt2, basically finds the line between them and then finds the point
	on that line corresponding to the coor in dir (dir = "x", "y", "z"). """
def makePsuedoPointBtweenEdge(pt1, pt2, dir, coor):

	p1Num = pt1.number
	p2Num = pt2.number

	xSlope = pt2.x - pt1.x
	ySlope = pt2.y - pt1.y
	zSlope = pt2.z - pt1.z

	""" lines of equations:
		x = a + xSlope * t
		y = b + ySlope * t
		z = c + zSlope * t """

	a = pt1.x
	b = pt1.y
	c = pt1.z

	t, x, y, z = 0, 0, 0, 0
	#set "n" to coor if dir == "n", n is either x, y, or z

	if (dir == 'x'): #solve for t using the x coordinate, find y and z
		if (xSlope == 0):
			return PsuedoPoint([pt1.x, pt1.y, pt1.z], p1Num, p2Num)
		x = coor
		t = (x - a) / xSlope
		y = b + ySlope * t
		z = c + zSlope * t
		return PsuedoPoint([coor, y, z], p1Num, p2Num)
	elif (dir == 'y'): #solve for t using the y coordinate, find x and z
		if (ySlope == 0):
			return PsuedoPoint([pt1.x, pt1.y, pt1.z], p1Num, p2Num)
		y = coor
		t = (y - b) / ySlope
		x = a + xSlope * t
		z = c + zSlope * t
		return PsuedoPoint([x, coor, z], p1Num, p2Num)
	elif (dir == 'z'):
		if (zSlope == 0):
			return PsuedoPoint([pt1.x, pt1.y, pt1.z], p1Num, p2Num)
		z = coor
		t = (z - c) / zSlope
		x = a + xSlope * t
		y = b + ySlope * t
		return PsuedoPoint([x, y, coor], p1Num, p2Num)
	else:
		print "Bad direction"

""" Iterates Through a list of BoundedBoxes until all  BoundedBoxes are given a Type >= 0, runs
	getPointsinBox on each box in boxLst. will also return a marker variable returns True if there
	are more boxes to examine through further iteration.
	Also returns a dictionary which maps psuedo points (represented by their coordinates in a tuple) mapped
	to a set of their associated polygons"""
def iterateBoxes(boxLst, ptsList, plyList, psdptDict):
	iterateMore = False
	curtype = -1
	newBoxLst = []
	for box in boxLst:
#		box.getPointsinBox(ptsList) #check run time
		box.setBoxType()
		curtype = box.boxType
#		curtype=box.getBoxType()
#		print 'boxtype=', curtype
		if (curtype == -1):
			iterateMore = True
			childBox1, childBox2 = box.spliceBoxLong()
			newBoxLst.append(childBox1)
			newBoxLst.append(childBox2)
			psdptDict = dict(psdptDict.items() + makeBarrierPseudoPts(childBox1, childBox2, ptsList, plyList).items())
			#for each of the new boxes, want them to inherit the appropriate pseudopoints from their parent

		else:#Then, the current box either has one point or no points inside of it
			#Now, a call to getFacesinBox() will be appropriate, and just looks at the psuedopts inside the box
			newBoxLst.append(box)
	return newBoxLst, psdptDict, iterateMore

""" void, runs getFacesinBox() on all the boxs in the BOXLST using the dictionary psdPtDict"""
def runGetFacesinBox(boxLst, psdPtDict):
	for box in boxLst:
		box.getFacesinBox(psdPtDict)

""" run all the iteration, using Point instances in POINTSLIST, and faces/Polygon instances in PLYLIST,
	constructs all bounding boxes (with initial jitter-buffer of b (default = 0.1), starting with an empty box list, return the final boxlist """
def runIteration(vertList,gridPoints, plyList, b=0.1):

	BiggestBox = makeInitialBox(vertList, gridPoints, b)

	psdPttoFace = {}

	moreToDo = True

	boxList = [BiggestBox]

	while (moreToDo):
		boxList, psdPttoFace, moreToDo = iterateBoxes(boxList, vertList, plyList, psdPttoFace)

	runGetFacesinBox(boxList, psdPttoFace)

	return BiggestBox, boxList, psdPttoFace

""" Returns the bounded box that PT is inside as a tuple of variable length.
	ASSUMES runIteration() WAS CALLED BEFORE THIS METHOD IS INVOKED.
	If a point is along the boundary of two bounded boxes, getBox(pt, startBox) will return the smaller of the two,
	or both if the boxes are the same size. """
def getBox(pt, startBox):
	if (not startBox.isPtinBox(pt) or startBox == None):
		#print "There is a problem"
		print "Pt must be in it's starting box"
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
			return goLeft #goRight #both of them
		else:
			return goRight
	if (leftBox.isPtinBox(pt)):
		return getBox(pt, leftBox)
	if (rightBox.isPtinBox(pt)):
		return getBox(pt, rightBox)

""" 2D test.
	Given a point instance PT, find the lowest level bounding box it belongs to and then return
	the distance to the nearest geometry (**not necessarily*** within the box) as well as the box it is inside itself. """
def get2DNearestDistToGeometryWithinBox(pt, startBox, psdptDict):
	lowestLevel = getBox(pt, startBox)

	if (type(lowestLevel) == tuple):
		min = sys.maxint
		for box in lowestLevel:
			dist = box.get2DGeometry(pt, psdptDict, startBox, True)
			if (dist < min):
				min = dist
		return dist, lowestLevel

	return lowestLevel.get2DGeometry(pt, psdptDict, startBox, True), lowestLevel

""" 3D test.
	Given a point instance PT, find the lowest level bounding box it belongs to and then return
	the distance to the nearest geometry (**not necessarily*** within the box) as well as the box it is inside itself. """
def get3DNearestDistToGeometryWithinBox(pt, startBox, psdptDict):
#	lowestLevel = getBox(pt, startBox)

	# if (type(lowestLevel) == tuple):
	# 	min = sys.maxint
	# 	for box in lowestLevel:
	# 		dist = box.get2DGeometry(pt, psdptDict, startBox, True)
	# 		if (dist < min):
	# 			min = dist
	# 	return dist, lowestLevel

	return startBox.get3DGeometry(pt, psdptDict, startBox, True), startBox

""" dictionary which maps Psuedo Points (in tuple form, to make an immutable key) to a set of their respective Polygons. """
psdPttoFace = {}
