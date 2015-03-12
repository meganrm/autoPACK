def binary_search(ptList, cutValue, axis):
    max_index=len(ptList)-1
    min_index=0
    while max_index > min_index:
      mid_index = min_index + ((max_index - min_index)/2)
#      print "new mid_index=",mid_index
      if ptList[mid_index].__dict__[axis] > cutValue:
        max_index = mid_index-1
#        print "ptList[mid_index] > cutValue, new max_index=",max_index
      else:  # so A[mid_index] <= x
        if ptList[mid_index+1].__dict__[axis] > cutValue: 
          return mid_index
        else:
          min_index = mid_index + 1
#          print "move mid_index up one, mid_index=",mid_index
    return min_index
