#!/usr/bin/env python
def get_boundary(Polygon_no,NROWS,NCOLS_EVEN,NCOLS_ODD):
    boundaries = [[1,2],[0,1] ,[-1,0], [4,-1], [3,4],[2,3]]
    #bottom side
    if Polygon_no in range(1,NCOLS_ODD):    
        boundaries = [[1,2],[0,1] ,[-1,0], [4,-1]]
    #left side odd
    if Polygon_no %(NCOLS_EVEN+NCOLS_ODD)==0:
        boundaries = [[1,2],[0,1],[2,3]]
    #even
    if Polygon_no %(NCOLS_EVEN+NCOLS_ODD)==NCOLS_ODD:
        boundaries = [[1,2],[0,1] ,[-1,0], [3,4],[2,3]]
    #right side odd
    if Polygon_no %(NCOLS_EVEN+NCOLS_ODD)==NCOLS_EVEN:
        boundaries = [[-1,0], [4,-1], [3,4]]
    #even
    if Polygon_no %(NCOLS_EVEN+NCOLS_ODD)==NCOLS_EVEN+NCOLS_ODD-1:
        boundaries = [[0,1] ,[-1,0], [4,-1], [3,4],[2,3]]
    
    # bottom left
    if Polygon_no == 0:
        boundaries = [[1,2],[0,1]]

    # bottom right
    if Polygon_no == NCOLS_EVEN:
        boundaries = [[-1,0], [4,-1]]
    if NROWS%2 == 1:
        if Polygon_no in range((NROWS/2)*(NCOLS_EVEN+NCOLS_ODD),(NROWS/2)*(NCOLS_EVEN+NCOLS_ODD)+NCOLS_EVEN):
            boundaries = [[1,2], [4,-1], [3,4],[2,3]]
        #top left
        if Polygon_no == NROWS/2*(NCOLS_EVEN+NCOLS_ODD):
            boundaries = [[1,2],[2,3]]
        #top right
        if Polygon_no == NROWS/2*(NCOLS_EVEN+NCOLS_ODD)+NCOLS_EVEN:
            boundaries = [[4,-1], [3,4]]
    else:
        if Polygon_no in range(int(NROWS/2*(NCOLS_EVEN+NCOLS_ODD)-NCOLS_EVEN),int(NROWS/2*(NCOLS_EVEN+NCOLS_ODD)-1)):
            boundaries = [[1,2], [4,-1], [3,4],[2,3]]
        #top left
        if Polygon_no == NROWS/2*(NCOLS_EVEN+NCOLS_ODD)-NCOLS_EVEN:
            boundaries = [[1,2], [3,4],[2,3]]
        #top right
        if Polygon_no == NROWS/2*(NCOLS_EVEN+NCOLS_ODD)-1:
            boundaries = [[4,-1], [3,4],[2,3]]

    #neighbors = [i for i in neighbor_direction]
    return boundaries
    
    