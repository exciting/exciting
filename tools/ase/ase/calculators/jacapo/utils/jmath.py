'''
mathematical utilities
'''
import numpy as np
import bisect

def vinterp3d(x,y,z,u,xi,yi,zi):

    p = np.array([xi,yi,zi])

    #1D arrays of cooridinates
    xv = x[:,0,0]
    yv = y[0,:,0]
    zv = z[0,0,:]

    # we subtract 1 because bisect tells us where to insert the
    # element to maintain an ordered list, so we want the index to the
    # left of that point
    i = bisect.bisect_right(xv,xi) - 1
    j = bisect.bisect_right(yv,yi) - 1
    k = bisect.bisect_right(zv,zi) - 1

    #occasionally we get the edge of the cell and we then go one index
    #back again
    if i == len(xv)-1: i-=1
    if j == len(yv)-1: j-=1
    if k == len(zv)-1: k-=1
    
    
    #points at edge of cell. We only need P1, P2, P3, and P5
    P1 = np.array([x[i,j,k],y[i,j,k],z[i,j,k]])
    P2 = np.array([x[i+1,j,k],y[i+1,j,k],z[i+1,j,k]])
    P3 = np.array([x[i,j+1,k],y[i,j+1,k],z[i,j+1,k]])
    P5 = np.array([x[i,j,k+1],y[i,j,k+1],z[i,j,k+1]])

    #values of u at edge of cell
    u1 = u[i,j,k]
    u2 = u[i+1,j,k]
    u3 = u[i,j+1,k]
    u4 = u[i+1,j+1,k]
    u5 = u[i,j,k+1]
    u6 = u[i+1,j,k+1]
    u7 = u[i,j+1,k+1]
    u8 = u[i+1,j+1,k+1]

    #cell basis vectors, not the unit cell, but the voxel cell containing the point
    cbasis = np.array([P2-P1,
                       P3-P1,
                       P5-P1])

    #now get interpolated point in terms of the cell basis
    s = np.dot(np.linalg.inv(cbasis.T),np.array([xi,yi,zi])-P1)

    #now s = (sa, sb, sc) which are fractional coordinates in the vector space
    #next we do the interpolations
    ui1 = u1 + s[0]*(u2-u1)
    ui2 = u3 + s[0]*(u4-u3)

    ui3 = u5 + s[0]*(u6-u5)
    ui4 = u7 + s[0]*(u8-u7)

    ui5 = ui1 + s[1]*(ui2-ui1)
    ui6 = ui3 + s[1]*(ui4-ui3)

    ui7 = ui5 + s[2]*(ui6-ui5)

    return ui7

