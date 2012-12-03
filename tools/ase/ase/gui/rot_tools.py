# Gives the rotation matrix which rotates theta degrees about
# vecU

#  Generates the rotation matrix that rotate theta degrees about the vecU
def rotate_about_vec(vecU, theta):
    import numpy as np
    vecU = np.array(vecU)
    vecU = vecU / (sum(vecU ** 2) ** 0.5)
    ux, uy, uz = vecU
    st = np.sin(theta)
    ct = np.cos(theta)
    mat = np.array([[ux ** 2 + ct * (1 - ux ** 2), 
                     ux * uy * (1 - ct) - uz * st, 
                     uz * ux * (1 - ct) + uy * st],
                    [ux * uy * (1 - ct) + uz * st, 
                     uy ** 2 + ct * (1 - uy ** 2),
                     uy * uz * (1 - ct) - ux * st],
                    [uz * ux * (1 - ct) - uy * st,
                     uy * uz * (1 - ct) + ux * st,
                     uz ** 2 + ct * (1 - uz **2)]])
    return (mat)

# Generates the rotation matrix which rotates aVec into intoVec
def rotate_vec_into_newvec(aVec, intoVec):
    def length(v):
        return((sum(v ** 2)) ** 0.5)

    import numpy as np
    from math import acos
    fac = 1.0
    aVec = np.array(aVec)
    intoVec = np.array(intoVec)
    nor = np.cross(aVec, intoVec)
    if length(nor) == 0:
        nor = np.array([1, 0, 0])
    nor = nor / length(nor)
    theta = acos(np.dot(aVec, intoVec) / (length(aVec) * length(intoVec)))
    if np.dot(aVec, intoVec) < 0:
        theta = theta + np.pi
        fac = -1
    return(fac * rotate_about_vec(nor, theta))

# Applies the rotation matrix to the vector and returns the rotated vector
def rotate_vec (rot_mat, vec):
    import numpy as np
    rot_vec = np.dot(rot_mat, vec)

    return (rot_vec)
