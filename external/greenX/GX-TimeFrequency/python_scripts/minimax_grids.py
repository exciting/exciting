import numpy as np
import matplotlib.pyplot as pl
from matplotlib import patches
from scipy.optimize import curve_fit, root, fsolve
from scipy.signal import argrelextrema
from numpy import dot, outer
from os import listdir
from os.path import isfile, join

def main():
    
    # Set parameters
    n_minimax = 16                     # Number of minimax points
    R_minimax = int(1E13)              # Range of the minimax approximation
    n_x       = 12000                   # total number of points on the x-axis for optimization
    eps_diff  = 10**(-10)

    R_increase_method = 'multiply'
#    R_increase_method = 'add'
    R_add = 1
    R_mult = 1.01
    R_fac_decrease = 1.2

    path = "alpha_beta_of_N_"+str(n_minimax)

    alphas_betas_init = np.loadtxt("../alpha_beta_of_N_"+str(n_minimax),dtype=np.float128)
#    alphas_betas_init = np.loadtxt("alpha_beta_of_N_30_R_0000000002840_E_1.997E-14",dtype=np.float128)
    converged_alphas_betas_E = alphas_betas_init

    ydata = np.zeros(n_x,dtype=np.float128)

    alphas_betas_E = np.append(alphas_betas_init,1)

    while True:

       E_old = alphas_betas_E[-1]*2
   
       i = 0
       while (alphas_betas_E[-1]/E_old < 1-eps_diff or alphas_betas_E[-1] > E_old):

           E_old = alphas_betas_E[-1]
           num_extrema = 0

           while num_extrema < 2*n_minimax+1:

              xdata = 10**(np.logspace(0,np.log10(np.log10(R_minimax)+1),n_x,dtype=np.float128))/10

              extrema_x = np.append(xdata[0], xdata[argrelextrema(eta_plotting(xdata,alphas_betas_E[0:np.size(alphas_betas_E)-1]), np.greater)[0]])
              if np.size(extrema_x) == n_minimax: 
                 extrema_x = np.append(extrema_x, xdata[-1])
              extrema_x = np.append(extrema_x, xdata[argrelextrema(eta_plotting(xdata,alphas_betas_E[0:np.size(alphas_betas_E)-1]), np.less)[0]])
              num_extrema = np.size(extrema_x)
              print("num_extrema", num_extrema)
              if(num_extrema < 2*n_minimax+1):
                  alphas_betas_E = converged_alphas_betas_E
                  if(R_increase_method == 'multiply'):
                      R_minimax = int(R_minimax*R_mult)
                  elif(R_increase_method == 'add'):
                      R_minimax += R_add

           alphas_betas_E[np.size(alphas_betas_E)-1] = np.average(np.abs(eta_plotting(extrema_x,alphas_betas_E[0:np.size(alphas_betas_E)-1])))
           i += 1
           print("iteration =", i, "E =",  alphas_betas_E[-1], "Range =", R_minimax)

           alphas_betas_E = my_fsolve(extrema_x, alphas_betas_E)

           break_i = False
           if(i>500): 
               break_i = True
               break
  
       converged_alphas_betas_E = alphas_betas_E

       sort_indices = np.argsort(alphas_betas_E[0:n_minimax])
       num_zeros = 13-len(str(R_minimax))
   
       np.savetxt("alpha_beta_of_N_"+str(n_minimax)+"_R_"+"0"*num_zeros+str(R_minimax)+"_E_"+\
               np.array2string(np.amax(np.abs(eta_plotting(extrema_x,alphas_betas_E))), formatter={'float_kind':lambda x: "%.3E" % x}), \
               np.append(alphas_betas_E[sort_indices],alphas_betas_E[sort_indices+n_minimax]) )
   
       fig1, (axis1) = pl.subplots(1,1)
       axis1.set_xlim((0.8,R_minimax))
       axis1.semilogx(xdata,eta_plotting(xdata,alphas_betas_E))
       axis1.semilogx([0.8,R_minimax], [alphas_betas_E[-1],alphas_betas_E[-1]])
       axis1.semilogx([0.8,R_minimax], [-alphas_betas_E[-1],-alphas_betas_E[-1]])
       pl.show()

       if(break_i):
           R_minimax += R_add
       else:
           R_minimax = int(R_minimax/R_fac_decrease)

def build_denominator_matrix(xdata,alphas):
    matrix = []
    for x in xdata:
      matrix.append((2*x/(x**2+np.array(alphas)**2))**2/np.pi)
    return matrix

def eta_plotting(x, *params):
    params_1d = np.transpose(params)[:,0]
    denominator_matrix = build_denominator_matrix(x, np.array(params_1d[0:np.size(params)//2]))
    return 1/x - np.array(denominator_matrix).dot(params_1d[np.size(params)//2:(np.size(params)//2)*2])

def my_fsolve(extrema_x, alphas_betas_E):
    size_problem = np.size(alphas_betas_E)
    n_minimax = (size_problem-1)//2

    E = np.empty(size_problem, dtype=np.float128)
    E[0:size_problem//2+1] = alphas_betas_E[-1]
    E[size_problem//2+1:] = -alphas_betas_E[-1]

    vec_f = eta_plotting(extrema_x, alphas_betas_E) - E

    mat_J = np.zeros((size_problem, size_problem+1),dtype=np.float128)

    for index_i in range(n_minimax):
        mat_J[:,index_i] = alphas_betas_E[index_i+n_minimax]*16/np.pi*alphas_betas_E[index_i]*extrema_x**2/(extrema_x**2+(alphas_betas_E[index_i])**2)**3
        mat_J[:,index_i+n_minimax] = -4/np.pi*extrema_x**2/(extrema_x**2+(alphas_betas_E[index_i])**2)**2

    mat_J[:,-2] = -np.sign(E[0:size_problem])
    mat_J[:,-1] = vec_f[0:size_problem]

    delta = gauss(mat_J)

    return alphas_betas_E - delta

def gauss(A):
    n = len(A)

    for i in range(0, n):
        # Search for maximum in this column
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k

        # Swap maximum row with current row (column by column)
        for k in range(i, n+1):
            tmp = A[maxRow][k]
            A[maxRow][k] = A[i][k]
            A[i][k] = tmp

        # Make all rows below this one 0 in current column
        for k in range(i+1, n):
            c = -A[k][i]/A[i][i]
            for j in range(i, n+1):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]

    # Solve equation Ax=b for an upper triangular matrix A
    x = [0 for i in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = A[i][n]/A[i][i]
        for k in range(i-1, -1, -1):
            A[k][n] -= A[k][i] * x[i]
    return x

if __name__ == "__main__":
    main()

