from sympy import Symbol
import numpy as np

# Sistema de ecuaciones a utilizar

def function_exercise(wxyz):
    w, x, y, z = wxyz
    return [(w**2)+x-(3*y)+(4*z)+(3/4),
            (3*w**2)+x-(y**2)+(z**2)+(13/4),
            (5*w)+(3*x**2)+y-(4*z**2)-(99/2),(8*w**2)-(14*x)+(6*y**2)-(7*z**2)+7]

#Derivada del sistema de ecuaciones a utilizar

def jacobian_exercise(wxyz):
    w, x, y, z = wxyz
    return [[2*w,1,-3,4],
            [6*w,1,-2*y,2*z],
            [5,6*x,1,-8*z],[16*w,-14,12*y,-14*z]]

def JF(x):
    
    J_save = np.array(jacobian_exercise(x))
    
    return J_save

def F(x):
    
    F_save =  np.array(function_exercise(x))
    
    return F_save

def solve(x_0, tol):

            w = []
            x = []
            y = []
            z = []
            err = []

            a = (1) 

            b = (-2)

            itera = 0

            ep = tol + 1

            while(ep > tol):  #Condición de parada

                Fx = F(x_0)         #F(x_k)

                Jx_kF = np.linalg.solve(JF(x_0),-Fx)        #Mult J-1*F(x)

                y_k = x_0 + 1/2*(Jx_kF)
                

                JF_sol = np.linalg.solve(JF(y_k),-F(x_0))
               
            
                z_k = x_0 + JF_sol

                M = (a *JF(x_0) + b *  JF(y_k))
                
               
                Fz = F(z_k)
                
                M_solve = np.linalg.solve(M,-Fz)
                
                
                
                x_k = z_k + M_solve
                
                

                x_0 = x_k

                itera += 1

                w.append(x_0[0])
    
                x.append(x_0[1])
                
                y.append(x_0[2])
                
                z.append(x_0[3])
                
                err.append(ep)
                
                ep = np.linalg.norm(F(x_0),2)
                print(ep)
                
            return x_k,itera

print(solve([2,1,2,1],0.1))
    
    
