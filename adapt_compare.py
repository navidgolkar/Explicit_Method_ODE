import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.integrate import odeint
import methods.explicit as expl

# Define a range of tolerances for time grids with adaptive (variable) step size.
#tol_v = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8] # long runtime
#tol_v = [1e-3, 1e-4, 1e-5, 1e-6, 1e-7] # long runtime
tol_v = [1e-2, 1e-3, 1e-4, 1e-5]

# Define a range step numbers for equidistant time grids.
#N_v = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8] # long runtime
#N_v = [1e3, 1e4, 1e5, 1e6, 1e7] # long runtime
N_v = [1e3, 1e4, 1e5, 1e6]

t0 = 0               # integration interval [t0,te] 
te = 25
y0 = np.array([0,8]) # initial condition

# time grid for odeint
h = 0.1
t_odeint = np.arange(t0,te+h,h)  

# various values for the factor mu appearing in the ODE function
# y'' -  mu*(1-y^2)y'+ y = 0.
mu = [0.01,1,15]
#mu_k = [15]

for k in range(len(mu)):
  mu_k = mu[k]
  # y'' =  mu*(1-y^2)y'+ y  
  # rewritten as a first order system:
  #   y0' = y1
  #   y1' = y0'' = mu*(1-y0^2)y1+ y0  
  vanderpol = lambda t,y: np.array([y[1], mu_k*(1-y[0]**2)*y[1] - y[0]])      
  errors = np.zeros((8,len(tol_v))) # initialize arrays for global errors
  times = np.zeros((8,len(tol_v)))  # initialize arrays for runtimes
  #y_odeint = odeint(vanderpol,y0,t_odeint,tfirst = True, rtol=1e-10, atol=1e-12)
  y_odeint = odeint(vanderpol,y0,t_odeint,tfirst = True)
      
  for n in range(len(tol_v)):
    
    # run methods with adaptive step size
    start = time.time()
    t_1, y_1 = expl.rk21(vanderpol,y0,t0,te,tol=tol_v[n])
    end = time.time()
    times[0,n] = end - start # runtime of improved Euler with adaptive step size
    errors[0,n] = max(abs(y_odeint[-1,:] - y_1[-1,:])) # global error at te
    # run Runge-Kutta with adaptive step size
    start = time.time()
    t_2, y_2 = expl.euler_adaptive(vanderpol,y0,t0,te,tol=tol_v[n])
    end = time.time()
    times[1,n] = end - start # runtime of improved Euler with adaptive step size
    errors[1,n] = max(abs(y_odeint[-1,:] - y_2[-1,:])) # global error at te
    
    
    # run methods with constant step size
    Nsteps = int(N_v[n])+1
    t_const = np.linspace(t0,te,Nsteps)
    
    #no.1: Explicit Euler
    start = time.time()
    y_const = expl.euler_exp(vanderpol, y0, t_const)
    end = time.time()
    times[2,n] = end - start  # runtime of improved Euler with constant step size
    errors[2,n] = max(abs(y_odeint[-1,:] - y_const[-1,:])) # global error at te
    
    #no.2: Improved Euler (RK-2)
    start = time.time()
    y_const = expl.euler_improved(vanderpol, y0, t_const)
    end = time.time()
    times[3,n] = end - start  # runtime of improved Euler with constant step size
    errors[3,n] = max(abs(y_odeint[-1,:] - y_const[-1,:])) # global error at te
    
    #no.3: Heun-2
    start = time.time()
    y_const = expl.heun2(vanderpol, y0, t_const)
    end = time.time()
    times[4,n] = end - start  # runtime of improved Euler with constant step size
    errors[4,n] = max(abs(y_odeint[-1,:] - y_const[-1,:])) # global error at te
    
    #no.4: Ralston
    start = time.time()
    y_const = expl.ralston(vanderpol, y0, t_const)
    end = time.time()
    times[5,n] = end - start  # runtime of improved Euler with constant step size
    errors[5,n] = max(abs(y_odeint[-1,:] - y_const[-1,:])) # global error at te
    
    #no.5: Heun-3
    start = time.time()
    y_const = expl.heun3(vanderpol, y0, t_const)
    end = time.time()
    times[6,n] = end - start  # runtime of improved Euler with constant step size
    errors[6,n] = max(abs(y_odeint[-1,:] - y_const[-1,:])) # global error at te
    
    #no.6: Classical RK-4
    start = time.time()
    y_const = expl.rk4(vanderpol, y0, t_const)
    end = time.time()
    times[7,n] = end - start  # runtime of improved Euler with constant step size
    errors[7,n] = max(abs(y_odeint[-1,:] - y_const[-1,:])) # global error at te
  
  
  plt.rcParams.update({'font.size': 6})
  # function plot of y_1(t)
  fig,ax = plt.subplots(3,1,sharex=True)
  
  ax[0].plot(t_odeint,y_odeint[:,0],'k',t_1, y_1[:,0],'-.g',t_2, y_2[:,0],'--b')
  ax[0].legend(['RK45','RK2 adaptive','Euler adaptive'], loc = 'best')
  ax[0].set_title('function plot for $\mu$ = ' + str(mu_k))
    
  ax[1].plot(t_1[0:-1],np.diff(t_1)) # Plot time vs stepsize 
  ax[1].set_title('stepsize adaptive RK-2 for $\mu$ = ' + str(mu_k))
  
  ax[2].plot(t_2[0:-1],np.diff(t_2)) # Plot time vs stepsize 
  ax[2].set_title('stepsize adaptive Euler for $\mu$ = ' + str(mu_k))
  
  plt.figure()
  # Plot runtimes vs global errors
  plt.loglog(times[0,:],errors[0,:],'-*b',times[1,:],errors[1,:],'-*g',times[2,:],errors[2,:],'-or',times[7,:],errors[7,:],'-oy')
  plt.xlabel('time [s]')
  plt.ylabel('error')
  plt.legend(['RK2 adaptive','Euler adaptive','Explicit Euler','Classical RK-4'], loc = 'best')   
  plt.title('runtime vs global error for $\mu$ = ' + str(mu_k))
  
  plt.figure()
  # Plot runtimes vs global errors
  plt.loglog(times[2,:],errors[2,:],'-or',times[3,:],errors[3,:],'-og',times[4,:],errors[4,:],'-ob',times[5,:],errors[5,:],'-oc',times[6,:],errors[6,:],'-om',times[7,:],errors[7,:],'-oy')
  plt.xlabel('time [s]')
  plt.ylabel('error')
  plt.legend(['Explicit Euler','Improved Euler','Heun-2','Ralston','Heun-3','Classical RK-4'], loc = 'best')   
  plt.title('runtime vs global error for $\mu$ = ' + str(mu_k))
