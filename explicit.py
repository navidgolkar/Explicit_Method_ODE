import numpy as np

'''
Explicit Euler Method
y_{n+1} = y_n + h*f(x_n, y_n)

Butcher tableau:
0   |   0
------------------
    |   1
'''
def euler_exp(f,y0,t):    
  y = np.zeros((len(t),len(y0)))         # initialize y 
  y[0,:] = y0                            # store the initial value y0 as the first row of y
  for n in range(len(t)-1):  
    h = t[n+1]-t[n]                      # current step size
    y[n+1,:] = y[n,:] + h*f(t[n],y[n,:])
  return y


'''
Improved Euler (Runge-Kutta 2-Step) Method
y_{n+1/2} = y_n + (h/2)*f(x_n, y_n)
y_{n+1} =  y_n + h*(f(x_{n+1/2}, y_{n+1/2}))

Butcher tableau:
0   |   0       0
1/2 |   1/2     0
------------------
    |   0       1
'''
def euler_improved(f,y0,t):    
  y = np.zeros((len(t),len(y0)));           # initialize y 
  y[0,:] = y0                               # store the initial value y0 as the first row of y  
  for n in range(len(t)-1):  
    h = t[n+1]-t[n]                         # current step size
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+h/2.0,y[n,:]+(h/2.0)*k1 )
    #y_half = y[n,:] + (h/2.0)*f(t[n],y[n,:])    # intermediate half step (y_{n+1/2})
    #y[n+1] = y[n,:] + h*f(t[n]+h/2.0,y_half)  # improved Euler method
    y[n+1] = y[n,:] + h*k2  # improved Euler method
  return y


'''
Heun 2-Step Method

Butcher tableau:
0   |   0       0
1   |   1       0
------------------
    |   1/2     1/2
'''
def heun2(f,y0,t):
  y = np.zeros((len(t),len(y0)))         # initialize y  
  y[0,:] = y0                            # store the initial condition in first row
  for n in range(len(t)-1):
    h = t[n+1]-t[n]  
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+h,y[n,:]+h*k1 )
    y[n+1,:] = y[n,:] + (h/2.0)*(k1 + k2)
  return y


'''
Ralston Method

Butcher tableau:
0   |   0       0
2/3 |   2/3     0
------------------
    |   1/4     3/4
'''
def ralston(f,y0,t):
  y = np.zeros((len(t),len(y0)))         # initialize y  
  y[0,:] = y0                            # store the initial condition in first row
  for n in range(len(t)-1):
    h = t[n+1]-t[n]  
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+(2.0/3.0)*h,y[n,:]+h*(2.0/3.0)*k1 )
    y[n+1,:] = y[n,:] + (h/4.0)*(k1 + 3*k2)
  return y


'''
General Runge-Kutta 2-Step method

Butcher tableau:
0   |   0           0
a   |   a           0
-------------------------
    |   1-1/(2a)    1/(2a)
    
which gives:
    a=1/2   ->  improved euler
    a=1     ->  Heun 2-step
    a=2/3   ->  Ralston
'''
def g_rk2(f,y0,t,alpha=0.5):
  y = np.zeros((len(t),len(y0)))         # initialize y  
  y[0,:] = y0                            # store the initial condition in first row
  for n in range(len(t)-1):
    h = t[n+1]-t[n]  
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+alpha*h,y[n,:]+h*alpha*k1 )
    y[n+1,:] = y[n,:] + (h/(2.0*alpha))*((2*alpha-1)*k1 + k2)
  return y


'''
Adams-Bashforth 2-Step Method
y_{n+2} = y_{n+1} + (3/2)*h*f(t_{n+1}, y_{n+1}) - (1/2)*h*f(t_n, y_n)
'''
def ab2(f, y0, y1, t):
    y = np.zeros((len(t),len(y0)))      # initialize y
    y[0,:] = y0                         # store the initial value y0 as the first row of y
    y[1,:] = y1
    for n in range(len(t)-2):
        h = t[n+1]-t[n]                 # current step size
        y[n+2,:] = y[n+1,:] + (3.0/2.0)*h*f(t[n+1], y[n+1,:]) - (1.0/2.0)*h*f(t[n], y[n,:])
    return y


'''
Adams-Bashforth 3-Step Method
y_{n+3} = y_{n+1} + (23/12)*h*f(t_{n+2}, y_{n+2}) - (16/12)*h*f(t_{n+1}, y_{n+1}) - (5/12)*h*f(t_n, y_n)
'''
def ab3(f, y0, y1, y2, t):
    y = np.zeros((len(t),len(y0)), dtype='complex')     # initialize y
    y[0,:] = y0                                         # store the initial value y0 as the first row of y
    y[1,:] = y1
    y[2,:] = y2
    for n in range(len(t)-3):
        h = t[n+1]-t[n]                                 # current step size
        y[n+3,:] = y[n+2,:] + h*((23.0/12.0)*f(t[n+2], y[n+2,:]) - (16.0/12.0)*f(t[n+1], y[n+1,:]) + (5.0/12.0)*f(t[n], y[n,:]) )
    return y


'''
Heun 3-Step Method

Butcher tableau:
0   |   0       0       0
1/3 |   1/3     0       0
2/3 |   0       2/3     0
--------------------------
    |   1/4     0       3/4
'''
def heun3(f,y0,t):
  y = np.zeros((len(t),len(y0)))         # initialize y  
  y[0,:] = y0                            # store the initial condition in first row
  for n in range(len(t)-1):
    h = t[n+1]-t[n]  
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+h/3.0,y[n,:]+(h/3.0)*k1 )
    k3 = f(t[n]+(2.0/3.0)*h,y[n,:]+(2.0/3.0)*h*k2 )
    y[n+1,:] = y[n,:] + (h/4.0)*(k1 + 3*k3)
  return y


'''
Classical Runge-Kutta 4-Step Method

Butcher tableau:
0   |   0       0       0       0
1/2 |   1/2     0       0       0
1/2 |   0       1/2     0       0
1   |   0       0       1       0
----------------------------------
    |   1/6     1/3     1/3     1/6
'''
def rk4(f,y0,t):
  y = np.zeros((len(t),len(y0)))         # initialize y  
  y[0,:] = y0                            # store the initial condition in first row
  for n in range(len(t)-1):
    h = t[n+1]-t[n]  
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+h/2.0,y[n,:]+(h/2.0)*k1)
    k3 = f(t[n]+h/2.0,y[n,:]+(h/2.0)*k2)
    k4 = f(t[n]+h,y[n,:]+(h*k3))
    y[n+1,:] = y[n,:] + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)
  return y


'''
(Adaptive) Runge-Ketta 2-Step Method with step size control

Butcher tableau:
0     |   0       0
1     |   1       0
------------------
p=2   |   1/2     1/2
p'=1  |   1       0
'''
def rk21(f, y0, t0, te, tol=1e-5, hstart=0.1, scale_factor=0.9, hmin=0.2, hmax=10):
  n = 0
  t = np.array([t0]) # initial condition 
  y = np.array([y0]) # y(t0) = y0
  htmp = hstart
  while htmp > 0:
    k1 = f(t[n],y[n,:])
    k2 = f(t[n]+htmp,y[n,:]+htmp*k1)
    y1 = y[n,:] + (htmp/2.0)*(k1 + k2)
    y2 = y[n,:] + htmp*k1
    phi = np.linalg.norm(y2-y1)
    # Note how the change in step size is limited to 0.2 and 10 (hmin and hmax)
    hneu = htmp * np.min([np.max([scale_factor*np.sqrt(tol/phi), hmin]), hmax])
    if phi > tol:
        htmp = hneu
    else:
        h = htmp
        t = np.append(t,[t[n]+h],axis = 0)   # Create a column vector.       
        y = np.append(y,[y1],axis = 0)  # Use the value provided by the more accurate scheme (which equals the improved Euler method). 
        htmp = np.min([te-t[n+1],hneu])      # Reduce the time stepsize if close to end time te
        n = n+1
  return t,y


'''
(Adaptive) Explicit Euler Method with step size control
'''
def euler_adaptive(f, y0, t0, te, tol=1e-5, hstart=0.1, scale_factor=0.9, hmin=0.2, hmax=10):
  n = 0
  t = np.array([t0]) # initial condition 
  y = np.array([y0]) # y(t0) = y0
  htmp = hstart
  while htmp > 0:
    fn = f(t[n],y[n,:]); 
    y1 = y[n,:] + htmp*fn
    y2 = y[n,:] + htmp/2.0*fn
    y2 = y2 + htmp/2.0 * f(t[n]+htmp/2.0,y2)
    phi = 2 * np.linalg.norm(y2-y1)
    # Note how the change in step size is limited to 0.2 and 10 (hmin and hmax)
    hneu = htmp * np.min([np.max([scale_factor*np.sqrt(tol/phi), hmin]), hmax])
    if phi > tol:
      htmp = hneu
    else:
      h = htmp
      t = np.append(t,[t[n]+h],axis = 0)   # Create a column vector.       
      y = np.append(y,[2*y2-y1],axis = 0)  # Use the value provided by the more accurate scheme (which equals the improved Euler method). 
      htmp = np.min([te-t[n+1],hneu])      # Reduce the time stepsize if close to end time te
      n = n+1
  return t,y
