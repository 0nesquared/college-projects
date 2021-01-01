import matplotlib.pyplot as plt 
import numpy
from scipy.integrate import odeint

full_output = 1

def f(y, x, params):
    a, m = y                                                            # unpack current values of y
    n1, n2, n3, n4, k1, k2, k3, alpha = params                          # unpack parameters
    derivs = [1.0 - a - m - k1*(a**n2)*((m+alpha*a)**n3)*((1-a-m)**(-n4)) +k3*(m**n1),           # list of dy/dt=f functions
             k1*(a**n2)*((m+alpha*a)**n3)*((1-a-m)**(-n4)) + k2*(m**n1)*(a+m-1.0) - k3*(m**n1)]
    return derivs

# Parameters
n1 = 1.8
n2 = 1.0
n3 = 3.0
n4 = 1.0    
k1 = 1.0     
k2 = 50.0
k3 = 0.0
alpha= 0.2     

# Initial values
a0 = 0.2     
m0 = 0.3     
s0 = 0.5

# Bundle parameters for ODE solver
params = [n1, n2, n3, n4, k1, k2, k3, alpha]

# Bundle initial conditions for ODE solver
y0 = [a0, m0]

# Make time array for solution
xStop = 20000.0
dx = 0.0002
x = numpy.arange(0., xStop, dx)

# Call the ODE solver
psoln = odeint(f, y0, x, args=(params,))

t = []
a_F = []
m_F = []
for i in numpy.arange(0.,xStop,10):
    t.append(i)
    a_F.append(float(psoln[i,0]))
    m_F.append(float(psoln[i,1]))
t = numpy.array(t)
a_F = numpy.array(a_F)
m_F = numpy.array(m_F)
s_F = 1.0 - a_F - m_F
 

# Plot results as functions of time
fig = plt.figure(1, figsize=(8,8))
# Plot a as a function of time
ax1 = fig.add_subplot(311)
ax1.plot(t, a_F)
ax1.set_xlabel('time')
ax1.set_ylabel('a')
# Plot m as a function of time
ax2 = fig.add_subplot(312)
ax2.plot(t, m_F)
ax2.set_xlabel('time')
ax2.set_ylabel('m')
# Plot s as a function of time
ax2 = fig.add_subplot(313)
ax2.plot(t, s_F)
ax2.set_xlabel('time')
ax2.set_ylabel('s')

plt.tight_layout()

figP = plt.figure(2, figsize=(8,8))
# Plot a as a function of time
ax1 = figP.add_subplot(311)
ax1.plot(m_F, s_F)
ax1.set_xlabel('MF')
ax1.set_ylabel('SF')
# Plot m as a function of time
ax2 = figP.add_subplot(312)
ax2.plot(m_F, a_F)
ax2.set_xlabel('MF')
ax2.set_ylabel('AF')
# Plot s as a function of time
ax2 = figP.add_subplot(313)
ax2.plot(s_F, a_F)
ax2.set_xlabel('SF')
ax2.set_ylabel('AF')

plt.tight_layout()

plt.show()