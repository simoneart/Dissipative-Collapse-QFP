import numpy as np

'''INITIALIZING THE PARAMETERS'''

m = 1.
R0 = 1.
gamma = 1.
hbar = 1. # 1.0545718e-34

beta = 1. #temperature of the noise
betac = 8.*m*R0**2/hbar**2 #critical value of the inverse temperature
                           #damping dominates as long as beta>betac (to verify)

D1 = (hbar**4)*(beta**2)*gamma/(32.*np.sqrt(np.pi)*(R0**3))
D2 = hbar**2*gamma*m/8./np.sqrt(np.pi)/R0**3 - 3*hbar**3*gamma*beta/64./R0**5 + \
    15*hbar**4*gamma*beta**2/2048./m/R0**7
D3 = 3.*hbar**4*beta**2*gamma/256./R0**5
D4 = hbar**2*D3
f = 3.*hbar**4*beta**2*gamma*np.sqrt(np.pi)/64./R0**5 + \
    hbar**2*gamma*beta*m*np.sqrt(np.pi)/4./R0**3


'''DYNAMICS SOLVER'''

def dyn(IC,Grd,timeAxis,dx,dt):
    
    #phase-space grid
    q_grid = Grd[0]
    p_grid = Grd[1]
    
    wt = np.zeros((len(timeAxis), len(q_grid), len(p_grid)))
    wt[0] = IC
    
    #the parameters for the 1D case are still missing so it wont work
    for t in range(1,len(timeAxis)):
        for i in range(2,len(q_grid)-2):
            for j in range(2,len(p_grid)-2): 
                wt[t,i,j] = wt[t-1,i,j]*(1.-dt/dx*(m*w**2*q_grid[i] + p_grid[j]*\
                (-1./m-f+4*D3+2.*D3*p_grid[j]/dx)-2*D3+2*(D1+D2)/dx+f*dx-D4/4./dx**3))+\
                wt[t-1,i,j+1]*dt/dx*(m*w**2*q_grid[i]+D2/dx+p_grid[j]*(4*D3-f+D3*p_grid[j]/dx))+\
                wt[t-1,i+1,j]*dt/dx*(D1/dx-p_grid[j]/m)+\
                wt[t-1,i,j-1]*dt/dx**2*(D2+D3*p_grid[j]**2)+\
                wt[t-1,i-1,j]*dt/dx**2*D1+\
                +dt*D4/16./dx**4*(wt[t-1,i+2,j+2]-2*wt[t-1,i,j+2]-2*wt[t-1,i,j-2]-\
                                  2*wt[t-1,i+2,j]-2*wt[t-1,i-2,j]+wt[t-1,i-2,j+2]+wt[t-1,i+2,j-2]+\
                                  wt[t-1,i-2,j-2])
    
    return wt

'''PARAMETERS OF THE HARMONIC OSCILLATOR AND GAUSSIAN WIGNER FUNCTION'''

w = 1. 
q_zpf = np.sqrt(hbar/(2*m*w))
p_zpf = np.sqrt(hbar*m*w/2)

def gauss_Wf(r,d,V):
    '''
    Parameters
    ----------
    r : (q,p) vector in phase-space
    d : (<q>, <p>) mean vector
    V : covariance matrix

    Returns
    -------
    w : gaussian Wigner function evaluated at point r of the phase-space
    '''
    
    #swtiching to quadratures 
    for i in range(int(len(r)/2)):
        r[i] /= q_zpf
        r[i+int(len(r)/2)] /= p_zpf
        d[i] /= q_zpf
        d[i+int(len(r)/2)] /= p_zpf
    
    #wrong normalization factor
    w = 1./(2*np.pi*np.sqrt(np.linalg.det(V)))*np.exp(-(r-d).T@np.linalg.inv(V)@(r-d))
    
    return w

'''WIGNER ENTROPY AND ENTROPY PRODUCTION RATE'''

def W_ent(W, domain):
    '''
    Parameters
    ----------
    W : Wigner function of the system
    domain : domain of integration (2D for now (q,p))

    Returns
    -------
    Wigner entropy of the system (approximated, the domain will be finite)
    '''
    
    integrand = np.zeros((len(domain[0]),len(domain[1])))
    
    I1 = np.zeros(len(domain[0]))
    integrand[:,:] = W[:,:] * np.log(W[:,:])
    
    I1[:] = np.trapz(integrand[:],domain[1])
    S = -1 * np.trapz(I1,domain[0])
    
    return S

def W_rel_ent(W1, W2, domain):
    '''
    Parameters
    ----------
    W1 : Wigner function of the system 1
    W2 : Wigner function of the system 2
    domain : domain of integration (2D for now (q,p))
    Returns
    -------
    Relative Wigner entropy between W1 and W2

    '''
    
    integrand = np.zeros((len(domain[0]),len(domain[1])))
    I1 = np.zeros(len(domain[0]))
    
    #integrand[:,:] = W1[:,:] * np.log(W1[:,:]/W2[:,:])
    for i in range(len(domain[0])):
        for j in range(len(domain[0])):
            if W1[i,j]/W2[i,j] <= 1e-14:
                integrand[i,j] = 0.
            else:
                integrand[i,j] = W1[i,j] * np.log(W1[i,j]/W2[i,j])
                
    I1[:] = np.trapz(integrand[:],domain[1])
    S = np.trapz(I1,domain[0])
    
    return S

def entropy_production(Wt, Wtarget, dt, Num_steps, domain):
    '''
    Parameters
    ----------
    Wt : Wigner functions of the system over time
    Wtarget : Wigner function of the target state
    Num_steps : number of time steps
    dt : time step
    domain : domain of integration

    Returns
    -------
    Entropy production over time

    '''
    
    Kt = np.array([W_rel_ent(Wt[n], Wtarget, domain) for n in range(Num_steps)])
    
    PI = -1./dt * np.diff(Kt)
    
    return PI
