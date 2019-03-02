"""M345SC Homework 2, part 2
Name: Ziyao Zhang
CID: 01181949
"""

def model1(G,x=0,params=(50,80,105,71,1,0),tf=6,Nt=400,display=False):
    """
    Question 2.1
    Simulate model with tau=0.

    Input:
    G: Networkx graph.
    x: node which is initially infected.
    params: Model parameters.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf.
    display: A plot of S(t) for the infected node is generated when true.

    Output:
    S: Array containing S(t) for infected node.
    """
    def RHS(y,t):
        # parameter setup
        a,theta0,theta1,g,k,tau = params
        theta = theta0 + theta1*(1-np.sin(2*np.pi*t))

        # system of ODEs
        S,I,V = y
        St = a*I - (g+k)*S
        It = theta*S*V - (k+a)*I
        Vt = k*(1-V) - theta*S*V
        return St,It,Vt

    # initial condition
    S0,I0,V0 = 0.05,0.05,0.1

    # solve
    tarray = np.linspace(0,tf,Nt+1)
    y = odeint(RHS,[S0,I0,V0],tarray)
    S = y[:,0]

    # plot S vs time
    if display:
        plt.figure()
        plt.plot(tarray,S)
        plt.xlabel(r'Time $t$')
        plt.ylabel(r'Fraction $S(t)$')
        plt.title('Fraction of Spreaders (by Ziyao)')

    return S



def modelN(G,x=0,params=(50,80,105,71,1,0.01),tf=6,Nt=400,display=False):
    """
    Question 2.2
    Simulate model with tau=0.

    Input:
    G: Networkx graph.
    x: node which is initially infected.
    params: Model parameters.
    tf,Nt: Solutions Nt time steps from t=0 to t=tf.
    display: A plot of S(t) for the infected node is generated when true.

    Output:
    Smean,Svar: Array containing mean and variance of S across network nodes at
                each time step.
    Imean,Ivar: Array containing mean and variance of I across network nodes at
                each time step.
    Vmean,Vvar: Array containing mean and variance of V across network nodes at
                each time step.
    """
    # flux F transpose
    def fluxT(A,tau):
        q = A.sum(axis=0)
        num = np.multiply(A,q)
        denom = num.sum(axis=1)
        return tau*np.divide(num,denom)

    F = fluxT(nx.adjacency_matrix(G).todense(),params[-1])
    N = nx.number_of_nodes(G)

    def RHS(y,t):
        """
        Compute RHS of model at time t.

        Input:
        y: A 3N x 1 array with y[:N],y[N:2*N],y[2*N:3*N] corresponding to S on
           nodes 0 to N-1, I on nodes 0 to N-1, and V on nodes 0 to N-1
           respectively.

        Output:
        dy: A 3N x 1 array corresponding to dy/dt.

        Discussion:
        The cost of RHS function is as follows:
        1) To assign values for S,I,V from the input, initialise and later
           assign values for the output dy, it requires A0*N operations in total
           where A0 is some constant.
        2) To compute dS/dt, the calculation involves four types of array matrix
           arithmetics - simple addition and subtraction between arrays,
           multipication of array and scalar number, elementwise multipication
           of arrays and matrix multiplication, which have computional costs of
           A1*N,A2*N,A3*N and B*N^2 respectively, where Ai and B are some
           constants.
        3) A few basic number arithmetics (eg, g+k+tau, k+a+tau, k+tau) brings
           a minor cost of C operations, where C is some constant.

        In conclusion, each call to function RHS requires B*N^2 + A*N + C
        operations.
        """
        # parameter setup
        a,theta0,theta1,g,k,tau = params
        theta = theta0 + theta1*(1-np.sin(2*np.pi*t))

        # system of ODEs
        S,I,V = y[0:N],y[N:2*N],y[2*N:3*N]
        dy = np.zeros(3*N)
        dy[0:N] = a*I - (g+k+tau)*S + S*F
        dy[N:2*N] = theta*S*V - (k+a+tau)*I + I*F
        dy[2*N:3*N] = k - (k+tau)*V - theta*S*V + V*F

        return dy

    # initial condition
    init = [0 for n in range(2*N)] + [1 for n in range(N)]
    init[x],init[x+N],init[x+2*N] = 0.05,0.05,0.1

    # solve
    tarray = np.linspace(0,tf,Nt+1)
    y = odeint(RHS,init,tarray)

    S,I,V = y[:,0:N],y[:,N:2*N],y[:,2*N:3*N]
    Smean = S.mean(1)
    Svar = S.var(1)
    Imean = I.mean(1)
    Ivar = I.var(1)
    Vmean = V.mean(1)
    Vvar = V.var(1)


    if display:
        plt.figure()
        plt.plot(tarray,Smean)
        plt.xlabel(r'Time $t$')
        plt.ylabel(r'Average Fraction $<S(t)>$')
        plt.title('Average Fraction of Spreaders (by Ziyao)')

        plt.figure()
        plt.plot(tarray,Svar)
        plt.xlabel(r'Time $t$')
        plt.ylabel(r'Variance of Fraction $<(S(t)−<S(t)>)^2>$')
        plt.title('Variance of Fraction of Spreaders (by Ziyao)')

    return Smean,Svar,Imean,Ivar,Vmean,Vvar



def diffusion(G,D,x=0,tf=6,Nt=400):
    """Analyze similarities and differences between simplified infection model
    and linear diffusion on Barabasi-Albert networks.

    Discussion:
        In the following context, S refers to Spreaders, I refers to Infected
    Cells and V refers to Vulnerable Cells.

    Fig1 & Fig2:
    The figure 1 suggests that the value of theta0 makes no difference on the
    dynamics of S in the infection model. This is as expected since S only
    transports between nodes it only depends on the flux (which is fixed over
    time).  The transportation based on linear diffusion performs the process
    much faster than the one based on flux, given the same transportation rate
    (tau=D). This is also understandable as in linear diffusion S simply moves
    from high density nodes to low density nodes, whereas in infection model the
    movement of S depends on the connectivity of nodes. Lastly tau controls the
    rate of transportation.

    Fig3 & Fig4:
    The first and third plots are strictly symmetric by x-axis. This makes sense
    because under current model conditions Vulnerable Cells are converted into
    Infected Cells only without recovering. The convertion is faster for larger
    theta0 as theta0 controls the rate of convertion.

    Similarity & Difference:
    In both infection model and linear diffusion, S tends to get evenly
    distributed over the nodes, which matches the nature of the two models.

    However, linear diffusion also averages I and V over the nodes, whereas in
    the infection model I and V tends to diverge and trend becomes faster as
    process goes on. Biologically speaking, Vulnerable Cells and Infected Cells
    do not convert into each other in the diffusion model and Vulnerable Cells
    keep getting infected in the infection model.
    """
    L = nx.laplacian_matrix(G).todense()
    N = nx.number_of_nodes(G)

    def RHS(y,t):
        # linear diffusion
        S,I,V = y[0:N],y[N:2*N],y[2*N:3*N]
        dy = np.zeros(3*N)
        dy[0:N] = -D*S*L
        dy[N:2*N] = -D*I*L
        dy[2*N:3*N] = -D*V*L

        return dy

    # initial condition
    init = [0 for n in range(2*N)] + [1 for n in range(N)]
    init[x],init[x+N],init[x+2*N] = 0.05,0.05,0.1

    # solve
    tarray = np.linspace(0,tf,Nt+1)
    y = odeint(RHS,init,tarray)

    S,I,V = y[:,0:N],y[:,N:2*N],y[:,2*N:3*N]
    Smean = S.mean(1)
    Svar = S.var(1)
    Imean = I.mean(1)
    Ivar = I.var(1)
    Vmean = V.mean(1)
    Vvar = V.var(1)

    return Smean,Svar,Imean,Ivar,Vmean,Vvar



if __name__ == '__main__':
    import numpy as np
    import networkx as nx
    from scipy.integrate import odeint
    import matplotlib.pyplot as plt
    from time import time

    # test
    Er = nx.erdos_renyi_graph(1000,0.2)
    model1(Er,display=True)
    print('launch modelN')
    t1 = time()
    modelN(Er,display=True)
    t2 = time()
    print('time=',t2-t1)

    # infection model vs linear diffusion
    G = nx.barabasi_albert_graph(100,5)

    # fix tau=D, change theta0
    tau = D = 0.05
    tf = 6
    Nt = 400

    dSmean,dSvar,dImean,dIvar,dVmean,dVvar = diffusion(G,D,x=0)
    tarray = np.linspace(0,tf,Nt+1)

    plt.figure(figsize=(12,10))
    plt.plot(tarray,dSvar,label='linear diffusion')
    for theta0 in (60,120,180):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Svar,label=r'infection $\theta_0$ = %.f' %theta0)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Variance of Fraction $<(S(t)−<S(t)>)^2>$')
    plt.title('Variance of Fraction of Spreaders (by Ziyao)')
    plt.legend()

    plt.figure(figsize=(14,14))
    plt.subplot(2,2,1)
    plt.plot(tarray,dImean,label='linear diffusion')
    for theta0 in (60,120,180):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Imean,label=r'infection $\theta_0$ = %.f' %theta0)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Average Fraction $<(I(t)>$')
    plt.title('Average Fraction of Infected Cells (by Ziyao)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(tarray,dIvar,label='linear diffusion')
    for theta0 in (60,120,180):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Ivar,label=r'infection $\theta_0$ = %.f' %theta0)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Variance of Fraction $<(I(t)−<I(t)>)^2>$')
    plt.title('Variance of Fraction of Infected Cells (by Ziyao)')
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(tarray,dVmean,label='linear diffusion')
    for theta0 in (60,120,180):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Vmean,label=r'infection $\theta_0$ = %.f' %theta0)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Average Fraction $<(V(t)>$')
    plt.title('Average Fraction of Vulnerable Cells (by Ziyao)')
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(tarray,dVvar,label='linear diffusion')
    for theta0 in (60,120,180):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Vvar,label=r'infection $\theta_0$ = %.f' %theta0)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Variance of Fraction $<(V(t)−<V(t)>)^2>$')
    plt.title('Variance of Fraction of Vulnerable Cells (by Ziyao)')
    plt.legend()

    # fix theta0=80, change tau
    tf = 6
    Nt = 400
    theta0 = 80

    plt.figure(figsize=(12,10))
    plt.plot(tarray,dSvar,label='linear diffusion')
    for tau in (0.05,0.1,0.2):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Svar,label=r'infection $\tau$ = %.2f' %tau)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Variance of Fraction $<(S(t)−<S(t)>)^2>$')
    plt.title('Variance of Fraction of Spreaders (by Ziyao)')
    plt.legend()

    plt.figure(figsize=(14,14))
    plt.subplot(2,2,1)
    plt.plot(tarray,dImean,label='linear diffusion')
    for tau in (0.05,0.1,0.2):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Imean,label=r'infection $\tau$ = %.2f' %tau)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Average Fraction $<(I(t)>$')
    plt.title('Average Fraction of Infected Cells (by Ziyao)')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(tarray,dIvar,label='linear diffusion')
    for tau in (0.05,0.1,0.2):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Ivar,label=r'infection $\tau$ = %.2f' %tau)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Variance of Fraction $<(I(t)−<I(t)>)^2>$')
    plt.title('Variance of Fraction of Infected Cells (by Ziyao)')
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(tarray,dVmean,label='linear diffusion')
    for tau in (0.05,0.1,0.2):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Vmean,label=r'infection $\tau$ = %.2f' %tau)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Average Fraction $<(V(t)>$')
    plt.title('Average Fraction of Vulnerable Cells (by Ziyao)')
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(tarray,dVvar,label='linear diffusion')
    for tau in (0.05,0.1,0.2):
        Smean,Svar,Imean,Ivar,Vmean,Vvar = modelN(G,0,(0,theta0,0,0,0,tau))
        plt.plot(tarray,Vvar,label=r'infection $\tau$ = %.2f' %tau)
    plt.xlabel(r'Time $t$')
    plt.ylabel(r'Variance of Fraction $<(V(t)−<V(t)>)^2>$')
    plt.title('Variance of Fraction of Vulnerable Cells (by Ziyao)')
    plt.legend()
