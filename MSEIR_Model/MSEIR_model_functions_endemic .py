### The Python script below was written and revised by Brandon Elford and Dr. Cedric Chauve to model the SIR model, and has been extended by Laura Gutierrez Funderburk to model the MSEIR model.

""" Importing Packages """

#Python packages to use, respectively the ODE solver function, Numpy package to handle arrays more easily, and the Plot package
from scipy.integrate import odeint 
import numpy as np                 
import matplotlib.pyplot as plt    

""" Setting up initial values """

## Rates of change from one class to the other for chickenpox
delta = 1/180
epsilon = 1/14
gamma = 1/7


## Set, respectively, rates of birth and death as 
b = 1/50
d = 1/60

## Recall that in this model, q  is defined as follows
q = b - d



#Initial infected population
I_0 = 200
R_0 = 200 
M_0 = 200
E_0 = 200
N_0 = 1000


#Initial susceptible population depends on the previous values
S_0 = N_0 - M_0 -  E_0 -I_0 - R_0 

#We let y_0 be the vector that holds the initial values for M,S, I,E, R
y_0 = [M_0,S_0, I_0, E_0, R_0,N_0]

#Defines the domain for t to be [0,160]
#and using 100 points to plot the curve
t = np.linspace(0, 700, 100)
maxt = 700

""" Function definition area """

def diffs(y,t,beta,delta,epsilon,gamma,b,d):
    M,S,E,I,R,N = y
    
    dM_dt = b*(N - S) - (delta + d)*M
    dS_dt = b*S + delta*M - beta*S*I/N - d*S
    dE_dt = beta*S*I/N - (epsilon + d)*E
    dI_dt = epsilon*E - (gamma  + d)*I
    dR_dt = gamma*I - d*R
    dN_dt = (b - d)*N
    
    return dM_dt,dS_dt,dE_dt,dI_dt,dR_dt,dN_dt

#Here we have defined a function plot_ODE that is used to plot our model
#Inputs - result from our ODE above for M,S, I,E,R 
#       - population N
#       - max time considered
#       - ticlen is the length between the tic marks on the x-axis
#       - fs is figure size when it prints out
#       - title will be the title of your plot
#Output - a visual representation of our model

def plot_ODE(M, S, E, I, R,N, maxt, ticlen, fs,title):
    ###This next block allows us to plot the results
    fig = plt.figure(facecolor='w', figsize= fs)
    ax = fig.add_subplot(111, facecolor='#dddddd')
    
    #This command sets the title of your graph
    ax.set_title(title)
    
    #Plots m(t) = M(t)/N over the domain of t and labels it 'Maternal Immunity'
    #'b' = colour line blue
    ax.plot(t, M, 'm',label='Maternal Immunity')
    #Plots s(t) = S(t)/N over the domain of t and labels it 'Susceptible'
    #'b' = colour line blue
    ax.plot(t, S, 'y',label='Susceptible')
    #Plots e(t) = E(t)/N over the domain of t and labels it 'Exposed'
    #'b' = colour line blue
    ax.plot(t, E, 'b',label='Exposed')
    #Plots i(t) = I(t)/N and labels it 'Infected'
    #'r' = colour line red
    ax.plot(t, I, 'r', label='Infected')
    #Plots r(t) = R(t)/N and labels is 'Recovered'
    #'g' = colour line green
    ax.plot(t, R, 'g', label='Recovered')

    #Labels the x axis
    h = ax.set_xlabel('Time /days')
    #Defined the number of tic's on the x-axis
    plt.xticks(np.arange(0,maxt+1,ticlen), rotation=270)
    
    #Labels the y axis
    ax.set_ylabel('Total Population')
    
    #sets the range of the y values
    ax.set_ylim(-0.1,1000)

    #Toggles the grid lines
    ax.grid(color = 'w', linestyle = '-')

    #Toggles legend, again alpha represents the transparency
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)

    #Toggles a border with 'True' and 'False'
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    
    #Prints the plot to your screen
    plt.show()


#This function takes M_0,E_0,I_0, R_0,N_0 as inputs and repeats the functions we used 
#above over some number of initial conditions with the same model
def phase_plane(M_0,E_0,I_0,R_0,N_0):
    
    #This is a check that both of your inputs are the same length
    #If they are not (one has too many or too little data points) then it will return an error message
    #if len(I_0) == len(R_0):
    start = 0
    end = len(I_0)
    #else:
    #    return "Input error"
    
    #Initial susceptible population depends on the previous values
    S_0 = [N_0[i] - M_0[i] - E_0[i] - I_0[i] - R_0[i] for i in range(start,end)]

    #We let y_0 be the vector that holds the initial values for M,S,E,I and R
    y_0 = []
    for i in range(start,end):
        y_0.append([M_0[i], S_0[i],E_0[i], I_0[i], R_0[i],N_0[i]])
    
    #This command solves the ODE using the differentials we defined above from diffs, giving the initial vector y_0 and domain of t
    #The "args" command inputs the arguments of N, beta, and gamma to the ODE solver
    #We let the solution be equal to "ret"
    reti = [odeint(diffs, y_0[i], t, args=(beta,delta,gamma,epsilon,b,d)) for i in range(start,end)]
    
    #ret is a 6 dimensional array that holds the results of all variables M(t),S(t),E(t),I(t) and R(t) 
    # for the given values of t
    #The ".T" after "ret" is a numpy command which gives us the transpose of the original array
    #We initialize three arrays of size 13 (Si,Ii,Ri) which will soon hold the data points we plot for our phase plane
    Mi = [0] * end
    Si = [0] * end
    Ei = [0] * end
    Ii = [0] * end
    Ri = [0] * end
    Ni = [0] * end
    for i in range(start,end):
        Mi[i],Si[i],Ei[i],Ii[i],Ri[i],Ni[i] = reti[i].T
    
    #after the calculations we will return Mi,Si, Ii, and the start and end length for the number of points
    return Mi, Si, Ei, Ii, Ri, Ni, start, end


#This function plots all of the resulting phase lines that we found with the function "Phase_plane"
def phase_plot(Si, Ii,Ni,start, end, title):

    #These commands are for formatting the resulting plot
    #figsize = (n,m) tells is to print out a picture with those same dimensions
    fig2 = plt.figure(facecolor='w', figsize=(8,8))
    ax2 = fig2.add_subplot(111, facecolor='#dddddd')

    #using the "start" and "end" from earlier, we plot each pair of susceptibles and infected 
    for i in range(start, end):
        ax2.plot(Si[i]/Ni[i], Ii[i]/Ni[i], 'black')

    #plotting the line from x=1,y=0 to x=0,y=1
    ax2.plot([0,1],[1,0],'black')

    #labelling the axis
    ax2.set_xlabel("s(t)")
    ax2.set_ylabel("i(t)")

    #setting the max and min values showin on the plots
    ax2.set_ylim(0,1)
    ax2.set_xlim(0,1)

    #Setting the title of the plot
    ax2.set_title(title)

    #More formatting for the plot
    ax2.grid(color='w', linestyle='-')

    #prints the plot to the screen
    plt.show()
   
    

""" Plotting input """
ticlen = 10
def equilibrium_choice(beta, title,ticlen):
    ret = odeint(diffs, y_0, t, args=(beta,delta,epsilon,gamma,b,d))
    M,S,E,I,R,N = ret.T
    plot_ODE(M,S,I,E,R,N,maxt,ticlen,(12,8),title)
