# -*- coding: utf-8 -*-
"""
Created on Wed Feb 20 10:16:22 2019

@author: Joshua
"""
### THIS CODE WAS WRITTEN TO SIMULATE 2D-HEAT TRANSFER ON A mxn GRID WHERE THE TEMPERATURES
### OF THE ELMEENTS ON THE EDEGE ARE (ASSUMED TO BE) FIXED. EACH ELEMENT CORRESPONDS TO ONE SUBDIVISION OF THE
### AREA.


import numpy as np
import math as mt
import matplotlib.pyplot as plt

def node_temp(W:float,H:float,th:float,Ta:float,h:float,hsink:float,k:float,m:int,n:int):
    
    #Tc contact temperature
    #Ta ambien temperature
    #W width of the plate
    #H height of the plate
    #m rows
    #n columns
    #hsink convection coefficient (overall) for blocks with heat sink
    #th thickness
    
    ##Since we didn't take this measurements, we will assume that the temperature of the block (at the surface, for convection) without
    #heat sink is 90C and that of the other two 70C.
    ##Note: in this simulation, thickness is set to 1 following the convention of the textbooks, this reduces the influence of h*As in the equations because As << Axy. Thus, convention might have been neglected. Heat transfer through the bottom face of the board was also neglected##
    
    T1=65.3 #Since we didn't take this measurements, we will assume that the temperature of the block (at the surface, for convection) without
    #heat sink is 85C and that of the other two 70C
    T2=80
    T3=60.2
    V=23.82 #Voltage
    R=85 #All blocks had R's between 84-86
    Qgen=V*V/R #Heat entering to the block
    Afin=(27.9**2 + 56*(1.4**2 + 4*1.4*15.2))/1000000
    ABh=((40*40-27.9**2)+ 4*40*20)/1000000 #Area of the block in contact with air
    Af=np.full((m,n),1) #Fraction of area exposed to air of most elements

    B1=np.zeros((m,n)) #Arrays containing one if there's a block on this element
    B1[2,5]=1
    B1[3,5]=1
    B1[4,5]=1
    B1[4,6]=1
    B1[4,7]=1
    B1[3,7]=1
    B1[2,7]=1
    B1[2,6]=1
    B1[3,6]=1
            
    Af[2,5]=3/4 #Array containing fraction area of the element exposed to air
    Af[3,5]=1/2
    Af[4,5]=3/4
    Af[4,6]=1/2
    Af[4,7]=3/4
    Af[3,7]=1/2
    Af[2,7]=3/4
    Af[2,6]=1/2
    Af[3,6]=0
    
    B2=np.zeros((m,n))
    Af[5,2]=3/4
    Af[6,2]=1/2
    Af[7,2]=3/4
    Af[7,3]=1/2
    Af[7,4]=3/4
    Af[6,4]=1/2
    Af[5,4]=3/4
    Af[5,3]=1/2
    Af[6,3]=0
    
    B2[5,2]=1
    B2[6,2]=1
    B2[7,2]=1
    B2[7,3]=1
    B2[7,4]=1
    B2[6,4]=1
    B2[5,4]=1
    B2[5,3]=1
    B2[6,3]=1
        
    B3=np.zeros((m,n))
    B3[8,6]=1
    B3[8,7]=1
    B3[8,8]=1
    B3[7,8]=1
    B3[6,8]=1
    B3[6,7]=1
    B3[6,6]=1
    B3[7,6]=1
    B3[7,7]=1
        
    Af[8,6]=3/4
    Af[8,7]=1/2
    Af[8,8]=3/4
    Af[7,8]=1/2
    Af[6,8]=3/4
    Af[6,7]=1/2
    Af[6,6]=3/4
    Af[7,6]=1/2
    Af[7,7]=0   
    ################################################################################
    hB1=hsink*B1 #USE THIS LINES TO CHOOSE WHAT BLOCKS HAVE A HEAT SINK
    hB2=h*B2
    hB3=hsink*B3
    ####################################################################################
    hs=hB1+hB2+hB3
    TB1=T1*B1
    TB2=T2*B2
    TB3=T3*B3
    TB=TB1+TB2+TB3
        
    T=np.zeros((m,n))
    Ti=np.zeros((m,n)) #Stores previous iteration
    its=0 #Number of iterations
    epsi=1 #Error
    dx=W/n
    dy=H/m
    As=dx*dy #Surface area of one side of one element
    Ay=dx*th
    Ax=dy*th

    
    #Sets values from experiment, aka boundary conditions
    T[0,0]=25.5
    T[0,1]=26
    T[0,2]=27
    T[0,3]=28.4
    T[0,4]=29.8
    T[0,5]=30.8
    T[0,6]=31.7
    T[0,7]=31
    T[0,8]=28.8
    T[0,9]=27.3
    T[0,10]=27.3
    
    T[1,0]=25.7
    T[2,0]=26.5
    T[3,0]=27.5
    T[4,0]=29.7
    T[5,0]=31.5
    T[6,0]=32.8
    T[7,0]=33.5
    T[8,0]=31.4
    T[9,0]=29.9
    T[10,0]=28.9
    
    T[m-1,1]=29.5
    T[m-1,2]=30.5
    T[m-1,3]=32.1
    T[m-1,4]=33.1
    T[m-1,5]=34
    T[m-1,6]=34.3
    T[m-1,7]=35.7
    T[m-1,8]=34.3
    T[m-1,9]=32.8
    T[m-1,10]=29.85
    
    T[1,n-1]=28.1
    T[2,n-1]=29.8
    T[3,n-1]=31.3
    T[4,n-1]=32.5
    T[5,n-1]=34.1
    T[6,n-1]=35.5
    T[7,n-1]=36.1
    T[8,n-1]=33.5
    T[9,n-1]=31.7
    
    #Sets Qin term in appropriate nodes
    
    #Block1 with heat sink
           
    while epsi>0.00005:#Stops when error reaches below
        its=its+1
        print('iteration:',its)
        for i in range(m):
            for j in range(n):
                if i!=0 and j!=0 and i!=m-1 and j!=n-1:
                    Tu=T[i-1,j] #Tupper
                    Tl=T[i,j-1] #Tleft
                    Tr=T[i,j+1] #Tright
                    Td=T[i+1,j] #Tdown
                    if hs[i,j]==hsink:
                        qin=(Qgen +(hs[i,j]*Afin+h*ABh)*(Ta-TB[i,j]))/(0.04*0.04) #Heat flux going into the circuit board
                    else:
                        qin=(Qgen+h*(0.04*0.04 + 0.02*0.04*4)*(Ta-TB[i,j]))/(0.04*0.04)
                    Qin=(1-Af[i,j])*As*qin
                    Ti[i,j]=T[i,j]
                    T[i,j]=(2*h*As*Af[i,j]*Ta + k*Ay/dy*(Td+Tu) + k*Ax/dx*(Tl+Tr) + Qin)/(2*h*As*Af[i,j] + 2*k*Ax/dx + 2*k*Ay/dy)     
                    
                    T[0,0]=25.5
                    T[0,1]=26
                    T[0,2]=27
                    T[0,3]=28.4
                    T[0,4]=29.8
                    T[0,5]=30.8
                    T[0,6]=31.7
                    T[0,7]=31
                    T[0,8]=28.8
                    T[0,9]=27.3
                    T[0,10]=27.3
                    
                    T[1,0]=25.7
                    T[2,0]=26.5
                    T[3,0]=27.5
                    T[4,0]=29.7
                    T[5,0]=31.5
                    T[6,0]=32.8
                    T[7,0]=33.5
                    T[8,0]=31.4
                    T[9,0]=29.9
                    T[10,0]=28.9
                    
                    T[m-1,1]=29.5
                    T[m-1,2]=30.5
                    T[m-1,3]=32.1
                    T[m-1,4]=33.1
                    T[m-1,5]=34
                    T[m-1,6]=34.3
                    T[m-1,7]=35.7
                    T[m-1,8]=34.3
                    T[m-1,9]=32.8
                    T[m-1,10]=29.85
                    
                    T[1,n-1]=28.1
                    T[2,n-1]=29.8
                    T[3,n-1]=31.3
                    T[4,n-1]=32.5
                    T[5,n-1]=34.1
                    T[6,n-1]=35.5
                    T[7,n-1]=36.1
                    T[8,n-1]=33.5
                    T[9,n-1]=31.7
                    
                    if T[i,j] > 10000: #Checks that values do not become too large or too small
                        print('Values became too high')
                        break
                    elif T[i,j] <-10000:
                        print('Values became too high')
                        break
        Tnew=T
        diffs=[]
        for i in range(m):
            for j in range(n):
                 if i!=0 and j!=0 and i!=m-1 and j!=n-1:
                     diffs.append(abs((Tnew[i,j]-Ti[i,j])/Tnew[i,j])) #Error on each element
        #Print('diffs list',diffs) ---- Used to check
        epsi=max(diffs) #Convergence based on maximum error       
    print('error:',epsi)
    print('iteration:',its)
    return T
############################ Function used eto fit h to data ###################
    
def node_temp_fit(W:float,H:float,th:float,Ta:float,h:float,hsink:float,k:float,m:int,n:int):
    
    #Tc contact temperature
    #Ta ambien temperature
    #W width of the plate
    #H height of the plate
    #m rows
    #n columns
    #hsink convection coefficient (overall) for blocks with heat sink
    #th thickness
    
    ##Since we didn't take this measurements, we will assume that the temperature of the block (at the surface, for convection) without
    #heat sink is 90C and that of the other two 70C.
    ##Note: in this simulation, thickness is set to 1 following the convention of the textbooks, this reduces the influence of h*As in the equations because As << Axy. Thus, convention might have been neglected. Heat transfer through the bottom face of the board was also neglected##
    
    T1=65.3 #Since we didn't take this measurements, we will assume that the temperature of the block (at the surface, for convection) without
    #heat sink is 85C and that of the other two 70C
    T2=60.2
    T3=80
    V=23.82 #Voltage
    R=85 #All blocks had R's between 84-86
    Qgen=V*V/R #Heat entering to the block
    Afin=(27.9**2 + 56*(1.4**2 + 4*1.4*15.2))/1000000
    ABh=((40*40-27.9**2)+ 4*40*20)/1000000 #Area of the block in contact with air
    Af=np.full((m,n),1) #Fraction of area exposed to air of most elements

    B1=np.zeros((m,n)) #Arrays containing one if there's a block on this element
    B1[2,5]=1
    B1[3,5]=1
    B1[4,5]=1
    B1[4,6]=1
    B1[4,7]=1
    B1[3,7]=1
    B1[2,7]=1
    B1[2,6]=1
    B1[3,6]=1
            
    Af[2,5]=3/4 #Array containing fraction area of the element exposed to air
    Af[3,5]=1/2
    Af[4,5]=3/4
    Af[4,6]=1/2
    Af[4,7]=3/4
    Af[3,7]=1/2
    Af[2,7]=3/4
    Af[2,6]=1/2
    Af[3,6]=0
    
    B2=np.zeros((m,n))
    Af[5,2]=3/4
    Af[6,2]=1/2
    Af[7,2]=3/4
    Af[7,3]=1/2
    Af[7,4]=3/4
    Af[6,4]=1/2
    Af[5,4]=3/4
    Af[5,3]=1/2
    Af[6,3]=0
    
    B2[5,2]=1
    B2[6,2]=1
    B2[7,2]=1
    B2[7,3]=1
    B2[7,4]=1
    B2[6,4]=1
    B2[5,4]=1
    B2[5,3]=1
    B2[6,3]=1
        
    B3=np.zeros((m,n))
    B3[8,6]=1
    B3[8,7]=1
    B3[8,8]=1
    B3[7,8]=1
    B3[6,8]=1
    B3[6,7]=1
    B3[6,6]=1
    B3[7,6]=1
    B3[7,7]=1
        
    Af[8,6]=3/4
    Af[8,7]=1/2
    Af[8,8]=3/4
    Af[7,8]=1/2
    Af[6,8]=3/4
    Af[6,7]=1/2
    Af[6,6]=3/4
    Af[7,6]=1/2
    Af[7,7]=0   
    ################################################################################
    hB1=hsink*B1 #USE THIS LINES TO CHOOSE WHAT BLOCKS HAVE A HEAT SINK
    hB2=hsink*B2
    hB3=h*B3
    ####################################################################################
    hs=hB1+hB2+hB3
    TB1=T1*B1
    TB2=T2*B2
    TB3=T3*B3
    TB=TB1+TB2+TB3
        
    T=np.zeros((m,n))
    Ti=np.zeros((m,n)) #Stores previous iteration
    its=0 #Number of iterations
    epsi=1 #Error
    dx=W/n
    dy=H/m
    As=dx*dy #Surface area of one side of one element
    Ay=dx*th
    Ax=dy*th

    
    #Sets values from experiment, aka boundary conditions
    T[0,0]=25.5
    T[0,1]=26
    T[0,2]=27
    T[0,3]=28.4
    T[0,4]=29.8
    T[0,5]=30.8
    T[0,6]=31.7
    T[0,7]=31
    T[0,8]=28.8
    T[0,9]=27.3
    T[0,10]=27.3
    
    T[1,0]=25.7
    T[2,0]=26.5
    T[3,0]=27.5
    T[4,0]=29.7
    T[5,0]=31.5
    T[6,0]=32.8
    T[7,0]=33.5
    T[8,0]=31.4
    T[9,0]=29.9
    T[10,0]=28.9
    
    T[m-1,1]=29.5
    T[m-1,2]=30.5
    T[m-1,3]=32.1
    T[m-1,4]=33.1
    T[m-1,5]=34
    T[m-1,6]=34.3
    T[m-1,7]=35.7
    T[m-1,8]=34.3
    T[m-1,9]=32.8
    T[m-1,10]=29.85
    
    T[1,n-1]=28.1
    T[2,n-1]=29.8
    T[3,n-1]=31.3
    T[4,n-1]=32.5
    T[5,n-1]=34.1
    T[6,n-1]=35.5
    T[7,n-1]=36.1
    T[8,n-1]=33.5
    T[9,n-1]=31.7
    
    #Sets Qin term in appropriate nodes
    
    #Block1 with heat sink
           
    while epsi>0.00005:#Stops when error reaches below
        its=its+1
        print('iteration:',its)
        for i in range(m):
            for j in range(n):
                if i!=0 and j!=0 and i!=m-1 and j!=n-1:
                    Tu=T[i-1,j] #Tupper
                    Tl=T[i,j-1] #Tleft
                    Tr=T[i,j+1] #Tright
                    Td=T[i+1,j] #Tdown
                    if hs[i,j]==hsink:
                        qin=(Qgen +(hs[i,j]*Afin+h*ABh)*(Ta-TB[i,j]))/(0.04*0.04) #Heat flux going into the circuit board
                    else:
                        qin=(Qgen+h*(0.04*0.04 + 0.02*0.04*4)*(Ta-TB[i,j]))/(0.04*0.04)
                    Qin=(1-Af[i,j])*As*qin
                    Ti[i,j]=T[i,j]
                    T[i,j]=(2*h*As*Af[i,j]*Ta + k*Ay/dy*(Td+Tu) + k*Ax/dx*(Tl+Tr) + Qin)/(2*h*As*Af[i,j] + 2*k*Ax/dx + 2*k*Ay/dy)     
                    
                    T[0,0]=25.5
                    T[0,1]=26
                    T[0,2]=27
                    T[0,3]=28.4
                    T[0,4]=29.8
                    T[0,5]=30.8
                    T[0,6]=31.7
                    T[0,7]=31
                    T[0,8]=28.8
                    T[0,9]=27.3
                    T[0,10]=27.3
                    
                    T[1,0]=25.7
                    T[2,0]=26.5
                    T[3,0]=27.5
                    T[4,0]=29.7
                    T[5,0]=31.5
                    T[6,0]=32.8
                    T[7,0]=33.5
                    T[8,0]=31.4
                    T[9,0]=29.9
                    T[10,0]=28.9
                    
                    T[m-1,1]=29.5
                    T[m-1,2]=30.5
                    T[m-1,3]=32.1
                    T[m-1,4]=33.1
                    T[m-1,5]=34
                    T[m-1,6]=34.3
                    T[m-1,7]=35.7
                    T[m-1,8]=34.3
                    T[m-1,9]=32.8
                    T[m-1,10]=29.85
                    
                    T[1,n-1]=28.1
                    T[2,n-1]=29.8
                    T[3,n-1]=31.3
                    T[4,n-1]=32.5
                    T[5,n-1]=34.1
                    T[6,n-1]=35.5
                    T[7,n-1]=36.1
                    T[8,n-1]=33.5
                    T[9,n-1]=31.7
                    
                    if T[i,j] > 10000: #Checks that values do not become too large or too small
                        print('Values became too high')
                        break
                    elif T[i,j] <-10000:
                        print('Values became too high')
                        break
        Tnew=T
        diffs=[]
        for i in range(m):
            for j in range(n):
                 if i!=0 and j!=0 and i!=m-1 and j!=n-1:
                     diffs.append(abs((Tnew[i,j]-Ti[i,j])/Tnew[i,j])) #Error on each element
        #Print('diffs list',diffs) ---- Used to check
        epsi=max(diffs) #Convergence based on maximum error       
    print('error:',epsi)
    print('iteration:',its)
    return T
           
#################################END OF FUNCTION######################################
    
def l_sq_method(A:np.ndarray,B:np.ndarray): 
    #Finds the RMS. Relies on both arrays being the same size. A is the sample data
    m=np.size(A[:,0])
    n=np.size(A[0])
    l_sqs=0 #Initiates sum of squares
    count=0 #Initializes count of valid sums
    for i in range(m):
        for j in range(n):
            if A[i,j]!=0: #Only takes into account data not equal to zero. This is because we do not have the temperatures across the entire board
                count=count+1
                l_sqs=l_sqs+((A[i,j]-B[i,j])*(A[i,j]-B[i,j])) #Adds square of difference of each element
    output=(mt.sqrt(l_sqs))/count #divides by number of elements and takes sqrt
    print('sqs',output)
    return output #Returns the RMS

##############################END OF FUNCTION####################################

def best_fit(hlow:float, hhigh:float, data:np.ndarray): 
    #Finds best-fit h from a list. 
    #The user may need to iterate based on the result. This could be automated.
    
    hrange=np.linspace(hlow,hhigh,200) #Creates a list with 200 elements. Can be changed for more accuracy.
    squares=[]
    for i in hrange:
        trial=node_temp_fit(W1,H1,th1,Ta1,h1,i,k1,m1,n1) #Used to compare with data
        sqs=l_sq_method(data,trial)
        squares.append(sqs)
    win=hrange[np.argmin(squares)] #Finds the h that generated the closest value
    print('h_best_fit:',win)
    
    #Plots and saves the image
    plt.figure(1)
    plt.plot(hrange,squares)
    plt.xticks(range(0,51,5))
    plt.plot(win,min(squares),'ro')
    plt.text(win+0.1,min(squares)+0.01,'Min')
    plt.ylabel('Root Mean Square Error')
    plt.xlabel('h [W/m^2.K]')
    #plt.savefig('RMS vs h.png',dpi=200)
    return win
######################### END OF FUNCTION ######################### 

#Important to define these before using best fit function
W1=0.220
H1=0.220
m1=11#Rows
n1=11#Columns
k1=8.730 #Calculated K
Ta1=20 #Ambient temperature
th1=0.00157 #From textbook Section 5-4, th=1 for 2-D simulation
h1=10
hsink1=10

hlow1=15
hhigh1=20

T=np.zeros((m1,n1))
T[2,2]=46
T[3,2]=51.3
T[4,2]=51.7
T[4,3]=59.4
T[4,4]=52
T[3,4]=52
T[2,4]=47.4
T[2,3]=51.3

T[5,2]=49.6
T[6,2]=54.5
T[7,2]=53.9
T[7,3]=56.4
T[7,4]=53.4
T[6,4]=62
T[5,4]=56.8
T[5,3]=55.2

T[8,6]=57.7
T[8,7]=60.7
T[8,8]=53.1
T[7,8]=59
T[6,8]=55
T[6,7]=62.4
T[6,6]=60.5
T[7,6]=62.5  

T[9,1]=32.2
T[2,8]=34

ExpData=T

hBestFit=best_fit(hlow1,hhigh1,ExpData) #Important to define parameters for node_temp function before running best fit

hBestFit=19.773
#RMS=6.636

M1=node_temp(W1,H1,th1,Ta1,h1,hBestFit,k1,m1,n1)
RMS=l_sq_method(ExpData,M1)

plt.figure(2)
plt.contourf(np.flip(M1,0))
plt.colorbar()
plt.show()
plt.savefig('Heatdistlab3-1.png',dpi=200)

print('Resulting Matrix of Temperatures')
print(M1)
print('MeanT=',np.mean(M1))
print('MaxT=',np.amax(M1))
print('Error=',RMS,'Deg')
