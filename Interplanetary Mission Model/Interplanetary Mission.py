# -*- coding: utf-8 -*-
"""
Created on Sun Oct 16 11:25:07 2016

@author: Aman
"""

import math
import numpy

AU=1.49597871e8
mewSun=1.327e11

mew_list = [2.2032e4 , 3.24859e5 , 3.986e5 , 4.282837e4 , 1.26686534e8 , 3.7931187e7 , 5.793939e6 , 6.836529e6 , 8.71e2]
radius_list = [2439.7 , 6051.8 , 6378.1 , 3396.2 , 71492 , 60268 , 25559 , 24764 , 1195]


y = []
m = []
d = []
h = []
mi = []
sec = [] 



#Function to input the time
def InputTime():
    
    year = input('year = ')
    y.append(float(year))

    month = input('month = ')
    m.append(float(month))

    day = input('day = ')
    d.append(float(day))
    
    hours = input('hours = ')
    h.append(float(hours))
    
    minutes = input('minutes = ')
    mi.append(float(minutes))
    
    seconds = input('seconds = ')
    sec.append(float(seconds))
    


#year  - range: 1901 - 2099
#month - range: 1 - 12
#day   - range: 1 - 31
#ut    - universal time
#        hour (0 - 23)
#        minute (0 - 60)
#        second (0 - 60)

#converting a date with universal time into the julian format
def julian(y,m,d,h,mi,sec):
    UT= h+(mi/60)+(sec/3600)      
    return ((367*y)-math.floor((7*(y+math.floor((m+9)/12)))/4)+math.floor((275*m)/9)+d+1721013.5+(UT/24))

print ('Enter the departure time (the time at which the spacecraft exits the sphere of influence of the planet)')
InputTime()
print ('Enter the arrival time (the time at which the spacecraft enters the sphere of influence of the (other) planet)')
InputTime()


JD1 = julian(y[0],m[0],d[0],h[0],mi[0],sec[0])
JD2 = julian(y[1],m[1],d[1],h[1],mi[1],sec[1])
print()
print('The corresponding Julian Days for the times you have entered are:')    
print('Departure time : ', JD1)
print('Arrival time : ', JD2)
print()
print ('Enter the codes for the planets between which you want the interplanetary trajectory.')
print ('0:Mercury')
print ('1:Venus')
print ('2:Earth')
print ('3:Mars')
print ('4:Jupiter')
print ('5:Saturn')
print ('6:Uranus')
print ('7:Neptune')
print ('8:Pluto')

n1 = input('code1 = ')
n1 = int(n1)
n2 = input('code2 = ')
n2 = int(n2)
print()


#Function to calculate the magnitude of a vector
def mag(v):
    return(math.sqrt((v[0]**2)+(v[1]**2)+(v[2]**2)))

#Solving for eccentric anomaly from kepler's equation
def Newton(M,e):

    if(M>math.pi):
        E=M-(e/2)
    else:
        E=M+(e/2)    

    f=E-e*math.sin(E)-M
    fD=1-e*math.cos(E)
    ratio=f/fD
 
    while(ratio>10**-8):
        E=E-ratio
        f=E-e*math.sin(E)-M
        fD=1-e*math.cos(E)
        ratio=f/fD
    
    return(E)



#Orbital Elements 
eleName = ['a' , 'e' , 'incl' , 'RA' , 'w_hat' , 'L']
#a = semi-major axis                            (AU)
#e = eccentricity
#incl = angle of inclination                    (degrees)
#RA = Right Ascension                           (degrees)
#w_hat = longitude of perihelion ( = RA + w)    (degrees)
#L = mean longitude ( = w_hat + M)              (degrees)


elements = [[[0.38709893,0.00000066],[0.72333199,0.00000092],[1.00000011,-0.00000005],[1.52366231,-0.00007221],[5.20336301,0.00060737],[9.53707032,-0.00301530],[19.19126393,0.00152025],[30.06896348,-0.00125196],[39.48168677,-0.00076912]], 
[[0.20563069,0.00002527],[0.00677323,-0.00004938],[0.01671022,-0.00003804],[0.09341233,0.00011902],[0.04839266,-0.00012880],[0.05415060,-0.00036762],[0.04716771,-0.00019150],[0.00858587,0.00002514],[0.24880766,0.00006465]],
[[7.00487,-(23.51/3600)],[3.39471,-(2.86/3600)],[0.00005,-(46.94/3600)],[1.85061,-(25.47/3600)],[1.30530,-(4.15/3600)],[2.48446,6.11/3600],[0.76986,-(2.09/3600)],[1.76917,-(3.64/3600)],[17.14175,11.07/3600]],
[[48.33167,-(446.30/3600)],[76.68069,-(996.89/3600)],[-11.26064,-(18228.25/3600)],[49.57854,-(1020.19/3600)],[100.55615,1217.17/3600],[113.71504,-(1591.05/3600)],[74.22988,-(1681.4/3600)],[131.72169,-(151.25/3600)],[110.30347,-(37.33/3600)]],
[[77.4545,573.57/3600],[131.53298,-(108.80/3600)],[102.94719,(1198.28/3600)],[336.04084,(1560.78/3600)],[14.75385,839.93/3600],[92.43194,-(1948.89/3600)],[170.96424,1312.56/3600],[44.97135,-(844.43/3600)],[224.06676,-(132.25/3600)]],
[[252.25084,(538101628.29/3600)],[181.97973,(210664136.06/3600)],[100.46435,(129597740.63/3600)],[355.45332,(68905103.78/3600)],[252.25084,538101628.29/3600],[181.97973,210664136.06/3600],[100.46435,129597740.63/3600],[355.45332,68905103.78/3600],[34.40438,10925078.35],[49.94432,4401052.95/3600],[313.23218,1542547.79/3600],[304.88003,786449.21/3600],[238.92881,522747.90/3600]]]   




#A function that returns the state vectors of the planet 
def elemUpdate(elem,T0,n):
    k=0
    while (k<6):
        elem.append(elements[k][n][0] + T0*elements[k][n][1])
        k=k+1
        
    if(elem[2]<0):
        elem[2]=abs(elem[2])
    if(elem[2]>180):
        elem[2]=elem[2]%180
    k=3
    while(k<6):
        if (elem[k]<0):
            elem[k]=360+elem[k]
        if(elem[k]>360):
            elem[k]=elem[k]%360
        k=k+1
    
    
    h=math.sqrt(mewSun*elem[0]*AU*(1-(elem[1]**2)))
    w=elem[4]-elem[3]
    M=elem[5]-elem[4]
    E=math.degrees(Newton(math.radians(M),elem[1]))
    if(E<0):
        E=360+E
    kP=math.sqrt((1+elem[1])/(1-elem[1]))
    theta=(numpy.arctan(kP*math.tan(math.radians(E/2))))*2
    
    r_hat=[((h**2)/(mewSun*(1+elem[1]*math.cos(theta))))*math.cos(theta) , ((h**2)/(mewSun*(1+elem[1]*math.cos(theta))))*math.sin(theta) , 0]
    
    v_hat=[(mewSun/h)*(-math.sin(theta)) , (mewSun/h)*(elem[1]+math.cos(theta)) , 0]
    
    
    p1 = math.cos(math.radians(elem[2]))
    p2 = math.sin(math.radians(elem[2]))
    p3 = math.cos(math.radians(elem[3]))
    p4 = math.sin(math.radians(elem[3]))
    p5 = math.cos(math.radians(w))
    p6 = math.sin(math.radians(w))

    #Calculating the transformation matrix    
    QT = [[(((-p4)*p1*p6)+(p3*p5)) , ((-p4*p1*p5)-(p3*p6)) , p4*p2] , 
       [((p3*p1*p6)+(p4*p5)) , ((p3*p1*p5)-(p4*p6)) , -p3*p2] , 
        [p2*p6 , p2*p5 , p1]]
        
           
    r=[0,0,0]
    v=[0,0,0]
        
    for i in range(3):
        for j in range(3):
            r[i]= r[i]+r_hat[j]*QT[i][j]        

    for i in range(3):
        for j in range(3):
            v[i]= v[i]+v_hat[j]*QT[i][j]
    
    
    
    Oelem=[0,0,0,0,0,0,0,0,0,0,0]
    
    Oelem[0]=h
    Oelem[1]=elem[1]
    Oelem[2]=elem[3]    
    Oelem[3]=elem[2]
    Oelem[4]=w
    if((math.degrees(theta))>0):
        Oelem[5]=math.degrees(theta)
    else:
        Oelem[5]=360+math.degrees(theta)    
    Oelem[6]=elem[0]
    Oelem[7]=elem[4]
    Oelem[8]=elem[5]
    Oelem[9]=M
    Oelem[10]=E
    
            
    k=2
    if(elem[3]<0):
        elem[3]=abs(elem[3])
    if(elem[3]>180):
        elem[3]=elem[3]%180
    while(k<11 and k!=3 and k!=6):
        if (elem[k]<0):
            elem[k]=360+elem[k]
        if(elem[k]>360):
            elem[k]=elem[k]%360
        k=k+1
    
        
    return r,v,Oelem
    
#Function to display the orbital elements of the planets   
def Output_oelem(elem):
             
              print('(h)angular momentum in(km^2/s) : %12.8E'%(elem[0])) 
              print('(e)eccentricity : %.8f'%(elem[1])) 
              print('(RA)right ascension in(deg) : %.4f'%(elem[2])) 
              print('(i)inclination in(deg) : %.4f'%(elem[3]))
              print('(w)argument of perihelion in(deg) : %.4f'%(elem[4]))
              print('(TA)true anomaly in(deg) : %.4f'%(elem[5]))
              print('(a)semimajor axis in(AU) : %.8f'%(elem[6]))
              print('(w_hat)longitude of perihelion ( = RA + w) in(deg) : %.4f'%(elem[7])) 
              print('(L)mean longitude ( = w_hat + M) in(deg) : %.4f'%(elem[8]))
              print('(M)mean anomaly in(deg) : %.4f'%(elem[9]))
              print('(E)eccentric anomaly in (deg) : %.4f '%(elem[10]))
              print()
 
   
#Calculating the number of centuries between J2000 and the given date
T1= (JD1-2451545)/36525    
T2= (JD2-2451545)/36525

elem1=[]
r1,v1,elem1 = elemUpdate(elem1,T1,n1) 
elem2=[]
r2,v2,elem2 = elemUpdate(elem2,T2,n2)

r1_mag=mag(r1)
r2_mag=mag(r2)
v1_mag=mag(v1)
v2_mag=mag(v2)

print('The State Vectors of planet 1 are : ')
print('r = (in AU) [%.8f , %.8f , %.8f]  '%(r1[0]/AU,r1[1]/AU,r1[2]/AU))
print('The radial distance from the sun(in AU) is %.8f.'%(r1_mag/AU))
print('v = (in km/s) [%.8f , %.8f , %.8f  '%(v1[0],v1[1],v1[2]))
print('The speed is %.8f km/s.'%(v1_mag))
print()
print('The orbital elements of planet 1 are : ')
Output_oelem(elem1)

print('The State Vectors of planet 2 are : ')
print('r = (in AU) [%.8f , %.8f , %.8f]  '%(r2[0]/AU,r2[1]/AU,r2[2]/AU))
print('The radial distance from the sun(in AU) is %.8f.'%(r2_mag/AU))
print('v = (in km/s) [%.8f , %.8f , %.8f]  '%(v2[0],v2[1],v2[2]))
print('The speed is %.8f km/s.'%(v2_mag))
print()
print('The orbital elements of planet 2 are : ')
Output_oelem(elem2)


del_t=abs(JD1-JD2)




r1_dot_r2=(r1[0]*r2[0])+(r1[1]*r2[1])+(r1[2]*r2[2])
r1_cross_r2=[((r1[1]*r2[2])-(r1[2]*r2[1])) , -((r1[0]*r2[2])-(r1[2]*r2[0])) , ((r1[0]*r2[1])-(r1[1]*r2[0]))]

print('Enter the number corresponding to the type of trajectory you want: ')
print('1 : Prograde trajectory')
print('2 : Retrograde trajectory')
num=int(input('n = '))
if(num==1):
    if(r1_cross_r2[2]>=0):
        del_theta=math.degrees(numpy.arccos((r1_dot_r2)/(r1_mag*r2_mag)))
    else:
        del_theta=360-math.degrees(numpy.arccos((r1_dot_r2)/(r1_mag*r2_mag)))
else:
    if(r1_cross_r2[2]<0):
        del_theta=math.degrees(numpy.arccos((r1_dot_r2)/(r1_mag*r2_mag)))
    else:
        del_theta=360-math.degrees(numpy.arccos((r1_dot_r2)/(r1_mag*r2_mag)))


A = (math.sin(math.radians(del_theta)))*(math.sqrt((r1_mag*r2_mag)/(1-math.cos(math.radians(del_theta)))))


#Function to calculate the value of z by Newton's method by solving the equations 
#of Universal Variable Formulation and Lagrange Coefficients
def Newton_prime(del_t):

    z=0
    
    C=0.5
    S=1.0/6
            
    Y = r1_mag+r2_mag+((A*((z*S)-1))/(math.sqrt(C)))
    
    f = (((Y/C)**(1.5))*S)+(A*math.sqrt(Y))-(del_t*math.sqrt(mewSun))
    
    fD = (((math.sqrt(2.0))/40)*(Y**(1.5))) + ((A/8)*((math.sqrt(Y))+(A*math.sqrt(0.5/Y))))
    
    ratio=f/fD
     
    while(abs(ratio)>10**-8):
        z=z-ratio
        
        if(z>0):
            C=(1-math.cos(math.sqrt(z)))/z
            S=(math.sqrt(z)-math.sin(math.sqrt(z)))/((math.sqrt(z))**3)
        if(z<0):
            C=((math.cosh(math.sqrt(-z)))-1)/(-z)
            S=((math.sinh(math.sqrt(-z)))-(math.sqrt(-z)))/((math.sqrt(-z))**3)
        if(z==0):
            C=0.5
            S=1.0/6
    
        Y = r1_mag+r2_mag+((A*((z*S)-1))/(math.sqrt(C)))
        
        f = (((Y/C)**(1.5))*S)+(A*math.sqrt(Y))-(del_t*math.sqrt(mewSun))
        if(z!=0):
            fD = (((Y/C)**(1.5))*(((1/(2*z))*(C-(1.5*S/C)))+((3*(S**2))/(4*C)))) + ((A/8)*((3*S*(math.sqrt(Y))/C)+(A*math.sqrt(C/Y))))
        else:
            fD = (((math.sqrt(2.0))/40)*(Y**(1.5))) + ((A/8)*((math.sqrt(Y))+(A*math.sqrt(0.5/Y))))
    
        ratio=f/fD
        
     
    return(z)



z = Newton_prime(del_t*86400)

if(z>0):
    C=(1-math.cos(math.sqrt(z)))/z
    S=(math.sqrt(z)-math.sin(math.sqrt(z)))/((math.sqrt(z))**3)
if(z<0):
    C=((math.cosh(math.sqrt(-z)))-1)/(-z)
    S=((math.sinh(math.sqrt(-z)))-(math.sqrt(-z)))/((math.sqrt(-z))**3)
if(z==0):
    C=0.5
    S=1.0/6
Y = r1_mag+r2_mag+((A*((z*S)-1))/(math.sqrt(C)))



#lagrange coefficients
f = 1 - (Y/r1_mag)
g = A*(math.sqrt(Y/mewSun))
g_dot = 1 - (Y/r2_mag)


vD=[0,0,0]
vA=[0,0,0]
for k in range(3):
    vD[k]=(1/g)*(r2[k]-(f*r1[k]))
    vA[k]=(1/g)*((g_dot*r2[k])-r1[k])
    

vD_mag=mag(vD) 
vA_mag=mag(vA) 


#Function to display the orbital elements for the interplanetary orbit that has been determined  
def Output_elem(elem):
             
                              
             
              print()
              print('(h)angular momentum in(km^2/s) : %12.8E '%(elem[0])) 
              print('(e)eccentricity : %.8f'%(elem[1])) 
              print('(RA)right ascension in(deg) : %.4f'%(elem[2])) 
              print('(i)inclination in(deg) : %.4f'%(elem[3]))
              print('(w)argument of perihelion in(deg) : %.4f'%(elem[4]))
              print('(TA)true anomaly in(deg) : %.4f'%(elem[5]))
              print()
              if(elem[1]>1):
                  print('The orbit is a HYPERBOLA with respect to the sun.')
              if(elem[1]<1):
                  print('The orbit is an ELLIPSE with respect to the sun.')
              if(elem[1]==1):
                  print('The orbit is a PARABOLA with respect to the sun.')
              print()    


#Function to return the orbital elements from a given state vector
def Orbital_elements(r,v):
    r_mag = mag(r)
    v_mag = mag(v)
    r_dot_v = (r[0]*v[0])+(r[1]*v[1])+(r[2]*v[2])  
    v_radial = r_dot_v/r_mag
    
    h=[((r[1]*v[2])-(r[2]*v[1])) , -((r[0]*v[2])-(r[2]*v[0])) , ((r[0]*v[1])-(r[1]*v[0]))]
    h_mag = mag(h)
    
    i = math.degrees(numpy.arccos(h[2]/h_mag))
    
    k=[0,0,1]
    N=[((k[1]*h[2])-(k[2]*h[1])) , -((k[0]*h[2])-(k[2]*h[0])) , ((k[0]*h[1])-(k[1]*h[0]))]
    N_mag = mag(N)
    
    if(N[1]>=0):
        RA = math.degrees(numpy.arccos(N[0]/N_mag))
    else:
        RA = 360-math.degrees(numpy.arccos(N[0]/N_mag))
        
    e=[0,0,0]    
    for k in range(3):
        e[k]=(1.0/mewSun)*((((v_mag**2)-(mewSun/r_mag))*r[k])-(r_mag*v_radial*v[k]))    
    
    e_mag = mag(e)   
    
    N_dot_e=(N[0]*e[0])+(N[1]*e[1])+(N[2]*e[2])
    
    if(e[2]>=0):
        w = math.degrees(numpy.arccos(N_dot_e/(N_mag*e_mag)))
    else:
        w = 360-math.degrees(numpy.arccos(N_dot_e/(N_mag*e_mag)))
    
    e_dot_r=(e[0]*r[0])+(e[1]*r[1])+(e[2]*r[2])    
    if(v_radial>=0):
        theta = math.degrees(numpy.arccos(e_dot_r/(e_mag*r_mag)))
    else:
        theta = 360-math.degrees(numpy.arccos(e_dot_r/(e_mag*r_mag)))
        
    elem=[0,0,0,0,0,0]
    elem[0]=h_mag
    elem[1]=e_mag
    elem[2]=RA
    elem[3]=i
    elem[4]=w
    elem[5]=theta    
    
    if(elem[3]<0):
        elem[3]=abs(elem[3])
    if(elem[3]>180):
        elem[3]=elem[3]%180
    k=4
    while(k<6):
        if (elem[k]<0):
            elem[k]=360+elem[k]
        if(elem[k]>360):
            elem[k]=elem[k]%360
        k=k+1
    return(elem)



elem_IO = Orbital_elements(r1,vD)
print()
print('The orbital elements for the interplanetary mission are: ' )
Output_elem(elem_IO)

v_inf_D=[0,0,0]
v_inf_A=[0,0,0]
for k in range(3):
    v_inf_D[k]=vD[k]-v1[k]
    v_inf_A[k]=vA[k]-v2[k]
    


v_inf_D_mag = mag(v_inf_D)
v_inf_A_mag = mag(v_inf_A)
print('The required departure velocity at the sphere of influence of planet1 in km/s :')
print('v = (in km/s) [%.8f , %.8f , %.8f]  '%(v_inf_D[0],v_inf_D[1],v_inf_D[2]))
print('Departure speed required : %.8f km/s.'%(v_inf_D_mag))
print('The required arrival velocity at the sphere of influence of planet2 in km/s :')
print('v = (in km/s) [%.8f , %.8f , %.8f]  '%(v_inf_A[0],v_inf_A[1],v_inf_A[2]))    
print('Entry speed required : %.8f km/s.'%(v_inf_A_mag))
print()




def departure(mew,radius):
    
    print('Enter the height of the initial circular parking orbit in km : ')
    height = float(input('h = '))

    e_hyp = 1.0+((height+radius)*(v_inf_D_mag**2)/mew)
    h_hyp = mew*(math.sqrt((e_hyp**2)-1))/v_inf_D_mag
    
    v_orbit = math.sqrt(mew/(height+radius))
    v_hyp = (h_hyp/(height+radius))
    
    delta_v_dep = v_hyp-v_orbit
    beta =abs(math.degrees(numpy.arccos(1.0/e_hyp)))
    
    
    print('The eccentricity of the required escape hyperbola is : %.8f'%(e_hyp))
    print('The delta_v required is : %.8f in km/s'%(delta_v_dep) , 'at an angle of %.4f'%(beta) , "with respect to the planet's heliocentric velocity vector.")
    
    print()
    
departure(mew_list[n1],radius_list[n1])

def arrival(mew,radius):
    print('Enter the height of the circular parking orbit where you want to place the spacecraft in km : ')
    height = float(input('h = '))

    vP_hyp=math.sqrt((v_inf_A_mag**2)+(2*mew/(height+radius)))
    e_hyp = 1.0+((height+radius)*(v_inf_A_mag**2)/mew)
    v_orbit = math.sqrt(mew/(height+radius))
    
    delta_v_arr = v_orbit-vP_hyp
    
    beta =abs(math.degrees(numpy.arccos(1.0/e_hyp)))
    
    print('The eccentricity of the required entry hyperbola is : %.8f'%(e_hyp))
    print('The delta_v required is : %.8f km/s'%(delta_v_arr) , 'at an angle of %.4f'%(beta) , "with respect to the planet's heliocentric velocity vector.")
    
    
arrival(mew_list[n2],radius_list[n2])    

print()
user = input('  ')

