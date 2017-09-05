import numpy as n
from numpy import random as r
import matplotlib as m
from matplotlib import pyplot as p
from pylab import *

def orbits(mstar,xstar,ystar):
    """
    Purpose: Calculates the closest approach of Earth to each planet 

    Arguments:
              mstar: the mass in jupiter masses of a second star
              xstar: the x-position of the star in AU from Sun
              ystar: the y-position of the star in AU from Sun
    Returns: closest approach of Earth and planets

    Written: Ryan A. Arneson, UCI, 5/2012
    """

    #f_e = (G*m_sun*m_e)/rse**2
    #a_e = f_e/m_e
    names =['Mercury','Venus','Earth','Mars','Jupiter','Saturn','Uranus','Neptune']
    #masses of sun & planets relative to Jupiter
    masses = [1047.0, 0.000174, 0.00256, 0.00315, 0.000338, 1.0, 0.299, 0.0457, 0.054, mstar]
    #initial position/velocity of mars [AU]/[AU/day] J2000 on 10/14/1986 from Sun center
    #data from JPL Horizons
    x_os = [0.0, 1.094050356328868E-1, 7.194440093584461E-1 , 9.340689685363645E-1 ,1.351605694231208, 4.911497561483949, -3.284886890712781,-2.657597541284523,2.776243778355217, xstar]
    y_os = [0.0, -4.383200854282770E-1,8.941943690872418E-2 , 3.500364313014421E-1, -2.940663739212329E-1, -7.732194781270249E-1, -9.448299632406481, -1.898439779002059E+1, -3.010866093908099E+1, ystar]
    z_os = [0.0, -4.584423177987438E-2, -4.031542348923110E-2, 1.395568488843737E-5,-3.941717590580762E-2, -1.067994388145757E-1,2.953664162377747E-1, -3.595484504191691E-2, 5.559822452731458E-1, 0.0]
    vx_o = [0.0, 2.165680820476263E-2, -2.573947699771772E-3 , -6.320062946246359E-3, 3.513131305001818E-3, 1.075638452457469E-3, 4.970374587731719E-3, 3.868292761512772E-3, 3.108940218824626E-03, 0.0]
    vy_o = [0.0, 8.254797988512067E-3 , 1.997954220673672E-2, 1.603987294908900E-2 , 1.486942594392115E-2,7.811612791701780E-3, -1.843926124056859E-3, -7.225554801124897E-4, 3.124172465610533E-4, 0.0]
    vz_o = [0.0,-1.314325716682900E-3 , 4.209228228452999E-4, 1.495322157966427E-8 ,2.249875571332907E-4 , -5.642166311988968E-5, -1.655256265475190E-4, -5.296781927226852E-5,-7.796160025424588E-5,0.0]
    ##Assume for now circular orbits all in plane
    #using Euler method
    dt = 0.2
    tf = 4000
    x = n.arange(0,tf,dt)
    y = n.arange(0,tf,dt)
    p.close()
    for k in range(1,5):
        x[0] = x_os[k]
        y[0] = y_os[k]
        vx = vx_o[k]
        vy = vy_o[k]
        for i in range(1,n.size(x)):
            r1 = n.sqrt(x[i-1]**2 + y[i-1]**2)
            r2 = n.sqrt((x[i-1]-x_os[9])**2 + (y[i-1]-y_os[9])**2)
            F_m1 = ((2.824E-7)*(masses[0]*masses[k]))/r1**2
            F_m2 = ((2.824E-7)*(masses[9]*masses[k]))/r2**2
            a_x1 = -(x[i-1]/r1)*(F_m1/masses[k])
            a_y1 = -(y[i-1]/r1)*(F_m1/masses[k])
            a_x2 = -((x[i-1]-x_os[9])/r2)*(F_m2/masses[k])
            a_y2 = -((y[i-1]-y_os[9])/r2)*(F_m2/masses[k])
            x[i] = x[i-1] + dt*(vx + (a_x1+a_x2)*dt)
            y[i] = y[i-1] + dt*(vy + (a_y1+a_y2)*dt)
            vx += (a_x1+a_x2)*dt
            vy += (a_y1+a_y2)*dt
        #n.save('x_'+str(k),x)
        #n.save('y_'+str(k),y)
        subplot(121)
        p.plot(x,y,'-',label=names[k-1])
        p.plot(x_os[9],y_os[9],'ro')
        p.plot(x_os[0],y_os[0],'yo')
        p.ylim((-2,2))
        p.xlim((-2,2))
        p.xlabel("AU")
        p.ylabel("AU")
        p.title("Inner Planets")
        p.legend()
        #outer planets (larger step size w/same resoltion)
    
    dt = 1.0
    tf = 80000
    x = n.arange(0,tf,dt)
    y = n.arange(0,tf,dt)
    for k in range(5,9):
        x[0] = x_os[k]
        y[0] = y_os[k]
        vx = vx_o[k]
        vy = vy_o[k]
        for i in range(1,n.size(x)):
            r1 = n.sqrt(x[i-1]**2 + y[i-1]**2)
            r2 = n.sqrt((x[i-1]-x_os[9])**2 + (y[i-1]-y_os[9])**2)
            F_m1 = ((2.824E-7)*(masses[0]*masses[k]))/r1**2
            F_m2 = ((2.824E-7)*(masses[9]*masses[k]))/r2**2
            a_x1 = -(x[i-1]/r1)*(F_m1/masses[k])
            a_y1 = -(y[i-1]/r1)*(F_m1/masses[k])
            a_x2 = -((x[i-1]-x_os[9])/r2)*(F_m2/masses[k])
            a_y2 = -((y[i-1]-x_os[9])/r2)*(F_m2/masses[k])
            #a_x2 = -((x_os[9]-x[i-1])/r2)*(F_m2/masses[k])
            #a_y2 = -((y_os[9]-y[i-1])/r2)*(F_m2/masses[k])
            x[i] = x[i-1] + dt*(vx + (a_x1+a_x2)*dt)
            y[i] = y[i-1] + dt*(vy + (a_y1+a_y2)*dt)
            vx += (a_x1+a_x2)*dt
            vy += (a_y1+a_y2)*dt
        #n.save('x_'+str(k),x)
        #n.save('y_'+str(k),y)
        subplot(122)
        p.plot(x,y,'-',label=names[k-1])
        p.plot(x_os[9],y_os[9],'ro')
        p.plot(x_os[0],y_os[0],'yo')
        p.ylim((-35,35))
        p.xlim((-35,35))
        p.xlabel("AU")
        p.ylabel("AU")
        p.title("Outer Planets")
        p.legend()
    p.show()
    return
