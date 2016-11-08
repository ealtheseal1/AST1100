# -*- coding: utf-8 -*-
"""
Created on Mon Nov 07 23:34:15 2016

@author: Eirik
"""
from AST1100SolarSystem import AST1100SolarSystem
import numpy as np
import time
import random
myseed = int(68074)

mySystem = AST1100SolarSystem(myseed)


AU = 149597870700 
G = 6.67408e-11 
mS = 1988500e24 
m_p = mySystem.mass[0]*mS 

m_sat = 1100 #kg 
k = 1.38064852e-23 

L = 0.000001
T = 10000 #1e4 #K
N = 100000 
m_h = 3.3476467e-27 

def printer(mySystem):
    print "Analyzing data homeplanet"
    print "------------------------------"
    time.sleep(2)
    print "Home planet mass: %g [solar mass]" %mySystem.starMass
    print "Home planet radius %g [km]" %mySystem.starRadius
    print "Number of planets in system %g" %mySystem.numberOfPlanets
    print "Temperature on Homeplanet %g [K]" %mySystem.temperature
    print 'Semimajor axis [AU] {0:g} \n Eccentricity {1:g} \n Angle {2:g} \n Angle semi-major {3:g} \n mass [Solar mass] {4:g} \n radius [km] {5:g} \n period {6:g} \n x-position {7:g},\n y-position {8:g},\n x-velocity {9:g},\n y-velocity {10:g},\n Atmosphere density {10:g}'.format(Homeplanet[0], Homeplanet[1], Homeplanet[2], Homeplanet[3], Homeplanet[4], Homeplanet[5], Homeplanet[6], Homeplanet[7],Homeplanet[8],Homeplanet[9],Homeplanet[10])
    print 'Calculated escape velocity is %g [m/s]' %(np.sqrt((2*G*m_p)/mySystem.radius[0]))
    print "------------------------------"    
    print "Analyzing data for destination planet" 
    
    time.sleep(2)
    print 'Semimajor axis [AU] {0:g} \n Eccentricity {1:g} \n Angle {2:g} \n Angle semi-major {3:g} \n mass [Solar mass] {4:g} \n radius [km] {5:g} \n period {6:g} \n x-position {7:g},\n y-position {8:g},\n x-velocity {9:g},\n y-velocity {10:g},\n Atmosphere density {10:g}'.format(Destinationplanet[0], Destinationplanet[1], Destinationplanet[2], Destinationplanet[3], Destinationplanet[4], Destinationplanet[5], Destinationplanet[6], Destinationplanet[7],Destinationplanet[8],Destinationplanet[9],Destinationplanet[10])    
    
    
    """ The other planets are not interesting at this point.
    for j in xrange(len(mySystem.a)):
    
        print "%g" %mySystem.omega[j]
        print "%g" %mySystem.psi[j]
        print "%g" %mySystem.mass[j]
        print "%g" %mySystem.radius[j]
        print "%g" %mySystem.period[j]
        print "%g" %mySystem.x0[j]
        print "%g" %mySystem.y0[j]
        print "%g" %mySystem.vx0[j]
        print "%g" %mySystem.vy0[j]
        print "%g" %mySystem.rho0[j]
    """
    return 0

def engine():
    random.seed(500) # seeding the generator to retain numbers between runs
    mu = 0
    o = np.sqrt((k*T)/m_h) # sigma, the standard deviation
    r = np.random.uniform(0,L,(N,3))
    v = np.random.normal(mu,o,(N,3))
    
    hist, bins = histogram(v, 'auto')
    width = 1.1*(bins[1] - bins[0])
    center = (bins[:-1]+bins[1:])/2
    bar(center, hist, align='center', width=width)
    title('velocity $v [m/s]$')
    show()
    
    import scipy.stats as stats
    test = stats.normaltest(r)
    test = np.array(test)
    if test[1,0] > 0.055:
        print "Not Gaussian normaldistr."
    else:
        print "Gaussian test successful"
    
    
    
    
if __name__ == '__main__':   
    Homeplanet = np.array([mySystem.a[0],mySystem.e[0],mySystem.omega[0],mySystem.psi[0],mySystem.radius[0],mySystem.mass[0],mySystem.period[0],mySystem.x0[0],mySystem.y0[0],mySystem.vx0[0],mySystem.vy0[0],mySystem.rho0[0]]) 
    Destinationplanet =np.array([mySystem.a[1],mySystem.e[1],mySystem.omega[1],mySystem.psi[1],mySystem.radius[1],mySystem.mass[1],mySystem.period[1],mySystem.x0[1],mySystem.y0[1],mySystem.vx0[1],mySystem.vy0[1],mySystem.rho0[1]]) 
    printer(mySystem)
    
    engine()
