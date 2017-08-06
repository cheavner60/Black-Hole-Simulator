from numpy import *
from matplotlib.pyplot import *
import scipy.integrate as integ
from astropy.io import ascii



def derivs(rv,t,M,m,c,lbar,dphist):
    r=rv[0]
    phi=rv[1]
    john=(1-(2*M)/r)*(1+(lbar**2/r**2)/c**2)
    ebarsq=(1-(2*M)/r)*(1+(lbar**2/r**2)/c**2)+0.0025
    drdtau=sqrt(ebarsq-john)
    #drdphi=sqrt((ebarsq-john)/(lbar**2/r**4))
    dphidtau=(1/r**2)*lbar
    #import pdb; pdb.set_trace()
    return ([drdtau,dphidtau])

def final(M2=1.99e30,m=6.0e24,phi=0.0,sh=1,lv=3.0e4,v=2.e4,n=1,tstep=1e-3,r=1.5e11):
    c=3e8
    M=6.673e-11*M2/(c**2)
    schw=2*M
    r=sh*schw
# lv is the velocity of the earth around the sun in its current orbit
    #lv=
    lbar=r*lv
    rv=(r,phi)
    x=rv[0]*cos(rv[1])                                                        
    y=rv[0]*sin(rv[1]) 
# p-seption
    p=2*pi*sqrt(r**3/(M*c**2))
# conversion from meters to seconds
    #p1=p/2.99e8
    #p=sqrt((4*pi**2*a**3)/(g*tm))
    dt=p*tstep
    nsteps=int(n/tstep)
    x9000=1/tstep/40
    dphi=6*pi*M/r
    dphist=dphi/250
    print ('lbar   =  '+str(lbar))
    print ('dphi   =  '+str(dphi))
    print ('r      =  '+str(r))
    print ('r^2    =  '+str(r**2))
    print ('M      =  '+str(M))
    print ('M2     =  '+str(M2))
    print ('p      =  '+str(p))
    #print ('p1     =  '+str(p1))
    print ('dt     =  '+str(dt))
    print ('x9000  =  '+str(x9000))
    print ('nsteps =  '+str(nsteps))
    print ('dphist =  '+str(dphist))
    #return('your mom')
    #dt=dt*1e4/2
    print dt
    #return('your mom again')
    figure()
    subplot(1,1,1)
    gca().set_aspect('equal')
    plot(x,y,'bo')
    plot(0,0,'go')
    xlim([-3*r,3*r])
    ylim([-3*r,3*r])
    t=(0,dt)
    #figure()
    r2d2=[]
    c3p0=[]
    gv=(0,0,0,0,0,0)
    for n in range(nsteps):
        #import pdb; pdb.set_trace()
        value=integ.odeint(derivs,rv,t,args=(M,m,c,lbar,dphist))
        rv=value[1]
        x=rv[0]*cos(rv[1]+dphist)
        y=rv[0]*sin(rv[1]+dphist)
        r2d2.append(rv[0])
        c3p0.append(rv[1])
        if n%x9000==0:
            #gv=(x,y,gv[0],gv[1],gv[2],gv[3])
            #cla()
            xlim([-3*r,3*r])
            ylim([-3*r,3*r])
            plot(x,y,'bo')
            #plot(gv[2],gv[3],'bo')
            #plot(gv[4],gv[5],'bo')
            plot(0,0,'go')
            draw()
    #print r2d2
    #print c3p0
