import math as m
import numpy as np
import random

class Particle:
    """
    Class to store data of particle
    """
    #def __init__(self, mass, x, y, z, vx, vy, vz):
    #    self.mass = mass
    #    self.x  = x
    #    self.y  = y
    #    self.z  = z
    #    self.vx = vx
    #    self.vy = vy
    #    self.vz = vz

    def __init__(self, mass, r_p, f, x, y, z, vx, vy, vz, ngb, flag):
        self.mass = mass
        self.r_p  = r_p
        self.f  = f
        self.x  = x
        self.y  = y
        self.z  = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
        self.ngb = ngb
        self.flag = flag
        
    def pos(self, target=None):
        if target is None:
            return [self.x, self.y, self.z]
        else:
            return [self.x-target.x, self.y-target.y, self.z-target.z]

    def pos2(self, target=None):
        x, y, z = self.pos(target)
        return x**2 + y**2 + z**2
        
    def vel(self, target=None):
        if target is None:
            return [self.vx, self.vy, self.vz]
        else:
            return [self.vx-target.vx, self.vy-target.vy, self.vz-target.vz]

    def vel2(self, target=None):
        vx, vy, vz = self.vel(target)
        return vx**2 + vy**2 + vz**2

    def getOrbitalElement(self, mu=1.):
        r = m.sqrt(self.pos2())
        v2 = self.vel2()
        rv = self.x*self.vx + self.y*self.vy + self.z*self.vz
        rxv = [self.y*self.vz - self.z*self.vy,\
               self.z*self.vx - self.x*self.vz,\
               self.x*self.vy - self.y*self.vx\
        ]
        
        self.ax  = 1./(2./r - v2/mu)
        self.ecc = m.sqrt( (1.-r/self.ax)**2 + (rv)**2/(mu*self.ax) )
        self.inc = m.atan2(m.sqrt(rxv[0]**2+rxv[1]**2), rxv[2])

        return [self.ax, self.ecc, self.inc]

def Get_Cartesian_from_Keplerian(mu,orb_Kep):
    axi=orb_Kep[0]
    ecc=orb_Kep[1]
    inc=orb_Kep[2]
    loa=orb_Kep[3]
    lop=orb_Kep[4]
    tra=orb_Kep[5]

    rn    = axi*(1.e0-ecc*ecc) / ( 1.e0 + ecc*m.cos(tra) )
    hn    = m.sqrt( mu*axi*( 1.e0-ecc*ecc ) )
    dr_dt = hn/(axi*(1.e0-ecc*ecc)) *ecc*m.sin(tra)
    df_dt = hn/(rn*rn)

    rl    = np.array( [ [rn*m.cos(tra)], \
                        [rn*m.sin(tra)], \
                        [0.e0] ] )
    vl    = np.array( [ [dr_dt*m.cos(tra)-rn*m.sin(tra)*df_dt], \
                        [dr_dt*m.sin(tra)+rn*m.cos(tra)*df_dt], \
                        [0.e0] ] )

    #
    A1=np.matrix( [ [m.cos(lop),   -m.sin(lop), 0], \
                    [m.sin(lop),    m.cos(lop), 0], \
                    [         0,             0, 1] \
                  ] )
    A2=np.matrix( [ [   1.e0,               0.e0,               0.e0], \
                    [   0.e0,  m.cos(inc), -m.sin(inc)], \
                    [   0.e0,  m.sin(inc),  m.cos(inc)], \
                  ] )
    A3=np.matrix( [ [m.cos(loa),   -m.sin(loa), 0.e0], \
                    [m.sin(loa),    m.cos(loa), 0.e0], \
                    [      0.e0,          0.e0, 1.e0] \
                  ] )
    rl= A1*rl
    rl= A2*rl
    rl= A3*rl
    vl= A1*vl
    vl= A2*vl
    vl= A3*vl
    orb_Car=[float(rl[0][0]),float(rl[1][0]),float(rl[2][0]),float(vl[0][0]),float(vl[1][0]),float(vl[2][0])]
    return orb_Car

def Create_InitialGPLUM(myfile,pp):
    with open(myfile, 'w') as fout:
        id=0
        for ptcl in pp.values():
            print("{:d} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:.15e} {:d} {:d}".format(id,ptcl.mass,ptcl.r_p,ptcl.f,ptcl.x,ptcl.y,ptcl.z,ptcl.vx,ptcl.vy,ptcl.vz,ptcl.ngb,ptcl.flag), file=fout)  
            id += 1
    return print("Finish create {:s}".format(myfile))


def Get_SolarSystem(pp,n0,n1):
    AU_mks  =   1.495e8    
    axi_SS  =   [0.387e0,0.718e0,1.e0,1.5e0,5.2e0,9.55e0]
    mass_SS =   [1.65e-7,2.54e-6,3.e-6,3.e-7,1.e-3,3.e-4]
    radi_SS =   [2439.4,6051.8,6378.137,3396.19,71492,60268]
    #
    z_SS=1.e-6
    f_radius = 1.e0
    idx=0
    for id in range(n0,n1):
        ptcl = Particle(float(mass_SS[id]),\
                        float(radi_SS[id]/AU_mks),\
                        float(f_radius),\
                        float(axi_SS[id]),\
                        float(0),\
                        float(z_SS),\
                        float(0),\
                        float(m.sqrt(1./axi_SS[id])),\
                        float(0),\
                        int(0),\
                        int(0) )
        pp[idx] = ptcl
        idx += 1
    return pp


def Add_Particle(pp,id,axi,ecc,inc,mp,rp):
    f_radius = 1.e0
    orb_Kep = [axi,ecc,inc,random.uniform(0.e0,2.e0*m.pi),random.uniform(0.e0,2.e0*m.pi),random.uniform(0.e0,2.e0*m.pi)]
    orb_Car = Get_Cartesian_from_Keplerian(1.,orb_Kep)
    ptcl = Particle(float(mp),\
                    float(rp),\
                    float(f_radius),\
                    float(orb_Car[0]),\
                    float(orb_Car[1]),\
                    float(orb_Car[2]),\
                    float(orb_Car[3]),\
                    float(orb_Car[4]),\
                    float(orb_Car[5]),\
                    int(0),\
                    int(0) )
    pp[id] = ptcl

    return pp


def Add_Particles_Ring(pp,n0,n1,a0,Sigma,mp):
    f_radius = 1.e0
    rp=4.266312374581940e-05*(mp/3.e-6)**(1.e0/3.e0)
    Nptc=n1-n0
    Mtot=mp*Nptc
    # set axi  
    daxi = Mtot/(2.*m.pi*a0*Sigma)

    idx=0
    for id in range(n0,n1):
        # axi
        axi = a0 -daxi/2. +daxi/Nptc*(id-n0)
        ecc = 1.e-3
        inc = 0.5e-3
        #
        orb_Kep = [axi,ecc,inc,random.uniform(0.e0,2.e0*m.pi),random.uniform(0.e0,2.e0*m.pi),random.uniform(0.e0,2.e0*m.pi)]
        orb_Car = Get_Cartesian_from_Keplerian(1.,orb_Kep)
        ptcl = Particle(float(mp),\
                        float(rp),\
                        float(f_radius),\
                        float(orb_Car[0]),\
                        float(orb_Car[1]),\
                        float(orb_Car[2]),\
                        float(orb_Car[3]),\
                        float(orb_Car[4]),\
                        float(orb_Car[5]),\
                        int(0),\
                        int(0) )
        pp[id] = ptcl
        idx += 1

    return pp


def Add_Particles_MMRs(pp,n0,n1,mp,p,q):
    f_radius = 1.e0
    Nptc=n1-n0
    rp=4.266312374581940e-05*(mp/3.e-6)**(1.e0/3.e0)
    # set axi
    a0  =  0.1e0
   
    idx=0
    for id in range(n0,n1):
        # axi
        axi = a0*pow(p/q,2.*id/3.)
        ecc = 1.e-3
        inc = 0.5e-3
        #
        orb_Kep = [axi,ecc,inc,random.uniform(0.e0,2.e0*m.pi),random.uniform(0.e0,2.e0*m.pi),random.uniform(0.e0,2.e0*m.pi)]
        orb_Car = Get_Cartesian_from_Keplerian(1.,orb_Kep)
        ptcl = Particle(float(mp),\
                        float(rp),\
                        float(f_radius),\
                        float(orb_Car[0]),\
                        float(orb_Car[1]),\
                        float(orb_Car[2]),\
                        float(orb_Car[3]),\
                        float(orb_Car[4]),\
                        float(orb_Car[5]),\
                        int(0),\
                        int(0) )
        pp[id] = ptcl
        idx += 1

    return pp



#>>>>>>
AU_cgs          = 1.49597870e13
AU_mks          = 1.49597870e11

M_ear_cgs       = 5.972e27
M_ear_mks       = 5.972e24
M_jup_cgs       = 1.898e30
M_jup_mks       = 1.898e27
M_sun_cgs       = 1.989e33
M_sun_mks       = 1.989e30

#
Nptc = 7999
ap=5.2
mp=3.e-6/Nptc #
Me=3.e-6
Sigma=20.e0 *(1./M_sun_cgs)/(1./AU_cgs)**2.

#
Nptc = 3
rho=5.*(1./M_sun_cgs)/(1./AU_cgs)**3.
mp=1.e-3
rp=(3.*mp/(4.*m.pi*rho))**(1./3.)

pp = {}
pp = Add_Particles(pp,0,10.,1.e-3,1.e-3,mp,rp)
pp = Add_Particles_Ring(pp,1,1+Nptc,ap*(2./3.)**(1./3.),Sigma,mp)
Create_InitialGPLUM("INIT00000.dat",pp)
