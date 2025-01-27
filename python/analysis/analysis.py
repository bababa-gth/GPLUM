from .classdef import *
from scipy import stats

PI = m.acos(-1.)
L = 14959787070000
M = 1.9884e33

def readSnap(filename, bHeader=True):
    """
    Read data of particles from a snapshot file

    Arguments
    ----------
    filename : character string. snapshot filename
    bHeader :  boolian. flag as to whether header exists in snapshot (default True)
    
    Returns
    ----------
    header : Header. header in snapshot
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    """
    
    i = 0
    pp = {}
    header = Header(0.,0,0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

    with open(filename) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]

            if i==0 and bHeader :
                header = Header(float(part[0]),\
                                int(part[1]),\
                                int(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]),\
                                float(part[12]) )
                
            else :
                idx = int(part[0])
                ptcl = Particle(float(part[1]),\
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                int(part[10]),\
                                int(part[11]) )
                pp[idx] = ptcl;

            i = i + 1

    return [header, pp]

def readSnapForParticles(filename, ids, bHeader=True):
    """
    Read data of paticular particles from a snapshot file

    Arguments
    ----------
    filename : character string. snapshot filename
    ids :      list of int. list of particle id of particles to be read
    bHeader :  boolian. flag as to whether header exists in snapshot
    
    Returns
    ----------
    header : Header. header in snapshot
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    """

    i = 0
    pp = {}
    
    with open(filename) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]

            if i==0 and bHeader :
                header = Header(float(part[0]),\
                                int(part[1]),\
                                int(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]),\
                                float(part[12]) )
                
            else :
                idx = int(part[0])
                if idx in ids :
                    ptcl = Particle(float(part[1]),\
                                    float(part[2]),\
                                    float(part[3]),\
                                    float(part[4]),\
                                    float(part[5]),\
                                    float(part[6]),\
                                    float(part[7]),\
                                    float(part[8]),\
                                    float(part[9]),\
                                    int(part[10]),\
                                    int(part[11]) )
                    pp[idx] = ptcl;

            i = i + 1

    return [header, pp]


def clacOrbitalElements(pp, m_sun=1.) :
    """
    Calculate orbital elements of paticles

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    m_sun :  float. mass of a central star (default 1.0)
    
    Returns
    ---------
    m_max:   float. largest particle mass
    m_ave:   float. mean particle mass
    ecc_rms: float. root mean square of eccentricities
    inc_rms: float. root mean square of inclinations
    """

    n_body = len(pp)

    m_max = 0.
    m_ave = 0.
    ecc_rms = 0.
    inc_rms = 0.
    for ptcl in pp.values():
        ax, ecc, inc = ptcl.getOrbitalElement(m_sun)
        ecc_rms = ecc_rms + ecc*ecc
        inc_rms = inc_rms + inc*inc
        m_ave = m_ave + ptcl.mass
        if ptcl.mass > m_max :
            m_max = ptcl.mass

    ecc_rms = m.sqrt(ecc_rms / n_body)
    inc_rms = m.sqrt(inc_rms / n_body)
    m_ave = m_ave / n_body

    return [m_max, m_ave, ecc_rms, inc_rms]


def calcEnergy(pp, m_sun=1., eps=0.):
    """
    Calculate energy of particle system

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    m_sun :  float. mass of a central star (default 1.0)
    eps      float. softening parameter (default 0.0)
    
    Returns
    ----------
    etot :        float. total energy
    ekin :        float. kinetic energy
    ephi_sun :    float. potential energy by central star
    ephi_planet : float. potential energy by interaction of particles
    """

    ekin = 0.
    ephi_sun = 0.
    ephi_planet = 0.
    
    for id, ptcl in pp.items():
        dr2 = ptcl.pos2() + eps**2
        dv2 = ptcl.vel2()
        
        ekin     += ptcl.mass * dv2
        ephi_sun += ptcl.mass / sqrt(dr2)

        for id2, ptcl2 in pp.items():
            if id < id2 :
                dr2 = ptcl.pos2(ptcl2) + eps**2               
                ephi_planet += ptcl.mass * ptcl2.mass / sqrt(dr2)
                
    ekin *= 0.5
    ephi_sun *= m_sun
    etot = ekin + ephi_sun + ephi_planet

    return [etot, ekin, ephi_sun, ephi_planet]


def calcCumulativeNumberDistribution(pp, m_sun=1.0, var="mass"):
    """
    Make set of the cumulative number distributuion.

    Arguments
    ----------
    pp :     dictionary {int : Particle}. particle system. key indicats particle id
    m_sun :  float. mass of the central star (default 1.0)
    var :    charactor string. variable of number distribution (default mass)
             "mass" -- mass distribution
             "ecc"  -- eccentricity distribution
             "inc"  -- inclination distribution
    
    Returns
    ----------
    N_c :    list of [float, int]. set of mass of particles and cumulative number
    """

    v_var = []
    
    for ptcl in pp.values():
        if var == "ecc" :
            if not hasattr(ptcl, "ecc"): ptcl.getOrbitalElement(m_sun)
            v_var.append(ptcl.ecc)
        elif var == "inc" :
            if not hasattr(ptcl, "inc"): ptcl.getOrbitalElement(m_sun)
            v_var.append(ptcl.inc)
        else :
            v_var.append(ptcl.mass)

    v_var.sort()
    v_var.reverse()

    N_c = []
    n = 1
    for v in v_var :
        N_c.append([v,n])
        n = n + 1

    return N_c

def calcEccentricityDistribution(pp, w_b=2., m_sun=1.0) :
    """
    Calculate distribution of eccentricity and inclination against mass

    Arguments
    ----------
    pp :      dictionary {int : Particle}. particle system. key indicats particle id
    w_b :     float. width of mass bin in log scale (default 2.0) 
    m_sun :   float. mass of a central star (default 1.0)
    
    Returns
    ----------
    dist :    list of [float, float, int, float, float, 
                       float, float, float, float, float, float]. 
              minimun mass of bin, maximum mass of bin, 
              the number of particle, total mass, mean mass, 
              rms of eccentricity, 
              bottom and top of 68% confidence interval of rms of eccentricity
              rms of inclination, 
              bottom and top of 68% confidence interval of rms of inclination
    """
    
    Ptcl = []
    for ptcl in pp.values():
        if not (hasattr(ptcl, "ecc") and hasattr(ptcl, "inc")):
            ptcl.getOrbitalElement(m_sun)
        Ptcl.append([ptcl.mass, ptcl.ecc, ptcl.inc])

    Ptcl.sort()
    
    m_max = Ptcl[-1][0]
    m0 = w_b**(int(m.log(min([i for i,j,k in Ptcl]))/m.log(w_b)))
    m1 = m0*w_b
    
    dist = []
    i = 0
    while m0 < m_max :
        mass_tot  = 0.
        mass_mean = 0.
        mu_ecc    = 0.
        mu_inc    = 0.
        s2_ecc = 0.
        s2_inc = 0.
        rms_ecc = 0.
        rms_inc = 0.
        btm_ecc = 0.
        top_ecc = 0.
        btm_inc = 0.
        top_inc = 0.
        n_p = 0
        while i < len(Ptcl) and Ptcl[i][0] < m1 :
            mass_tot  = mass_tot  + Ptcl[i][0]
            mass_mean = mass_mean + m.log(Ptcl[i][0])
            mu_ecc    = mu_ecc    + Ptcl[i][1]
            mu_inc    = mu_inc    + Ptcl[i][2]
            s2_ecc    = s2_ecc + Ptcl[i][1]*Ptcl[i][1]
            s2_inc    = s2_inc + Ptcl[i][2]*Ptcl[i][2]
            n_p  = n_p + 1
            i = i + 1

        if n_p > 0 :
            mass_mean = m.exp(mass_mean/n_p)
            mu_ecc    = mu_ecc/n_p
            mu_inc    = mu_inc/n_p
        else :
            mass_tot  = 0.
            mass_mean = None
            mu_ecc    = None
            mu_inc    = None
        if n_p > 1 :
            rms_ecc = m.sqrt(s2_ecc/n_p)
            rms_inc = m.sqrt(s2_inc/n_p)
            
            s2_ecc = s2_ecc - n_p*mu_ecc*mu_ecc
            s2_inc = s2_inc - n_p*mu_inc*mu_inc

            sigma2_mu_ecc = s2_ecc/(n_p*(n_p-1))
            sigma2_mu_inc = s2_inc/(n_p*(n_p-1))

            k = n_p-1
            chi2_30, chi2_70 = stats.chi2.interval(alpha=0.68, df=k)

            btm_ecc = m.sqrt(mu_ecc*mu_ecc + s2_ecc/k
                             - m.sqrt((s2_ecc*(1./k-1./chi2_70))**2. + sigma2_mu_ecc**2))
            top_ecc = m.sqrt(mu_ecc*mu_ecc + s2_ecc/k
                             + m.sqrt((s2_ecc*(1./chi2_30-1./k))**2. + sigma2_mu_ecc**2))
            btm_inc = m.sqrt(mu_inc*mu_inc + s2_inc/k
                             - m.sqrt((s2_inc*(1./k-1./chi2_70))**2. + sigma2_mu_inc**2))
            top_inc = m.sqrt(mu_inc*mu_inc + s2_inc/k
                             + m.sqrt((s2_inc*(1./chi2_30-1./k))**2. + sigma2_mu_inc**2))
            
        else :
            rms_ecc = m.sqrt(s2_ecc/n_p)
            rms_inc = m.sqrt(s2_inc/n_p)
            
            btm_ecc = None
            top_ecc = None
            btm_inc = None
            top_inc = None



            
        dist.append([m0, m1, n_p, mass_tot, mass_mean, \
                     rms_ecc, btm_ecc, top_ecc, \
                     rms_inc, btm_inc, top_inc])

        m0 = m1
        m1 = m0*w_b

    return dist


def outputOrbitalElements(filename0, filename1, bHeader=True, m_sun=1.0) :
    """
    Make snapshot file with orbital elements

    Arguments
    ----------
    filename0 : character string. input snapshot filename
    filename1 : character string. output snapshot filename
    m_sun :     float. mass of the central star (default 1.0)
    bHeader :   boolian. flag as to whether header exists in snapshot (default True)
    
    Returns
    ----------
    """

    header, pp = readSnap(filename0, bHeader)
    m_max, m_ave, ecc_rms, inc_rms = clacOrbitalElements(pp, m_sun)

    f = open(filename1, mode='w')
    
    line = "{:.8f}".format(header.time) + "\t" \
           + "{:d}".format(header.n_body) + "\t" \
           + "{:.15e}".format(m_max) + "\t" \
           + "{:.15e}".format(m_ave) + "\t" \
           + "{:.15e}".format(ecc_rms) + "\t" \
           + "{:.15e}".format(inc_rms) + "\n"
    f.write(line)

    for i, ptcl in pp.items():
        line = "{:d}".format(i) + "\t" \
               + "{:.15e}".format(ptcl.mass) + "\t" \
               + "{:.15e}".format(ptcl.x) + "\t" \
               + "{:.15e}".format(ptcl.y) + "\t" \
               + "{:.15e}".format(ptcl.z) + "\t" \
               + "{:.15e}".format(ptcl.vx) + "\t" \
               + "{:.15e}".format(ptcl.vy) + "\t" \
               + "{:.15e}".format(ptcl.vz) + "\t" \
               + "{:.15e}".format(ptcl.ax) + "\t" \
               + "{:.15e}".format(ptcl.ecc) + "\t" \
               + "{:.15e}".format(ptcl.inc) + "\n"
        f.write(line)

    f.close()
    

def convertSnap(filename0, filename1, bHeader=True) :
    """
    Convert snapshot file in version 2.0 or earlier format to in version 2.1 format

    Arguments
    ----------
    filename0 : character string. input snapshot filename
    filename1 : character string. output snapshot filename
    
    Returns
    ----------
    """

    i = 0
    pp = {}
    header = Header(0.,0,0,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.)

    with open(filename0) as f:
        for line in f:
            part = [p.strip() for p in line.split("\t")]
            
            if i==0 and bHeader :
                header = Header(float(part[0]),\
                                int(part[1]),\
                                int(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                float(part[8]),\
                                float(part[9]),\
                                float(part[10]),\
                                float(part[11]),\
                                float(part[12]) )
                
            else :
                idx = int(part[0])
                ptcl = Particle(float(part[1]),\
                                0., 0.,
                                float(part[2]),\
                                float(part[3]),\
                                float(part[4]),\
                                float(part[5]),\
                                float(part[6]),\
                                float(part[7]),\
                                int(part[8]),\
                                0 )
                pp[idx] = ptcl;
                
            i = i + 1

    f = open(filename1, mode='w')
    
    line = "{:.8f}".format(header.time) + "\t" \
           + "{:d}".format(header.n_body) + "\t" \
           + "{:d}".format(header.id_next) + "\t" \
           + "{:.15e}".format(header.etot0) + "\t" \
           + "{:.15e}".format(header.ekin0) + "\t" \
           + "{:.15e}".format(header.ephi_sun0) + "\t" \
           + "{:.15e}".format(header.ephi_planet0) + "\t" \
           + "{:.15e}".format(header.edisp0) + "\t" \
           + "{:.15e}".format(header.etot1) + "\t" \
           + "{:.15e}".format(header.ekin1) + "\t" \
           + "{:.15e}".format(header.ephi_sun1) + "\t" \
           + "{:.15e}".format(header.ephi_planet1) + "\t" \
           + "{:.15e}".format(header.edisp1) + "\n" 
    f.write(line)

    for i, ptcl in pp.items():
        line = "{:d}".format(i) + "\t" \
               + "{:.15e}".format(ptcl.mass) + "\t" \
               + "{:.15e}".format(ptcl.r_p) + "\t" \
               + "{:.15e}".format(ptcl.f) + "\t" \
               + "{:.15e}".format(ptcl.x) + "\t" \
               + "{:.15e}".format(ptcl.y) + "\t" \
               + "{:.15e}".format(ptcl.z) + "\t" \
               + "{:.15e}".format(ptcl.vx) + "\t" \
               + "{:.15e}".format(ptcl.vy) + "\t" \
               + "{:.15e}".format(ptcl.vz) + "\t" \
               + "{:d}".format(ptcl.ngb) + "\t" \
               + "{:d}".format(ptcl.flag) + "\n"
        f.write(line)

    f.close()
        
        
    

    

    
