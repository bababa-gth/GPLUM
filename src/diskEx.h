#pragma once

class DiskProperty {
public:
    DiskProperty();
    ~DiskProperty();

    void Information();
    void calcSigmaGas(PS::F64 rxy);
    void calcSigmaGas_up(PS::F64 rxy);
    void calcSigmaGas_fgap(PS::F64 mp_ms, PS::F64 r_h, PS::F64 alpha_vis);
    void calcMidplaneTemperture(PS::F64 rxy);
    void calcSoundSpeed(PS::F64 Tmid);
    void calcScaleHeight(PS::F64 csound, PS::F64 rxy, PS::F64 Msun);
    void calcrhoGas(PS::F64 SigmaGas, PS::F64 hscale, PS::F64 z2);
    void calcUvel(PS::F64 hscale_rxy_2,PS::F64 rz_hscale_2,PS::F64 alpha,PS::F64 beta,PS::F64vec pvel,PS::F64vec vkep);

    static PS::S32 DiskStg;
    static PS::F64 SigmaGas0_cgs;
    static PS::F64 SigmaGas0;
    static PS::F64 Tmid0_K;
    static PS::F64 Tmid0;
    static PS::F64 alpha_gas;
    static PS::F64 beta_gas;
    static PS::F64 Cd;
    static PS::F64 mu;
    static PS::F64 alpha_vis;

    // Disk parameters
    PS::F64 SigmaGas;
    PS::F64 SigmaGas_up;
    PS::F64 SigmaGas_fgap;

    PS::F64 rhoGas;
    PS::F64 Tmid;
    PS::F64 hscale;
    PS::F64 csound;
    PS::F64vec uvel;
    PS::F64vec vgas;
    PS::F64 alpha_gas_c;
    PS::F64 beta_gas_c;
    PS::F64 eta_gas_c;

private:
    static constexpr PS::F64 m_pi = 3.141592653589793238462643383279502884L;
    static constexpr PS::F64 L_MKS = 149597870700;
    static constexpr PS::F64 L_CGS = 14959787070000;
    static constexpr PS::F64 M_MKS = 1.9884e30;
    static constexpr PS::F64 M_CGS = 1.9884e33;
    static constexpr PS::F64 T     = 365.25*24.*60.*60./(2.*m_pi);

    static constexpr PS::F64 k_B = 1.380649e-16 /(M_CGS*L_CGS*L_CGS)*T*T;
    static constexpr PS::F64 N_A = 6.022140857e23;
    static constexpr PS::F64 m_H = 1./N_A /M_CGS;

};

DiskProperty::DiskProperty() {
}

DiskProperty::~DiskProperty() {
}

void DiskProperty::Information() {

    std::cout << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << "SigmaGas0_cgs = " << SigmaGas0_cgs << std::endl;
    std::cout << "Tmid0_K = " << Tmid0_K << std::endl;
    std::cout << "alpha_gas = " << alpha_gas << std::endl;
    std::cout << "beta_gas = " << beta_gas << std::endl;
    std::cout << "Cd = " << Cd << std::endl;
    std::cout << "mu = " << mu << std::endl;
    std::cout << "--------------------------------" << std::endl;
    std::cout << std::endl;
    std::cout << std::endl;

}

void DiskProperty::calcSigmaGas(PS::F64 rxy) {
    SigmaGas = SigmaGas0 *pow(rxy, -alpha_gas);
}

void DiskProperty::calcSigmaGas_up(PS::F64 rxy) {
    SigmaGas_up = SigmaGas0 *pow(rxy, -alpha_gas);
}

void DiskProperty::calcSigmaGas_fgap(PS::F64 mp_ms, PS::F64 r_h, PS::F64 alpha_vis) {
    PS::F64 K=mp_ms*mp_ms *r_h*r_h*r_h*r_h*r_h /alpha_vis;
    SigmaGas_fgap = 1./(1.+0.04*K);
}

void DiskProperty::calcMidplaneTemperture(PS::F64 rxy) {
    Tmid = Tmid0 *pow(rxy, -beta_gas);
}

void DiskProperty::calcSoundSpeed(PS::F64 Tmid) {
    csound = sqrt(k_B * Tmid / (mu * m_H));
}

void DiskProperty::calcScaleHeight(PS::F64 csound, PS::F64 rxy, PS::F64 Msun) {
    hscale = csound*sqrt(rxy*rxy*rxy/Msun);
}

void DiskProperty::calcrhoGas(
    PS::F64 SigmaGas,
    PS::F64 hscale,
    PS::F64 z2
    ) {
    rhoGas = SigmaGas/(sqrt(2.*m_pi)*hscale) *exp(-z2/(2.*hscale*hscale));
}

void DiskProperty::calcUvel(
    PS::F64 hscale_rxy_2,
    PS::F64 rz_hscale_2,
    PS::F64 alpha,
    PS::F64 beta,
    PS::F64vec pvel,
    PS::F64vec vkep
    ) {
    eta_gas_c = 0.5e0 *hscale_rxy_2 *(1.5e0*(1.-rz_hscale_2) +alpha +beta *(1. +rz_hscale_2));
    vgas = (1.0 - eta_gas_c)*vkep;
    uvel = pvel - vgas;
}


PS::S32 DiskProperty::DiskStg = 0;
PS::F64 DiskProperty::SigmaGas0_cgs = 1.7e3;
PS::F64 DiskProperty::SigmaGas0     = 1.7e3*(1./DiskProperty::M_CGS)/pow(1./DiskProperty::L_CGS, 2.);
PS::F64 DiskProperty::Tmid0_K   = 280.;
PS::F64 DiskProperty::Tmid0     = Tmid0_K;
PS::F64 DiskProperty::alpha_gas = 1.5e0;
PS::F64 DiskProperty::beta_gas = 0.5e0;
PS::F64 DiskProperty::Cd = 1.e0;
PS::F64 DiskProperty::mu = 1.e0;
PS::F64 DiskProperty::alpha_vis = 1.e-3;








class DragForce {
public:
    DragForce();
    ~DragForce();

    // Aerodynamic Drags
    PS::F64vec calcAerodynamicGasDrag_Adachi1976(PS::F64vec u, PS::F64 mass, PS::F64 radius, PS::F64 rho_gas, PS::F64 coeff);

    // Tidal Drags
    PS::F64vec calcTidalDragAxi_usingTimescale(PS::F64vec vp, PS::F64 tau_axi);
    PS::F64vec calcTidalDragEcc_usingTimescale(PS::F64vec rp, PS::F64vec vp, PS::F64 tau_ecc);
    PS::F64vec calcTidalDragInc_usingTimescale(PS::F64 vz, PS::F64 tau_inc);
    // Timescales
    PS::F64 calc_TauTideKanagawa2018(PS::F64 m, PS::F64 ms, PS::F64 r, PS::F64 Sigma, PS::F64 hs);
    PS::F64 calc_TauWave_Tanaka2002(PS::F64 m, PS::F64 ms, PS::F64 r, PS::F64 Sigma, PS::F64 hs);
    PS::F64 calc_TauTideAxi_Cresswell2008(PS::F64 twave, PS::F64 r_h, PS::F64 e, PS::F64 i, PS::F64 alpha);
    PS::F64 calc_TauTideEcc_Cresswell2008(PS::F64 twave, PS::F64 r_h, PS::F64 e, PS::F64 i);
    PS::F64 calc_TauTideInc_Cresswell2008(PS::F64 twave, PS::F64 r_h, PS::F64 e, PS::F64 i);

    template <class Tpsys>
    void calcGasDrag(Tpsys & pp,
                     PS::F64 time);

    static PS::S32 setting;

private:
    static constexpr PS::F64 m_pi = 3.141592653589793238462643383279502884L;
    
    static constexpr PS::F64 L_MKS = 149597870700;
    static constexpr PS::F64 L_CGS = 14959787070000;
    static constexpr PS::F64 M_MKS = 1.9884e30;
    static constexpr PS::F64 M_CGS = 1.9884e33;
    static constexpr PS::F64 T     = 365.25*24.*60.*60./(2.*m_pi);

    static constexpr PS::F64 c_mig = 2.;
    

};

DragForce::DragForce() 
{
    if ( PS::Comm::getRank() == 0 ) {
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "--------------------------------" << std::endl;
        std::cout   << "Use DragForce class: " << std::endl
                    << "Modle id setting is set as " << setting << std::endl
                    << std::endl;                        
        std::cout << "--------------------------------" << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
    }
}

DragForce::~DragForce() 
{
}


template <class Tpsys>
void DragForce::calcGasDrag(Tpsys & pp,
                        PS::F64 time)
{
    const PS::S32 n_loc = pp.getNumberOfParticleLocal(); // Number of Particles in one MPI process
    #pragma omp parallel for
    for(PS::S32 i=0; i<n_loc; i++){
        //
        PS::F64vec fdrag_aero=0.;
        PS::F64vec fdrag_tide=0.;
        PS::F64vec fdrag=0.;

        // Calculate position in the protoplanetary disk
        PS::F64 rxy2 = pp[i].pos.x*pp[i].pos.x + pp[i].pos.y*pp[i].pos.y;
        PS::F64 rxy_inv = 1./sqrt(rxy2);
        PS::F64 rxy = rxy2 * rxy_inv;
        PS::F64 rz2 = pp[i].pos.z*pp[i].pos.z;
        PS::F64 rz  = sqrt(rz2);

        // generate DiskProperty object
        DiskProperty disk;
        disk.calcMidplaneTemperture(rxy);
        disk.calcSoundSpeed(disk.Tmid);
        disk.calcScaleHeight(disk.csound,rxy,FP_t::m_sun); // FP_t::m_sun is the mass of central star
        disk.calcSigmaGas_up(rxy);

        //
        switch (setting) {
            case 0: // aerodynamic gas drag
            {
                disk.SigmaGas   = disk.SigmaGas_up;
                disk.calcrhoGas(disk.SigmaGas,disk.hscale,rz2);
                // get disk gas velocity
                disk.alpha_gas_c = disk.alpha_gas;
                disk.beta_gas_c = disk.beta_gas;
                PS::F64 hscae2 = disk.hscale*disk.hscale;
                PS::F64 hscale_rxy_2 = hscae2 *rxy_inv *rxy_inv;
                PS::F64 rz_hscale_2  = rz2/hscae2;
                PS::F64vec ev(-pp[i].pos.y*rxy_inv, pp[i].pos.x*rxy_inv, 0.0);
                PS::F64vec vkep = sqrt(FP_t::m_sun * rxy_inv) * ev; // FP_t::m_sun is the mass of central star
                disk.calcUvel(hscale_rxy_2,rz_hscale_2,disk.alpha_gas_c,disk.beta_gas_c,pp[i].vel,vkep);
                //
                fdrag_aero  = calcAerodynamicGasDrag_Adachi1976(disk.uvel,pp[i].mass,pp[i].r_planet,disk.rhoGas,0.5e0*disk.Cd*m_pi);
                fdrag   =   fdrag_aero;
                break;
            }
            case 1: // tidal gas drag by Kanagawa+2018
            {
                // surface density at gap bottom
                PS::F64 m_ms = pp[i].mass/FP_t::m_sun;
                PS::F64 r_hscale = rxy/disk.hscale;
                disk.calcSigmaGas_fgap(m_ms, r_hscale, disk.alpha_vis);
                disk.SigmaGas   = disk.SigmaGas_fgap *disk.SigmaGas_up;
                //
                PS::F64 tau_tide_axi = calc_TauTideKanagawa2018(pp[i].mass, FP_t::m_sun, rxy, disk.SigmaGas, disk.hscale);
                // std::cout << tau_tide_axi << std::endl;
                fdrag_aero  = 0.;
                fdrag_tide  += calcTidalDragAxi_usingTimescale(pp[i].vel,tau_tide_axi);
                fdrag   =   fdrag_tide;
                break;
            }
            case 2: // tidal gas drag by Tanaka+2002, Tanaka & Ward 2004, Cresswell & Nelson 2008 with reduction by Kanagawa+2018
            {
                // surface density at gap bottom
                PS::F64 m_ms = pp[i].mass/FP_t::m_sun;
                PS::F64 r_hscale = rxy/disk.hscale;
                disk.calcSigmaGas_fgap(m_ms, r_hscale, disk.alpha_vis);
                disk.SigmaGas   = disk.SigmaGas_fgap *disk.SigmaGas_up;
                disk.alpha_gas_c = disk.alpha_gas;
                //
                PS::F64 tau_wave     = calc_TauWave_Tanaka2002(pp[i].mass, FP_t::m_sun, rxy, disk.SigmaGas, disk.hscale);
                PS::F64 tau_tide_axi = calc_TauTideAxi_Cresswell2008(tau_wave, r_hscale, pp[i].getEccentricity(), pp[i].getInclination(), disk.alpha_gas_c);
                PS::F64 tau_tide_ecc = calc_TauTideEcc_Cresswell2008(tau_wave, r_hscale, pp[i].getEccentricity(), pp[i].getInclination());
                PS::F64 tau_tide_inc = calc_TauTideInc_Cresswell2008(tau_wave, r_hscale, pp[i].getEccentricity(), pp[i].getInclination());

                fdrag_tide  += calcTidalDragAxi_usingTimescale(pp[i].vel,tau_tide_axi);
                fdrag_tide  += calcTidalDragEcc_usingTimescale(pp[i].pos, pp[i].vel, tau_tide_ecc);
                fdrag_tide  += calcTidalDragInc_usingTimescale(pp[i].vel.z, tau_tide_inc);
                fdrag   =   fdrag_tide;
                break;
            }
            default:
            {
                break;
            }
        }
        pp[i].acc_gd = 0.;
        pp[i].acc_gd += fdrag;
        pp[i].acc += pp[i].acc_gd;
    }
}

PS::F64vec DragForce::calcAerodynamicGasDrag_Adachi1976
    (   PS::F64vec u,
        PS::F64 mass,
        PS::F64 radius,
        PS::F64 rho_gas,
        PS::F64 coeff // pi * Cd / 2
    ) 
{
    return -coeff *radius *radius *rho_gas *sqrt(u*u) *u /mass;
}


PS::F64vec DragForce::calcTidalDragAxi_usingTimescale
    (   PS::F64vec vp,
        PS::F64 tau_axi
    ) 
{
    return -0.5*vp/tau_axi;
}


PS::F64vec DragForce::calcTidalDragEcc_usingTimescale
    (   PS::F64vec rp,
        PS::F64vec vp,
        PS::F64 tau_ecc
    ) 
{
    return -2.*(vp*rp)*rp/((rp*rp)*tau_ecc);
}

PS::F64vec DragForce::calcTidalDragInc_usingTimescale
    (   PS::F64 vz,
        PS::F64 tau_inc
    ) 
{
    PS::F64vec veck =0.;
    veck.z=1.;
    return -vz/tau_inc *veck ;
}


PS::F64 DragForce::calc_TauTideKanagawa2018(PS::F64 m, PS::F64 ms, PS::F64 r, PS::F64 Sigma, PS::F64 hs)
{
    return ( pow(ms,1.5)*hs*hs )/(2.*c_mig*m*Sigma*pow(r,2.5));
}

PS::F64 DragForce::calc_TauWave_Tanaka2002(PS::F64 m, PS::F64 ms, PS::F64 r, PS::F64 Sigma, PS::F64 hs)
{
    return ( pow(ms,1.5)*hs*hs*hs*hs )/(m*Sigma*pow(r,4.5));
}

PS::F64 DragForce::calc_TauTideAxi_Cresswell2008(PS::F64 twave, PS::F64 r_h, PS::F64 e, PS::F64 i, PS::F64 alpha)
{
    PS::F64 e_rh = e*r_h;
    PS::F64 i_rh = i*r_h;
    PS::F64 Pe = (1.+pow(0.444*e_rh,1.2) +pow(0.352*e_rh,6.))/(1.-pow(0.495*e_rh,4.));
    // Note: factor 2 is not in this equation. In calcDrag
    return twave/(2.7+1.1*alpha) *r_h*r_h 
            *( Pe+Pe/abs(Pe)*  
                (
                    0.070*i_rh +0.085*i_rh*i_rh*i_rh*i_rh -0.080*e_rh*i_rh*i_rh
                )
            );
}

PS::F64 DragForce::calc_TauTideEcc_Cresswell2008(PS::F64 twave, PS::F64 r_h, PS::F64 e, PS::F64 i)
{
    PS::F64 e_rh = e*r_h;
    PS::F64 i_rh = i*r_h;
    return 1.28*twave
            *(
                1.-0.14*e_rh*e_rh +0.06*e_rh*e_rh*e_rh +0.18*e_rh*i_rh*i_rh
            );
}

PS::F64 DragForce::calc_TauTideInc_Cresswell2008(PS::F64 twave, PS::F64 r_h, PS::F64 e, PS::F64 i)
{
    PS::F64 e_rh = e*r_h;
    PS::F64 i_rh = i*r_h;
    return 1.84*twave
            *(
                1.-0.30*i_rh*i_rh +0.24*i_rh*i_rh*i_rh +0.14*e_rh*e_rh*i_rh
            );
}

PS::S32 DragForce::setting = 0;






