#include<iostream>
#include<fstream>
#include<unistd.h>
#include<sys/stat.h>
#include<particle_simulator.hpp>
#ifdef ENABLE_PHANTOM_GRAPE_X86
#include <gp5util.h>
#endif
#ifdef ENABLE_GPU_CUDA
#define MULTI_WALK
#include"force_gpu_cuda.hpp"
#endif
#include "user-defined.hpp"
#include"force_sunway_impl.hpp"
#include"ic.hpp"
#include"collision.hpp"
#include"nbody.hpp"


PS::F64 Epi::eps      = 1.0/32.0;
PS::F64 Epj::r_search = 1.0/32.0;
PS::F64 Epj::r_coll   = 1.0/32.0;
PS::F64 FPGrav::r_coll   = Epj::r_coll;

int main(int argc, char *argv[]) {
    std::cout<<std::setprecision(15);
    std::cerr<<std::setprecision(15);

    PS::Initialize(argc, argv);

    const PS::F64 PI = 4.0 * atan(1.0);
    PS::F32 theta = 0.5;
    PS::S32 n_leaf_limit = 8;
    PS::S32 n_group_limit = 64;
    //PS::F32 dt = 1.0 / 32.0;
    PS::F32 dt = 1.0 / 64.0;
    PS::F32 time_end = dt*128;
    PS::F32 dt_diag = 1.0 / 8.0;
    PS::F32 dt_snap = 1.0;
    char dir_name[1024];
    PS::S64 n_glb = 1024;
    PS::F64 ax_in = 0.45;
    PS::F64 ax_out = 0.55;
    PS::F64 ecc_rms = 2.0;
    PS::F64 inc_rms = 1.0;
    PS::F64 tau = 3.0;
    PS::S32 c;
    sprintf(dir_name,"./result");
    opterr = 0;
    while((c=getopt(argc,argv,"i:o:d:D:t:T:l:n:N:hs:")) != -1){
        switch(c){
        case 'o':
            sprintf(dir_name,optarg);
            break;
        case 't':
            theta = atof(optarg);
            std::cerr << "theta =" << theta << std::endl;
            break;
        case 'T':
            time_end = atof(optarg);
            std::cerr << "time_end = " << time_end << std::endl;
            break;
        case 's':
            dt = atof(optarg);
            std::cerr << "time_step = " << dt << std::endl;
            break;
        case 'd':
            dt_diag = atof(optarg);
            std::cerr << "dt_diag = " << dt_diag << std::endl;
            break;
        case 'D':
            dt_snap = atof(optarg);
            std::cerr << "dt_snap = " << dt_snap << std::endl;
            break;
        case 'l':
            n_leaf_limit = atoi(optarg);
            std::cerr << "n_leaf_limit = " << n_leaf_limit << std::endl;
            break;
        case 'n':
            n_group_limit = atoi(optarg);
            std::cerr << "n_group_limit = " << n_group_limit << std::endl;
            break;
        case 'N':
            n_glb = atoi(optarg);
            std::cerr << "n_glb = " << n_glb << std::endl;
            break;
        case 'h':
            if(PS::Comm::getRank() == 0) {
                printHelp();
            }
            PS::Finalize();
            return 0;
        default:
            if(PS::Comm::getRank() == 0) {
                std::cerr<<"No such option! Available options are here."<<std::endl;
                printHelp();
            }
            PS::Abort();
        }
    }
    makeOutputDirectory(dir_name);
    std::ofstream fout_eng;
    if(PS::Comm::getRank() == 0) {
        char sout_de[1024];
        sprintf(sout_de, "%s/t-de.dat", dir_name);
        fout_eng.open(sout_de);
        fprintf(stdout, "This is a sample program of N-body simulation on FDPS!\n");
        fprintf(stdout, "Number of processes: %d\n", PS::Comm::getNumberOfProc());
        fprintf(stdout, "Number of threads per process: %d\n", PS::Comm::getNumberOfThread());
    }

    Planet planet(1.0, PS::F64vec(0.0), PS::F64vec(0.0));
    const int n_sat = 1;
    Satellite * satellite = new Satellite[n_sat];
    satellite[0].mass = 0.0;
    satellite[0].pos = PS::F64vec(1.0, 0.0, 0.0);
    satellite[0].vel = PS::F64vec(0.0, 1.0, 0.0);

    PS::ParticleSystem<FPGrav> system_grav;
    system_grav.initialize();

    //tau  = n*PI*r_phy^2/(PI*(ax_out^2-ax_in^2))
    PS::F64 r_phy = sqrt(tau*(ax_out*ax_out - ax_in*ax_in) / n_glb);
    PS::F64 ratio_r_phy_r_hill = 1.0;
    PS::F64 r_hill = r_phy / ratio_r_phy_r_hill;
    // r_hill = ax*qbrt(2*mass_dust / (3*mass_pla));
    PS::F64 mass_dust = (r_hill/ax_out)*(r_hill/ax_out)*(r_hill/ax_out)*3.0*planet.mass*0.5;
    PS::F32 time_sys = 0.0;
    PS::F64 power = 0.0;
    PS::S32 seed = 0;
    PS::F64 dens = mass_dust * n_glb / (PI*(ax_out*ax_out - ax_in*ax_in));
    SetParticleKeplerDisk(system_grav, n_glb, ax_in, ax_out, ecc_rms, inc_rms, 
                          dens, planet.mass, 0.0, 1.0, power, seed);
    PS::S32 n_loc    = system_grav.getNumberOfParticleLocal();
    //for(PS::S32 i=0; i<n_loc; i++) system_grav[i].r_search = r_hill * 10.0;
    Epj::r_search = r_hill * 10.0;
    Epj::r_coll = r_phy;
    std::cerr<<"r_hill= "<<r_hill
             <<" mass_dust_total= "<<mass_dust * n_glb<<std::endl;

    satellite[0].r_coll = Epj::r_coll * pow((satellite[0].mass/mass_dust), 1.0/3.0);



    const PS::F32 coef_ema = 0.3;
    PS::DomainInfo dinfo;
    dinfo.initialize(coef_ema);
    dinfo.decomposeDomainAll(system_grav);
    system_grav.exchangeParticle(dinfo);
    n_loc = system_grav.getNumberOfParticleLocal();
    for(PS::S32 i=0; i<n_loc; i++) system_grav[i].rank_org = PS::Comm::getRank();


#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_open();
    g5_set_eps_to_all(FPGrav::eps);
#endif


    ///////////////
    /// chose TREE TYPE    
    //PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::Monopole tree_grav;
    //typedef PS::SPJMonopole SPJ;
    //PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::MonopoleWithScatterSearch tree_grav;
    typedef Epi EPI;
    typedef Epj EPJ;
    typedef Force FORCE;
    typedef PS::SPJMonopoleScatter SPJ;
    PS::TreeForForceLong<FORCE, EPI, EPJ>::MonopoleWithScatterSearch tree_grav;
    //PS::TreeForForceLong<FPGrav, FPGrav, FPGrav>::MonopoleWithSymmetrySearch2 tree_grav;
    //typedef PS::SPJMonopoleInAndOut SPJ;

    tree_grav.initialize(n_glb, theta, n_leaf_limit, n_group_limit);
    std::cerr<<"check a"<<std::endl;
#if 1
    tree_grav.calcForceAllAndWriteBackReuseListMultiWalk(DispatchKernelWithSP<EPI, EPJ, SPJ>,
                                                         RetrieveKernel<FORCE>,
                                                         0,
                                                         system_grav,
                                                         dinfo,
                                                         N_WALK_LIMIT,
                                                         true,
                                                         false);
#else
    tree_grav.calcForceAllAndWriteBackReuseList(CalcGravity<FPGrav>,
                                                CalcGravity<SPJ>,
                                                system_grav,
                                                dinfo,
                                                true,
                                                false);
#endif
    std::cerr<<"check b"<<std::endl;
    PS::ReallocatableArray<Force> force_dust_from_satellite;
    force_dust_from_satellite.resizeNoInitialize(n_loc);
    for(PS::S32 i=0; i<n_sat; i++) satellite[i].clear();
    CalcForceFromPlanet(system_grav, satellite, n_sat, planet);
    CalcForceBetweenDustAndSatellite(system_grav, satellite, force_dust_from_satellite.getPointer(), n_sat);
    bool flag_regular_ic = true;
    if(flag_regular_ic){
        PS::S64 n_coll_loc = 0;
        for(PS::S32 i=0; i<n_loc; i++){
            n_coll_loc += tree_grav.getForce(i).n_coll;
        }
        PS::S64 n_coll_glb = PS::Comm::getSum(n_coll_loc);
        std::cerr<<"n_coll_glb= "<<n_coll_glb<<std::endl;
        /*
        while(n_coll_glb != 0){
            Collision::correction(system_grav, tree_grav, 1.0, 1.0, 1.0+1e-10);
            tree_grav.calcForceAllAndWriteBackReuseListMultiWalk(DispatchKernelWithSP<EPI, EPJ, SPJ>,
                                                                 RetrieveKernel<FORCE>,
                                                                 0,
                                                                 system_grav,
                                                                 dinfo,
                                                                 N_WALK_LIMIT,
                                                                 true,
                                                                 true);
            for(PS::S32 i=0; i<n_sat; i++) satellite[i].clear();
            force_dust_from_satellite.resizeNoInitialize(n_loc);
            CalcForceFromPlanet(system_grav, satellite, n_sat, planet);
            CalcForceBetweenDustAndSatellite(system_grav, satellite, force_dust_from_satellite.getPointer(), n_sat);
            n_coll_loc = 0;
            for(PS::S32 i=0; i<n_loc; i++){
                n_coll_loc += tree_grav.getForce(i).n_coll;
            }
            n_coll_glb = PS::Comm::getSum(n_coll_loc);
            std::cerr<<"n_coll_glb= "<<n_coll_glb<<std::endl;
        }
        */

    }


    PS::F64 Epot0, Ekin0, Etot0, Epot1, Ekin1, Etot1;
    calcEnergy(system_grav, satellite, n_sat, Etot0, Ekin0, Epot0);
    std::cerr<<"Ekin0= "<<Ekin0<<" Epot0= "<<Epot0<<" Etot0= "<<Etot0<<std::endl;
    PS::F64 time_diag = 0.0;
    PS::F64 time_snap = 0.0;
    PS::S64 n_loop = 0;
    PS::S32 id_snap = 0;


#if 0
    while(time_sys < time_end){
        if( (time_sys >= time_snap) || ( (time_sys + dt) - time_snap ) > (time_snap - time_sys) ){
            char filename[256];
            sprintf(filename, "%s/%04d.dat", dir_name, id_snap++);
            FileHeader header;
            header.time   = time_sys;
            header.n_body = system_grav.getNumberOfParticleGlobal();
            system_grav.writeParticleAscii(filename, header);
            time_snap += dt_snap;
        }

        calcEnergy(system_grav, satellite, n_sat, Etot1, Ekin1, Epot1);
        
        if(PS::Comm::getRank() == 0){
            if( (time_sys >= time_diag) || ( (time_sys + dt) - time_diag ) > (time_diag - time_sys) ){
                fout_eng << time_sys << "   " << (Etot1 - Etot0) / Etot0 << std::endl;
                fprintf(stdout, "time: %10.7f energy error: %+e\n",
                        time_sys, (Etot1 - Etot0) / Etot0);
                time_diag += dt_diag;
            }            
        }


        time_sys += dt;
#if 0

        if(PS::Comm::getRank() == 0){
            for(PS::S32 i=0; i<10; i++){
                std::cerr<<"i= "<<i<<" system_grav[i].acc= "<<system_grav[i].acc<<std::endl;
            }
        }
        for(PS::S32 i=0; i<system_grav.getNumberOfParticleLocal(); i++){
            system_grav[i].mass *= 2.0;
        }

        tree_grav.calcForceAllAndWriteBackReuseList(CalcGravity<FPGrav>,
                                                    CalcGravity<SPJ>,
                                                    system_grav,
                                                    dinfo,
                                                    true,
                                                    true);
#elif 0
        kick(system_grav, dt * 0.5);
        drift(system_grav, dt);
        dinfo.decomposeDomainAll(system_grav);
        system_grav.exchangeParticle(dinfo);
        tree_grav.calcForceAllAndWriteBackReuseList(CalcGravity<FPGrav>,
                                                    CalcGravity<SPJ>,
                                                    system_grav,
                                                    dinfo,
                                                    true,
                                                    false);
        kick(system_grav, dt * 0.5);
#elif 0
        kick(system_grav, dt * 0.5);
        drift(system_grav, dt);
        tree_grav.calcForceAllAndWriteBackReuseList(CalcGravity<FPGrav>,
                                                    CalcGravity<SPJ>,
                                                    system_grav,
                                                    dinfo,
                                                    true,
                                                    true);
        kick(system_grav, dt * 0.5);
#elif 0
        kick(system_grav, dt * 0.5);
        kick(satellite, n_sat, dt*0.5);
        drift(system_grav, dt);
        drift(satellite, n_sat, dt*0.5);
        dinfo.decomposeDomainAll(system_grav);
        system_grav.exchangeParticle(dinfo);
        tree_grav.calcForceAllAndWriteBackReuseListMultiWalk(DispatchKernelWithSP<EPI, EPJ, SPJ>,
                                                             RetrieveKernel<FORCE>,
                                                             0,
                                                             system_grav,
                                                             dinfo,
                                                             N_WALK_LIMIT,
                                                             true,
                                                             false);
        //CalcForceFromPlanet(system_grav, planet);
        for(PS::S32 i=0; i<n_sat; i++) satellite[i].clear();
        CalcForceFromPlanet(system_grav, satellite, n_sat, planet);
        CalcForceBetweenDustAndSatellite(system_grav, satellite, n_sat);
        kick(system_grav, dt * 0.5);
        kick(satellite, n_sat, dt*0.5);
#elif 1
        kick(system_grav, dt * 0.5);
        kick(satellite, n_sat, dt*0.5);
        drift(system_grav, dt);
        drift(satellite, n_sat, dt);

        bool reuse_flag = true;
        if((n_loop+1) % 10 == 0 ){
            dinfo.decomposeDomainAll(system_grav);
            system_grav.exchangeParticle(dinfo);
            n_loc = system_grav.getNumberOfParticleLocal();
            reuse_flag = false;
        }
        tree_grav.calcForceAllAndWriteBackReuseListMultiWalk(DispatchKernelWithSP<EPI, EPJ, SPJ>,
                                                             RetrieveKernel<FORCE>,
                                                             0,
                                                             system_grav,
                                                             dinfo,
                                                             N_WALK_LIMIT,
                                                             true,
                                                             reuse_flag);
        //CalcForceFromPlanet(system_grav, planet);
        for(PS::S32 i=0; i<n_sat; i++) satellite[i].clear();
        force_dust_from_satellite.resizeNoInitialize(n_loc);
        CalcForceFromPlanet(system_grav, satellite, n_sat, planet);
        CalcForceBetweenDustAndSatellite(system_grav, satellite, force_dust_from_satellite.getPointer(), n_sat);
        if(PS::Comm::getRank() == 0){
            for(PS::S32 i=0; i<n_loc; i++){
                if(tree_grav.getForce(i).n_coll > 0){
                    std::cout<<"id= "<<system_grav[i].id<<" n_coll= "<<tree_grav.getForce(i).n_coll<<std::endl;
                }
            }
        }
        kick(system_grav, dt * 0.5);
        kick(satellite, n_sat, dt*0.5);
        Collision::correction(system_grav, tree_grav);
#else
        // IT WORKS WELL
        // just rotate along z axis
        PS::F64 rotation_angle = 4.3789; // some random number
        PS::F64 cos_phi = cos(rotation_angle);
        PS::F64 sin_phi = sin(rotation_angle);
        for(PS::S32 i=0; i<n_loc; i++){
            PS::F64 x  = system_grav[i].pos.x;
            PS::F64 y  = system_grav[i].pos.y;
            PS::F64 vx = system_grav[i].vel.x;
            PS::F64 vy = system_grav[i].vel.y;
            system_grav[i].pos.x = x*cos_phi  - y*sin_phi;
            system_grav[i].pos.y = x*sin_phi  + y*cos_phi;
            system_grav[i].vel.x = vx*cos_phi - vy*sin_phi;
            system_grav[i].vel.y = vx*sin_phi + vy*cos_phi;
        }
        for(PS::S32 i=0; i<n_sat; i++){
            PS::F64 x  = satellite[i].pos.x;
            PS::F64 y  = satellite[i].pos.y;
            PS::F64 vx = satellite[i].vel.x;
            PS::F64 vy = satellite[i].vel.y;
            satellite[i].pos.x = x*cos_phi  - y*sin_phi;
            satellite[i].pos.y = x*sin_phi  + y*cos_phi;
            satellite[i].vel.x = vx*cos_phi - vy*sin_phi;
            satellite[i].vel.y = vx*sin_phi + vy*cos_phi;
        }
        tree_grav.calcForceAllAndWriteBackReuseListMultiWalk(DispatchKernelWithSP<EPI, EPJ, SPJ>,
                                                             RetrieveKernel<FORCE>,
                                                             0,
                                                             system_grav,
                                                             dinfo,
                                                             N_WALK_LIMIT,
                                                             true,
                                                             true);
        for(PS::S32 i=0; i<n_sat; i++) satellite[i].clear();
        CalcForceFromPlanet(system_grav, satellite, n_sat, planet);
        CalcForceBetweenDustAndSatellite(system_grav, satellite, n_sat);
#endif
        n_loop++;
    }
    
#ifdef ENABLE_PHANTOM_GRAPE_X86
    g5_close();
#endif

#endif

    PS::Finalize();
    return 0;
}

