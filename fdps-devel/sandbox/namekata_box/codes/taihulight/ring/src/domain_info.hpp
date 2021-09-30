#pragma once

#include<iostream>
#include<functional>
#include<algorithm>
#include<cmath> // To detect NaN
#include<string>
#include<sstream>

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include<mpi.h>
#endif

namespace ParticleSimulator{
    struct CompYDir{
        bool operator() (const F64vec & left, const F64vec & right) const {
            return left.y < right.y;
        }
        bool operator() (const F64vec & left, const F64 right) const {
            return left.y < right;
        }
        bool operator() (const F64 left, const F64vec & right) const {
            return left < right.y;
        }
    };
    
    template<S32 DIM>
    void SetNumberOfDomainMultiDimension(S32 np[], S32 rank[]){
        for(S32 i=0; i<DIMENSION_LIMIT; i++){
            np[i] = 1;
            rank[i] = 1;
        }
        std::vector<S32> npv;
        npv.resize(DIM);
        S32 np_tmp = Comm::getNumberOfProc();
        for(S32 d=DIM, cid=0; cid<DIM-1; d--, cid++){
            S32 tmp = (S32)pow(np_tmp+0.000001, (1.0/d)*1.000001 );
            while(np_tmp%tmp){
                tmp--;
            }
            npv[cid] = tmp;
            np_tmp /= npv[cid];
        }
        npv[DIM-1] = np_tmp;
        S32 rank_tmp = Comm::getRank();
        std::sort(npv.begin(), npv.end(), std::greater<S32>());
        for(S32 i=DIM-1; i>=0; i--){
            np[i] = npv[i];
            rank[i] = rank_tmp % np[i];
            rank_tmp /= np[i];
        }
    }

    class DomainInfo{
    public:
        F64ort * pos_domain_;
    private:

        TimeProfile time_profile_;

        F64vec * pos_sample_tot_;
        F64vec * pos_sample_loc_;
	
        //F64ort * pos_domain_;
        F64ort * pos_domain_temp_;

        S32 * n_smp_array_;
        S32 * n_smp_disp_array_;

        F32 coef_ema_;
        S32 target_number_of_sample_particle_;
        S32 number_of_sample_particle_tot_;
        S32 number_of_sample_particle_loc_;
        S32 n_domain_[DIMENSION_LIMIT]; // in 2-dim, n_domain_[2] is always 1.

        F64ort pos_root_domain_;
                
        bool first_call_by_initialize;
        bool first_call_by_decomposeDomain;

        S32 boundary_condition_;
        bool periodic_axis_[DIMENSION_LIMIT]; // in 2-dim, periodic_axis_[2] is always false.

#ifdef USE_SUPER_DOMAIN
        ReallocatableArray<F64ort> pos_super_domain_; // size is n_domain in phi direction
        ReallocatableArray<F64ort> pos_sub_domain_; // size is n_domain in phi direction
#endif

        
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
        // NEW
        MPI_Comm comm_1d_[DIMENSION_LIMIT];
        MPI_Comm comm_sub_[DIMENSION_LIMIT];
        int rank_1d_[DIMENSION_LIMIT];
        int rank_sub_[DIMENSION_LIMIT];
        int n_proc_sub_[DIMENSION_LIMIT];
#endif
	/*	
        void calculateBoundaryOfDomain(const S32 np,
                                       const F64vec pos_sample[],
                                       const S32 cid,
                                       const S32 istart,
                                       const S32 iend,
                                       F64 & xlow,
                                       F64 & xhigh) {
            if(istart == 0) {
                xlow  = pos_root_domain_.low_[cid];
            }
	    else {
                xlow  = 0.5 * (pos_sample[istart-1][cid] + pos_sample[istart][cid]);
            }
            if(iend == np - 1) {
                xhigh = pos_root_domain_.high_[cid];
            }
	    else {
                xhigh = 0.5 * (pos_sample[iend][cid] + pos_sample[iend+1][cid]);
            }
        }
	*/
	
        void calculateBoundaryOfDomainX(const S32 np,
					const F64vec pos_sample[],
					const S32 istart,
					const S32 iend,
					F64 & xlow,
					F64 & xhigh) {
            if(istart == 0) xlow  = pos_root_domain_.low_.x;
	    else xlow  = 0.5 * (pos_sample[istart-1].x + pos_sample[istart].x);
            if(iend == np - 1) xhigh = pos_root_domain_.high_.x;
	    else xhigh = 0.5 * (pos_sample[iend].x + pos_sample[iend+1].x);
        }
        void calculateBoundaryOfDomainY(const S32 np,
					const F64vec pos_sample[],
					const S32 istart,
					const S32 iend,
					F64 & xlow,
					F64 & xhigh) {
            if(istart == 0) xlow  = pos_root_domain_.low_.y;
	    else xlow  = 0.5 * (pos_sample[istart-1].y + pos_sample[istart].y);
            if(iend == np - 1) xhigh = pos_root_domain_.high_.y;
	    else xhigh = 0.5 * (pos_sample[iend].y + pos_sample[iend+1].y);
        }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION	
        void calculateBoundaryOfDomainZ(const S32 np,
					const F64vec pos_sample[],
					const S32 istart,
					const S32 iend,
					F64 & xlow,
					F64 & xhigh) {
            if(istart == 0) xlow  = pos_root_domain_.low_.z;
	    else xlow  = 0.5 * (pos_sample[istart-1].z + pos_sample[istart].z);
            if(iend == np - 1) xhigh = pos_root_domain_.high_.z;
	    else xhigh = 0.5 * (pos_sample[iend].z + pos_sample[iend+1].z);
        }	
#endif
    public:

        TimeProfile getTimeProfile() const {
            return time_profile_;
        }
        void clearTimeProfile(){
            time_profile_.clear();
        }
        DomainInfo() {
            first_call_by_initialize = true;
            first_call_by_decomposeDomain = true;
	    periodic_axis_[0] = periodic_axis_[1] = false;
	    pos_root_domain_.low_.x  = -LARGE_FLOAT;
	    pos_root_domain_.high_.x = LARGE_FLOAT;
	    pos_root_domain_.low_.y  = -LARGE_FLOAT;
	    pos_root_domain_.high_.y = LARGE_FLOAT;
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    periodic_axis_[2] = false;
	    pos_root_domain_.low_.z  = -LARGE_FLOAT;
	    pos_root_domain_.high_.z = LARGE_FLOAT;
#endif
            boundary_condition_ = BOUNDARY_CONDITION_OPEN;
        }

        void initialize(const F32 coef_ema = 1.0){
            if( coef_ema < 0.0 || coef_ema > 1.0){
                PARTICLE_SIMULATOR_PRINT_ERROR("The smoothing factor of an exponential moving average is must between 0 and 1.");
                std::cerr<<"The smoothing factor of an exponential moving average is must between 0 and 1."<<std::endl;
                Abort(-1);
            }
            assert(first_call_by_initialize);
            first_call_by_initialize = false;
            pos_sample_tot_ = NULL;
            pos_sample_loc_ = NULL;
	    
            n_smp_array_ = new S32[Comm::getNumberOfProc()];
            n_smp_disp_array_ = new S32[Comm::getNumberOfProc() + 1];

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_      = new F64ort[Comm::getNumberOfProc()];
            pos_domain_temp_ = new F64ort[Comm::getNumberOfProc()];
#else
            pos_domain_     = new F64ort[1];
            pos_domain_temp_= new F64ort[1];
#endif

            coef_ema_ = coef_ema;
            target_number_of_sample_particle_ = 0;
            number_of_sample_particle_tot_ = 0;
            number_of_sample_particle_loc_ = 0;

            //S32 rank_tmp[DIMENSION];
	    S32 rank_tmp[DIMENSION_LIMIT];
            SetNumberOfDomainMultiDimension<DIMENSION>(n_domain_, rank_tmp);

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            // NEW
            int rank_glb = Comm::getRank();
            for(S32 d=DIMENSION-1; d>=0; d--){
                rank_1d_[d] = rank_glb % n_domain_[d];
                rank_glb /= n_domain_[d];
                MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_+d);
                MPI_Comm_rank(comm_sub_[d], rank_sub_+d);
                MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_+d);
                MPI_Comm_size(comm_sub_[d], n_proc_sub_+d);
            }
	    for(S32 d=DIMENSION-1; d>=0; d--){
		Comm::setRankMultiDim(d, rank_tmp[d]);
		Comm::setNumberOfProcMultiDim(d, n_domain_[d]);
	    }
#endif
        }

        void setNumberOfDomainMultiDimension(const S32 nx, const S32 ny, const S32 nz=1){
            S32 n_proc = Comm::getNumberOfProc();
            if(n_proc != nx*ny*nz){
                PARTICLE_SIMULATOR_PRINT_ERROR("devided number of domains is not consistent with total processe number");
                Abort(-1);
            }
            n_domain_[0] = nx;
            n_domain_[1] = ny;
#ifdef PARTICLE_SIMULATOR_TWO_DIMENSION
            n_domain_[2] = 1;
#else
            n_domain_[2] = nz;
#endif
	    
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
            int rank_glb = Comm::getRank();
	    // make comunicator
            for(S32 d=DIMENSION-1; d>=0; d--){
                rank_1d_[d] = rank_glb % n_domain_[d];
                rank_glb /= n_domain_[d];
                MPI_Comm_split(MPI_COMM_WORLD, rank_1d_[d], rank_glb, comm_sub_+d);
                MPI_Comm_rank(comm_sub_[d], rank_sub_+d);
                MPI_Comm_split(MPI_COMM_WORLD, rank_sub_[d], rank_glb, comm_1d_+d);
                MPI_Comm_size(comm_sub_[d], n_proc_sub_+d);
            }
	    int rank_tmp = Comm::getRank();
	    for(S32 d=DIMENSION-1; d>=0; d--){
		Comm::setRankMultiDim(d, rank_tmp%n_domain_[d]);
		rank_tmp /= n_domain_[d];
	    }
	    Comm::setNumberOfProcMultiDim(0, nx);
	    Comm::setNumberOfProcMultiDim(1, ny);
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    Comm::setNumberOfProcMultiDim(2, ny);
#endif //PARTICLE_SIMULATOR_TWO_DIMENSION
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
        }

        void setDomain(const S32 nx, const S32 ny, const S32 nz=1){
            setNumberOfDomainMultiDimension(nx, ny, nz);
        }


        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const bool clear,
                                   const F32 weight) {
            F64 time_offset = GetWtime();
            Comm::barrier();
            if(psys.getFirstCallByDomainInfoCollectSampleParticle()) {
                F64vec *temp_loc = new F64vec[target_number_of_sample_particle_];
                Comm::barrier();
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
                    temp_loc[i] = pos_sample_loc_[i];                
                Comm::barrier();
                target_number_of_sample_particle_ += psys.getTargetNumberOfSampleParticle();
                delete [] pos_sample_tot_;
                delete [] pos_sample_loc_;
                Comm::barrier();
                pos_sample_tot_ = new F64vec[target_number_of_sample_particle_];
                Comm::barrier();
                pos_sample_loc_ = new F64vec[target_number_of_sample_particle_];
                for(S32 i = 0; i < number_of_sample_particle_loc_; i++)
                    pos_sample_loc_[i] = temp_loc[i];
                delete [] temp_loc;
            }
            if(clear) {
                number_of_sample_particle_loc_ = 0;
            }
            S32 number_of_sample_particle = 0;
            psys.getSampleParticle(number_of_sample_particle, &pos_sample_loc_[number_of_sample_particle_loc_], weight);
            number_of_sample_particle_loc_ += number_of_sample_particle;
            time_profile_.collect_sample_particle += GetWtime() - time_offset;
            return;
        }


        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys,
                                   const bool clear){
            const F32 wgh = psys.getNumberOfParticleLocal();
            collectSampleParticle(psys, clear, wgh);
        }

        template<class Tpsys>
        void collectSampleParticle(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
        }

	/*
	void decomposeDomainHierarchy(){
	    F64 time_offset = GetWtime();
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#else //PARTICLE_SIMULATOR_MPI_PARALLEL
	    static bool first = true;
	    static F64vec * pos_sample_loc_tmp;
            if(first){
		pos_sample_loc_tmp = new F64vec[target_number_of_sample_particle_];
	    }
	    const S32 nynz = n_domain_[1] * n_domain_[2];
	    S32 n_smp = number_of_sample_particle_loc_ / nynz;
	    for(S32 i=0; i<n_smp; i++){
		pos_sample_loc_tmp[i] = pos_sample_loc_[i*nynz];
	    }
#ifdef __HPC_ACE__
            Comm::allGather(&n_smp, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            Comm::allGatherV(pos_sample_loc_tmp, n_smp, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#else //__HPC_ACE__
            Comm::gather(&n_smp, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            Comm::gatherV(pos_sample_loc_tmp, n_smp, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#endif //__HPC_ACE__
            if(myrank == 0) {
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
		std::sort(pos_sample_loc_tmp, pos_sample_loc_tmp+number_of_sample_particle_tot_, LessOPX());
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
		    F64 x0, x1;
                    calculateBoundaryOfDomain(number_of_sample_particle_tot_, pos_sample_tot_, 0, istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        //pos_domain_temp_[i].low_[0]  = x0;
                        //pos_domain_temp_[i].high_[0] = x1;
                        pos_domain_temp_[i].low_.x  = x0;
                        pos_domain_temp_[i].high_.x = x1;			
                    }
                }
		delete [] istart;
		delete [] iend;
	    }
	    n_smp = number_of_sample_particle_loc_ / n_domain_[2];
            if(rank_sub_[0] == 0) {
#ifdef __HPC_ACE__
		Comm::allGather(&n_smp, 1, n_smp_array_);
		n_smp_disp_array_[0] = 0;
		for(S32 i=0; i<nproc; i++){
		    n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
		}
		Comm::allGatherV(pos_sample_loc_tmp, n_smp, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
		number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#else //__HPC_ACE__
		Comm::gather(&n_smp, 1, n_smp_array_);
		n_smp_disp_array_[0] = 0;
		for(S32 i=0; i<nproc; i++){
		    n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
		}
		Comm::gatherV(pos_sample_loc_tmp, n_smp, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
		number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#endif //__HPC_ACE__
	    }
#endif //PARTICLE_SIMULATOR_MPI_PARALLEL
	    time_profile_.decompose_domain = GetWtime() - time_offset;
	}
	*/
	

#if 1 //UNDER_CONSTRUCTION
	// new version multi-dimensional gathering
        void decomposeDomainMultiStep() {
            F64 time_offset = GetWtime();
	    //assert(!first_call_by_decomposeDomain);
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL // ifNdef
            pos_domain_[0] = pos_root_domain_;
#else
            static bool first = true;
            static S32 * n_send;
            static S32 * n_recv;
            static S32 * n_send_disp;
            static S32 * n_recv_disp;
            static S32 * i_head;
            static S32 * i_tail;
            static F64vec * pos_sample_buf;
            static F64 *  coord_buf;
            static F64 *  coord_tot;
            static F64 * x_coord;
            static F64 * y_coord;
            static F64ort * pos_domain_temp_buf;
            S32 n_proc_glb = Comm::getNumberOfProc();
            //S32 rank_glb = Comm::getRank();

            if(first){
                n_send = new S32[n_proc_glb];
                n_recv = new S32[n_proc_glb];
                n_send_disp = new S32[n_proc_glb + 1];
                n_recv_disp = new S32[n_proc_glb + 1];
                i_head = new S32[n_proc_glb];
                i_tail = new S32[n_proc_glb];
                pos_sample_buf = new F64vec[target_number_of_sample_particle_];
                coord_buf = new F64[n_proc_glb * 2];
                coord_tot = new F64[n_proc_glb * 2];
                x_coord = new F64[n_proc_glb + 1];
                y_coord = new F64[n_proc_glb + 1];
                pos_domain_temp_buf = new F64ort[n_proc_glb];
                first = false;
            }

            //std::cout<<"rank_glb="<<rank_glb<<" number_of_sample_particle_loc_="<<number_of_sample_particle_loc_<<std::endl;
            ///////////// sort particles along x direction
            std::sort(pos_sample_loc_, pos_sample_loc_+number_of_sample_particle_loc_, LessOPX());


            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
            S32 id_domain_3d = 0;
            S32 id_domain_x = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                while( pos_domain_[id_domain_3d].high_.x <= pos_sample_loc_[i].x ){
                    id_domain_3d += n_proc_sub_[0];
                    id_domain_x++;
                    //std::cout<<"rank_glb="<<rank_glb<<" id_domain_3d="<<id_domain_3d<<std::endl;
                }
                n_send[id_domain_x]++;
            }
            //for(S32 i=0; i<n_domain_[0]; i++)std::cout<<"rank_glb="<<rank_glb<<" n_send[i]="<<n_send[i]<<std::endl;
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            MPI_Alltoallv(pos_sample_loc_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
            //std::cout<<"rank_glb="<<rank_glb<<" n_recv_disp[n_domain_[0]]="<<n_recv_disp[n_domain_[0]]<<std::endl;

	    
            ///////////// allgather particles in Y-Z plane
            S32 n_send_tmp = n_recv_disp[ n_domain_[0] ]; // # of particle in own cell.
            MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_sub_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            S32 n_par_slab = n_recv_disp[ n_proc_sub_[0] ];
            //std::cout<<"rank_glb="<<rank_glb<<" n_par_slab="<<n_par_slab<<std::endl;
            MPI_Allgatherv(pos_sample_buf,  n_send_tmp, GetDataType<F64vec>(),
                           pos_sample_tot_, n_recv, n_recv_disp, GetDataType<F64vec>(), comm_sub_[0]);
	    
            ///////////// sort particles along x direction again
            std::sort(pos_sample_tot_, pos_sample_tot_+n_par_slab, LessOPX());
	    /*
	    if(Comm::getRank() == 1){
		std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
		for(S32 i=0; i<n_par_slab; i++){
		    std::cout<<"pos_sample_tot_[i]="<<pos_sample_tot_[i]<<std::endl;
		}
	    }
	    */
	    ///////////// determine X coord
            MPI_Allgather(&n_par_slab, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            number_of_sample_particle_tot_ = n_recv_disp[n_domain_[0]];
	    //std::cout<<"rank_glb="<<rank_glb<<" number_of_sample_particle_tot_="<<number_of_sample_particle_tot_<<std::endl;

	    // get index of 
            S32 n_ave = number_of_sample_particle_tot_ / n_domain_[0];
            for(S32 i=0; i<n_domain_[0]; i++){
                i_head[i] = n_ave * i;
                if( i < number_of_sample_particle_tot_ % n_domain_[0]){
                    i_head[i] += i;
                }
                else{
                    i_head[i] += number_of_sample_particle_tot_ % n_domain_[0];
                }
                if(i > 0) i_tail[i-1] = i_head[i] - 1;
            }
            i_tail[n_domain_[0]-1] = number_of_sample_particle_tot_ - 1;
/*
            if(Comm::getRank() == 1){
                for(S32 i=0; i<n_domain_[0]; i++) std::cout<<"i_head[i]="<<i_head[i]<<" i_tail[i]="<<i_tail[i]<<std::endl;
            }
*/
            n_send_tmp = 0; // temporally used
            for(S32 i=0; i<n_domain_[0]; i++){
                if( n_recv_disp[rank_1d_[0]] <= i_head[i] &&  i_head[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_head[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
                if( n_recv_disp[rank_1d_[0]] <= i_tail[i] &&  i_tail[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_tail[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
            }
	    //for(S32 i=0; i<n_send_tmp; i++) std::cout<<"rank_glb="<<rank_glb<<" coord_buf[i]="<<coord_buf[i]<<std::endl;

	    MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(),
                      n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
	    n_recv_disp[0] = 0;
	    for(S32 i=0; i<n_domain_[0]; i++){
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
	    }
	    //std::cout<<"rank_glb="<<rank_glb<<" n_recv_disp[n_domain_[0]]="<<n_recv_disp[n_domain_[0]]<<std::endl;

	    MPI_Allgatherv(coord_buf, n_send_tmp, GetDataType<>(coord_buf[0]),
                       coord_tot, n_recv, n_recv_disp, GetDataType<>(coord_buf[0]), comm_1d_[0]);
	    
	    assert( n_recv_disp[n_domain_[0]] == n_domain_[0]*2);

/*
	    if(rank_glb == 1){
            //for(S32 i=0; i<n_send_tmp; i++) std::cout<<"coord_buf[i]="<<coord_buf[i]<<std::endl;
            for(S32 i=0; i<n_domain_[0]*2; i++) std::cout<<"coord_tot[i]="<<coord_tot[i]<<std::endl;
	    }
*/
	    
	    // size of x_coord_buf is n_domain_[0]+1
	    x_coord[0] = pos_root_domain_.low_.x;
	    x_coord[n_domain_[0]] = pos_root_domain_.high_.x;

	    for(S32 i=1; i<n_domain_[0]; i++){
            x_coord[i] = (coord_tot[i*2] + coord_tot[i*2-1]) * 0.5;
/*
            if(rank_glb == 1){
                std::cout<<"i="<<i<<std::endl;
                std::cout<<"coord_tot[i*2-1]="<<coord_tot[i*2-1]<<std::endl;
                std::cout<<"coord_tot[i*2]="<<coord_tot[i*2]<<std::endl;
                std::cout<<"x_coord[i]="<<x_coord[i]<<std::endl;
            }
*/
	    }

/*
	    if(rank_glb == 1){
            for(S32 i=0; i<n_domain_[0]+1; i++) std::cout<<"x_coord[i]="<<x_coord[i]<<std::endl;
	    }
*/
	    
        ///////////// migrate particles along x direction
        for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
	    id_domain_x = 0;
	    for(S32 i=0; i<n_par_slab; i++){
            while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
            n_send[id_domain_x]++;
        }
	    //for(S32 i=0; i<n_domain_[0]; i++)std::cout<<"rank_glb="<<rank_glb<<" n_send[i]="<<n_send[i]<<std::endl;

            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
	    MPI_Alltoallv(pos_sample_tot_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
	    //std::cout<<"rank_glb="<<rank_glb<<" n_par_slab="<<n_par_slab<<" n_recv_disp[n_domain_[0]]="<<n_recv_disp[n_domain_[0]]<<std::endl;
	    n_par_slab = n_recv_disp[n_domain_[0]];
	    // OK
	    
	    ////////////////////////////////////
            ///////////// determine y corrdinate
	    std::sort(pos_sample_buf, pos_sample_buf+n_par_slab, LessOPY());

            // get index of
	    n_ave = n_par_slab / n_domain_[1];
            for(S32 i=0; i<n_domain_[1]; i++){
                i_head[i] = n_ave * i;
                if( i < n_par_slab % n_domain_[1]) i_head[i] += i;
                else i_head[i] += n_par_slab % n_domain_[1];
                if(i > 0) i_tail[i-1] = i_head[i] - 1;
            }
            i_tail[n_domain_[1]-1] = n_par_slab - 1;
	    /*
	    if(Comm::getRank() == 0){
		for(S32 i=0; i<n_domain_[1]; i++) std::cout<<"i_head[i]="<<i_head[i]<<" i_tail[i]="<<i_tail[i]<<std::endl;
	    }
	    */
	    
	    // size of y_coord is n_domain_[1]+1
	    y_coord[0] = pos_root_domain_.low_.y;
	    y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
        for(S32 i=1; i<n_domain_[1]; i++) y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;

	    //if(Comm::getRank() == 0) for(S32 i=0; i<n_domain_[1]+1; i++) std::cout<<"y_coord[i]="<<y_coord[i]<<std::endl;


#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    ////////////////////////////////////
            ///////////// determine z corrdinate
	    //#pragma omp parallel for
	    for(S32 iy=0; iy<n_domain_[1]; iy++){
            const S32 iy_ptcl_head = i_head[iy];
            const S32 iy_ptcl_tail = i_tail[iy];
            const S32 nz_tot = iy_ptcl_tail - iy_ptcl_head + 1;
            std::sort(pos_sample_buf+iy_ptcl_head, pos_sample_buf+iy_ptcl_head+nz_tot, LessOPZ());
/*
		if(Comm::getRank() == 0){
		    for(S32 i=0; i<nz_tot; i++){
			std::cout<<"pos_sample_buf[iy_ptcl_head+i]="<<pos_sample_buf[iy_ptcl_head+i]<<std::endl;
		    }
		}
*/
            S32 nz_ave_tmp = nz_tot / n_domain_[2];
            pos_domain_temp_buf[iy*n_domain_[2]].low_.z = pos_root_domain_.low_.z;
            pos_domain_temp_buf[(iy+1)*n_domain_[2] - 1].high_.z = pos_root_domain_.high_.z;
            pos_domain_temp_buf[iy*n_domain_[2]].low_.y = y_coord[iy];
            pos_domain_temp_buf[iy*n_domain_[2]].high_.y = y_coord[iy+1];
            for(S32 iz=1; iz<n_domain_[2]; iz++){
                pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.y = y_coord[iy];
                pos_domain_temp_buf[iy*n_domain_[2]+iz].high_.y = y_coord[iy+1];
                S32 iz_tmp = nz_ave_tmp * iz;
                if(iz < nz_tot % n_domain_[2]) iz_tmp += iz;
                else iz_tmp += nz_tot % n_domain_[2];
                F64 z_coord_tmp = (pos_sample_buf[iy_ptcl_head+iz_tmp].z + pos_sample_buf[iy_ptcl_head+iz_tmp-1].z) * 0.5;
/*
		    if(Comm::getRank() == 0){
			std::cout<<"iz_tmp="<<iz_tmp<<" z_coord_tmp="<<z_coord_tmp<<std::endl;
		    }
*/
                pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.z = z_coord_tmp;
                pos_domain_temp_buf[iy*n_domain_[2]+iz-1].high_.z = z_coord_tmp;
            }
	    }
#endif
	    for(S32 i=0; i<n_proc_sub_[0]; i++){
                pos_domain_temp_buf[i].low_.x = x_coord[rank_1d_[0]];
                pos_domain_temp_buf[i].high_.x = x_coord[rank_1d_[0]+1];
	    }

/*
	    if(Comm::getRank() == 4){
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                std::cout<<"pos_domain_temp_buf[i]="<<pos_domain_temp_buf[i]<<std::endl;
            }
	    }
*/
	    //////////////////////////////////////////////
        ///////////// exchange pos_domain_tmp
	    MPI_Allgather(pos_domain_temp_buf, n_proc_sub_[0], GetDataType<F64ort>(),
                      pos_domain_temp_, n_proc_sub_[0], GetDataType<F64ort>(), comm_1d_[0]);

/*
	    if(Comm::getRank() == 0){
            for(S32 i=0; i<Comm::getNumberOfProc(); i++){
                std::cout<<"pos_domain_temp_[i]="<<pos_domain_temp_[i]<<std::endl;
            }
	    }
*/
/*
	    for(S32 i = 0; i < n_proc_glb; i++) {
            pos_domain_[i].low_  = pos_domain_temp_[i].low_;
            pos_domain_[i].high_ = pos_domain_temp_[i].high_;
	    }
*/

	    if(first_call_by_decomposeDomain) {
            first_call_by_decomposeDomain = false;
            for(S32 i = 0; i < n_proc_glb; i++) {
                //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                pos_domain_[i].high_ = pos_domain_temp_[i].high_;
            }
	    } else {
            for(S32 i = 0; i < n_proc_glb; i++) {
                //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                    + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                    + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                
            }
	    }
        //time_profile_.decompose_domain += GetWtime() - time_offset;
#endif // PARTICLE_SIMULATOR_MPI_PARALLEL
	    time_profile_.decompose_domain = GetWtime() - time_offset;
        }
#endif // UNDER_CONSTRUCTION

        void decomposeDomain() {
            static int number_of_dd=0; // for debug
            number_of_dd++;
            Comm::barrier();
            if (Comm::getRank() == 0) 
                std::cout << "number_of_dd = " << number_of_dd << std::endl;
            if (Comm::getRank() == 0) std::cerr<<"check 0"<<std::endl;            
            F64 time_offset = GetWtime();
            // ****** collect sample particles to process 0. ****** 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    S32 nproc  = Comm::getNumberOfProc();
            S32 myrank = Comm::getRank();
//#ifdef __HPC_ACE__
#if 1
            Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            // 2sec 16384 process
            Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];

            if (number_of_dd == 7) {
                std::ofstream fp;
                fp.open("check_loc.dat",std::ios::trunc);
                for (S32 i=0; i<number_of_sample_particle_loc_; i++) {
                    fp << pos_sample_loc_[i].x << "  "
                       << pos_sample_loc_[i].y << "  "
                       << pos_sample_loc_[i].z << std::endl;
                }
                fp.close();
            }
#else //__HPC_ACE__
            Comm::barrier();
            if(Comm::getRank()==0) std::cerr<<"check 1"<<std::endl;
            Comm::gather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            Comm::barrier();
            if(Comm::getRank()==0) std::cerr<<"check 2"<<std::endl;
            // 510sec 16384 process
            Comm::gatherV(pos_sample_loc_, number_of_sample_particle_loc_, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
#endif //__HPC_ACE__
            // ****************************************************
            // *** decompose domain *******************************
            Comm::barrier();
            if(myrank == 0) {
                if(Comm::getRank()==0) std::cerr<<"check 3"<<std::endl;
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
                std::cerr << "check 3-1" << std::endl;
                std::cerr << "nsmpltot = " << number_of_sample_particle_tot_<<std::endl;
                if (number_of_dd == 7) {
                    std::ofstream fp;
                    fp.open("check.dat",std::ios::trunc);
                    for (S32 i=0; i<number_of_sample_particle_tot_; i++) {
                        //if (isnan(pos_sample_tot_->x) ||
                        //    isnan(pos_sample_tot_->y) ||
                        //    isnan(pos_sample_tot_->z)) {
                        //   std::cout << "[NaN] i = 0 " << i << std::endl;
                        //}
                        fp << pos_sample_tot_[i].x << "  "
                           << pos_sample_tot_[i].y << "  "
                           << pos_sample_tot_[i].z << std::endl;
                    }
                    fp.close();
                }
                std::sort(pos_sample_tot_, pos_sample_tot_+number_of_sample_particle_tot_, Cmpvec(&F64vec::x));
                std::cerr<<"check 3-2"<<std::endl;
                
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
                if(Comm::getRank()==0) std::cerr<<"check 4"<<std::endl;
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    F64 x0 = 0.0;
		    F64 x1 = 0.0;
                    //calculateBoundaryOfDomain(number_of_sample_particle_tot_, pos_sample_tot_, 0, istart[ix0], iend[ix1-1], x0, x1);
		    calculateBoundaryOfDomainX(number_of_sample_particle_tot_, pos_sample_tot_, istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        pos_domain_temp_[i].low_.x  = x0;
                        pos_domain_temp_[i].high_.x = x1;
                    }
                }
                if(Comm::getRank()==0) std::cerr<<"check 5"<<std::endl;
                // ------------------------------------------
                // --- y direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[ix0], iend[ix1-1], 1);
		    std::sort(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), Cmpvec(&F64vec::y));
/*
                    for(S32 i=istart[ix0]; i<iend[ix1-1]+1; i++){
                        if(pos_sample_tot_[i+1].y < pos_sample_tot_[i].y){
                            std::cout<<"y sort is wrong: i="<<i<<std::endl;
                            std::cout<<"pos_sample_tot_[i+1]="<<pos_sample_tot_[i+1]<<std::endl;
                            std::cout<<"pos_sample_tot_[i]="<<pos_sample_tot_[i]<<std::endl;
                        }
                    }
*/
                    //if(Comm::getRank()==0) std::cerr<<"check 6"<<std::endl;
                    S32 number_of_sample_particle_tot_y = iend[ix1-1] - istart[ix0] + 1;
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //F64 y0, y1;
			F64 y0 = 0.0;
			F64 y1 = 0.0;
                        //std::cout<<"ix="<<ix<<" ix0="<<ix0<<" ix1="<<ix1<<" iy="<<iy<<" iy0="<<iy0<<" iy1="<<iy1<<" y0="<<y0<<" y1="<<y1<<std::endl;
                        //calculateBoundaryOfDomain(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], 1, istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
			calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                        for(S32 i = iy0; i < iy1; i++) {
                            pos_domain_temp_[i].low_.y  = y0;
                            pos_domain_temp_[i].high_.y = y1;
                        }
                    }
                }
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- z direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    //std::cerr<<"ix= "<<ix<<" n_domain_[0]= "<<n_domain_[0]<<std::endl;
                    S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //if(Comm::getRank()==0) std::cerr<<"check 7"<<std::endl;
                        //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
			std::sort(pos_sample_tot_+istart[iy0], pos_sample_tot_+(iend[iy1-1]+1), Cmpvec(&F64vec::z));
                        S32 number_of_sample_particle_tot_z = iend[iy1-1] - istart[iy0] + 1;
                        for(S32 iz = 0; iz < n_domain_[2]; iz++) {
                            S32 iz0 = iy0 + iz;
                            //F64 z0, z1;
			    F64 z0 = 0.0;
			    F64 z1 = 0.0;						    
                            //calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
			    calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            pos_domain_temp_[iz0].low_.z  = z0;
                            pos_domain_temp_[iz0].high_.z = z1;
                        }
                    }
                }
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- process first ------------------------
                if(first_call_by_decomposeDomain) {
                    first_call_by_decomposeDomain = false;
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                    }
                } else {
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                        
                    }
                }
                // ------------------------------------------
                delete [] istart;
                delete [] iend;
            }
            Comm::barrier();
            // ****************************************************
            // *** broad cast pos_domain_ *************************
            MPI::COMM_WORLD.Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0);
            //std::cout<<"end of bcast: "<<"time: "<<GetWtime() - Tbegin<<std::endl;
            //Comm::broadcast(pos_domain_, nproc);
            // ****************************************************
#else       // PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#endif     // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"pos_root_domain_="<<pos_root_domain_<<std::endl;
            std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
#endif
            time_profile_.decompose_domain += GetWtime() - time_offset;

            if (number_of_dd == 2) {
                if (Comm::getRank() == 0) {
                    std::ofstream fp;
                    std::string filename;
                    std::ostringstream filenum;
                    filenum << std::setfill('0') << std::setw(5) << Comm::getRank();
                    filename = "./dinfo/pos_domain" + filenum.str() + ".dat";
                    fp.open(filename.c_str(),std::ios::trunc);
                        for (S32 i=0; i<nproc; i++) {
                           fp << pos_domain_[i].low_.x << "  "
                              << pos_domain_[i].low_.y << "  "
                              << pos_domain_[i].low_.z << "  "
                              << i << std::endl;
                           fp << pos_domain_[i].high_.x << "  "
                              << pos_domain_[i].high_.y << "  "
                              << pos_domain_[i].high_.z << "  "
                              << i << std::endl;
                        }
                    fp.close();
                    std::cout << "**** pos_domain*.dat output ****." << std::endl;
                }
            }
            Comm::barrier();
        }

        void decomposeDomain2() {
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
            Comm::barrier();
            if(Comm::getRank()==0) std::cerr<<"check 0"<<std::endl;
#endif
            F64 time_offset = GetWtime();
            // ****** collect sample particles to process 0. ****** 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    S32 nproc  = Comm::getNumberOfProc();
            S32 myrank = Comm::getRank();
            assert(n_domain_[1] % 2 == 0);
            Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            // 2sec 16384 process
            Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
            // ****************************************************
            // *** decompose domain *******************************
            if(myrank == 0) {
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 3"<<std::endl;
#endif
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
		std::sort(pos_sample_tot_, pos_sample_tot_+number_of_sample_particle_tot_, Cmpvec(&F64vec::x));
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 4"<<std::endl;
#endif
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    F64 x0 = 0.0;
		    F64 x1 = 0.0;
		    calculateBoundaryOfDomainX(number_of_sample_particle_tot_, pos_sample_tot_, istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        pos_domain_temp_[i].low_.x  = x0;
                        pos_domain_temp_[i].high_.x = x1;
                    }
                }
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 5"<<std::endl;
#endif
                // ------------------------------------------
                // --- y direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
		    std::sort(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), Cmpvec(&F64vec::y));
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                    if(Comm::getRank()==0) std::cerr<<"check 6"<<std::endl;
#endif
                    S32 number_of_sample_particle_tot_y = iend[ix1-1] - istart[ix0] + 1;
                    S32 number_of_sample_particle_sep_y_lower = 0; // index of samples
                    for(S32 i=istart[ix0]; i<iend[ix1-1]+1; i++){
                        if(pos_sample_tot_[i].y < 0.0){
                            number_of_sample_particle_sep_y_lower++;
                        }
                    }
                    S32 number_of_sample_particle_sep_y_upper = number_of_sample_particle_tot_y - number_of_sample_particle_sep_y_lower;
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                    if(Comm::getRank()==0){
                        std::cerr<<"number_of_sample_particle_tot_y= "<<number_of_sample_particle_tot_y<<std::endl;
                        std::cerr<<"number_of_sample_particle_sep_y_lower= "<<number_of_sample_particle_sep_y_lower<<std::endl;
                        std::cerr<<"number_of_sample_particle_sep_y_upper= "<<number_of_sample_particle_sep_y_upper<<std::endl;
                    }
#endif

                    // lower part
                    S32 dn = number_of_sample_particle_sep_y_lower / (n_domain_[1]/2);
                    for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        F64 y0 = 0.0;
                        F64 y1 = 0.0;
                        if(iy==0){
                            y0 = pos_root_domain_.low_.y;
                        }
                        else{
                            //y0 = (pos_sample_tot_[istart[ix0]+dn*(iy-1)].y + pos_sample_tot_[istart[ix0]+dn*iy].y) * 0.5;
                            y0 = (pos_sample_tot_[istart[ix0]+dn*iy].y + pos_sample_tot_[istart[ix0]+dn*iy+1].y) * 0.5;
                        }
                        if(iy==(n_domain_[1]/2-1)){
                            y1 = 0.0;
                        }
                        else{
                            y1 = (pos_sample_tot_[istart[ix0]+dn*(iy+1)].y + pos_sample_tot_[istart[ix0]+dn*(iy+1)+1].y) * 0.5;
                        }
                        for(S32 i = iy0; i < iy1; i++) {
                            pos_domain_temp_[i].low_.y  = y0;
                            pos_domain_temp_[i].high_.y  = y1;
                        }
                        //std::cerr<<"iy= "<<iy<<" y0= "<<y0<<" y1= "<<y1<<std::endl;
                    }
                    //upper part
                    //number_of_sample_particle_sep_y = number_of_sample_particle_tot_y - number_of_sample_particle_sep_y;
                    dn = number_of_sample_particle_sep_y_upper / (n_domain_[1]/2);
                    for(S32 iy = n_domain_[1]/2; iy<n_domain_[1]; iy++) {
                        S32 iy_tmp = iy - n_domain_[1]/2;
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        F64 y0 = 0.0;
                        F64 y1 = 0.0;
                        if(iy==n_domain_[1]/2){
                            y0 = 0.0;
                        }
                        else{
                            y0 = (pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*iy_tmp].y
                                  + pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*iy_tmp+1].y) * 0.5;
                        }
                        if(iy==n_domain_[1]-1){
                            y1 = pos_root_domain_.high_.y;
                        }
                        else{
                            y1 = (pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*(iy_tmp+1)].y
                                  + pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*(iy_tmp+1)+1].y) * 0.5;                            
                        }
                        for(S32 i = iy0; i < iy1; i++) {
                            pos_domain_temp_[i].low_.y  = y0;
                            pos_domain_temp_[i].high_.y = y1;
                        }
                        //std::cerr<<"iy= "<<iy<<" y0= "<<y0<<" y1= "<<y1<<std::endl;
                    }
                    
                    ///////////////
                    // FOR DEBUG
                    for(S32 iy = 0; iy<n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        for(S32 i = iy0; i < iy1; i++) {
                            if(pos_domain_temp_[i].high_.y <= pos_domain_temp_[i].low_.y){
                                std::cerr<<"i= "<<i
                                         <<" pos_domain_temp_[i].high_.y= "<<pos_domain_temp_[i].high_.y
                                         <<" pos_domain_temp_[i].low_.y= "<<pos_domain_temp_[i].low_.y
                                         <<std::endl;
                            }
                            assert(pos_domain_temp_[i].high_.y > pos_domain_temp_[i].low_.y);
                            if(iy > 0){
                                if(fabs(pos_domain_temp_[i-1].high_.y - pos_domain_temp_[i].low_.y) >= 1e-13){
                                    std::cerr<<"i= "<<i
                                             <<"pos_domain_temp_[i-1].high_.y= "<<pos_domain_temp_[i-1].high_.y
                                             <<"pos_domain_temp_[i].low_.y= "<<pos_domain_temp_[i].low_.y
                                             <<std::endl;
                                }
                                assert( fabs(pos_domain_temp_[i-1].high_.y - pos_domain_temp_[i].low_.y) < 1e-13);
                            }
                            if(iy < n_domain_[1]-1){
                                if(fabs(pos_domain_temp_[i+1].low_.y - pos_domain_temp_[i].high_.y) >= 1e-13 ){
                                    std::cerr<<"i= "<<i
                                             <<"pos_domain_temp_[i+1].low_.y= "<<pos_domain_temp_[i+1].low_.y
                                             <<"pos_domain_temp_[i].high_.y= "<<pos_domain_temp_[i].high_.y
                                             <<std::endl;
                                }
                                assert( fabs(pos_domain_temp_[i+1].low_.y - pos_domain_temp_[i].high_.y) < 1e-13);
                            }
                        }
                    }
                    // FOR DEBUG
                    ///////////////
                }
                
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- z direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    //std::cerr<<"ix= "<<ix<<" n_domain_[0]= "<<n_domain_[0]<<std::endl;
                    S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //if(Comm::getRank()==0) std::cerr<<"check 7"<<std::endl;
                        //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
			std::sort(pos_sample_tot_+istart[iy0], pos_sample_tot_+(iend[iy1-1]+1), Cmpvec(&F64vec::z));
                        S32 number_of_sample_particle_tot_z = iend[iy1-1] - istart[iy0] + 1;
                        for(S32 iz = 0; iz < n_domain_[2]; iz++) {
                            S32 iz0 = iy0 + iz;
                            //F64 z0, z1;
			    F64 z0 = 0.0;
			    F64 z1 = 0.0;						    
                            //calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
			    calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            pos_domain_temp_[iz0].low_.z  = z0;
                            pos_domain_temp_[iz0].high_.z = z1;
                        }
                    }
                }
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- process first ------------------------
                if(first_call_by_decomposeDomain) {
                    first_call_by_decomposeDomain = false;
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                    }
                } else {
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                        
                    }
                }
                // ------------------------------------------
                delete [] istart;
                delete [] iend;
            }
            // ****************************************************
            // *** broad cast pos_domain_ *************************
            MPI::COMM_WORLD.Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0);
            //std::cout<<"end of bcast: "<<"time: "<<GetWtime() - Tbegin<<std::endl;
            //Comm::broadcast(pos_domain_, nproc);
            // ****************************************************
#else       // PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#endif     // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"pos_root_domain_="<<pos_root_domain_<<std::endl;
            std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
#endif
            time_profile_.decompose_domain += GetWtime() - time_offset;
        }


        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys,
                                const F32 wgh){
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }

        template<class Tpsys>
        void decomposeDomainAll(Tpsys & psys){
            const F32 wgh = psys.getNumberOfParticleLocal();
            const bool clear = true;
            collectSampleParticle(psys, clear, wgh);
            decomposeDomain();
        }

        void getRootDomain(FILE *fp) {
            fprintf(fp, "%+e %+e %+e\n",
                pos_root_domain_.low_[0],
                pos_root_domain_.low_[1],
                pos_root_domain_.low_[2]);
            fprintf(fp, "%+e %+e %+e\n",
                pos_root_domain_.high_[0],
                pos_root_domain_.high_[1],
                pos_root_domain_.high_[2]);

            return;
        }


#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
        void getSampleParticleLocal(FILE *fp) {
            for(S32 i = 0; i < number_of_sample_particle_loc_; i++) {
                fprintf(fp, "%+e %+e %+e\n", pos_sample_loc_[i].x, pos_sample_loc_[i].y, pos_sample_loc_[i].z);
            }            
            return;
        }

        void getSampleParticleTotal(FILE *fp) {
            for(S32 i = 0; i < number_of_sample_particle_tot_; i++) {
                fprintf(fp, "%+e %+e %+e\n", pos_sample_tot_[i].x, pos_sample_tot_[i].y, pos_sample_tot_[i].z);
            }            
            return;
        }
#endif
	
        void getPosDomainTotal(FILE *fp) {
            //S32 nproc = MPI::COMM_WORLD.Get_size();
            S32 nproc = Comm::getNumberOfProc();
            for(S32 i = 0; i < nproc; i++) {
                for(S32 k = 0; k < DIMENSION; k++)
                    fprintf(fp, "%+e ", pos_domain_[i].low_[k]);
                for(S32 k = 0; k < DIMENSION; k++)
                    fprintf(fp, "%+e ", pos_domain_[i].high_[k]);
                fprintf(fp, "\n");
            }
        }

        S32 * getPointerOfNDomain(){return n_domain_;};

        /* AT_DEBUG
        F32ort * getPointerOfPosDomain(){return pos_domain_;};
        */
        F64ort * getPointerOfPosDomain(){return pos_domain_;};

        // A. Tanikawa need this method for Particle Mesh...
        S32 getNDomain(const S32 dim) const {return n_domain_[dim];};

        /* AT_DEBUG
        // for DEBUG
        F32vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
        // for DEBUG // A. Tanikawa would not like to delete this method...
        F32ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
        // for DEBUG
        void setPosDomain(const S32 id, const F32ort & pos){ pos_domain_[id] = pos;}
        */
        // for DEBUG
        F64vec & getPosSample(const S32 id=0){return pos_sample_tot_[id];}
        // for DEBUG // A. Tanikawa would not like to delete this method...
        F64ort & getPosDomain(const S32 id=0) const {return pos_domain_[id];}
        // for DEBUG
        void setPosDomain(const S32 id, const F64ort & pos){ pos_domain_[id] = pos;}
        // added by DN
        void initializePosDomain(F64ort & pos_domain){
            S32 myrank = Comm::getRank();
            S32 nprocs = Comm::getNumberOfProc();
            pos_domain_temp_[myrank].low_  = pos_domain.low_;
            pos_domain_temp_[myrank].high_ = pos_domain.high_;
            // Broadcast
            MPI::COMM_WORLD.Allgather(&pos_domain,1,GetDataType<F64ort>(),
                                      &pos_domain_temp_[0],1,GetDataType<F64ort>());
            // copy
            for (S32 irank=0; irank<nprocs; irank++) {
                pos_domain_[irank].low_ = pos_domain_temp_[irank].low_;
                pos_domain_[irank].high_ = pos_domain_temp_[irank].high_;
            }
            // output to check
            //std::cout << "pos_domain.low_.x = " << pos_domain.low_.x << std::endl;
            //std::cout << "pos_domain.low_.y = " << pos_domain.low_.y << std::endl;
            //std::cout << "pos_domain.low_.z = " << pos_domain.low_.z << std::endl;
            //std::cout << "pos_domain.high_.x = " << pos_domain.high_.x << std::endl;
            //std::cout << "pos_domain.high_.y = " << pos_domain.high_.y << std::endl;
            //std::cout << "pos_domain.high_.z = " << pos_domain.high_.z << std::endl;
            //std::cout << "pos_domain_[myrank].low_.x = " << pos_domain_[myrank].low_.x << std::endl;
            //std::cout << "pos_domain_[myrank].low_.y = " << pos_domain_[myrank].low_.y << std::endl;
            //std::cout << "pos_domain_[myrank].low_.z = " << pos_domain_[myrank].low_.z << std::endl;
            //std::cout << "pos_domain_[myrank].high_.x = " << pos_domain_[myrank].high_.x << std::endl;
            //std::cout << "pos_domain_[myrank].high_.y = " << pos_domain_[myrank].high_.y << std::endl;
            //std::cout << "pos_domain_[myrank].high_.z = " << pos_domain_[myrank].high_.z << std::endl;
        }

        void setBoundaryCondition(enum BOUNDARY_CONDITION bc){
            boundary_condition_ = bc;
            if(DIMENSION == 2 && 
               (bc == BOUNDARY_CONDITION_PERIODIC_XYZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_XZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_YZ ||
                bc == BOUNDARY_CONDITION_PERIODIC_Z ) ){
                throw "PS_ERROR: in setBoundaryCondition(enum BOUNDARY_CONDITION) \n boundary condition is incompatible with DIMENSION";
            }
            if(bc == BOUNDARY_CONDITION_PERIODIC_X) periodic_axis_[0] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_Y) periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_Z) periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XY) periodic_axis_[0] = periodic_axis_[1] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XZ) periodic_axis_[0] = periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_YZ) periodic_axis_[1] = periodic_axis_[2] = true;
            else if(bc == BOUNDARY_CONDITION_PERIODIC_XYZ) periodic_axis_[0] = periodic_axis_[1] = periodic_axis_[2] = true;
        }

        S32 getBoundaryCondition() const { return boundary_condition_; }

        /* AT_DEBUG
        void setPosRootDomain(const F32vec & low, const F32vec & high){
        */
        void setPosRootDomain(const F64vec & low, const F64vec & high){
            for(S32 k=0; k<DIMENSION; k++){
                if(low[k] > high[k]){
                    PARTICLE_SIMULATOR_PRINT_ERROR("The coodinate of the root domain is inconsistent.");
                    std::cerr<<"The coordinate of the low vertex of the rood domain="<<low<<std::endl;
                    std::cerr<<"The coordinate of the high vertex of the rood domain="<<high<<std::endl;
                    Abort(-1);
                }
            }

            for(S32 i=0; i<DIMENSION; i++){
                //std::cerr<<"low[i]="<<low[i]<<std::endl;
                //std::cerr<<"high[i]="<<high[i]<<std::endl;
                //std::cerr<<"periodic_axis_[i]="<<periodic_axis_[i]<<std::endl;
#ifndef PHI_R_TREE
                if( periodic_axis_[i] == false ) continue;
#endif
                if(low[i] < high[i]){
                    pos_root_domain_.low_[i] = low[i];
                    pos_root_domain_.high_[i] = high[i];
                }
            }
        }

        /* AT_DEBUG
        F32ort getPosRootDomain() const { return pos_root_domain_; }
        */
        const F64ort getPosRootDomain() const { return pos_root_domain_; }

        void getPeriodicAxis(bool pa[]) const {
            for(S32 i=0; i<DIMENSION; i++) pa[i] = periodic_axis_[i];
        }

        template<class Tpsys>
        bool checkCollectSampleParticleSubset(Tpsys & psys);
        template<class Tpsys>
        bool checkCollectSampleParticleAverage(Tpsys & psys);
        template<class Tpsys>
        bool checkDecomposeDomain(Tpsys & psys);

        S32 getRank1d(const S32 dim) const {return rank_1d_[dim];}


        // FOR RING_SYSTEM
	// new version multi-dimensional gathering
        void decomposeDomainMultiStep2(const bool flag_smaple_sort=false,
                                       const bool split_y_at_the_center=false,
                                       const F64  y_coord_split=0.0) {
            F64 time_offset = GetWtime();
            if(flag_smaple_sort){
                const F32 freq_smp_sort = 0.1;
                decomposeDomain3(freq_smp_sort, split_y_at_the_center, y_coord_split);
                first_call_by_decomposeDomain = true;
            }
            Comm::barrier(); if (Comm::getRank() == 0) std::cout << "etime2(flag_smaple_sort) = " << GetWtime() - time_offset << std::endl;
	    //assert(!first_call_by_decomposeDomain);
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL 
            pos_domain_[0] = pos_root_domain_;
#else //PARTICLE_SIMULATOR_MPI_PARALLEL 
            static bool first = true;
            static S32 * n_send;
            static S32 * n_recv;
            static S32 * n_send_disp;
            static S32 * n_recv_disp;
            static S32 * i_head;
            static S32 * i_tail;
            static F64vec * pos_sample_buf;
            static F64 *  coord_buf;
            static F64 *  coord_tot;
            static F64 * x_coord;
            static F64 * y_coord;
            static F64ort * pos_domain_temp_buf;
            static MPI_Request * req_send;
            static MPI_Request * req_recv;
            static MPI_Status * stat_send;
            static MPI_Status * stat_recv;
            const S32 n_proc_glb = Comm::getNumberOfProc();
            const S32 my_rank_glb = Comm::getRank();
            if(first){
                n_send = new S32[n_proc_glb];
                n_recv = new S32[n_proc_glb];
                n_send_disp = new S32[n_proc_glb + 1];
                n_recv_disp = new S32[n_proc_glb + 1];
                i_head = new S32[n_proc_glb];
                i_tail = new S32[n_proc_glb];
                pos_sample_buf = new F64vec[target_number_of_sample_particle_];
                coord_buf = new F64[n_proc_glb * 2];
                coord_tot = new F64[n_proc_glb * 2];
                x_coord = new F64[n_proc_glb + 1];
                y_coord = new F64[n_proc_glb + 1];
                pos_domain_temp_buf = new F64ort[n_proc_glb];
                req_send = new MPI_Request[n_proc_glb];
                req_recv = new MPI_Request[n_proc_glb];
                stat_send = new MPI_Status[n_proc_glb];
                stat_recv = new MPI_Status[n_proc_glb];
                first = false;
            }
            //std::cout<<"rank_glb="<<rank_glb<<" number_of_sample_particle_loc_="<<number_of_sample_particle_loc_<<std::endl;
            ///////////// sort particles along x direction
            std::sort(pos_sample_loc_, pos_sample_loc_+number_of_sample_particle_loc_, LessOPX());

            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
            S32 id_domain_3d = 0;
            S32 id_domain_x = 0;
    #if 1
        #ifdef DEBUG_DD_MULTI_STEP2
            S32 n_cnt_err = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                if(!pos_root_domain_.contained(pos_sample_loc_[i])){
                    if(my_rank_glb == n_proc_glb/2){
                        std::cerr<<"pos_sample_loc_[i]= "<<pos_sample_loc_[i]<<std::endl;
                    }
                    n_cnt_err++;
                }
            }
            std::cerr<<"n_cnt_err= "<<n_cnt_err<<std::endl;
        #endif //DEBUG_DD_MULTI_STEP2
            //S32 disp_domain_x_max_loc = (abs(id_domain_x-rank_1d_[0]) < (n_domain_[0]-abs(id_domain_x-rank_1d_[0]))) ? abs(id_domain_x-rank_1d_[0]) : n_domain_[0]-abs(id_domain_x-rank_1d_[0]);
            S32 disp_domain_x_max_loc = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                while( pos_domain_[id_domain_3d].high_.x <= pos_sample_loc_[i].x ){
                    id_domain_3d += n_proc_sub_[0];
                    id_domain_x++;
                }
                if(n_send[id_domain_x] == 0){
                    S32 disp_domain_tmp = abs(id_domain_x-rank_1d_[0]);
                    disp_domain_tmp = (disp_domain_tmp < (n_domain_[0]-disp_domain_tmp) ) ? disp_domain_tmp : (n_domain_[0]-disp_domain_tmp);
                    if(disp_domain_x_max_loc < disp_domain_tmp) disp_domain_x_max_loc = disp_domain_tmp;
                }
                n_send[id_domain_x]++;
            }
            S32 disp_domain_x_max_glb = 0;
            MPI_Allreduce(&disp_domain_x_max_loc, &disp_domain_x_max_glb, 1, MPI_INT, MPI_MAX, comm_1d_[0]);
        #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 1 @decomposeDomainMultiStep2"<<std::endl;
        #endif

            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            S32 my_rank_1d = 0;
            MPI_Comm_rank(comm_1d_[0], &my_rank_1d);
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                MPI_Isend(n_send+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_send+n_proc_send);
                MPI_Irecv(n_recv+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_recv+n_proc_recv);
                n_proc_send++;
                n_proc_recv++;
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            
        #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 2 @decomposeDomainMultiStep2"<<std::endl;
        #endif
            
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            
            n_proc_send = n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                if(n_send[id_tmp] > 0){
                    MPI_Isend(pos_sample_loc_+n_send_disp[id_tmp], n_send[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_send+n_proc_send);
                    n_proc_send++;
                }
                if(n_recv[id_tmp] > 0){
                    MPI_Irecv(pos_sample_buf+n_recv_disp[id_tmp], n_recv[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_recv+n_proc_recv);
                    n_proc_recv++;
                }
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);

        #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 3 @decomposeDomainMultiStep2"<<std::endl;
        #endif
    #else //RING_SYSTEM
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            MPI_Alltoallv(pos_sample_loc_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
    #endif //RING_SYSTEM

            //////////////////
            ///////////// allgather particles in Y-Z plane
            // this part is the same as original
            S32 n_send_tmp = n_recv_disp[ n_domain_[0] ]; // # of particle in own cell.
            MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_sub_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            S32 n_par_slab = n_recv_disp[ n_proc_sub_[0] ];
            MPI_Allgatherv(pos_sample_buf,  n_send_tmp, GetDataType<F64vec>(),
                           pos_sample_tot_, n_recv, n_recv_disp, GetDataType<F64vec>(), comm_sub_[0]);

            ///////////// sort particles along x direction again.
            ///////////// after sort, particles are sorted in global.
            // this part is the same as original
            std::sort(pos_sample_tot_, pos_sample_tot_+n_par_slab, LessOPX());

            
            ///////////////////////////////
	    ///////////// determine X coord
            // to determine global rank in x
            MPI_Allgather(&n_par_slab, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            number_of_sample_particle_tot_ = n_recv_disp[n_domain_[0]];

	    // get index of 
            S32 n_ave = number_of_sample_particle_tot_ / n_domain_[0];
            for(S32 i=0; i<n_domain_[0]; i++){
                i_head[i] = n_ave * i;
                if( i < number_of_sample_particle_tot_ % n_domain_[0]){
                    i_head[i] += i;
                }
                else{
                    i_head[i] += number_of_sample_particle_tot_ % n_domain_[0];
                }
                if(i > 0) i_tail[i-1] = i_head[i] - 1;
            }
            i_tail[n_domain_[0]-1] = number_of_sample_particle_tot_ - 1;
            n_send_tmp = 0; // temporally used
            for(S32 i=0; i<n_domain_[0]; i++){
                if( n_recv_disp[rank_1d_[0]] <= i_head[i] &&  i_head[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_head[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
                if( n_recv_disp[rank_1d_[0]] <= i_tail[i] &&  i_tail[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_tail[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
            }

	    MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(),
                          n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
	    n_recv_disp[0] = 0;
	    for(S32 i=0; i<n_domain_[0]; i++){
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
	    }

	    MPI_Allgatherv(coord_buf, n_send_tmp, GetDataType<>(coord_buf[0]),
                       coord_tot, n_recv, n_recv_disp, GetDataType<>(coord_buf[0]), comm_1d_[0]);
	    
	    assert( n_recv_disp[n_domain_[0]] == n_domain_[0]*2);

	    // size of x_coord_buf is n_domain_[0]+1
	    x_coord[0] = pos_root_domain_.low_.x;
	    x_coord[n_domain_[0]] = pos_root_domain_.high_.x;

	    for(S32 i=1; i<n_domain_[0]; i++){
                x_coord[i] = (coord_tot[i*2] + coord_tot[i*2-1]) * 0.5;
	    }
            ////////////////////////////////////////////
            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
	    id_domain_x = 0;
        #ifdef RING_SYSTEM
            disp_domain_x_max_loc = 0;
            for(S32 i=0; i<n_par_slab; i++){
                while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
                if(n_send[id_domain_x] == 0){
                    S32 disp_domain_tmp = abs(id_domain_x-rank_1d_[0]);
                    disp_domain_tmp = (disp_domain_tmp < (n_domain_[0]-disp_domain_tmp) ) ? disp_domain_tmp : (n_domain_[0]-disp_domain_tmp);
                    if(disp_domain_x_max_loc < disp_domain_tmp) disp_domain_x_max_loc = disp_domain_tmp;
                }
                n_send[id_domain_x]++;
            }
            disp_domain_x_max_glb = 0;
            MPI_Allreduce(&disp_domain_x_max_loc, &disp_domain_x_max_glb, 1, MPI_INT, MPI_MAX, comm_1d_[0]);

            #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 4 @decomposeDomainMultiStep2"<<std::endl;
            #endif
            n_proc_send = 0;
            n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                MPI_Isend(n_send+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_send+n_proc_send);
                MPI_Irecv(n_recv+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_recv+n_proc_recv);
                n_proc_send++;
                n_proc_recv++;                
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 5 @decomposeDomainMultiStep2"<<std::endl;
            #endif
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            n_proc_send = n_proc_recv = 0;

            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                if(n_send[id_tmp] > 0){
                    MPI_Isend(pos_sample_tot_+n_send_disp[id_tmp], n_send[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_send+n_proc_send);
                    n_proc_send++;
                }
                if(n_recv[id_tmp] > 0){
                    MPI_Irecv(pos_sample_buf+n_recv_disp[id_tmp], n_recv[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_recv+n_proc_recv);
                    n_proc_recv++;
                }
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 6 @decomposeDomainMultiStep2"<<std::endl;
            #endif            
        #else //RING_SYSTEM
	    for(S32 i=0; i<n_par_slab; i++){
                while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
                n_send[id_domain_x]++;
            }
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
	    MPI_Alltoallv(pos_sample_tot_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
        #endif //RING_SYSTEM
	    n_par_slab = n_recv_disp[n_domain_[0]];
	    // OK

    #if 1
	    ////////////////////////////////////
            ///////////// determine y corrdinate
	    std::sort(pos_sample_buf, pos_sample_buf+n_par_slab, LessOPY());
            if(split_y_at_the_center){
                assert(n_domain_[1]%2==0);
                F64vec * y_coord_mid = std::lower_bound(pos_sample_buf, pos_sample_buf+n_par_slab, y_coord_split, CompYDir());
                S32 number_of_sample_particle_sep_y_low = y_coord_mid - pos_sample_buf; // i.e. pos_sample_buf[number_of_sample_particle_sep_y_low-1] <= y_coord_split

                // lower part
                S32 n_ave_low = number_of_sample_particle_sep_y_low / (n_domain_[1]/2);
                for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                    i_head[iy] = n_ave_low * iy;
                    if( iy < number_of_sample_particle_sep_y_low % (n_domain_[1]/2) ){ i_head[iy] += iy;}
                    else{ i_head[iy] += number_of_sample_particle_sep_y_low % (n_domain_[1]/2);}
                    if(iy > 0) i_tail[iy-1] = i_head[iy] - 1;
                }
                i_tail[n_domain_[1]/2-1] = number_of_sample_particle_sep_y_low - 1;

                // higher part
                S32 number_of_sample_particle_sep_y_high = n_par_slab - number_of_sample_particle_sep_y_low;
                S32 n_ave_high = number_of_sample_particle_sep_y_high / (n_domain_[1]/2);
                for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                    S32 iy_glb = iy+n_domain_[1]/2;
                    i_head[iy_glb] = n_ave_high * iy + number_of_sample_particle_sep_y_low;
                    if( iy < number_of_sample_particle_sep_y_high % (n_domain_[1]/2) ){ i_head[iy_glb] += iy;}
                    else{ i_head[iy_glb] += number_of_sample_particle_sep_y_high % (n_domain_[1]/2);}
                    if(iy > 0) i_tail[iy_glb-1] = i_head[iy_glb] - 1;
                }
                i_tail[n_domain_[1]-1] = n_par_slab - 1;
                
                y_coord[0] = pos_root_domain_.low_.y;
                y_coord[n_domain_[1]/2] = y_coord_split;
                y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
                for(S32 i=1; i<n_domain_[1]/2; i++){ y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;}
                for(S32 i=n_domain_[1]/2+1; i<n_domain_[1]; i++){ y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;}
                if(Comm::getRank() == 0){
                    for(S32 i=0; i<n_domain_[1]; i++){
                        std::cerr<<"i_head[i]= "<<i_head[i]<<" i_tail[i]= "<<i_tail[i]<<std::endl;
                    }
                    for(S32 i=0; i<n_domain_[1]+1; i++){
                        std::cerr<<"y_coord[i]= "<<y_coord[i]<<std::endl;
                    }
                }
            }
            else{
                // get index of
                n_ave = n_par_slab / n_domain_[1];
                for(S32 i=0; i<n_domain_[1]; i++){
                    i_head[i] = n_ave * i;
                    if( i < n_par_slab % n_domain_[1]) i_head[i] += i;
                    else i_head[i] += n_par_slab % n_domain_[1];
                    if(i > 0) i_tail[i-1] = i_head[i] - 1;
                }
                i_tail[n_domain_[1]-1] = n_par_slab - 1;
                // size of y_coord is n_domain_[1]+1
                y_coord[0] = pos_root_domain_.low_.y;
                y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
                for(S32 i=1; i<n_domain_[1]; i++) y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;
            }


        #ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    ////////////////////////////////////
            ///////////// determine z corrdinate
	    //#pragma omp parallel for
	    for(S32 iy=0; iy<n_domain_[1]; iy++){
            const S32 iy_ptcl_head = i_head[iy];
            const S32 iy_ptcl_tail = i_tail[iy];
            const S32 nz_tot = iy_ptcl_tail - iy_ptcl_head + 1;
            std::sort(pos_sample_buf+iy_ptcl_head, pos_sample_buf+iy_ptcl_head+nz_tot, LessOPZ());
            S32 nz_ave_tmp = nz_tot / n_domain_[2];
            pos_domain_temp_buf[iy*n_domain_[2]].low_.z = pos_root_domain_.low_.z;
            pos_domain_temp_buf[(iy+1)*n_domain_[2] - 1].high_.z = pos_root_domain_.high_.z;
            pos_domain_temp_buf[iy*n_domain_[2]].low_.y = y_coord[iy];
            pos_domain_temp_buf[iy*n_domain_[2]].high_.y = y_coord[iy+1];
            for(S32 iz=1; iz<n_domain_[2]; iz++){
                pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.y = y_coord[iy];
                pos_domain_temp_buf[iy*n_domain_[2]+iz].high_.y = y_coord[iy+1];
                S32 iz_tmp = nz_ave_tmp * iz;
                if(iz < nz_tot % n_domain_[2]) iz_tmp += iz;
                else iz_tmp += nz_tot % n_domain_[2];
                F64 z_coord_tmp = (pos_sample_buf[iy_ptcl_head+iz_tmp].z + pos_sample_buf[iy_ptcl_head+iz_tmp-1].z) * 0.5;

                pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.z = z_coord_tmp;
                pos_domain_temp_buf[iy*n_domain_[2]+iz-1].high_.z = z_coord_tmp;
            }
	    }
        #endif //PARTICLE_SIMULATOR_TWO_DIMENSION
	    for(S32 i=0; i<n_proc_sub_[0]; i++){
                pos_domain_temp_buf[i].low_.x = x_coord[rank_1d_[0]];
                pos_domain_temp_buf[i].high_.x = x_coord[rank_1d_[0]+1];
	    }

	    //////////////////////////////////////////////
            ///////////// exchange pos_domain_tmp
	    MPI_Allgather(pos_domain_temp_buf, n_proc_sub_[0], GetDataType<F64ort>(),
                      pos_domain_temp_, n_proc_sub_[0], GetDataType<F64ort>(), comm_1d_[0]);

	    if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;
                for(S32 i = 0; i < n_proc_glb; i++) {
                    //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                    pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                    pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                }
	    } else {
                for(S32 i = 0; i < n_proc_glb; i++) {
                    //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                    pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                    pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                }
	    }
            //time_profile_.decompose_domain += GetWtime() - time_offset;
     #endif // PARTICLE_SIMULATOR_MPI_PARALLEL
	    time_profile_.decompose_domain = GetWtime() - time_offset;

#endif //#if 1
            
        }

        void decomposeDomainMultiStep3(const bool flag_smaple_sort=false,
                                       const bool split_y_at_the_center=false,
                                       const F64  y_coord_split=0.0) {
            F64 time_offset = GetWtime();
            if(flag_smaple_sort){
                const F32 freq_smp_sort = 0.1;
                decomposeDomain3(freq_smp_sort, split_y_at_the_center, y_coord_split);
                first_call_by_decomposeDomain = true;
            }
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#else //PARTICLE_SIMULATOR_MPI_PARALLEL 
            static bool first = true;
            static S32 * n_send;
            static S32 * n_recv;
            static S32 * n_send_disp;
            static S32 * n_recv_disp;
            static S32 * i_head;
            static S32 * i_tail;
            static F64vec * pos_sample_buf;
            static F64 *  coord_buf;
            static F64 *  coord_tot;
            static F64 * x_coord;
            static F64 * y_coord;
            static F64ort * pos_domain_temp_buf;
            static MPI_Request * req_send;
            static MPI_Request * req_recv;
            static MPI_Status * stat_send;
            static MPI_Status * stat_recv;
            const S32 n_proc_glb = Comm::getNumberOfProc();
            const S32 my_rank_glb = Comm::getRank();
            if(first){
                n_send = new S32[n_proc_glb];
                n_recv = new S32[n_proc_glb];
                n_send_disp = new S32[n_proc_glb + 1];
                n_recv_disp = new S32[n_proc_glb + 1];
                i_head = new S32[n_proc_glb];
                i_tail = new S32[n_proc_glb];
                pos_sample_buf = new F64vec[target_number_of_sample_particle_];
                coord_buf = new F64[n_proc_glb * 2];
                coord_tot = new F64[n_proc_glb * 2];
                x_coord = new F64[n_proc_glb + 1];
                y_coord = new F64[n_proc_glb + 1];
                pos_domain_temp_buf = new F64ort[n_proc_glb];
                req_send = new MPI_Request[n_proc_glb];
                req_recv = new MPI_Request[n_proc_glb];
                stat_send = new MPI_Status[n_proc_glb];
                stat_recv = new MPI_Status[n_proc_glb];
                first = false;
            }
            ///////////// sort particles along x direction
            std::sort(pos_sample_loc_, pos_sample_loc_+number_of_sample_particle_loc_, LessOPX());
            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
            S32 id_domain_3d = 0;
            S32 id_domain_x = 0;
    #if 1
        #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            S32 n_cnt_err = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                if(!pos_root_domain_.contained(pos_sample_loc_[i])){
                    if(my_rank_glb == n_proc_glb/2){
                        std::cerr<<"pos_sample_loc_[i]= "<<pos_sample_loc_[i]<<std::endl;
                    }
                    n_cnt_err++;
                }
            }
            std::cerr<<"n_cnt_err= "<<n_cnt_err<<std::endl;
        #endif //DEBUG_PRINT_DD_MULTI_STEP3
            //S32 disp_domain_x_max_loc = (abs(id_domain_x-rank_1d_[0]) < (n_domain_[0]-abs(id_domain_x-rank_1d_[0]))) ? abs(id_domain_x-rank_1d_[0]) : n_domain_[0]-abs(id_domain_x-rank_1d_[0]);
            S32 disp_domain_x_max_loc = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                while( pos_domain_[id_domain_3d].high_.x <= pos_sample_loc_[i].x ){
                    id_domain_3d += n_proc_sub_[0];
                    id_domain_x++;
                }
                if(n_send[id_domain_x] == 0){
                    S32 disp_domain_tmp = abs(id_domain_x-rank_1d_[0]);
                    disp_domain_tmp = (disp_domain_tmp < (n_domain_[0]-disp_domain_tmp) ) ? disp_domain_tmp : (n_domain_[0]-disp_domain_tmp);
                    if(disp_domain_x_max_loc < disp_domain_tmp) disp_domain_x_max_loc = disp_domain_tmp;
                }
                n_send[id_domain_x]++;
            }
            S32 disp_domain_x_max_glb = 0;
            MPI_Allreduce(&disp_domain_x_max_loc, &disp_domain_x_max_glb, 1, MPI_INT, MPI_MAX, comm_1d_[0]);
        #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            if(Comm::getRank() == 0) std::cerr<<"OK 01 @decomposeDomainMultiStep3"<<std::endl;
        #endif

            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            S32 my_rank_1d = 0;
            MPI_Comm_rank(comm_1d_[0], &my_rank_1d);
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                MPI_Isend(n_send+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_send+n_proc_send);
                MPI_Irecv(n_recv+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_recv+n_proc_recv);
                n_proc_send++;
                n_proc_recv++;
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            
        #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            if(Comm::getRank() == 0) std::cerr<<"OK 02 @decomposeDomainMultiStep3"<<std::endl;
        #endif
            
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            
            n_proc_send = n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                if(n_send[id_tmp] > 0){
                    MPI_Isend(pos_sample_loc_+n_send_disp[id_tmp], n_send[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_send+n_proc_send);
                    n_proc_send++;
                }
                if(n_recv[id_tmp] > 0){
                    MPI_Irecv(pos_sample_buf+n_recv_disp[id_tmp], n_recv[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_recv+n_proc_recv);
                    n_proc_recv++;
                }
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);

        #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            if(Comm::getRank() == 0) std::cerr<<"OK 03 @decomposeDomainMultiStep3"<<std::endl;
        #endif
    #else // #if 1
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            MPI_Alltoallv(pos_sample_loc_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
    #endif // #if 1

            //////////////////
            ///////////// allgather particles in Y-Z plane
            // this part is the same as original
            S32 n_send_tmp = n_recv_disp[ n_domain_[0] ]; // # of particle in own cell.
            MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_sub_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            S32 n_par_slab = n_recv_disp[ n_proc_sub_[0] ];
            MPI_Allgatherv(pos_sample_buf,  n_send_tmp, GetDataType<F64vec>(),
                           pos_sample_tot_, n_recv, n_recv_disp, GetDataType<F64vec>(), comm_sub_[0]);

            ///////////// sort particles along x direction again.
            ///////////// after sort, particles are sorted in global.
            // this part is the same as original
            std::sort(pos_sample_tot_, pos_sample_tot_+n_par_slab, LessOPX());

            
            ///////////////////////////////
            ///////////// determine X coord
            // to determine global rank in x
            MPI_Allgather(&n_par_slab, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            number_of_sample_particle_tot_ = n_recv_disp[n_domain_[0]];

            // get index of 
            S32 n_ave = number_of_sample_particle_tot_ / n_domain_[0];
            for(S32 i=0; i<n_domain_[0]; i++){
                i_head[i] = n_ave * i;
                if( i < number_of_sample_particle_tot_ % n_domain_[0]){
                    i_head[i] += i;
                }
                else{
                    i_head[i] += number_of_sample_particle_tot_ % n_domain_[0];
                }
                if(i > 0) i_tail[i-1] = i_head[i] - 1;
            }
            i_tail[n_domain_[0]-1] = number_of_sample_particle_tot_ - 1;
            n_send_tmp = 0; // temporally used
            for(S32 i=0; i<n_domain_[0]; i++){
                if( n_recv_disp[rank_1d_[0]] <= i_head[i] &&  i_head[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_head[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
                if( n_recv_disp[rank_1d_[0]] <= i_tail[i] &&  i_tail[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_tail[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
            }

            MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(),
                          n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
            }

            MPI_Allgatherv(coord_buf, n_send_tmp, GetDataType<>(coord_buf[0]),
                           coord_tot, n_recv, n_recv_disp, GetDataType<>(coord_buf[0]), comm_1d_[0]);
	    
            assert( n_recv_disp[n_domain_[0]] == n_domain_[0]*2);

            // size of x_coord_buf is n_domain_[0]+1
            x_coord[0] = pos_root_domain_.low_.x;
            x_coord[n_domain_[0]] = pos_root_domain_.high_.x;
            for(S32 i=1; i<n_domain_[0]; i++){
                x_coord[i] = (coord_tot[i*2] + coord_tot[i*2-1]) * 0.5;
            }

            #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            Comm::barrier();
            if(Comm::getRank() == 0){
                for(S32 i=0; i<n_domain_[0]+1; i++){
                    std::cerr<<"x_coord[i]= "<<x_coord[i]<<std::endl;
                }
            }
            #endif //DEBUG_PRINT_DD_MULTI_STEP3

            ////////////////////////////////////////////
            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
            id_domain_x = 0;
        #if 1
            disp_domain_x_max_loc = 0;
            for(S32 i=0; i<n_par_slab; i++){
                while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
                if(n_send[id_domain_x] == 0){
                    S32 disp_domain_tmp = abs(id_domain_x-rank_1d_[0]);
                    disp_domain_tmp = (disp_domain_tmp < (n_domain_[0]-disp_domain_tmp) ) ? disp_domain_tmp : (n_domain_[0]-disp_domain_tmp);
                    if(disp_domain_x_max_loc < disp_domain_tmp) disp_domain_x_max_loc = disp_domain_tmp;
                }
                n_send[id_domain_x]++;
            }
            disp_domain_x_max_glb = 0;
            MPI_Allreduce(&disp_domain_x_max_loc, &disp_domain_x_max_glb, 1, MPI_INT, MPI_MAX, comm_1d_[0]);

            #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"OK 04 @decomposeDomainMultiStep3"<<std::endl;
                std::cerr<<"disp_domain_x_max_glb= "<<disp_domain_x_max_glb<<std::endl;
            }
            {
                std::vector<S32> n_send_tmp(n_domain_[0]);
                for(S32 i=0; i<n_domain_[0]; i++) n_send_tmp[i] = 0;
                S32 rank_domain_x = 0;
                for(S32 i=0; i<n_par_slab; i++){
                    for(S32 j=0; j<n_domain_[0]; j++){
                        if( x_coord[j] <= pos_sample_tot_[i].x && pos_sample_tot_[i].x < x_coord[j+1] ){
                            n_send_tmp[j]++;
                        }
                    }
                }
                if(Comm::getRank() == 0){
                    for(S32 j=0; j<n_domain_[0]; j++){
                        std::cerr<<"n_send[j]= "<<n_send[j]<<" n_send_tmp[j]= "<<n_send_tmp[j]<<std::endl;
                    }
                }
                for(S32 j=0; j<n_domain_[0]; j++)  assert(n_send[j] == n_send_tmp[j]);
            }
            #endif

            if(n_domain_[0] % 2 == 1){ assert(disp_domain_x_max_glb*2+1 < n_domain_[0]); }
            else{ assert(disp_domain_x_max_glb*2 <= n_domain_[0]); } // TODO: must be fixed
            n_proc_send = 0;
            n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                MPI_Isend(n_send+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_send+n_proc_send);
                MPI_Irecv(n_recv+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_recv+n_proc_recv);
                n_proc_send++;
                n_proc_recv++;                
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);

            #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"OK 05 @decomposeDomainMultiStep3"<<std::endl;
                for(S32 i=0; i<n_domain_[0]; i++){
                    std::cerr<<"i= "<<i<<" n_send[i]= "<<n_send[i]<<" n_recv[i]= "<<n_recv[i]<<std::endl;
                }
            }
            //exit(1);
            /// upto here
            #endif
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            n_proc_send = n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                //if(id_tmp == rank_1d_[0] ) continue;
                if(n_send[id_tmp] > 0){
                    MPI_Isend(pos_sample_tot_+n_send_disp[id_tmp], n_send[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_send+n_proc_send);
                    n_proc_send++;
                }
                if(n_recv[id_tmp] > 0){
                    MPI_Irecv(pos_sample_buf+n_recv_disp[id_tmp], n_recv[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_recv+n_proc_recv);
                    n_proc_recv++;
                }
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            #ifdef DEBUG_PRINT_DD_MULTI_STEP3
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"OK 06 @decomposeDomainMultiStep3"<<std::endl;
            }
            for(S32 i=0; i<n_recv_disp[n_domain_[0]]; i++){
                if( x_coord[my_rank_1d] > pos_sample_buf[i].x || pos_sample_buf[i].x >= x_coord[my_rank_1d+1] ){
                    std::cerr<<"my_rank_1d= "<<my_rank_1d<<std::endl;
                    std::cerr<<"x_coord[my_rank_1d]= "<<x_coord[my_rank_1d]<<" x_coord[my_rank_1d+1]= "<<x_coord[my_rank_1d+1]<<std::endl;
                    std::cerr<<"pos_sample_buf[i]= "<<pos_sample_buf[i]<<std::endl;
                }
                assert( x_coord[my_rank_1d] <= pos_sample_buf[i].x && pos_sample_buf[i].x < x_coord[my_rank_1d+1] );
            }
            //exit(1);
            #endif
        #else //#if 1
            for(S32 i=0; i<n_par_slab; i++){
                while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
                n_send[id_domain_x]++;
            }
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            MPI_Alltoallv(pos_sample_tot_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
        #endif //#if 1
            n_par_slab = n_recv_disp[n_domain_[0]];
	    // OK

    #if 1
            ////////////////////////////////////
            ///////////// determine y corrdinate
            std::sort(pos_sample_buf, pos_sample_buf+n_par_slab, LessOPY());
            if(split_y_at_the_center){
                assert(n_domain_[1]%2==0);
                F64vec * y_coord_mid = std::lower_bound(pos_sample_buf, pos_sample_buf+n_par_slab, y_coord_split, CompYDir());
                S32 number_of_sample_particle_sep_y_low = y_coord_mid - pos_sample_buf; // i.e. pos_sample_buf[number_of_sample_particle_sep_y_low-1] <= y_coord_split

                // lower part
                S32 n_ave_low = number_of_sample_particle_sep_y_low / (n_domain_[1]/2);
                for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                    i_head[iy] = n_ave_low * iy;
                    if( iy < number_of_sample_particle_sep_y_low % (n_domain_[1]/2) ){ i_head[iy] += iy;}
                    else{ i_head[iy] += number_of_sample_particle_sep_y_low % (n_domain_[1]/2);}
                    if(iy > 0) i_tail[iy-1] = i_head[iy] - 1;
                }
                i_tail[n_domain_[1]/2-1] = number_of_sample_particle_sep_y_low - 1;

                // higher part
                S32 number_of_sample_particle_sep_y_high = n_par_slab - number_of_sample_particle_sep_y_low;
                S32 n_ave_high = number_of_sample_particle_sep_y_high / (n_domain_[1]/2);
                for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                    S32 iy_glb = iy+n_domain_[1]/2;
                    i_head[iy_glb] = n_ave_high * iy + number_of_sample_particle_sep_y_low;
                    if( iy < number_of_sample_particle_sep_y_high % (n_domain_[1]/2) ){ i_head[iy_glb] += iy;}
                    else{ i_head[iy_glb] += number_of_sample_particle_sep_y_high % (n_domain_[1]/2);}
                    if(iy > 0) i_tail[iy_glb-1] = i_head[iy_glb] - 1;
                }
                i_tail[n_domain_[1]-1] = n_par_slab - 1;
                
                y_coord[0] = pos_root_domain_.low_.y;
                y_coord[n_domain_[1]/2] = y_coord_split;
                y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
                for(S32 i=1; i<n_domain_[1]/2; i++){ y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;}
                for(S32 i=n_domain_[1]/2+1; i<n_domain_[1]; i++){ y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;}
                if(Comm::getRank() == 0){
                    for(S32 i=0; i<n_domain_[1]; i++){
                        std::cerr<<"i_head[i]= "<<i_head[i]<<" i_tail[i]= "<<i_tail[i]<<std::endl;
                    }
                    for(S32 i=0; i<n_domain_[1]+1; i++){
                        std::cerr<<"y_coord[i]= "<<y_coord[i]<<std::endl;
                    }
                }
            }
            else{
                // get index of
                n_ave = n_par_slab / n_domain_[1];
                for(S32 i=0; i<n_domain_[1]; i++){
                    i_head[i] = n_ave * i;
                    if( i < n_par_slab % n_domain_[1]) i_head[i] += i;
                    else i_head[i] += n_par_slab % n_domain_[1];
                    if(i > 0) i_tail[i-1] = i_head[i] - 1;
                }
                i_tail[n_domain_[1]-1] = n_par_slab - 1;
                // size of y_coord is n_domain_[1]+1
                y_coord[0] = pos_root_domain_.low_.y;
                y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
                for(S32 i=1; i<n_domain_[1]; i++) y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;
            }

        #ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
            ////////////////////////////////////
            ///////////// determine z corrdinate
            //#pragma omp parallel for
            for(S32 iy=0; iy<n_domain_[1]; iy++){
                const S32 iy_ptcl_head = i_head[iy];
                const S32 iy_ptcl_tail = i_tail[iy];
                const S32 nz_tot = iy_ptcl_tail - iy_ptcl_head + 1;
                std::sort(pos_sample_buf+iy_ptcl_head, pos_sample_buf+iy_ptcl_head+nz_tot, LessOPZ());
                S32 nz_ave_tmp = nz_tot / n_domain_[2];
                pos_domain_temp_buf[iy*n_domain_[2]].low_.z = pos_root_domain_.low_.z;
                pos_domain_temp_buf[(iy+1)*n_domain_[2] - 1].high_.z = pos_root_domain_.high_.z;
                pos_domain_temp_buf[iy*n_domain_[2]].low_.y = y_coord[iy];
                pos_domain_temp_buf[iy*n_domain_[2]].high_.y = y_coord[iy+1];
                for(S32 iz=1; iz<n_domain_[2]; iz++){
                    pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.y = y_coord[iy];
                    pos_domain_temp_buf[iy*n_domain_[2]+iz].high_.y = y_coord[iy+1];
                    S32 iz_tmp = nz_ave_tmp * iz;
                    if(iz < nz_tot % n_domain_[2]) iz_tmp += iz;
                    else iz_tmp += nz_tot % n_domain_[2];
                    F64 z_coord_tmp = (pos_sample_buf[iy_ptcl_head+iz_tmp].z + pos_sample_buf[iy_ptcl_head+iz_tmp-1].z) * 0.5;
                    pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.z = z_coord_tmp;
                    pos_domain_temp_buf[iy*n_domain_[2]+iz-1].high_.z = z_coord_tmp;
                }
            }
        #endif //PARTICLE_SIMULATOR_TWO_DIMENSION
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                pos_domain_temp_buf[i].low_.x = x_coord[rank_1d_[0]];
                pos_domain_temp_buf[i].high_.x = x_coord[rank_1d_[0]+1];
            }

            //////////////////////////////////////////////
            ///////////// exchange pos_domain_tmp
            MPI_Allgather(pos_domain_temp_buf, n_proc_sub_[0], GetDataType<F64ort>(),
                          pos_domain_temp_, n_proc_sub_[0], GetDataType<F64ort>(), comm_1d_[0]);

#ifdef DEBUG_PRINT_DD_MULTI_STEP3
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"OK 07 @decomposeDomainMultiStep3"<<std::endl;
                for(S32 i=0; i<n_proc_glb; i++){
                    std::cerr<<"pos_domain_temp_[i]= "<<pos_domain_temp_[i]<<std::endl;
                }
            }
#endif
            if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;
                for(S32 i = 0; i < n_proc_glb; i++) {
                    //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                    pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                    pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                }
            } 
            else {
                for(S32 i = 0; i < n_proc_glb; i++) {
                    //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                    pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                    pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                }
            }
            //time_profile_.decompose_domain += GetWtime() - time_offset;
     #endif // PARTICLE_SIMULATOR_MPI_PARALLEL
            time_profile_.decompose_domain += GetWtime() - time_offset;
#endif //#if 1

#ifdef DEBUG_PRINT_DD_MULTI_STEP3
            Comm::barrier();
            if(Comm::getRank() == 0){
                std::cerr<<"OK 08 @decomposeDomainMultiStep3"<<std::endl;
            }
            //exit(1);
#endif
        }






        
        void decomposeDomain3(const F32 freq_smp_sort=1.0,
                              const bool split_y_at_the_center=false,
                              const F64 y_coord_split=0.0) {
            F64 time_offset = GetWtime();
            //F64 y_coord_split = 1.0;
            // ****** collect sample particles to process 0. ****** 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    S32 nproc  = Comm::getNumberOfProc();
            S32 myrank = Comm::getRank();
            //S32 number_of_sample_particle_loc_tmp = number_of_sample_particle_loc_;
            S32 number_of_sample_particle_loc_tmp = (S32)((F64)number_of_sample_particle_loc_*freq_smp_sort);
            if(number_of_sample_particle_loc_tmp < 5) number_of_sample_particle_loc_tmp = 5;
            if(number_of_sample_particle_loc_tmp < n_domain_[1]) number_of_sample_particle_loc_tmp = n_domain_[1];
            //Comm::allGather(&number_of_sample_particle_loc_, 1, n_smp_array_);
            Comm::allGather(&number_of_sample_particle_loc_tmp, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            //Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_tmp, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
            // ****************************************************
            // *** decompose domain *******************************
            if(myrank == 0) {
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 3"<<std::endl;
#endif
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
		std::sort(pos_sample_tot_, pos_sample_tot_+number_of_sample_particle_tot_, Cmpvec(&F64vec::x));
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 4"<<std::endl;
#endif
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    F64 x0 = 0.0;
		    F64 x1 = 0.0;
		    calculateBoundaryOfDomainX(number_of_sample_particle_tot_, pos_sample_tot_, istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        pos_domain_temp_[i].low_.x  = x0;
                        pos_domain_temp_[i].high_.x = x1;
                    }
                }
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 5"<<std::endl;
#endif
                // ------------------------------------------
                // --- y direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
		    std::sort(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), Cmpvec(&F64vec::y));
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                    if(Comm::getRank()==0) std::cerr<<"check 6"<<std::endl;
#endif
                    S32 number_of_sample_particle_tot_y = iend[ix1-1] - istart[ix0] + 1;
                    if(split_y_at_the_center){
                        F64vec * y_coord_mid = std::lower_bound(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), y_coord_split, CompYDir());
                        S32 number_of_sample_particle_sep_y_lower = y_coord_mid - (pos_sample_tot_+istart[ix0]); // i.e. pos_sample_buf[number_of_sample_particle_sep_y_low-1] <= y_coord_split
                        S32 number_of_sample_particle_sep_y_upper = number_of_sample_particle_tot_y - number_of_sample_particle_sep_y_lower;
                        // lower part
                        S32 dn = number_of_sample_particle_sep_y_lower / (n_domain_[1]/2);
                        for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                            S32 iy0 = ix0 +  iy      * n_domain_[2];
                            S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                            F64 y0 = 0.0;
                            F64 y1 = 0.0;
                            if(iy==0){
                                y0 = pos_root_domain_.low_.y;
                            }
                            else{
                                //y0 = (pos_sample_tot_[istart[ix0]+dn*(iy-1)].y + pos_sample_tot_[istart[ix0]+dn*iy].y) * 0.5;
                                y0 = (pos_sample_tot_[istart[ix0]+dn*iy].y + pos_sample_tot_[istart[ix0]+dn*iy+1].y) * 0.5;
                            }
                            if(iy==(n_domain_[1]/2-1)){
                                //y1 = 0.0;
                                y1 = y_coord_split;
                            }
                            else{
                                y1 = (pos_sample_tot_[istart[ix0]+dn*(iy+1)].y + pos_sample_tot_[istart[ix0]+dn*(iy+1)+1].y) * 0.5;
                            }
                            for(S32 i = iy0; i < iy1; i++) {
                                pos_domain_temp_[i].low_.y  = y0;
                                pos_domain_temp_[i].high_.y  = y1;
                            }
                            //std::cerr<<"iy= "<<iy<<" y0= "<<y0<<" y1= "<<y1<<std::endl;
                        }
                        //upper part
                        //number_of_sample_particle_sep_y = number_of_sample_particle_tot_y - number_of_sample_particle_sep_y;
                        dn = number_of_sample_particle_sep_y_upper / (n_domain_[1]/2);
                        for(S32 iy = n_domain_[1]/2; iy<n_domain_[1]; iy++) {
                            S32 iy_tmp = iy - n_domain_[1]/2;
                            S32 iy0 = ix0 +  iy      * n_domain_[2];
                            S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                            F64 y0 = 0.0;
                            F64 y1 = 0.0;
                            if(iy==n_domain_[1]/2){
                                //y0 = 0.0;
                                y0 = y_coord_split;
                            }
                            else{
                                y0 = (pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*iy_tmp].y
                                      + pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*iy_tmp+1].y) * 0.5;
                            }
                            if(iy==n_domain_[1]-1){
                                y1 = pos_root_domain_.high_.y;
                            }
                            else{
                                y1 = (pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*(iy_tmp+1)].y
                                      + pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*(iy_tmp+1)+1].y) * 0.5;                            
                            }
                            for(S32 i = iy0; i < iy1; i++) {
                                pos_domain_temp_[i].low_.y  = y0;
                                pos_domain_temp_[i].high_.y = y1;
                            }
                        }
                    }
                    else{
                        for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                            S32 iy0 = ix0 +  iy      * n_domain_[2];
                            S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                            //F64 y0, y1;
                            F64 y0 = 0.0;
                            F64 y1 = 0.0;
                            calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                            for(S32 i = iy0; i < iy1; i++) {
                                pos_domain_temp_[i].low_.y  = y0;
                                pos_domain_temp_[i].high_.y = y1;
                            }
                        }
                    }
                    ///////////////
                    // FOR DEBUG
                    for(S32 iy = 0; iy<n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        for(S32 i = iy0; i < iy1; i++) {
                            if(pos_domain_temp_[i].high_.y <= pos_domain_temp_[i].low_.y){
                                std::cerr<<"i= "<<i
                                         <<" pos_domain_temp_[i].high_.y= "<<pos_domain_temp_[i].high_.y
                                         <<" pos_domain_temp_[i].low_.y= "<<pos_domain_temp_[i].low_.y
                                         <<std::endl;
                            }
                            assert(pos_domain_temp_[i].high_.y > pos_domain_temp_[i].low_.y);
                            if(iy > 0){
                                if(fabs(pos_domain_temp_[i-1].high_.y - pos_domain_temp_[i].low_.y) >= 1e-13){
                                    std::cerr<<"i= "<<i
                                             <<"pos_domain_temp_[i-1].high_.y= "<<pos_domain_temp_[i-1].high_.y
                                             <<"pos_domain_temp_[i].low_.y= "<<pos_domain_temp_[i].low_.y
                                             <<std::endl;
                                }
                                assert( fabs(pos_domain_temp_[i-1].high_.y - pos_domain_temp_[i].low_.y) < 1e-13);
                            }
                            if(iy < n_domain_[1]-1){
                                if(fabs(pos_domain_temp_[i+1].low_.y - pos_domain_temp_[i].high_.y) >= 1e-13 ){
                                    std::cerr<<"i= "<<i
                                             <<"pos_domain_temp_[i+1].low_.y= "<<pos_domain_temp_[i+1].low_.y
                                             <<"pos_domain_temp_[i].high_.y= "<<pos_domain_temp_[i].high_.y
                                             <<std::endl;
                                }
                                assert( fabs(pos_domain_temp_[i+1].low_.y - pos_domain_temp_[i].high_.y) < 1e-13);
                            }
                        }
                    }
                    // FOR DEBUG
                    ///////////////
                }
                
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- z direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    //std::cerr<<"ix= "<<ix<<" n_domain_[0]= "<<n_domain_[0]<<std::endl;
                    S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //if(Comm::getRank()==0) std::cerr<<"check 7"<<std::endl;
                        //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
			std::sort(pos_sample_tot_+istart[iy0], pos_sample_tot_+(iend[iy1-1]+1), Cmpvec(&F64vec::z));
                        S32 number_of_sample_particle_tot_z = iend[iy1-1] - istart[iy0] + 1;
                        for(S32 iz = 0; iz < n_domain_[2]; iz++) {
                            S32 iz0 = iy0 + iz;
                            //F64 z0, z1;
			    F64 z0 = 0.0;
			    F64 z1 = 0.0;						    
                            //calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
			    calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            pos_domain_temp_[iz0].low_.z  = z0;
                            pos_domain_temp_[iz0].high_.z = z1;
                        }
                    }
                }
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- process first ------------------------
                if(first_call_by_decomposeDomain) {
                    //first_call_by_decomposeDomain = false;
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                    }
                } else {
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                        
                    }
                }
                // ------------------------------------------
                delete [] istart;
                delete [] iend;
            }
            // ****************************************************
            // *** broad cast pos_domain_ *************************
            MPI::COMM_WORLD.Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0);
            if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;
                MPI::COMM_WORLD.Bcast(&first_call_by_decomposeDomain, 1, GetDataType<bool>(), 0);
            }
            //std::cout<<"end of bcast: "<<"time: "<<GetWtime() - Tbegin<<std::endl;
            //Comm::broadcast(pos_domain_, nproc);
            // ****************************************************
#else       // PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#endif     // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"pos_root_domain_="<<pos_root_domain_<<std::endl;
            std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
#endif
            time_profile_.decompose_domain += GetWtime() - time_offset;
        }


        
#ifdef USE_SUPER_DOMAIN
        S32 getRankSub(const S32 dim) const {return rank_sub_[dim];}
        MPI_Comm getComm1d(const S32 dim=0) const { return comm_1d_[dim];}
        MPI_Comm getCommSub(const S32 dim=0) const { return comm_sub_[dim];}
        
        void initializeSuperDomain(){
            pos_super_domain_.resizeNoInitialize(n_domain_[0]);
            pos_sub_domain_.resizeNoInitialize(n_proc_sub_[0]);
            S32 my_rank_1d  = rank_1d_[0];
            S32 my_rank_sub = rank_sub_[0];
            S32 my_rank = Comm::getRank();
            MPI_Allgather(pos_domain_+my_rank, 1, GetDataType<F64ort>(),
                          pos_sub_domain_.getPointer(), 1, GetDataType<F64ort>(), comm_sub_[0]);
            /*
            if(my_rank_sub == 0 && my_rank_1d == 0){
                std::cerr<<"my_rank= "<<my_rank<<std::endl;
                std::cerr<<"pos_sub_domain_.size()= "<<pos_sub_domain_.size()<<std::endl;
                for(S32 i=0; i<pos_sub_domain_.size(); i++){
                    std::cerr<<"i= "<<i
                             <<" pos_sub_domain_[i]= "<<pos_sub_domain_[i]
                             <<std::endl;
                }
            }
            */
            F64ort pos_super_domain_tmp;
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                pos_super_domain_tmp.merge(pos_sub_domain_[i]);
            }
            if(my_rank_sub == 0){
                MPI_Allgather(&pos_super_domain_tmp, 1, GetDataType<F64ort>(),
                              pos_super_domain_.getPointer(), 1, GetDataType<F64ort>(),
                              comm_1d_[0]);
            }
            /*
            if(my_rank_sub == 0 && my_rank_1d == 0){
                std::cerr<<"my_rank= "<<my_rank<<std::endl;
                for(S32 i=0; i<pos_super_domain_.size(); i++){
                    std::cerr<<"i= "<<i
                             <<" pos_super_domain_[i]= "<<pos_super_domain_[i]
                             <<std::endl;
                }
            }
            */
        }

    #if 0
        void decomposeDomainUsingSuperDomain(const F32 freq_smp_sort=1.0,
                                             const bool split_y_at_the_center=false,
                                             const F64 y_coord_split=0.0) {
            F64 time_offset = GetWtime();
            // ****** collect sample particles to process 0 in super domain. ****** 
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	    S32 nproc  = Comm::getNumberOfProc();
            S32 myrank = Comm::getRank();
            S32 number_of_sample_particle_loc_tmp = number_of_sample_particle_loc_;
            if(freq_smp_sort < 1.0){
                number_of_sample_particle_loc_tmp = (S32)((F64)number_of_sample_particle_loc_*freq_smp_sort);
            }
            if(number_of_sample_particle_loc_tmp < 5) number_of_sample_particle_loc_tmp = 5;
            if(number_of_sample_particle_loc_tmp < n_domain_[1]) number_of_sample_particle_loc_tmp = n_domain_[1];
            Comm::allGather(&number_of_sample_particle_loc_tmp, 1, n_smp_array_);
            n_smp_disp_array_[0] = 0;
            for(S32 i=0; i<nproc; i++){
                n_smp_disp_array_[i+1] = n_smp_disp_array_[i] + n_smp_array_[i];
            }
            Comm::allGatherV(pos_sample_loc_, number_of_sample_particle_loc_tmp, pos_sample_tot_, n_smp_array_, n_smp_disp_array_);
            number_of_sample_particle_tot_ = n_smp_disp_array_[nproc];
            // ****************************************************
            // *** decompose domain *******************************
            if(myrank == 0) {
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 3"<<std::endl;
#endif
                S32 * istart = new S32[nproc];
                S32 * iend   = new S32[nproc];
                // --- x direction --------------------------
		std::sort(pos_sample_tot_, pos_sample_tot_+number_of_sample_particle_tot_, Cmpvec(&F64vec::x));
                for(S32 i = 0; i < nproc; i++) {
                    istart[i] = ((S64)(i) * (S64)(number_of_sample_particle_tot_)) / (S64)(nproc);
                    if(i > 0)
                        iend[i-1] = istart[i] - 1;
                }
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 4"<<std::endl;
#endif
                iend[nproc-1] = number_of_sample_particle_tot_ - 1;
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
                    F64 x0 = 0.0;
		    F64 x1 = 0.0;
		    calculateBoundaryOfDomainX(number_of_sample_particle_tot_, pos_sample_tot_, istart[ix0], iend[ix1-1], x0, x1);
                    for(S32 i = ix0; i < ix1; i++) {
                        pos_domain_temp_[i].low_.x  = x0;
                        pos_domain_temp_[i].high_.x = x1;
                    }
                }
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                if(Comm::getRank()==0) std::cerr<<"check 5"<<std::endl;
#endif
                // ------------------------------------------
                // --- y direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    S32 ix0 =  ix      * n_domain_[1] * n_domain_[2];
                    S32 ix1 = (ix + 1) * n_domain_[1] * n_domain_[2];
		    std::sort(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), Cmpvec(&F64vec::y));
#ifdef PARTICLE_SIMULATOR_DINFO_DEBUG_PRINT
                    if(Comm::getRank()==0) std::cerr<<"check 6"<<std::endl;
#endif
                    S32 number_of_sample_particle_tot_y = iend[ix1-1] - istart[ix0] + 1;
                    if(split_y_at_the_center){
                        F64vec * y_coord_mid = std::lower_bound(pos_sample_tot_+istart[ix0], pos_sample_tot_+(iend[ix1-1]+1), y_coord_split, CompYDir());
                        S32 number_of_sample_particle_sep_y_lower = y_coord_mid - (pos_sample_tot_+istart[ix0]); // i.e. pos_sample_buf[number_of_sample_particle_sep_y_low-1] <= y_coord_split
                        S32 number_of_sample_particle_sep_y_upper = number_of_sample_particle_tot_y - number_of_sample_particle_sep_y_lower;
                        // lower part
                        S32 dn = number_of_sample_particle_sep_y_lower / (n_domain_[1]/2);
                        for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                            S32 iy0 = ix0 +  iy      * n_domain_[2];
                            S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                            F64 y0 = 0.0;
                            F64 y1 = 0.0;
                            if(iy==0){
                                y0 = pos_root_domain_.low_.y;
                            }
                            else{
                                //y0 = (pos_sample_tot_[istart[ix0]+dn*(iy-1)].y + pos_sample_tot_[istart[ix0]+dn*iy].y) * 0.5;
                                y0 = (pos_sample_tot_[istart[ix0]+dn*iy].y + pos_sample_tot_[istart[ix0]+dn*iy+1].y) * 0.5;
                            }
                            if(iy==(n_domain_[1]/2-1)){
                                //y1 = 0.0;
                                y1 = y_coord_split;
                            }
                            else{
                                y1 = (pos_sample_tot_[istart[ix0]+dn*(iy+1)].y + pos_sample_tot_[istart[ix0]+dn*(iy+1)+1].y) * 0.5;
                            }
                            for(S32 i = iy0; i < iy1; i++) {
                                pos_domain_temp_[i].low_.y  = y0;
                                pos_domain_temp_[i].high_.y  = y1;
                            }
                            //std::cerr<<"iy= "<<iy<<" y0= "<<y0<<" y1= "<<y1<<std::endl;
                        }
                        //upper part
                        //number_of_sample_particle_sep_y = number_of_sample_particle_tot_y - number_of_sample_particle_sep_y;
                        dn = number_of_sample_particle_sep_y_upper / (n_domain_[1]/2);
                        for(S32 iy = n_domain_[1]/2; iy<n_domain_[1]; iy++) {
                            S32 iy_tmp = iy - n_domain_[1]/2;
                            S32 iy0 = ix0 +  iy      * n_domain_[2];
                            S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                            F64 y0 = 0.0;
                            F64 y1 = 0.0;
                            if(iy==n_domain_[1]/2){
                                //y0 = 0.0;
                                y0 = y_coord_split;
                            }
                            else{
                                y0 = (pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*iy_tmp].y
                                      + pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*iy_tmp+1].y) * 0.5;
                            }
                            if(iy==n_domain_[1]-1){
                                y1 = pos_root_domain_.high_.y;
                            }
                            else{
                                y1 = (pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*(iy_tmp+1)].y
                                      + pos_sample_tot_[istart[ix0]+number_of_sample_particle_sep_y_lower+dn*(iy_tmp+1)+1].y) * 0.5;                            
                            }
                            for(S32 i = iy0; i < iy1; i++) {
                                pos_domain_temp_[i].low_.y  = y0;
                                pos_domain_temp_[i].high_.y = y1;
                            }
                        }
                    }
                    else{
                        for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                            S32 iy0 = ix0 +  iy      * n_domain_[2];
                            S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                            //F64 y0, y1;
                            F64 y0 = 0.0;
                            F64 y1 = 0.0;
                            calculateBoundaryOfDomainY(number_of_sample_particle_tot_y, pos_sample_tot_+istart[ix0], istart[iy0]-istart[ix0], iend[iy1-1]-istart[ix0], y0, y1);
                            for(S32 i = iy0; i < iy1; i++) {
                                pos_domain_temp_[i].low_.y  = y0;
                                pos_domain_temp_[i].high_.y = y1;
                            }
                        }
                    }
                    ///////////////
                    // FOR DEBUG
                    for(S32 iy = 0; iy<n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        for(S32 i = iy0; i < iy1; i++) {
                            if(pos_domain_temp_[i].high_.y <= pos_domain_temp_[i].low_.y){
                                std::cerr<<"i= "<<i
                                         <<" pos_domain_temp_[i].high_.y= "<<pos_domain_temp_[i].high_.y
                                         <<" pos_domain_temp_[i].low_.y= "<<pos_domain_temp_[i].low_.y
                                         <<std::endl;
                            }
                            assert(pos_domain_temp_[i].high_.y > pos_domain_temp_[i].low_.y);
                            if(iy > 0){
                                if(fabs(pos_domain_temp_[i-1].high_.y - pos_domain_temp_[i].low_.y) >= 1e-13){
                                    std::cerr<<"i= "<<i
                                             <<"pos_domain_temp_[i-1].high_.y= "<<pos_domain_temp_[i-1].high_.y
                                             <<"pos_domain_temp_[i].low_.y= "<<pos_domain_temp_[i].low_.y
                                             <<std::endl;
                                }
                                assert( fabs(pos_domain_temp_[i-1].high_.y - pos_domain_temp_[i].low_.y) < 1e-13);
                            }
                            if(iy < n_domain_[1]-1){
                                if(fabs(pos_domain_temp_[i+1].low_.y - pos_domain_temp_[i].high_.y) >= 1e-13 ){
                                    std::cerr<<"i= "<<i
                                             <<"pos_domain_temp_[i+1].low_.y= "<<pos_domain_temp_[i+1].low_.y
                                             <<"pos_domain_temp_[i].high_.y= "<<pos_domain_temp_[i].high_.y
                                             <<std::endl;
                                }
                                assert( fabs(pos_domain_temp_[i+1].low_.y - pos_domain_temp_[i].high_.y) < 1e-13);
                            }
                        }
                    }
                    // FOR DEBUG
                    ///////////////
                }
                
#ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- z direction --------------------------
                for(S32 ix = 0; ix < n_domain_[0]; ix++) {
                    //std::cerr<<"ix= "<<ix<<" n_domain_[0]= "<<n_domain_[0]<<std::endl;
                    S32 ix0 = ix * n_domain_[1] * n_domain_[2];
                    for(S32 iy = 0; iy < n_domain_[1]; iy++) {
                        S32 iy0 = ix0 +  iy      * n_domain_[2];
                        S32 iy1 = ix0 + (iy + 1) * n_domain_[2];
                        //if(Comm::getRank()==0) std::cerr<<"check 7"<<std::endl;
                        //sortCoordinateOfSampleParticle(pos_sample_tot_, istart[iy0], iend[iy1-1], 2);
			std::sort(pos_sample_tot_+istart[iy0], pos_sample_tot_+(iend[iy1-1]+1), Cmpvec(&F64vec::z));
                        S32 number_of_sample_particle_tot_z = iend[iy1-1] - istart[iy0] + 1;
                        for(S32 iz = 0; iz < n_domain_[2]; iz++) {
                            S32 iz0 = iy0 + iz;
                            //F64 z0, z1;
			    F64 z0 = 0.0;
			    F64 z1 = 0.0;						    
                            //calculateBoundaryOfDomain(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], 2, istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
			    calculateBoundaryOfDomainZ(number_of_sample_particle_tot_z, pos_sample_tot_+istart[iy0], istart[iz0]-istart[iy0], iend[iz0]-istart[iy0], z0, z1);
                            pos_domain_temp_[iz0].low_.z  = z0;
                            pos_domain_temp_[iz0].high_.z = z1;
                        }
                    }
                }
#endif // PARTICLE_SIMULATOR_TWO_DIMENSION
                // ------------------------------------------
                // --- process first ------------------------
                if(first_call_by_decomposeDomain) {
                    //first_call_by_decomposeDomain = false;
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                        pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                    }
                } else {
                    for(S32 i = 0; i < nproc; i++) {
                        //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                        pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                        pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                            + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                        
                    }
                }
                // ------------------------------------------
                delete [] istart;
                delete [] iend;
            }
            // ****************************************************
            // *** broad cast pos_domain_ *************************
            MPI::COMM_WORLD.Bcast(pos_domain_, nproc, GetDataType<F64ort>(), 0);
            if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;
                MPI::COMM_WORLD.Bcast(&first_call_by_decomposeDomain, 1, GetDataType<bool>(), 0);
            }
            //std::cout<<"end of bcast: "<<"time: "<<GetWtime() - Tbegin<<std::endl;
            //Comm::broadcast(pos_domain_, nproc);
            // ****************************************************
#else       // PARTICLE_SIMULATOR_MPI_PARALLEL
            pos_domain_[0] = pos_root_domain_;
#endif     // PARTICLE_SIMULATOR_MPI_PARALLEL
#ifdef PARTICLE_SIMULATOR_DEBUG_PRINT
            PARTICLE_SIMULATOR_PRINT_LINE_INFO();
            std::cout<<"pos_root_domain_="<<pos_root_domain_<<std::endl;
            std::cout<<"pos_domain_[Comm::getRank()]="<<pos_domain_[Comm::getRank()]<<std::endl;
#endif
            time_profile_.decompose_domain += GetWtime() - time_offset;
        }        

        
        // FOR RING_SYSTEM
	// new version multi-dimensional gathering
        void decomposeDomainMultiStepUsingSyperDomain(const bool flag_smaple_sort=false,
                                                      const bool split_y_at_the_center=false,
                                                      const F64  y_coord_split=0.0) {
            F64 time_offset = GetWtime();
            if(flag_smaple_sort){
                const F32 freq_smp_sort = 0.1;
                decomposeDomain3(freq_smp_sort, split_y_at_the_center, y_coord_split);
                first_call_by_decomposeDomain = true;
            }
#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL 
            pos_domain_[0] = pos_root_domain_;
#else //PARTICLE_SIMULATOR_MPI_PARALLEL 
            static bool first = true;
            static S32 * n_send;
            static S32 * n_recv;
            static S32 * n_send_disp;
            static S32 * n_recv_disp;
            static S32 * i_head;
            static S32 * i_tail;
            static F64vec * pos_sample_buf;
            static F64 *  coord_buf;
            static F64 *  coord_tot;
            static F64 * x_coord;
            static F64 * y_coord;
            static F64ort * pos_domain_temp_buf;
            static MPI_Request * req_send;
            static MPI_Request * req_recv;
            static MPI_Status * stat_send;
            static MPI_Status * stat_recv;
            const S32 n_proc_glb = Comm::getNumberOfProc();
            const S32 my_rank_glb = Comm::getRank();
            if(first){
                n_send = new S32[n_proc_glb];
                n_recv = new S32[n_proc_glb];
                n_send_disp = new S32[n_proc_glb + 1];
                n_recv_disp = new S32[n_proc_glb + 1];
                i_head = new S32[n_proc_glb];
                i_tail = new S32[n_proc_glb];
                pos_sample_buf = new F64vec[target_number_of_sample_particle_];
                coord_buf = new F64[n_proc_glb * 2];
                coord_tot = new F64[n_proc_glb * 2];
                x_coord = new F64[n_proc_glb + 1];
                y_coord = new F64[n_proc_glb + 1];
                pos_domain_temp_buf = new F64ort[n_proc_glb];
                req_send = new MPI_Request[n_proc_glb];
                req_recv = new MPI_Request[n_proc_glb];
                stat_send = new MPI_Status[n_proc_glb];
                stat_recv = new MPI_Status[n_proc_glb];
                first = false;
            }
            //std::cout<<"rank_glb="<<rank_glb<<" number_of_sample_particle_loc_="<<number_of_sample_particle_loc_<<std::endl;
            ///////////// sort particles along x direction
            std::sort(pos_sample_loc_, pos_sample_loc_+number_of_sample_particle_loc_, LessOPX());

            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
            S32 id_domain_3d = 0;
            S32 id_domain_x = 0;
    #if 1
        #ifdef DEBUG_DD_MULTI_STEP_SUPER_DOMAIN
            S32 n_cnt_err = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                if(!pos_root_domain_.contained(pos_sample_loc_[i])){
                    if(my_rank_glb == n_proc_glb/2){
                        std::cerr<<"pos_sample_loc_[i]= "<<pos_sample_loc_[i]<<std::endl;
                    }
                    n_cnt_err++;
                }
            }
            std::cerr<<"n_cnt_err= "<<n_cnt_err<<std::endl;
        #endif //DEBUG_DD_MULTI_STEP_SUPER_DOMAIN
            S32 disp_domain_x_max_loc = 0;
            for(S32 i=0; i<number_of_sample_particle_loc_; i++){
                while( pos_domain_[id_domain_3d].high_.x <= pos_sample_loc_[i].x ){
                    id_domain_3d += n_proc_sub_[0];
                    id_domain_x++;
                }
                if(n_send[id_domain_x] == 0){
                    S32 disp_domain_tmp = abs(id_domain_x-rank_1d_[0]);
                    disp_domain_tmp = (disp_domain_tmp < (n_domain_[0]-disp_domain_tmp) ) ? disp_domain_tmp : (n_domain_[0]-disp_domain_tmp);
                    if(disp_domain_x_max_loc < disp_domain_tmp) disp_domain_x_max_loc = disp_domain_tmp;
                }
                n_send[id_domain_x]++;
            }
            S32 disp_domain_x_max_glb = 0;
            MPI_Allreduce(&disp_domain_x_max_loc, &disp_domain_x_max_glb, 1, MPI_INT, MPI_MAX, comm_1d_[0]);
        #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 1 @decomposeDomainMultiStep2"<<std::endl;
        #endif

            S32 n_proc_send = 0;
            S32 n_proc_recv = 0;
            S32 my_rank_1d = 0;
            MPI_Comm_rank(comm_1d_[0], &my_rank_1d);
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                MPI_Isend(n_send+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_send+n_proc_send);
                MPI_Irecv(n_recv+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_recv+n_proc_recv);
                n_proc_send++;
                n_proc_recv++;
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            
        #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 2 @decomposeDomainMultiStep2"<<std::endl;
        #endif
            
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            
            n_proc_send = n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                if(n_send[id_tmp] > 0){
                    MPI_Isend(pos_sample_loc_+n_send_disp[id_tmp], n_send[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_send+n_proc_send);
                    n_proc_send++;
                }
                if(n_recv[id_tmp] > 0){
                    MPI_Irecv(pos_sample_buf+n_recv_disp[id_tmp], n_recv[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_recv+n_proc_recv);
                    n_proc_recv++;
                }
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);

        #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 3 @decomposeDomainMultiStep2"<<std::endl;
        #endif
    #else //RING_SYSTEM
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            MPI_Alltoallv(pos_sample_loc_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
    #endif //RING_SYSTEM

            //////////////////
            ///////////// allgather particles in Y-Z plane
            // this part is the same as original
            S32 n_send_tmp = n_recv_disp[ n_domain_[0] ]; // # of particle in own cell.
            MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_sub_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_proc_sub_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            S32 n_par_slab = n_recv_disp[ n_proc_sub_[0] ];
            MPI_Allgatherv(pos_sample_buf,  n_send_tmp, GetDataType<F64vec>(),
                           pos_sample_tot_, n_recv, n_recv_disp, GetDataType<F64vec>(), comm_sub_[0]);

            ///////////// sort particles along x direction again.
            ///////////// after sort, particles are sorted in global.
            // this part is the same as original
            std::sort(pos_sample_tot_, pos_sample_tot_+n_par_slab, LessOPX());

            
            ///////////////////////////////
	    ///////////// determine X coord
            // to determine global rank in x
            MPI_Allgather(&n_par_slab, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            number_of_sample_particle_tot_ = n_recv_disp[n_domain_[0]];

	    // get index of 
            S32 n_ave = number_of_sample_particle_tot_ / n_domain_[0];
            for(S32 i=0; i<n_domain_[0]; i++){
                i_head[i] = n_ave * i;
                if( i < number_of_sample_particle_tot_ % n_domain_[0]){
                    i_head[i] += i;
                }
                else{
                    i_head[i] += number_of_sample_particle_tot_ % n_domain_[0];
                }
                if(i > 0) i_tail[i-1] = i_head[i] - 1;
            }
            i_tail[n_domain_[0]-1] = number_of_sample_particle_tot_ - 1;
            n_send_tmp = 0; // temporally used
            for(S32 i=0; i<n_domain_[0]; i++){
                if( n_recv_disp[rank_1d_[0]] <= i_head[i] &&  i_head[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_head[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
                if( n_recv_disp[rank_1d_[0]] <= i_tail[i] &&  i_tail[i] < n_recv_disp[rank_1d_[0]]+n_par_slab){
                    S32 i_tmp = i_tail[i] - n_recv_disp[rank_1d_[0]];
                    coord_buf[n_send_tmp++] = pos_sample_tot_[i_tmp].x;
                }
            }

	    MPI_Allgather(&n_send_tmp, 1, GetDataType<S32>(),
                          n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
	    n_recv_disp[0] = 0;
	    for(S32 i=0; i<n_domain_[0]; i++){
            n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i]; 
	    }

	    MPI_Allgatherv(coord_buf, n_send_tmp, GetDataType<>(coord_buf[0]),
                       coord_tot, n_recv, n_recv_disp, GetDataType<>(coord_buf[0]), comm_1d_[0]);
	    
	    assert( n_recv_disp[n_domain_[0]] == n_domain_[0]*2);

	    // size of x_coord_buf is n_domain_[0]+1
	    x_coord[0] = pos_root_domain_.low_.x;
	    x_coord[n_domain_[0]] = pos_root_domain_.high_.x;

	    for(S32 i=1; i<n_domain_[0]; i++){
                x_coord[i] = (coord_tot[i*2] + coord_tot[i*2-1]) * 0.5;
	    }
            ////////////////////////////////////////////
            ///////////// migrate particles along x direction
            for(S32 i=0; i<n_domain_[0]; i++) n_send[i] = n_recv[i] = 0;
	    id_domain_x = 0;
        #ifdef RING_SYSTEM
            disp_domain_x_max_loc = 0;
            for(S32 i=0; i<n_par_slab; i++){
                while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
                if(n_send[id_domain_x] == 0){
                    S32 disp_domain_tmp = abs(id_domain_x-rank_1d_[0]);
                    disp_domain_tmp = (disp_domain_tmp < (n_domain_[0]-disp_domain_tmp) ) ? disp_domain_tmp : (n_domain_[0]-disp_domain_tmp);
                    if(disp_domain_x_max_loc < disp_domain_tmp) disp_domain_x_max_loc = disp_domain_tmp;
                }
                n_send[id_domain_x]++;
            }
            disp_domain_x_max_glb = 0;
            MPI_Allreduce(&disp_domain_x_max_loc, &disp_domain_x_max_glb, 1, MPI_INT, MPI_MAX, comm_1d_[0]);

            #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 4 @decomposeDomainMultiStep2"<<std::endl;
            #endif
            n_proc_send = 0;
            n_proc_recv = 0;
            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                MPI_Isend(n_send+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_send+n_proc_send);
                MPI_Irecv(n_recv+id_tmp, 1, MPI_INT, id_tmp, 0, comm_1d_[0], req_recv+n_proc_recv);
                n_proc_send++;
                n_proc_recv++;                
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 5 @decomposeDomainMultiStep2"<<std::endl;
            #endif
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
            n_proc_send = n_proc_recv = 0;

            for(S32 i=-disp_domain_x_max_glb; i<=disp_domain_x_max_glb; i++){
                S32 id_tmp = my_rank_1d+i;
                if(id_tmp < 0) id_tmp += n_domain_[0];
                else if(id_tmp >= n_domain_[0]) id_tmp -= n_domain_[0];
                if(n_send[id_tmp] > 0){
                    MPI_Isend(pos_sample_tot_+n_send_disp[id_tmp], n_send[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_send+n_proc_send);
                    n_proc_send++;
                }
                if(n_recv[id_tmp] > 0){
                    MPI_Irecv(pos_sample_buf+n_recv_disp[id_tmp], n_recv[id_tmp],
                              GetDataType<F64vec>(), id_tmp, 0, comm_1d_[0],
                              req_recv+n_proc_recv);
                    n_proc_recv++;
                }
            }
            MPI_Waitall(n_proc_send, req_send, stat_send);
            MPI_Waitall(n_proc_recv, req_recv, stat_recv);
            #ifdef DEBUG_DD_MULTI_STEP2
            if(Comm::getRank() == 0) std::cerr<<"OK 6 @decomposeDomainMultiStep2"<<std::endl;
            #endif            
        #else //RING_SYSTEM
	    for(S32 i=0; i<n_par_slab; i++){
                while( x_coord[id_domain_x+1] <= pos_sample_tot_[i].x ) id_domain_x++;
                n_send[id_domain_x]++;
            }
            MPI_Alltoall(n_send, 1, GetDataType<S32>(), n_recv, 1, GetDataType<S32>(), comm_1d_[0]);
            n_send_disp[0] = n_recv_disp[0] = 0;
            for(S32 i=0; i<n_domain_[0]; i++){
                n_send_disp[i+1] = n_send_disp[i] + n_send[i];
                n_recv_disp[i+1] = n_recv_disp[i] + n_recv[i];
            }
	    MPI_Alltoallv(pos_sample_tot_, n_send, n_send_disp, GetDataType<F64vec>(),
                          pos_sample_buf,  n_recv, n_recv_disp, GetDataType<F64vec>(), comm_1d_[0]);
        #endif //RING_SYSTEM
	    n_par_slab = n_recv_disp[n_domain_[0]];
	    // OK

    #if 1
	    ////////////////////////////////////
            ///////////// determine y corrdinate
	    std::sort(pos_sample_buf, pos_sample_buf+n_par_slab, LessOPY());
            if(split_y_at_the_center){
                assert(n_domain_[1]%2==0);
                F64vec * y_coord_mid = std::lower_bound(pos_sample_buf, pos_sample_buf+n_par_slab, y_coord_split, CompYDir());
                S32 number_of_sample_particle_sep_y_low = y_coord_mid - pos_sample_buf; // i.e. pos_sample_buf[number_of_sample_particle_sep_y_low-1] <= y_coord_split

                // lower part
                S32 n_ave_low = number_of_sample_particle_sep_y_low / (n_domain_[1]/2);
                for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                    i_head[iy] = n_ave_low * iy;
                    if( iy < number_of_sample_particle_sep_y_low % (n_domain_[1]/2) ){ i_head[iy] += iy;}
                    else{ i_head[iy] += number_of_sample_particle_sep_y_low % (n_domain_[1]/2);}
                    if(iy > 0) i_tail[iy-1] = i_head[iy] - 1;
                }
                i_tail[n_domain_[1]/2-1] = number_of_sample_particle_sep_y_low - 1;

                // higher part
                S32 number_of_sample_particle_sep_y_high = n_par_slab - number_of_sample_particle_sep_y_low;
                S32 n_ave_high = number_of_sample_particle_sep_y_high / (n_domain_[1]/2);
                for(S32 iy = 0; iy<n_domain_[1]/2; iy++) {
                    S32 iy_glb = iy+n_domain_[1]/2;
                    i_head[iy_glb] = n_ave_high * iy + number_of_sample_particle_sep_y_low;
                    if( iy < number_of_sample_particle_sep_y_high % (n_domain_[1]/2) ){ i_head[iy_glb] += iy;}
                    else{ i_head[iy_glb] += number_of_sample_particle_sep_y_high % (n_domain_[1]/2);}
                    if(iy > 0) i_tail[iy_glb-1] = i_head[iy_glb] - 1;
                }
                i_tail[n_domain_[1]-1] = n_par_slab - 1;
                
                y_coord[0] = pos_root_domain_.low_.y;
                y_coord[n_domain_[1]/2] = y_coord_split;
                y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
                for(S32 i=1; i<n_domain_[1]/2; i++){ y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;}
                for(S32 i=n_domain_[1]/2+1; i<n_domain_[1]; i++){ y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;}
                if(Comm::getRank() == 0){
                    for(S32 i=0; i<n_domain_[1]; i++){
                        std::cerr<<"i_head[i]= "<<i_head[i]<<" i_tail[i]= "<<i_tail[i]<<std::endl;
                    }
                    for(S32 i=0; i<n_domain_[1]+1; i++){
                        std::cerr<<"y_coord[i]= "<<y_coord[i]<<std::endl;
                    }
                }
            }
            else{
                // get index of
                n_ave = n_par_slab / n_domain_[1];
                for(S32 i=0; i<n_domain_[1]; i++){
                    i_head[i] = n_ave * i;
                    if( i < n_par_slab % n_domain_[1]) i_head[i] += i;
                    else i_head[i] += n_par_slab % n_domain_[1];
                    if(i > 0) i_tail[i-1] = i_head[i] - 1;
                }
                i_tail[n_domain_[1]-1] = n_par_slab - 1;
                // size of y_coord is n_domain_[1]+1
                y_coord[0] = pos_root_domain_.low_.y;
                y_coord[n_domain_[1]] = pos_root_domain_.high_.y;
                for(S32 i=1; i<n_domain_[1]; i++) y_coord[i] = (pos_sample_buf[i_head[i]].y + pos_sample_buf[i_tail[i-1]].y) * 0.5;
            }


        #ifndef PARTICLE_SIMULATOR_TWO_DIMENSION
	    ////////////////////////////////////
            ///////////// determine z corrdinate
	    //#pragma omp parallel for
	    for(S32 iy=0; iy<n_domain_[1]; iy++){
            const S32 iy_ptcl_head = i_head[iy];
            const S32 iy_ptcl_tail = i_tail[iy];
            const S32 nz_tot = iy_ptcl_tail - iy_ptcl_head + 1;
            std::sort(pos_sample_buf+iy_ptcl_head, pos_sample_buf+iy_ptcl_head+nz_tot, LessOPZ());
            S32 nz_ave_tmp = nz_tot / n_domain_[2];
            pos_domain_temp_buf[iy*n_domain_[2]].low_.z = pos_root_domain_.low_.z;
            pos_domain_temp_buf[(iy+1)*n_domain_[2] - 1].high_.z = pos_root_domain_.high_.z;
            pos_domain_temp_buf[iy*n_domain_[2]].low_.y = y_coord[iy];
            pos_domain_temp_buf[iy*n_domain_[2]].high_.y = y_coord[iy+1];
            for(S32 iz=1; iz<n_domain_[2]; iz++){
                pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.y = y_coord[iy];
                pos_domain_temp_buf[iy*n_domain_[2]+iz].high_.y = y_coord[iy+1];
                S32 iz_tmp = nz_ave_tmp * iz;
                if(iz < nz_tot % n_domain_[2]) iz_tmp += iz;
                else iz_tmp += nz_tot % n_domain_[2];
                F64 z_coord_tmp = (pos_sample_buf[iy_ptcl_head+iz_tmp].z + pos_sample_buf[iy_ptcl_head+iz_tmp-1].z) * 0.5;

                pos_domain_temp_buf[iy*n_domain_[2]+iz].low_.z = z_coord_tmp;
                pos_domain_temp_buf[iy*n_domain_[2]+iz-1].high_.z = z_coord_tmp;
            }
	    }
        #endif //PARTICLE_SIMULATOR_TWO_DIMENSION
	    for(S32 i=0; i<n_proc_sub_[0]; i++){
                pos_domain_temp_buf[i].low_.x = x_coord[rank_1d_[0]];
                pos_domain_temp_buf[i].high_.x = x_coord[rank_1d_[0]+1];
	    }

	    //////////////////////////////////////////////
            ///////////// exchange pos_domain_tmp
	    MPI_Allgather(pos_domain_temp_buf, n_proc_sub_[0], GetDataType<F64ort>(),
                      pos_domain_temp_, n_proc_sub_[0], GetDataType<F64ort>(), comm_1d_[0]);

	    if(first_call_by_decomposeDomain) {
                first_call_by_decomposeDomain = false;
                for(S32 i = 0; i < n_proc_glb; i++) {
                    //std::cout<<"pos_domain_temp_[i](first)= "<<pos_domain_temp_[i]<<std::endl;
                    pos_domain_[i].low_  = pos_domain_temp_[i].low_;
                    pos_domain_[i].high_ = pos_domain_temp_[i].high_;
                }
	    } else {
                for(S32 i = 0; i < n_proc_glb; i++) {
                    //std::cout<<"pos_domain_temp_[i](other)= "<<pos_domain_temp_[i]<<std::endl;
                    pos_domain_[i].low_  = (F64)coef_ema_ * pos_domain_temp_[i].low_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].low_;
                    pos_domain_[i].high_ = (F64)coef_ema_ * pos_domain_temp_[i].high_ 
                        + (F64)(1. - coef_ema_) * pos_domain_[i].high_;
                }
	    }
            //time_profile_.decompose_domain += GetWtime() - time_offset;
     #endif // PARTICLE_SIMULATOR_MPI_PARALLEL
	    time_profile_.decompose_domain = GetWtime() - time_offset;

#endif //#if 1
            
        }

    #endif



        
#endif
        
    };
}
