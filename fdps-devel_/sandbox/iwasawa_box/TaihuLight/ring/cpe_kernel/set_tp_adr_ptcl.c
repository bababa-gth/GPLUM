#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include"slave.h"
#include"dma.h"
#include"cpe_func.h"
#include"cpe_prof.h"

/* Local macros */
//#define CHECK_SET_PTCL_SORTED_GLB_CPE
//#define CHECK_SET_ADR_GLB_1ST
//#define CHECK_SET_ADR_GLB_2ND

#define VERSION_OF_SET_PTCL_SORTED_GLB_CPE (2)
#define VERSION_OF_SET_ADR_GLB_1ST (2)
#define VERSION_OF_SET_ADR_GLB_2ND (2)

extern int MY_RANK_MPI;

static inline U32_ GetMsb(U32_ val){
    return (val>>31) & 0x1;
}

static inline U32_ SetMsb(U32_ val){
    return val | 0x80000000;
}

static inline U32_ ClearMsb(U32_ val){
    return val & 0x7fffffff;
}

static inline void CopyKeyOnly(tpLM * src, tpLM * dst){
    dst->key_hi_ = src->key_hi_;
#ifdef USE_96BIT_KEY
    dst->key_lo_ = src->key_lo_;
#endif
}

/* CPE communication functions */
static inline void cpe_bcast_int32(const int root_cpe_id, int* data) {
    int my_id = athread_get_id(-1);
    int my_col_id, my_row_id;
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    int root_col_id, root_row_id;
    root_row_id = root_cpe_id/8;
    root_col_id = root_cpe_id%8;

    if (my_id == root_cpe_id) {
        REG_PUTR(*data,8);
        REG_PUTC(*data,8);
    }
    else {
        if (my_row_id == root_row_id) {
            REG_GETR(*data);
            REG_PUTC(*data,8);
        }
        else {
            REG_GETC(*data);
        }
    }
}

void SetTpAdrPtcl(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_tot = (int)(((unsigned long*)args)[0]);
    void * adr_tp = (void *)((unsigned long*)args)[1];
    int my_n = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp = sizeof(tpLM);
    volatile tpLM  * tp     = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    volatile tpLM  * tp_new = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    int my_n_ep = 0;
    int my_n_sp = 0;
    // just count number to determine offset to write
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            if(GetMsb(tp[j].adr_ptcl_)==0){
                my_n_ep++;
            }
            else{
                my_n_sp++;
            }
        }
    }
    /*
    sync_array_();
    for(i=0;i<my_id*10000000; i++) NOP();
    if(MY_RANK_MPI==0){
        printf("A) my_id=%d, my_n_ep=%d , my_n_sp=%d \n",
               my_id, my_n_ep, my_n_sp);
    }
    sync_array_();
    */
    
    int n_ep_ar[NUMBER_OF_CPE];
    int n_sp_ar[NUMBER_OF_CPE];
    for (i=0; i<NUMBER_OF_CPE; i++) {
        n_ep_ar[i] = n_sp_ar[i] = 0;
    }
    n_ep_ar[my_id] = my_n_ep;
    n_sp_ar[my_id] = my_n_sp;
    for (i=0; i<NUMBER_OF_CPE; i++) {
        cpe_bcast_int32(i, &n_ep_ar[i]);
        cpe_bcast_int32(i, &n_sp_ar[i]);
    }
    int my_n_disp_ep = 0;
    int my_n_disp_sp = 0;
    for(i=0; i<my_id; i++){
        my_n_disp_ep += n_ep_ar[i];
        my_n_disp_sp += n_sp_ar[i];
    }
    /*
    sync_array_();
    for(i=0;i<my_id*10000000; i++) NOP();
    if(MY_RANK_MPI==0){
        printf("B) my_id=%d, my_n_disp_ep=%d , my_n_disp_sp=%d \n",
               my_id, my_n_disp_ep, my_n_disp_sp);
    }
    sync_array_();
    */
    my_n_ep = my_n_sp = 0;
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            CopyKeyOnly(&(tp[j]), &(tp_new[j]));
            if(GetMsb(tp[j].adr_ptcl_)==0){
                tp_new[j].adr_ptcl_ = my_n_disp_ep+my_n_ep;
                my_n_ep++;
            }
            else{
                tp_new[j].adr_ptcl_ = SetMsb(my_n_disp_sp+my_n_sp);
                my_n_sp++;
            }
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp_new));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(tp_new, bsize_tp*CHUNK_SIZE);

}


void SetTpAdrPtclEpOnly(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_tot = (int)(((unsigned long*)args)[0]);
    void * adr_tp = (void *)((unsigned long*)args)[1];
    int my_n = n_tot/NUMBER_OF_CPE + ( (my_id < n_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_tot % NUMBER_OF_CPE) ? my_id : n_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp = sizeof(tpLM);
    volatile tpLM  * tp     = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    volatile tpLM  * tp_new = (tpLM *) ldm_malloc(bsize_tp*CHUNK_SIZE);
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            CopyKeyOnly(&(tp[j]), &(tp_new[j]));
            tp_new[j].adr_ptcl_ = my_offset+i+j;
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp_new));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(tp_new, bsize_tp*CHUNK_SIZE);
}

void slave_SetTpLocAdrPtclFromTpCpe(void * args){
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n = (int)(((unsigned long*)args)[0]);
    void * adr_tp       = (void *)((unsigned long*)args)[1];
    void * adr_adr_ptcl = (void *)((unsigned long*)args)[2];
    int my_n = n/NUMBER_OF_CPE + ( (my_id < n % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n/NUMBER_OF_CPE)*my_id + ( (my_id < n % NUMBER_OF_CPE) ? my_id : n % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_adr = sizeof(U64_);
    volatile tpLM  * tp      = (tpLM *) ldm_malloc(bsize_tp  * CHUNK_SIZE);
    volatile U64_  * adr_ptcl = (U64_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            adr_ptcl[j] = tp[j].adr_ptcl_;
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr * nn);
        dma(dma_put, (long*)((U64_ *)adr_adr_ptcl+my_offset+i), (long*)(adr_ptcl));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(adr_ptcl, bsize_adr*CHUNK_SIZE);
}

void SetAdrGlb1stCpe(void * args){
#if VERSION_OF_SET_ADR_GLB_1ST == 2
    volatile int my_id, my_col_id, my_row_id;
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    //-(for general DMA)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //-(for DMA get of adr_epj_org2glb)
    volatile dma_desc dma_get_adr;
    volatile int reply_get_adr = 0;
    dma_descriptor_init(&dma_get_adr, &reply_get_adr);
    //* Process the arguments
    int n                      = (int)(((unsigned long*)args)[0]);
    void * adr_tp              = (void *)((unsigned long*)args)[1];
    void * adr_adr_epj_org2glb = (void *)((unsigned long*)args)[2];
    void * adr_adr_epj_loc2glb = (void *)((unsigned long*)args)[3];
#ifdef CHECK_SET_ADR_GLB_1ST
    if (my_id == 0) {
        printf("n                   = %d\n",n);
        printf("adr_tp              = %lu\n",(unsigned long)adr_tp);
        printf("adr_adr_epj_org2glb = %lu\n",(unsigned long)adr_adr_epj_org2glb);
        printf("adr_adr_epj_loc2glb = %lu\n",(unsigned long)adr_adr_epj_loc2glb);
    }
#endif
    //* Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n/NUMBER_OF_CPE + ( (my_id < n % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n/NUMBER_OF_CPE)*my_id + ( (my_id < n % NUMBER_OF_CPE) ? my_id : n % NUMBER_OF_CPE );
#ifdef CHECK_SET_ADR_GLB_1ST
    if (my_id == 0) {
        printf("n_loc     = %d\n",n_loc);
        printf("my_offset = %d\n",my_offset);
    }
#endif
    //* Local variables
    //-(loop counters)
    int i,j,k,id;
    //-(local buffers)
    enum {
        CHUNK_SIZE = 64,
    };
#ifdef REMOVE_TP_LOC
    size_t bsize_tp = sizeof(U64_);
    size_t bsize_tp_chunk = bsize_tp * CHUNK_SIZE;
    volatile U64_ * tp = (U64_ *)  ldm_malloc(bsize_tp_chunk);
#else
    size_t bsize_tp = sizeof(tpLM);
    size_t bsize_tp_chunk = bsize_tp * CHUNK_SIZE;
    volatile tpLM * tp = (tpLM *)  ldm_malloc(bsize_tp_chunk);
#endif
    size_t bsize_adr = sizeof(U32_);
    size_t bsize_adr_chunk = bsize_adr * CHUNK_SIZE;
    volatile U32_ * adr_epj_loc2glb = (U32_ *) ldm_malloc(bsize_adr_chunk);
    //-(cache)
    volatile S32_ adr_adr_cache, n_adr_cache, n_adr_get;
    volatile U32_ * adr_epj_org2glb_cache = (U32_ *) ldm_malloc(bsize_adr_chunk);

    //* Update adr_epj_loc2glb
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
#ifdef CHECK_SET_ADR_GLB_1ST
        sync_array_();
        if (my_id == 0) {
            printf("i = %d\n",i);
        }
        sync_array_();
#endif
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        if (nn > 0) {
            //** Get tp[]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tp * nn);
#ifdef REMOVE_TP_LOC
            dma(dma_get, (long*)((U64_ *)adr_tp + my_offset + i), (long*)(tp));
#else
            dma(dma_get, (long*)((tpLM *)adr_tp + my_offset + i), (long*)(tp));
#endif
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Make adr_epj_org2glb_cache[]
#ifdef REMOVE_TP_LOC
            adr_adr_cache = tp[0];
#else
            adr_adr_cache = tp[0].adr_ptcl_;
#endif
            n_adr_cache   = nn;
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_adr * n_adr_cache);
            dma(dma_get, 
                (long*)((U32_ *)adr_adr_epj_org2glb + adr_adr_cache), 
                (long*)(adr_epj_org2glb_cache));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            //** Calculate adr_epj_loc2glb[]
            reply_get_adr = 0;
            dma_set_op(&dma_get_adr, DMA_GET);
            dma_set_mode(&dma_get_adr, PE_MODE);
            dma_set_reply(&dma_get_adr, &reply_get_adr);
            dma_set_size(&dma_get_adr, bsize_adr);
            n_adr_get = 0;
            for (j=0; j<nn; j++) {
#ifdef REMOVE_TP_LOC
                U32_ adr_org = tp[j]; // effectively, tp_loc_adr_ptcl_
#else
                U32_ adr_org = tp[j].adr_ptcl_;                            
#endif
                if ((adr_adr_cache <= adr_org) && (adr_org < adr_adr_cache + n_adr_cache)) {
                    adr_epj_loc2glb[j] = adr_epj_org2glb_cache[adr_org - adr_adr_cache];
                } else {
                    dma(dma_get_adr,
                        (long*)((U32_ *)adr_adr_epj_org2glb + adr_org),
                        (long*)((U32_ *)adr_epj_loc2glb + j));
                    n_adr_get++;
                }
            }
            if (n_adr_get > 0) {
                dma_wait(&reply_get, n_adr_get);
                while (reply_get != n_adr_get) {}
            }
            //** Put adr_epj_loc2glb[]
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  bsize_adr * nn);
            dma(dma_put, (long*)((U32_ *)adr_adr_epj_loc2glb + my_offset + i), (long*)(adr_epj_loc2glb));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }

    }

    //* Release memory
    ldm_free(tp, bsize_tp_chunk);
    ldm_free(adr_epj_loc2glb, bsize_adr_chunk);
    ldm_free(adr_epj_org2glb_cache, bsize_adr_chunk);
#elif VERSION_OF_SET_ADR_GLB_1ST == 1
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n = (int)(((unsigned long*)args)[0]);
    void * adr_tp       = (void *)((unsigned long*)args)[1];
    void * adr_adr_epj_org2glb = (void *)((unsigned long*)args)[2];
    void * adr_adr_epj_loc2glb = (void *)((unsigned long*)args)[3];
    int my_n = n / NUMBER_OF_CPE + ( (my_id < n % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset = (n / NUMBER_OF_CPE)*my_id + ( (my_id < n % NUMBER_OF_CPE) ? my_id : n % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
#ifdef REMOVE_TP_LOC
    size_t bsize_tp = sizeof(U64_);
    volatile U64_ * tp = (U64_ *)  ldm_malloc(bsize_tp * CHUNK_SIZE);
#else
    size_t bsize_tp = sizeof(tpLM);
    volatile tpLM * tp = (tpLM *)  ldm_malloc(bsize_tp * CHUNK_SIZE);
#endif

    size_t bsize_adr = sizeof(S32_);
    volatile S32_  * adr_epj_org2glb = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    volatile S32_  * adr_epj_loc2glb = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    volatile S32_ adr_org_tmp[CHUNK_SIZE];
    volatile S32_ adr_glb_tmp[CHUNK_SIZE];
    for (i=0; i<my_n; i+=CHUNK_SIZE) {
        int nrem = my_n - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        
        // get tp
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
#ifdef REMOVE_TP_LOC
        dma(dma_get, (long*)((U64_ *)adr_tp+my_offset+i), (long*)(tp));
#else
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset+i), (long*)(tp));
#endif
        while (reply_get != 1) {}
        dma_wait(&reply_get, 1);

        // prepare for getting adr_epj_org2glb
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_adr);
        for(j=0; j<nn; j++){
#ifdef REMOVE_TP_LOC
            adr_org_tmp[j] = (S32_)tp[j];
#else
            adr_org_tmp[j] = (S32_)tp[j].adr_ptcl_;
#endif
            dma(dma_get, (long*)((S32_ *)adr_adr_epj_org2glb+adr_org_tmp[j]),
                (long*)(&adr_glb_tmp[j]));
        }
        while (reply_get != nn) {}
        dma_wait(&reply_get, nn);

        // prepare for putting adr_epj_loc2glb
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr);
        for(j=0; j<nn; j++){
            dma(dma_put, (long*)((S32_ *)adr_adr_epj_loc2glb+my_offset+i+j),
                (long*)(&adr_glb_tmp[j]));
        }
        while (reply_put != nn) {}
        dma_wait(&reply_put, nn);
    }
    //* Release memory
    ldm_free(tp, bsize_tp*CHUNK_SIZE);
    ldm_free(adr_epj_org2glb, bsize_adr*CHUNK_SIZE);
    ldm_free(adr_epj_loc2glb, bsize_adr*CHUNK_SIZE);
#else
#error The value of `VERSION_OF_SET_ADR_GLB_1ST` is invalid.
#endif
}

void SetAdrGlb2ndCpe(void * args){
#if VERSION_OF_SET_ADR_GLB_2ND == 2
    volatile int my_id, my_col_id, my_row_id;
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    //-(for general DMA)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //-(for DMA get of adr_epj_org2glb)
    volatile dma_desc dma_get_adr;
    volatile int reply_get_adr = 0;
    dma_descriptor_init(&dma_get_adr, &reply_get_adr);
    //* Process the arguments
    int n_epj_sorted           = (int)(((unsigned long*)args)[0]);
    int n_adr_epj_loc2glb      = (int)(((unsigned long*)args)[1]);
    void * adr_adr_epj_glb2loc = (void *)((unsigned long*)args)[2];
    void * adr_adr_epj_loc2glb = (void *)((unsigned long*)args)[3];
#ifdef CHECK_SET_ADR_GLB_2ND
    if (my_id == 0) {
        printf("n_epj_sorted        = %d\n",n_epj_sorted);
        printf("n_adr_epj_loc2glb   = %d\n",n_adr_epj_loc2glb);
        printf("adr_adr_epj_glb2loc = %lu\n",(unsigned long)adr_adr_epj_glb2loc);
        printf("adr_adr_epj_loc2glb = %lu\n",(unsigned long)adr_adr_epj_loc2glb);
    }
#endif
    //* Local variables
    //-(loop counters)
    int i,j,k,id;
    //-(local buffers)
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_adr = sizeof(U32_);
    size_t bsize_adr_chunk = bsize_adr * CHUNK_SIZE;
    volatile U32_ * adr_epj_glb2loc = (U32_ *) ldm_malloc(bsize_adr_chunk);
    volatile U32_ * adr_epj_loc2glb = (U32_ *) ldm_malloc(bsize_adr_chunk);
    volatile S32_ * adr_adr_epj_glb2loc_ind = (S32_ *) ldm_malloc(bsize_adr_chunk);
    //-(group info.)
    volatile S32_  n_groups;
    volatile S32_ * adr_group_head  = (S32_ *) ldm_malloc(bsize_adr_chunk);
    volatile S32_ * n_group_members = (S32_ *) ldm_malloc(bsize_adr_chunk);

    //* Initialize adr_epj_glb2loc[] with -1 
    //** Compute the task of each CPE
    int n_loc,my_offset;
    n_loc = n_epj_sorted/NUMBER_OF_CPE + ( (my_id < n_epj_sorted % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_epj_sorted/NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted % NUMBER_OF_CPE) ? my_id : n_epj_sorted % NUMBER_OF_CPE );
#ifdef CHECK_SET_ADR_GLB_2ND
    if (my_id == 0) {
        printf("n_loc     = %d\n",n_loc);
        printf("my_offset = %d\n",my_offset);
    }
#endif
    //** Put -1 into  adr_epj_glb2loc[]
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
#ifdef CHECK_SET_ADR_GLB_2ND
        //sync_array_();
        //if (my_id == 0) {
        //    printf("i(1) = %d\n",i);
        //}
        //sync_array_();
#endif
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        if (nn > 0) {
            // Make a put buffer
            for (j=0; j<nn; j++)
                adr_epj_glb2loc[j] = -1;
            //** Put adr_epj_glb2loc[]
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  bsize_adr * nn);
            dma(dma_put,
                (long*)((U32_ *)adr_adr_epj_glb2loc + my_offset + i),
                (long*)(adr_epj_glb2loc));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
        }
    }
     
    //* Update adr_epj_glb2loc[]
    //** Compute the task of each CPE
    n_loc = n_adr_epj_loc2glb/NUMBER_OF_CPE + ( (my_id < n_adr_epj_loc2glb % NUMBER_OF_CPE) ? 1 : 0 );
    my_offset = (n_adr_epj_loc2glb/NUMBER_OF_CPE)*my_id + ( (my_id < n_adr_epj_loc2glb % NUMBER_OF_CPE) ? my_id : n_adr_epj_loc2glb % NUMBER_OF_CPE );
#ifdef CHECK_SET_ADR_GLB_2ND
    if (my_id == 0) {
        printf("n_loc     = %d\n",n_loc);
        printf("my_offset = %d\n",my_offset);
    }
#endif
    //** Update
    for (i=0; i<n_loc; i+=CHUNK_SIZE) {
#ifdef CHECK_SET_ADR_GLB_2ND
        //sync_array_();
        //if (my_id == 0) {
        //    printf("i(2) = %d\n",i);
        //}
        //sync_array_();
#endif
        int nrem = n_loc - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        if (nn > 0) {
            // Get adr_epj_loc2glb[]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_adr * nn);
            dma(dma_get,
                (long*)((U32_ *)adr_adr_epj_loc2glb + my_offset + i),
                (long*)(adr_epj_loc2glb));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            // Calculate adr_epj_glb2loc[]
            for (j=0; j<nn; j++) {
                U32_ adr = adr_epj_loc2glb[j];
                adr_adr_epj_glb2loc_ind[j] = adr;
                adr_epj_glb2loc[j] = my_offset + i + j;
            }
            // Grouping
            n_groups = 0;
            S32_ adr_prev = -1;
            for (j=0; j<nn; j++) {
                S32_ adr = adr_adr_epj_glb2loc_ind[j];
                if ((adr_prev == -1) ||
                    ((adr_prev != -1) && (adr != adr_prev + 1))) {
                    // a new group is found.
                    n_groups++;
                    adr_group_head[n_groups-1] = j;
                    n_group_members[n_groups-1] = 1;
                } else {
                    n_group_members[n_groups-1]++;
                }
                adr_prev = adr;
            }
            // Put adr_epj_glb2loc[]
#if 1
            for (j=0; j<n_groups; j++) {
                U32_ jloc = adr_group_head[j];
                U32_ adr = adr_adr_epj_glb2loc_ind[jloc];
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  bsize_adr * n_group_members[j]);
                dma(dma_put,
                    (long*)((U32_ *)adr_adr_epj_glb2loc + adr),
                    (long*)(&adr_epj_glb2loc[jloc]));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
            }
#else
            //[Put individually]
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  bsize_adr);
            for (j=0; j<nn; j++) {
                U32_ adr = adr_adr_epj_glb2loc_ind[j];
                dma(dma_put,
                    (long*)((U32_ *)adr_adr_epj_glb2loc + adr),
                    (long*)(&adr_epj_glb2loc[j]));
            }
            dma_wait(&reply_put, nn);
            while (reply_put != nn) {}
#endif
        }

    }

    //* Release memory
    ldm_free(adr_epj_glb2loc, bsize_adr_chunk);
    ldm_free(adr_epj_loc2glb, bsize_adr_chunk);
    ldm_free(adr_adr_epj_glb2loc_ind, bsize_adr_chunk);
    ldm_free(adr_group_head, bsize_adr_chunk);
    ldm_free(n_group_members, bsize_adr_chunk);
#elif VERSION_OF_SET_ADR_GLB_2ND == 1
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    int my_id = athread_get_id(-1);
    int n_epj_sorted = (int)(((unsigned long*)args)[0]);
    int n_adr = (int)(((unsigned long*)args)[1]);
    void * adr_adr_epj_glb2loc = (void *)((unsigned long*)args)[2];
    void * adr_adr_epj_loc2glb = (void *)((unsigned long*)args)[3];
    int my_n_epj_sorted = n_epj_sorted/NUMBER_OF_CPE
        + ( (my_id < n_epj_sorted % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset_epj_sorted = (n_epj_sorted / NUMBER_OF_CPE)*my_id + ( (my_id < n_epj_sorted % NUMBER_OF_CPE) ? my_id : n_epj_sorted % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_adr = sizeof(S32_);
    volatile S32_  * adr_epj_glb2loc = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    volatile S32_  * adr_epj_loc2glb = (S32_ *)  ldm_malloc(bsize_adr * CHUNK_SIZE);
    for (i=0; i<my_n_epj_sorted; i+=CHUNK_SIZE) {
        int nrem = my_n_epj_sorted - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_adr * nn);
        dma(dma_get, (long*)((S32_ *)adr_adr_epj_glb2loc+my_offset_epj_sorted+i),
            (long*)(adr_epj_glb2loc));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            adr_epj_glb2loc[j] = -1;
        }
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr * nn);
        dma(dma_put, (long*)((S32_ *)adr_adr_epj_glb2loc+my_offset_epj_sorted+i),
            (long*)(adr_epj_glb2loc));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
    }
    sync_array_();
    // buffer for DMA
    volatile S32_ adr_tmp[CHUNK_SIZE]; 
    volatile S32_ val_tmp[CHUNK_SIZE];
    int my_n_adr = n_adr / NUMBER_OF_CPE
        + ( (my_id < n_adr % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset_adr = (n_adr / NUMBER_OF_CPE)*my_id + ( (my_id < n_adr % NUMBER_OF_CPE) ? my_id : n_adr % NUMBER_OF_CPE );
    for (i=0; i<my_n_adr; i+=CHUNK_SIZE) {
        int nrem = my_n_adr - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        // get adr_epj_loc2glb
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_adr * nn);
        dma(dma_get, (long*)((S32_ *)adr_adr_epj_loc2glb+my_offset_adr+i),
            (long*)(adr_epj_loc2glb));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}

        // prepare for puting adr_epj_glb2loc
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_adr);
        for (j=0; j<nn; j++) {
            adr_tmp[j] = adr_epj_loc2glb[j];
            val_tmp[j] = my_offset_adr + i + j;
            dma(dma_put, (long*)((S32_ *)adr_adr_epj_glb2loc+adr_tmp[j]),
                (long*)(&val_tmp[j]));
        }

        dma_wait(&reply_put, nn);
        while (reply_put != nn) {}
    }

    //* Release memory
    ldm_free(adr_epj_glb2loc, bsize_adr*CHUNK_SIZE);
    ldm_free(adr_epj_loc2glb, bsize_adr*CHUNK_SIZE);
#else
#error The value of `VERSION_OF_SET_ADR_GLB_2ND` is invalid.
#endif
}


void SetPtclSortedGlbCpe(void * args){
#if VERSION_OF_SET_PTCL_SORTED_GLB_CPE == 2
    volatile int my_id, my_col_id, my_row_id;
    my_id = athread_get_id(-1);
    get_col_id_(&my_col_id);
    get_row_id_(&my_row_id);
    //* Initialize local vars. for DMA comm.
    //-(for general DMA)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);
    //-(for DMA get of epj_org)
    volatile dma_desc dma_get_epj_org;
    volatile int reply_get_epj_org = 0;
    dma_descriptor_init(&dma_get_epj_org, &reply_get_epj_org);
    //-(for DMA get of spj_org)
    volatile dma_desc dma_get_spj_org;
    volatile int reply_get_spj_org = 0;
    dma_descriptor_init(&dma_get_spj_org, &reply_get_spj_org);
    //* Process the arguments
    int n_glb_tot = (int)(((unsigned long*)args)[0]);
    int n_loc_tot = (int)(((unsigned long*)args)[1]);
    void * adr_tp = (void *)((unsigned long*)args)[2];
    void * adr_epj_sorted = (void *)((unsigned long*)args)[3];
    void * adr_epj_org = (void *)((unsigned long*)args)[4];
    void * adr_spj_sorted = (void *)((unsigned long*)args)[5];
    void * adr_spj_org = (void *)((unsigned long*)args)[6];
    void * adr_adr_epj_org2glb = (void *)((unsigned long*)args)[7];
    void * adr_adr_epj_buf2glb = (void *)((unsigned long*)args)[8];
    void * adr_adr_spj_buf2glb = (void *)((unsigned long*)args)[9];
    int cap_epj_org = (int)(((unsigned long*)args)[10]);
    int cap_spj_org = (int)(((unsigned long*)args)[11]);
  
    //* Local variables
    //-(Loop counters)
    int i,j,k;
    //-(prefix sum calc.)
    volatile int beg,end,psum;
    volatile int n_cnt_ep_tot, n_cnt_sp_tot;
    volatile int n_cnt_ep_loc, n_cnt_sp_loc;
    volatile int n_disp_ep_loc, n_disp_sp_loc;
    //-(local buffers)
    enum {
        CHUNK_SIZE = 64,
        CYCLE_SIZE = CHUNK_SIZE * NUMBER_OF_CPE,
    };
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_adr = sizeof(S32_);
    size_t bsize_epj = sizeof(epjLM);
    size_t bsize_spj = sizeof(spjLM);
    size_t bsize_tp_chunk  = bsize_tp  * CHUNK_SIZE;
    size_t bsize_adr_chunk = bsize_adr * CHUNK_SIZE;
    size_t bsize_epj_chunk = bsize_epj * CHUNK_SIZE;
    size_t bsize_spj_chunk = bsize_spj * CHUNK_SIZE;
    volatile tpLM * tp          = (tpLM *) ldm_malloc( bsize_tp_chunk );
    volatile epjLM * epj_sorted = (epjLM *) ldm_malloc( bsize_epj_chunk );
    volatile spjLM * spj_sorted = (spjLM *) ldm_malloc( bsize_spj_chunk );
    volatile U32_ * adr_epj_org2glb = (U32_ *) ldm_malloc( bsize_adr_chunk );
    volatile U32_ * adr_epj_buf2glb = (U32_ *) ldm_malloc( bsize_adr_chunk );
    volatile U32_ * adr_spj_buf2glb = (U32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ * adr_adr_epj_org2glb_ind = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ * adr_adr_epj_buf2glb_ind = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ * adr_adr_spj_buf2glb_ind = (S32_ *) ldm_malloc( bsize_adr_chunk );
    //-(cache)
    volatile S32_ adr_ep_cache, n_ep_cache, n_ep_get;
    volatile S32_ adr_sp_cache, n_sp_cache, n_sp_get;
    volatile epjLM * epj_org_cache = (epjLM *) ldm_malloc( bsize_epj_chunk );
    volatile spjLM * spj_org_cache = (spjLM *) ldm_malloc( bsize_spj_chunk );
    //-(group info.)
    volatile S32_ n_groups_epj_org2glb;
    volatile S32_ * adr_group_head_epj_org2glb  = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ * n_group_members_epj_org2glb = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ n_groups_epj_buf2glb;
    volatile S32_ * adr_group_head_epj_buf2glb  = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ * n_group_members_epj_buf2glb = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ n_groups_spj_buf2glb;
    volatile S32_ * adr_group_head_spj_buf2glb  = (S32_ *) ldm_malloc( bsize_adr_chunk );
    volatile S32_ * n_group_members_spj_buf2glb = (S32_ *) ldm_malloc( bsize_adr_chunk );
   
    //* Copy *_org to *_sorted and update tp[] and adr_*[]
    n_cnt_ep_tot = n_cnt_sp_tot = 0;
    for (i=0; i<n_glb_tot; i+=CYCLE_SIZE) {
#ifdef CHECK_SET_PTCL_SORTED_GLB_CPE
        sync_array_();
        if (my_id == 0) {
            printf("i = %d\n",i);
        }
        sync_array_();
#endif
        int nrem = n_glb_tot - i;
        int nn   = nrem < CYCLE_SIZE ? nrem : CYCLE_SIZE;
        //** Compute the task of each CPE
        int n_loc,my_offset;
        n_loc = nn/NUMBER_OF_CPE + ( (my_id < nn % NUMBER_OF_CPE) ? 1 : 0 );
        my_offset = i + (nn/NUMBER_OF_CPE)*my_id + ( (my_id < nn % NUMBER_OF_CPE) ? my_id : nn % NUMBER_OF_CPE );
        //** Get tp[] and count the numbers of EPJ and SPJ
        n_cnt_ep_loc = n_cnt_sp_loc = 0;
        adr_ep_cache = adr_sp_cache = -1;
        if (n_loc > 0) {
            // Get tp[]
            reply_get = 0;
            dma_set_op(&dma_get, DMA_GET);
            dma_set_mode(&dma_get, PE_MODE);
            dma_set_reply(&dma_get, &reply_get);
            dma_set_size(&dma_get, bsize_tp * n_loc);
            dma(dma_get, (long*)((tpLM *)adr_tp + my_offset), (long*)(tp));
            dma_wait(&reply_get, 1);
            while (reply_get != 1) {}
            // Count the numbers of EPJ and SPJ
            for (j=0; j<n_loc; j++) {
                if (GetMsb(tp[j].adr_ptcl_) == 0) {
                    if (adr_ep_cache == -1)
                        adr_ep_cache = tp[j].adr_ptcl_;
                    n_cnt_ep_loc++;
                } else {
                    if (adr_sp_cache == -1)
                        adr_sp_cache = ClearMsb(tp[j].adr_ptcl_);
                    n_cnt_sp_loc++;
                }
            }
            assert(n_cnt_ep_loc + n_cnt_sp_loc == n_loc);
        }
#ifdef CHECK_SET_PTCL_SORTED_GLB_CPE
        sync_array_();
        if (my_id == 0) {
            printf("Counts # of EPJ and SPJ completed.\n");
        }
        sync_array_();
#endif
        //** Compute the prefix sum (Collective CPE communication)
        prefix_sum(n_cnt_ep_loc,&beg,&end,&psum);
        n_disp_ep_loc = n_cnt_ep_tot + beg;
        n_cnt_ep_tot += psum; // for the next cycle
        prefix_sum(n_cnt_sp_loc,&beg,&end,&psum);
        n_disp_sp_loc = n_cnt_sp_tot + beg;
        n_cnt_sp_tot += psum; // for the next cycle
#ifdef CHECK_SET_PTCL_SORTED_GLB_CPE
        sync_array_();
        if (my_id == 0) {
            printf("Prefix sum calculation complected.\n");
        }
        sync_array_();
#endif
        //** Pack *_sorted[] & adr_*[] and update tp[]
        n_cnt_ep_loc = n_cnt_sp_loc = 0; // reset
        reply_get_epj_org = reply_get_spj_org = 0;
        dma_set_op(&dma_get_epj_org, DMA_GET);
        dma_set_mode(&dma_get_epj_org, PE_MODE);
        dma_set_reply(&dma_get_epj_org, &reply_get_epj_org);
        dma_set_size(&dma_get_epj_org, bsize_epj);
        dma_set_op(&dma_get_spj_org, DMA_GET);
        dma_set_mode(&dma_get_spj_org, PE_MODE);
        dma_set_reply(&dma_get_spj_org, &reply_get_spj_org);
        dma_set_size(&dma_get_spj_org, bsize_spj);
        if (n_loc > 0) {
#ifdef CHECK_SET_PTCL_SORTED_GLB_CPE
            if (my_id == 0) {
                printf("Update tp[] etc. started.\n");
            }
#endif
            // Make epj_org_cache & spj_org_cache
            if (adr_ep_cache != -1) {
                n_ep_cache = n_loc;
                if (n_loc > cap_epj_org - adr_ep_cache)
                    n_ep_cache = cap_epj_org - adr_ep_cache; 
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_epj * n_ep_cache);
                dma(dma_get, (long*)((epjLM *)adr_epj_org + adr_ep_cache), (long*)(epj_org_cache));
                dma_wait(&reply_get, 1);
                while (reply_get != 1) {}
            }
            if (adr_sp_cache != -1) {
                n_sp_cache = n_loc;
                if (n_loc > cap_spj_org - adr_sp_cache)
                    n_sp_cache = cap_spj_org - adr_sp_cache;
                reply_get = 0;
                dma_set_op(&dma_get, DMA_GET);
                dma_set_mode(&dma_get, PE_MODE);
                dma_set_reply(&dma_get, &reply_get);
                dma_set_size(&dma_get, bsize_spj * n_sp_cache);
                dma(dma_get, (long*)((spjLM *)adr_spj_org + adr_sp_cache), (long*)(spj_org_cache));
                dma_wait(&reply_get, 1);
                while (reply_get != 1) {}
            }
            
            n_ep_get = n_sp_get = 0;
            for (j=0; j<n_loc; j++) {
                U32_ adr = tp[j].adr_ptcl_;
                U32_ n_cnt_ep = n_disp_ep_loc + n_cnt_ep_loc;
                U32_ n_cnt_sp = n_disp_sp_loc + n_cnt_sp_loc;
                adr_adr_epj_org2glb_ind[j] = -1;
                adr_adr_epj_buf2glb_ind[j] = -1;
                adr_adr_spj_buf2glb_ind[j] = -1;
                if (GetMsb(adr) == 0) {
                    // Get epj_org 
                    if ((adr_ep_cache <= adr) && (adr < adr_ep_cache + n_ep_cache)) {
                        epj_sorted[n_cnt_ep_loc] = epj_org_cache[adr - adr_ep_cache];
                    } else {
                        dma(dma_get_epj_org,
                            (long*)((epjLM*)adr_epj_org + adr),
                            (long*)((epjLM*)epj_sorted + n_cnt_ep_loc));
                        n_ep_get++;
                    }
                    // Update tp[]
                    tp[j].adr_ptcl_ = n_cnt_ep;
                    // Set adr*[]
                    if (adr < n_loc_tot) { // particles are in own process
                        adr_adr_epj_org2glb_ind[j] = adr;
                        adr_epj_org2glb[j]         = n_cnt_ep;
                    } else {
                        adr_adr_epj_buf2glb_ind[j] = adr - n_loc_tot;
                        adr_epj_buf2glb[j]         = n_cnt_ep;
                    }
                    n_cnt_ep_loc++;
                } else {
                    U32_ adr_new = ClearMsb(adr);
                    // Get spj_org
                    if ((adr_sp_cache <= adr_new) && (adr_new < adr_sp_cache + n_sp_cache)) {
                        spj_sorted[n_cnt_sp_loc] = spj_org_cache[adr_new - adr_sp_cache];
                    } else {
                        dma(dma_get_spj_org,
                            (long*)((spjLM*)adr_spj_org + adr_new),
                            (long*)((spjLM*)spj_sorted + n_cnt_sp_loc));
                        n_sp_get++;
                    }
                    // Set adr*[]
                    adr_adr_spj_buf2glb_ind[j] = adr_new;
                    adr_spj_buf2glb[j]         = n_cnt_sp;
                    // Update tp[]
                    tp[j].adr_ptcl_ = SetMsb(n_cnt_sp);
                    n_cnt_sp_loc++;
                }
            }
            assert(n_cnt_ep_loc + n_cnt_sp_loc == n_loc);
            if (n_ep_get > 0) {
                dma_wait(&reply_get_epj_org, n_ep_get);
                while (reply_get_epj_org != n_ep_get) {}
            }
            if (n_sp_get > 0) {
                dma_wait(&reply_get_spj_org, n_sp_get);
                while (reply_get_spj_org != n_sp_get) {}
            }
#ifdef CHECK_SET_PTCL_SORTED_GLB_CPE
            if (my_id == 0) {
                printf("DMA put tp[] etc. started.\n");
            }
#endif
            // Put tp[], epj_sorted[], spj_sorted[]
            //-(tp)
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  sizeof(tpLM) * n_loc);
            dma(dma_put, (long*)((tpLM *)adr_tp + my_offset), (long*)(tp));
            dma_wait(&reply_put, 1);
            while (reply_put != 1) {}
            //-(epj_sorted)
            if (n_cnt_ep_loc > 0) {
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  sizeof(epjLM) * n_cnt_ep_loc);
                dma(dma_put, (long*)((epjLM *)adr_epj_sorted + n_disp_ep_loc), (long*)(epj_sorted));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
            }
            //-(spj_sorted)
            if (n_cnt_sp_loc > 0) {
                reply_put = 0;
                dma_set_op(&dma_put,    DMA_PUT);
                dma_set_mode(&dma_put,  PE_MODE);
                dma_set_reply(&dma_put, &reply_put);
                dma_set_size(&dma_put,  sizeof(spjLM) * n_cnt_sp_loc);
                dma(dma_put, (long*)((spjLM *)adr_spj_sorted + n_disp_sp_loc), (long*)(spj_sorted));
                dma_wait(&reply_put, 1);
                while (reply_put != 1) {}
            }
#if 1
            // Grouping
            S32_ j_prev_epj_org2glb = -10;
            S32_ j_prev_epj_buf2glb = -10;
            S32_ j_prev_spj_buf2glb = -10;
            S32_ adr_prev_epj_org2glb = -1;
            S32_ adr_prev_epj_buf2glb = -1;
            S32_ adr_prev_spj_buf2glb = -1;
            n_groups_epj_org2glb = 0;
            n_groups_epj_buf2glb = 0;
            n_groups_spj_buf2glb = 0;
            for (j=0; j<n_loc; j++) {
                // grouping adr_epj_org2glb[]
                volatile S32_ adr = adr_adr_epj_org2glb_ind[j];
                if (adr != -1) {
                    if (adr_prev_epj_org2glb == -1) {
                        // a new group is found.
                        n_groups_epj_org2glb++;
                        adr_group_head_epj_org2glb[n_groups_epj_org2glb-1] = j;
                        n_group_members_epj_org2glb[n_groups_epj_org2glb-1] = 1;
                    } else {
                        if (!((adr == adr_prev_epj_org2glb + 1) &&
                            (j == j_prev_epj_org2glb + 1))) {
                            // a new group is found.
                            n_groups_epj_org2glb++;
                            adr_group_head_epj_org2glb[n_groups_epj_org2glb-1] = j;
                            n_group_members_epj_org2glb[n_groups_epj_org2glb-1] = 1;
                        } else {
                            n_group_members_epj_org2glb[n_groups_epj_org2glb-1]++;  
                        }
                    }
                    j_prev_epj_org2glb = j;
                    adr_prev_epj_org2glb = adr;
                }
                // grouping adr_epj_buf2glb[]
                adr = adr_adr_epj_buf2glb_ind[j];
                if (adr != -1) {
                    if (adr_prev_epj_buf2glb == -1) {
                        // a new group is found.
                        n_groups_epj_buf2glb++;
                        adr_group_head_epj_buf2glb[n_groups_epj_buf2glb-1] = j;
                        n_group_members_epj_buf2glb[n_groups_epj_buf2glb-1] = 1;
                    } else {
                        if (!((adr == adr_prev_epj_buf2glb + 1) &&
                            (j == j_prev_epj_buf2glb + 1))) {
                            // a new group is found.
                            n_groups_epj_buf2glb++;
                            adr_group_head_epj_buf2glb[n_groups_epj_buf2glb-1] = j;
                            n_group_members_epj_buf2glb[n_groups_epj_buf2glb-1] = 1;
                        } else {
                            n_group_members_epj_buf2glb[n_groups_epj_buf2glb-1]++;  
                        }
                    }
                    j_prev_epj_buf2glb = j;
                    adr_prev_epj_buf2glb = adr;
                }
                // grouping adr_spj_buf2glb[]
                adr = adr_adr_spj_buf2glb_ind[j];
                if (adr != -1) {
                    if (adr_prev_spj_buf2glb == -1) {
                        // a new group is found.
                        n_groups_spj_buf2glb++;
                        adr_group_head_spj_buf2glb[n_groups_spj_buf2glb-1] = j;
                        n_group_members_spj_buf2glb[n_groups_spj_buf2glb-1] = 1;
                    } else {
                        if (!((adr == adr_prev_spj_buf2glb + 1) &&
                            (j == j_prev_spj_buf2glb + 1))) {
                            // a new group is found.
                            n_groups_spj_buf2glb++;
                            adr_group_head_spj_buf2glb[n_groups_spj_buf2glb-1] = j;
                            n_group_members_spj_buf2glb[n_groups_spj_buf2glb-1] = 1;
                        } else {
                            n_group_members_spj_buf2glb[n_groups_spj_buf2glb-1]++;  
                        }
                    }
                    j_prev_spj_buf2glb = j;
                    adr_prev_spj_buf2glb = adr;
                }
            }
            // Put adr*[]
            if (n_groups_epj_org2glb > 0) {
                for (j=0; j<n_groups_epj_org2glb; j++) {
                    U32_ jloc = adr_group_head_epj_org2glb[j];
                    U32_ adr  = adr_adr_epj_org2glb_ind[jloc];
                    reply_put = 0;
                    dma_set_op(&dma_put,    DMA_PUT);
                    dma_set_mode(&dma_put,  PE_MODE);
                    dma_set_reply(&dma_put, &reply_put);
                    dma_set_size(&dma_put,  bsize_adr * n_group_members_epj_org2glb[j]);
                    dma(dma_put,
                        (long*)((U32_ *)adr_adr_epj_org2glb + adr),
                        (long*)(&adr_epj_org2glb[jloc]));
                    dma_wait(&reply_put, 1);
                    while (reply_put != 1) {}
                }
            }
            if (n_groups_epj_buf2glb > 0) {
                for (j=0; j<n_groups_epj_buf2glb; j++) {
                    U32_ jloc = adr_group_head_epj_buf2glb[j];
                    U32_ adr  = adr_adr_epj_buf2glb_ind[jloc];
                    reply_put = 0;
                    dma_set_op(&dma_put,    DMA_PUT);
                    dma_set_mode(&dma_put,  PE_MODE);
                    dma_set_reply(&dma_put, &reply_put);
                    dma_set_size(&dma_put,  bsize_adr * n_group_members_epj_buf2glb[j]);
                    dma(dma_put,
                        (long*)((U32_ *)adr_adr_epj_buf2glb + adr),
                        (long*)(&adr_epj_buf2glb[jloc]));
                    dma_wait(&reply_put, 1);
                    while (reply_put != 1) {}
                }
            }
            if (n_groups_spj_buf2glb > 0) {
                for (j=0; j<n_groups_spj_buf2glb; j++) {
                    U32_ jloc = adr_group_head_spj_buf2glb[j];
                    U32_ adr  = adr_adr_spj_buf2glb_ind[jloc];
                    reply_put = 0;
                    dma_set_op(&dma_put,    DMA_PUT);
                    dma_set_mode(&dma_put,  PE_MODE);
                    dma_set_reply(&dma_put, &reply_put);
                    dma_set_size(&dma_put,  bsize_adr * n_group_members_spj_buf2glb[j]);
                    dma(dma_put,
                        (long*)((U32_ *)adr_adr_spj_buf2glb + adr),
                        (long*)(&adr_spj_buf2glb[jloc]));
                    dma_wait(&reply_put, 1);
                    while (reply_put != 1) {}
                }
            }
#else
            // Put adr*[]
            reply_put = 0;
            dma_set_op(&dma_put,    DMA_PUT);
            dma_set_mode(&dma_put,  PE_MODE);
            dma_set_reply(&dma_put, &reply_put);
            dma_set_size(&dma_put,  sizeof(U32_));
            volatile S32_ adr;
            for (j=0; j<n_loc; j++) {
                adr = adr_adr_epj_org2glb_ind[j];
                if (adr != -1) {
                    dma(dma_put,
                        (long*)((U32_ *)adr_adr_epj_org2glb + adr),
                        (long*)(&adr_epj_org2glb[j]));
                }
                adr = adr_adr_epj_buf2glb_ind[j];
                if (adr != -1) {
                    dma(dma_put,
                        (long*)((U32_ *)adr_adr_epj_buf2glb + adr),
                        (long*)(&adr_epj_buf2glb[j]));
                }
                adr = adr_adr_spj_buf2glb_ind[j];
                if (adr != -1) {
                    dma(dma_put,
                        (long*)((U32_ *)adr_adr_spj_buf2glb + adr),
                        (long*)(&adr_spj_buf2glb[j]));
                }
            }
            dma_wait(&reply_put, n_loc);
            while (reply_put != n_loc) {}
#endif
#ifdef CHECK_SET_PTCL_SORTED_GLB_CPE
            if (my_id == 0) {
                printf("DMA put tp[] etc. completed.\n");
            }
#endif
        }
   }

   //* Release memory
   ldm_free(tp, bsize_tp_chunk);
   ldm_free(epj_sorted, bsize_epj_chunk);
   ldm_free(spj_sorted, bsize_spj_chunk);
   ldm_free(adr_epj_org2glb, bsize_adr_chunk);
   ldm_free(adr_epj_buf2glb, bsize_adr_chunk);
   ldm_free(adr_spj_buf2glb, bsize_adr_chunk);
   ldm_free(adr_adr_epj_org2glb_ind, bsize_adr_chunk);
   ldm_free(adr_adr_epj_buf2glb_ind, bsize_adr_chunk);
   ldm_free(adr_adr_spj_buf2glb_ind, bsize_adr_chunk);
 
   ldm_free(epj_org_cache, bsize_epj_chunk);
   ldm_free(spj_org_cache, bsize_spj_chunk);

   ldm_free(adr_group_head_epj_org2glb, bsize_adr_chunk);
   ldm_free(n_group_members_epj_org2glb, bsize_adr_chunk);
   ldm_free(adr_group_head_epj_buf2glb, bsize_adr_chunk);
   ldm_free(n_group_members_epj_buf2glb, bsize_adr_chunk);
   ldm_free(adr_group_head_spj_buf2glb, bsize_adr_chunk);
   ldm_free(n_group_members_spj_buf2glb, bsize_adr_chunk);
#elif VERSION_OF_SET_PTCL_SORTED_GLB_CPE == 1
    U32_ i,j,k; 
    // (local vars. for DMA comm.)
    volatile dma_desc dma_get, dma_put;
    volatile int reply_get = 0;
    volatile int reply_put = 0;
    dma_descriptor_init(&dma_get, &reply_get);
    dma_descriptor_init(&dma_put, &reply_put);

    volatile dma_desc dma_put_adr_epj_org, dma_put_adr_epj_buf, dma_put_adr_spj_buf;
    volatile int reply_put_adr_epj_org = 0;
    volatile int reply_put_adr_epj_buf = 0;
    volatile int reply_put_adr_spj_buf = 0;
    dma_descriptor_init(&dma_put_adr_epj_org, &reply_put_adr_epj_org);
    dma_descriptor_init(&dma_put_adr_epj_buf, &reply_put_adr_epj_buf);
    dma_descriptor_init(&dma_put_adr_spj_buf, &reply_put_adr_spj_buf);
    
    volatile dma_desc dma_get_epj_org, dma_get_spj_org;
    volatile int reply_get_epj_org = 0;
    volatile int reply_get_spj_org = 0;
    dma_descriptor_init(&dma_get_epj_org, &reply_get_epj_org);
    dma_descriptor_init(&dma_get_spj_org, &reply_get_spj_org);

    volatile dma_desc dma_put_epj_sorted, dma_put_spj_sorted;
    volatile int reply_put_epj_sorted = 0;
    volatile int reply_put_spj_sorted = 0;
    dma_descriptor_init(&dma_put_epj_sorted, &reply_put_epj_sorted);
    dma_descriptor_init(&dma_put_spj_sorted, &reply_put_spj_sorted);
    
    int my_id = athread_get_id(-1);
    int n_glb_tot = (int)(((unsigned long*)args)[0]);
    int n_loc_tot = (int)(((unsigned long*)args)[1]);
    void * adr_tp = (void *)((unsigned long*)args)[2];
    void * adr_epj_sorted = (void *)((unsigned long*)args)[3];
    void * adr_epj_org = (void *)((unsigned long*)args)[4];
    void * adr_spj_sorted = (void *)((unsigned long*)args)[5];
    void * adr_spj_org = (void *)((unsigned long*)args)[6];
    void * adr_adr_epj_org2glb = (void *)((unsigned long*)args)[7];
    void * adr_adr_epj_buf2glb = (void *)((unsigned long*)args)[8];
    void * adr_adr_spj_buf2glb = (void *)((unsigned long*)args)[9];
    int my_n_glb = n_glb_tot/NUMBER_OF_CPE + ( (my_id < n_glb_tot % NUMBER_OF_CPE) ? 1 : 0 );
    int my_offset_glb = (n_glb_tot/NUMBER_OF_CPE)*my_id + ( (my_id < n_glb_tot % NUMBER_OF_CPE) ? my_id : n_glb_tot % NUMBER_OF_CPE );
    enum {
        CHUNK_SIZE = 64,
    };
    size_t bsize_tp  = sizeof(tpLM);
    size_t bsize_adr = sizeof(S32_);
    size_t bsize_ep  = sizeof(epjLM);
    size_t bsize_sp  = sizeof(spjLM);
    volatile tpLM tp[CHUNK_SIZE];
    volatile tpLM tp_new[CHUNK_SIZE];
    volatile int my_n_epj = 0;
    volatile int my_n_spj = 0;
    // just count number to determine offset to write
    for (i=0; i<my_n_glb; i+=CHUNK_SIZE) {
        int nrem = my_n_glb - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset_glb+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        for (j=0; j<nn; j++) {
            if(GetMsb(tp[j].adr_ptcl_)==0){
                my_n_epj++;
            }
            else{
                my_n_spj++;
            }
        }
    }
    sync_array_();

    // get offset to write by scan
    //int n_epj_ar[NUMBER_OF_CPE];
    //int n_spj_ar[NUMBER_OF_CPE];
    volatile int * n_epj_ar = (int *) ldm_malloc(sizeof(int)*NUMBER_OF_CPE);
    volatile int * n_spj_ar = (int *) ldm_malloc(sizeof(int)*NUMBER_OF_CPE);
    
    for (i=0; i<NUMBER_OF_CPE; i++) {
        n_epj_ar[i] = n_spj_ar[i] = 0;
    }
    n_epj_ar[my_id] = my_n_epj;
    n_spj_ar[my_id] = my_n_spj;
    sync_array_();
    for (i=0; i<NUMBER_OF_CPE; i++) {
        cpe_bcast_int32(i, &n_epj_ar[i]);
        cpe_bcast_int32(i, &n_spj_ar[i]);
    }
    sync_array_();
    volatile int my_n_disp_epj = 0;
    volatile int my_n_disp_spj = 0;
    for(i=0; i<my_id; i++){
        my_n_disp_epj += n_epj_ar[i];
        my_n_disp_spj += n_spj_ar[i];
    }
    ldm_free(n_epj_ar, sizeof(int)*NUMBER_OF_CPE);
    ldm_free(n_spj_ar, sizeof(int)*NUMBER_OF_CPE);
    sync_array_();

    /*
    EpjLM epj_org[CHUNK_SIZE];
    SpjLM spj_org[CHUNK_SIZE];
    S32_ my_n_epj_ar[CHUNK_SIZE];
    S32_ my_n_spj_ar[CHUNK_SIZE];
    S32_ adr_epj_ar[CHUNK_SIZE];
    S32_ adr_spj_ar[CHUNK_SIZE];
    */
    volatile epjLM * epj_org = (epjLM *) ldm_malloc(bsize_ep * CHUNK_SIZE);
    volatile spjLM * spj_org = (spjLM *) ldm_malloc(bsize_sp * CHUNK_SIZE);
    volatile S32_ * my_n_epj_ar = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    volatile S32_ * my_n_spj_ar = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    volatile S32_ * adr_epj_ar  = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    volatile S32_ * adr_spj_ar  = (S32_*) ldm_malloc(sizeof(S32_)*CHUNK_SIZE);
    
    sync_array_();
    my_n_epj = my_n_spj = 0;
    for (i=0; i<my_n_glb; i+=CHUNK_SIZE) {
        int nrem = my_n_glb - i;
        int nn   = nrem < CHUNK_SIZE ? nrem : CHUNK_SIZE;
        // get tp
        reply_get = 0;
        dma_set_op(&dma_get, DMA_GET);
        dma_set_mode(&dma_get, PE_MODE);
        dma_set_reply(&dma_get, &reply_get);
        dma_set_size(&dma_get, bsize_tp * nn);
        dma(dma_get, (long*)((tpLM *)adr_tp+my_offset_glb+i), (long*)(tp));
        dma_wait(&reply_get, 1);
        while (reply_get != 1) {}
        
        ////////////
        // to put adr
        // preparing for puting adr_epj_org
        reply_put_adr_epj_org = 0;
        dma_set_op(&dma_put_adr_epj_org, DMA_PUT);
        dma_set_mode(&dma_put_adr_epj_org, PE_MODE);
        dma_set_reply(&dma_put_adr_epj_org, &reply_put_adr_epj_org);
        dma_set_size(&dma_put_adr_epj_org, bsize_adr);

        // preparing for puting adr_epj_buf
        reply_put_adr_epj_buf = 0;
        dma_set_op(&dma_put_adr_epj_buf, DMA_PUT);
        dma_set_mode(&dma_put_adr_epj_buf, PE_MODE);
        dma_set_reply(&dma_put_adr_epj_buf, &reply_put_adr_epj_buf);
        dma_set_size(&dma_put_adr_epj_buf, bsize_adr);
        
        // preparing for puting adr_spj_buf
        reply_put_adr_spj_buf = 0;
        dma_set_op(&dma_put_adr_spj_buf, DMA_PUT);
        dma_set_mode(&dma_put_adr_spj_buf, PE_MODE);
        dma_set_reply(&dma_put_adr_spj_buf, &reply_put_adr_spj_buf);
        dma_set_size(&dma_put_adr_spj_buf, bsize_adr);
        // to put adr
        ////////////
        
        ////////////
        // to get epj_org and spj_org to put epj_sorted and spj_sorted
        // prepare for gettin epj_org
        reply_get_epj_org = 0;
        dma_set_op(&dma_get_epj_org, DMA_GET);
        dma_set_mode(&dma_get_epj_org, PE_MODE);
        dma_set_reply(&dma_get_epj_org, &reply_get_epj_org);
        dma_set_size(&dma_get_epj_org, bsize_ep);

        // prepare for gettin spj_org
        reply_get_spj_org = 0;
        dma_set_op(&dma_get_spj_org, DMA_GET);
        dma_set_mode(&dma_get_spj_org, PE_MODE);
        dma_set_reply(&dma_get_spj_org, &reply_get_spj_org);
        dma_set_size(&dma_get_spj_org, bsize_sp);
        // to get epj_org and spj_org to put epj_sorted and spj_sorted
        ////////////
        
        volatile int my_n_epj_tmp = 0;
        volatile int my_n_epj_org_tmp = 0;
        volatile int my_n_epj_buf_tmp = 0;
        volatile int my_n_spj_tmp = 0;
        for (j=0; j<nn; j++) {
            U32_ adr = tp[j].adr_ptcl_;
            CopyKeyOnly(&(tp[j]), &(tp_new[j]));
            if(GetMsb(adr)==0){
                tp_new[j].adr_ptcl_ = my_n_disp_epj + my_n_epj;
                adr_epj_ar[my_n_epj_tmp] = adr;
                my_n_epj_ar[my_n_epj_tmp] = my_n_disp_epj + my_n_epj;
                if(adr < n_loc_tot){ // particles are in own process
                    dma(dma_put_adr_epj_org,
                        (long*)((S32_*)adr_adr_epj_org2glb+adr_epj_ar[my_n_epj_tmp]),
                        (long*)((S32_*)my_n_epj_ar+my_n_epj_tmp));
                    my_n_epj_org_tmp++;
                }
                else{
                    dma(dma_put_adr_epj_buf,
                        (long*)((S32_*)adr_adr_epj_buf2glb+adr_epj_ar[my_n_epj_tmp]-n_loc_tot),
                        (long*)((S32_*)my_n_epj_ar+my_n_epj_tmp));
                    my_n_epj_buf_tmp++;
                }
                dma(dma_get_epj_org,
                    (long*)((epjLM*)adr_epj_org+adr_epj_ar[my_n_epj_tmp]),
                    (long*)((epjLM*)epj_org+my_n_epj_tmp));
                my_n_epj++;
                my_n_epj_tmp++;
            }
            else{
                tp_new[j].adr_ptcl_ = SetMsb(my_n_disp_spj + my_n_spj); // set tp
                adr_spj_ar[my_n_spj_tmp]  = ClearMsb(adr); // temporally store
                my_n_spj_ar[my_n_spj_tmp] = my_n_disp_spj + my_n_spj; // temporally store
                dma(dma_put_adr_spj_buf,
                    (long*)((S32_*)adr_adr_spj_buf2glb+adr_spj_ar[my_n_spj_tmp]),
                    (long*)((S32_*)my_n_spj_ar+my_n_spj_tmp)); //put adr_spj_buf2glb
                dma(dma_get_spj_org,
                    (long*)((spjLM*)adr_spj_org+adr_spj_ar[my_n_spj_tmp]),
                    (long*)((spjLM*)spj_org+my_n_spj_tmp));
                my_n_spj++;
                my_n_spj_tmp++;
            }
        }
        ///////////
        // put tp
        reply_put = 0;
        dma_set_op(&dma_put, DMA_PUT);
        dma_set_mode(&dma_put, PE_MODE);
        dma_set_reply(&dma_put, &reply_put);
        dma_set_size(&dma_put, bsize_tp * nn);
        dma(dma_put, (long*)((tpLM *)adr_tp+my_offset_glb+i), (long*)(tp_new));
        dma_wait(&reply_put, 1);
        while (reply_put != 1) {}
        // put tp
        ///////////
        
        ///////////
        // put tp, adr_epj_org2glb, adr_epj_buf2glb, adr_spj_buf2glb
        // wait puting adr_epj_org2glb
        dma_wait(&reply_put_adr_epj_org, my_n_epj_org_tmp);
        while (reply_put_adr_epj_org != my_n_epj_org_tmp) {}
        
        // wait puting adr_epj_buf2glb
        dma_wait(&reply_put_adr_epj_buf, my_n_epj_buf_tmp);
        while (reply_put_adr_epj_buf !=  my_n_epj_buf_tmp){}
        
        // wait puting adr_spj_buf2glb
        dma_wait(&reply_put_adr_spj_buf, my_n_spj_tmp);
        while (reply_put_adr_spj_buf != my_n_spj_tmp) {}
        // put tp, adr_epj_org2glb, adr_epj_buf2glb, adr_spj_buf2glb
        ///////////

        /////////
        // get epj_org, spj_org to put epj_sorted and spj_sorted
        // wait getting epj_org
        dma_wait(&reply_get_epj_org, my_n_epj_tmp);
        while (reply_get_epj_org != my_n_epj_tmp) {}
        
        // wait getting spj_org
        dma_wait(&reply_get_spj_org, my_n_spj_tmp);
        while (reply_get_spj_org != my_n_spj_tmp) {}
        // get epj_org, spj_org
        /////////

        ///////////
        // put epj_sorted, spj_sorted
        // put epj_sorted
        if( my_n_epj_tmp > 0) {
            reply_put_epj_sorted = 0;
            dma_set_op(&dma_put_epj_sorted, DMA_PUT);
            dma_set_mode(&dma_put_epj_sorted, PE_MODE);
            dma_set_reply(&dma_put_epj_sorted, &reply_put_epj_sorted);
            dma_set_size(&dma_put_epj_sorted, bsize_ep*my_n_epj_tmp);
            dma(dma_put_epj_sorted,
                (long*)((epjLM*)adr_epj_sorted + (my_n_disp_epj + my_n_epj - my_n_epj_tmp)),
                (long*)(epj_org));
            dma_wait(&reply_put_epj_sorted, 1);
            while (reply_put_epj_sorted != 1) {}
        }

        // put spj_sorted
        if( my_n_spj_tmp > 0) {
            reply_put_spj_sorted = 0;
            dma_set_op(&dma_put_spj_sorted, DMA_PUT);
            dma_set_mode(&dma_put_spj_sorted, PE_MODE);
            dma_set_reply(&dma_put_spj_sorted, &reply_put_spj_sorted);
            dma_set_size(&dma_put_spj_sorted, bsize_sp*my_n_spj_tmp);
            dma(dma_put_spj_sorted,
                (long*)((spjLM*)adr_spj_sorted + (my_n_disp_spj + my_n_spj - my_n_spj_tmp)),
                (long*)((spjLM*)spj_org));
            dma_wait(&reply_put_spj_sorted, 1);
            while (reply_put_spj_sorted != 1) {}
        }
        // put epj_sorted, spj_sorted
        ///////////

    }
    
    ldm_free(epj_org, bsize_ep*CHUNK_SIZE);
    ldm_free(spj_org, bsize_sp*CHUNK_SIZE);
    ldm_free(my_n_epj_ar, sizeof(S32_)*CHUNK_SIZE);
    ldm_free(my_n_spj_ar, sizeof(S32_)*CHUNK_SIZE);
    ldm_free(adr_epj_ar, sizeof(S32_)*CHUNK_SIZE);
    ldm_free(adr_spj_ar, sizeof(S32_)*CHUNK_SIZE);
#else
#error The value of `VERSION_OF_SET_PTCL_SORTED_GLB_CPE` is invalid.
#endif
}

