// Flag conflict check.

// General part.
#ifdef USE_GRAD_N //{
#ifndef USE_GRAD_H //{
#define USE_GRAD_H
#endif //}
#endif //}

#ifdef USE_DISPH //{
#ifdef USE_SPSPH //{
#error Both USE_DISPH and USE_SPSPH are used.
#endif // USE_SPSPH //}
#endif // USE_DISPH //}

#ifdef USE_SPSPH //{
#ifdef USE_DISPH //{
#error Both USE_DISPH and USE_SPSPH are used.
#endif // USE_DISPH //}
#endif // USE_SPSPH //}

#ifdef SET_SNII_TEMPERATURE_GAS //{
#undef SET_SNII_TEMPERATURE
#endif //}

#ifdef USE_SYMMETRIZED_SOFTENING //{
#define SYMMETRIC
#endif // USE_SYMMETRIZED_SOFTENING //}

// Specific part for individual run.

#ifdef TASK_HYDROSTATIC //{
#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (2)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (2)

#ifndef USE_PARTICLE_TAG
#define USE_PARTICLE_TAG
#endif //USE_PARTICLE_TAG

#ifndef PERIODIC_RUN
#define PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef FARULTRAVIOLET_HEATING
#undef FARULTRAVIOLET_HEATING 
#endif //FARULTRAVIOLET_HEATING
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif // TASK_HYDROSTATIC //}


#ifdef TASK_COLD_COLLAPSE //{
#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //TASK_COLD_COLLAPSE //}

#ifdef TASK_TEST_STROMGRENSPHERE
#ifndef USE_HIIREGION_MODEL
#define USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#endif //TASK_TEST_STROMGRENSPHERE

#ifdef TASK_1D_SHOCKE_TUBE //{
#undef DIMENSION 
#define DIMENSION           (1)
#undef BisectionDimension 
#define BisectionDimension  (1)

#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE
#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT //{
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT //}
#endif //TASK_1D_SHOCKE_TUBE //}


#ifdef TASK_BLAST_WAVE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)


#ifndef PERIODIC_RUN
#define PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif // TASK_BLAST_WAVE //}


#ifdef TASK_SINUSOIDAL_WAVE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (1)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (1)


#ifndef PERIODIC_RUN
#define PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_FARULTRAVIOLET_HEATING
#undef USE_FARULTRAVIOLET_HEATING 
#endif //USE_FARULTRAVIOLET_HEATING
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE
#endif // TASK_SINUSOIDAL_WAVE //}


#ifdef TASK_KELVINHELMHOLTZ_INSTABILITY //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (2)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (2)


#ifndef PERIODIC_RUN
#define PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifndef USE_PARTICLE_TAG
#define USE_PARTICLE_TAG
#endif //USE_PARTICLE_TAG

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif // TASK_KELVINHELMHOLTZ_INSTABILITY //}



#ifdef TASK_3D_COLLAPSE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT


#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#ifdef USE_CELIB
#undef USE_CELIB
#endif // USE_CELIB
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK
#endif // DELAYED_FEEDBACK
#ifdef USE_METAL_DIFFUSION
#undef USE_METAL_DIFFUSION
#endif // USE_METAL_DIFFUSION

#endif //TASK_3D_COLLAPSE //}


#ifdef TASK_SANTABARBARA //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifndef COSMOLOGICAL_RUN 
#define COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif //TASK_SANTABARBARA //}


#ifdef TASK_NFW //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN

#endif //TASK_NFW //}


#ifdef TASK_MERGER //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN

#ifndef STARFORMATION
#define STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif //TASK_MERGER //}


#ifdef TASK_DICE_RUN //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN

#ifndef STARFORMATION
#define STARFORMATION
#endif //STARFORMATION
#ifdef USE_FORCIBLE_STARFORMATION 
#undef USE_FORCIBLE_STARFORMATION 
#endif //USE_FORCIBLE_STARFORMATION 
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif //TASK_DICE_RUN //}


#ifdef TASK_MW //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN

#ifndef STARFORMATION
#define STARFORMATION
#endif //STARFORMATION
#ifdef USE_FORCIBLE_STARFORMATION 
#undef USE_FORCIBLE_STARFORMATION 
#endif //USE_FORCIBLE_STARFORMATION 
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#ifndef USE_BOUNDARY_CONDITION //{
#define USE_BOUNDARY_CONDITION
#undef BOUNDARY_CONDITION_SPHERICAL_SHELL_EDGE
#define BOUNDARY_CONDITION_SPHERICAL_SHELL_EDGE (0.1*MPC_CGS)
#endif // USE_BOUNDARY_CONDITION //}

#endif //TASK_MW //}


#ifdef TASK_GALACTIC_CENTER //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN

#ifndef COOLING_RUN
#define COOLING_RUN 
#endif //COOLING_RUN
#ifndef USE_SPAANS2008_COOLING_FUNCTIONS //{
#define USE_SPAANS2008_COOLING_FUNCTIONS 
#endif // USE_SPAANS2008_COOLING_FUNCTIONS //}
#ifdef USE_CLOUDY_COOLING_FUNCTIONS //{
#undef USE_CLOUDY_COOLING_FUNCTIONS 
#endif // USE_CLOUDY_COOLING_FUNCTIONS //}
#ifdef USE_SPAANS1997_COOLING_FUNCTIONS //{
#undef USE_SPAANS1997_COOLING_FUNCTIONS 
#endif // USE_SPAANS1997_COOLING_FUNCTIONS //}
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK
#endif //DELAYED_FEEDBACK
#ifdef USE_CELIB //{
#undef USE_CELIB
#endif //USE_CELIB //}
#ifndef USE_SINK_PARTICLE
#define USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif //TASK_GALACTIC_CENTER //}

#ifdef TASK_TURBULENCE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
// #ifdef COOLING_RUN
// #undef COOLING_RUN 
// #endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
// #ifdef USE_SINK_PARTICLE
// #undef USE_SINK_PARTICLE
// #endif //USE_SINK_PARTICLE

#endif //TASK_TURBULENCE //}

#ifdef TASK_M2_COLLAPSE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
// #ifdef USE_SINK_PARTICLE
// #undef USE_SINK_PARTICLE
// #endif //USE_SINK_PARTICLE

#endif //TASK_M2_COLLAPSE //}

#ifdef TASK_TEST_EQUILIBRIUM_TEMPERATURE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

// #ifndef COSMOLOGICAL_RUN 
// #define COSMOLOGICAL_RUN 
// #endif //COSMOLOGICAL_RUN
// #ifndef GRAVITY_RUN 
// #define GRAVITY_RUN 
// #endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifndef COOLING_RUN
#define COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif // TASK_TEST_EQUILIBRIUM_TEMPERATURE //}

#ifdef TASK_TEST_EQUILIBRIUM_TEMPERATURE //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

// #ifndef COSMOLOGICAL_RUN 
// #define COSMOLOGICAL_RUN 
// #endif //COSMOLOGICAL_RUN
// #ifndef GRAVITY_RUN 
// #define GRAVITY_RUN 
// #endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifndef COOLING_RUN
#define COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif // TASK_TEST_EQUILIBRIUM_TEMPERATURE //}

#if (defined(TASK_TEST_1D_THERMAL_CONDUCTIVITY)||defined(TASK_TEST_3D_THERMAL_CONDUCTIVITY)) //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#ifdef TASK_TEST_1D_THERMAL_CONDUCTIVITY //{
#define DIMENSION           (1)
#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY //}
#ifdef TASK_TEST_3D_THERMAL_CONDUCTIVITY //{
#define DIMENSION           (3)
#endif // TASK_TEST_3D_THERMAL_CONDUCTIVITY //}
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef USE_CELIB //{
#undef USE_CELIB
#endif // USE_CELIB //}
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#ifndef PERIODIC_RUN
#define PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifdef THERMAL_CONDUCTIVITY_COEF_TYPE //{
#undef THERMAL_CONDUCTIVITY_COEF_TYPE
#endif // THERMAL_CONDUCTIVITY_COEF_TYPE //}
#define THERMAL_CONDUCTIVITY_COEF_TYPE (-1)

// #undef MAXIMUM_TIME_HIERARCHY 
// #define MAXIMUM_TIME_HIERARCHY (0 )

// #ifdef TFactorCourant //{
// #undef TFactorCourant
// #endif // TFactorCourant //}
// #define TFactorCourant  (0.01)

#endif // TASK_TEST_1D_THERMAL_CONDUCTIVITY//TASK_TEST_3D_THERMAL_CONDUCTIVITY //}

#ifdef TASK_MAKE_SMOOTHDISK //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifndef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#define PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN

#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK
#endif //DELAYED_FEEDBACK
#ifdef USE_CELIB //{
#undef USE_CELIB
#endif //USE_CELIB //}
#ifndef USE_SINK_PARTICLE
#define USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif // TASK_MAKE_SMOOTHDISK //}


#ifdef TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef USE_VARIABLE_TIMESTEP 
#undef USE_VARIABLE_TIMESTEP
#endif //USE_VARIABLE_TIMESTEP


#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif //TASK_TEST_SYMMETRIZED_POTENTIAL_ERROR //}


#ifdef TASK_1D_TWOFLUIDS //{
#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (1)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (1)

#ifndef USE_PARTICLE_TAG
#define USE_PARTICLE_TAG
#endif //USE_PARTICLE_TAG

#ifndef PERIODIC_RUN
#define PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef GRAVITY_RUN 
#undef GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef FARULTRAVIOLET_HEATING
#undef FARULTRAVIOLET_HEATING 
#endif //FARULTRAVIOLET_HEATING
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#ifdef USE_CELIB
#undef USE_CELIB
#endif // USE_CELIB
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK
#endif // DELAYED_FEEDBACK
#ifdef USE_METAL_DIFFUSION
#undef USE_METAL_DIFFUSION
#endif // USE_METAL_DIFFUSION

#endif // TASK_1D_TWOFLUIDS //}

#ifdef TASK_KEPLER //{
#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (2)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (2)

#ifndef USE_PARTICLE_TAG
#define USE_PARTICLE_TAG
#endif //USE_PARTICLE_TAG


#undef ASCIIDATA_DUMP_INTERVAL
#define ASCIIDATA_DUMP_INTERVAL (1)
//#endif //ASCIIDATA_DUMP_INTERVAL

#ifdef PERIODIC_RUN
#undef PERIODIC_RUN 
#endif // PERIODIC_RUN

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifdef COSMOLOGICAL_RUN 
#undef COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_FAST_SCHEME
#undef USE_FAST_SCHEME
#endif //USE_FAST_SCHEME
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK 
#endif //DELAYED_FEEDBACK
#ifdef COOLING_RUN
#undef COOLING_RUN 
#endif //COOLING_RUN
#ifdef FARULTRAVIOLET_HEATING
#undef FARULTRAVIOLET_HEATING 
#endif //FARULTRAVIOLET_HEATING
#ifdef USE_HIIREGION_MODEL
#undef USE_HIIREGION_MODEL
#endif //USE_HIIREGION_MODEL
#ifdef STARFORMATION
#undef STARFORMATION
#endif //STARFORMATION
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#ifdef USE_CELIB
#undef USE_CELIB
#endif // USE_CELIB
#ifdef DELAYED_FEEDBACK
#undef DELAYED_FEEDBACK
#endif // DELAYED_FEEDBACK
#ifdef USE_METAL_DIFFUSION
#undef USE_METAL_DIFFUSION
#endif // USE_METAL_DIFFUSION

#undef TFactor
#define TFactor    (0.1)
#undef TFactorCourant  
#define TFactorCourant  (0.1)
#undef TFactorU
#define TFactorU        (0.1)
#undef TFactorVel      
#define TFactorVel      (0.1)
#undef TFactorAcc      
#define TFactorAcc      (0.1)

#endif // TASK_KEPLER //}



#ifdef TASK_GALAXY_FORMATION //{

#ifdef DIMENSION 
#undef DIMENSION 
#endif //DIMENSION
#define DIMENSION           (3)
#ifdef BisectionDimension 
#undef BisectionDimension 
#endif //BisectionDimension
#define BisectionDimension  (3)

#ifdef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#undef PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT
#endif //PRINT_LOG_TIMESTEP_IN_PHYSICAL_UNIT

#ifndef COSMOLOGICAL_RUN 
#define COSMOLOGICAL_RUN 
#endif //COSMOLOGICAL_RUN
#ifndef GRAVITY_RUN 
#define GRAVITY_RUN 
#endif //GRAVITY_RUN
#ifdef USE_SINK_PARTICLE
#undef USE_SINK_PARTICLE
#endif //USE_SINK_PARTICLE

#endif //TASK_GALAXYFORMATION //}
