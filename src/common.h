/*
 * File: common.h
 * 
 * File Contents: Preprocessor directives for adjusting how the model executes.
 *
 * Created on: Jan 28, 2015
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef COMMON_H_
#define COMMON_H_

#include "enums.h"

#define quote(x) #x

#include <float.h>
#include <vector>

#ifdef DBL_DECIMAL_DIG
  #define OP_DBL_Digs (DBL_DECIMAL_DIG)
#else
  #ifdef DECIMAL_DIG
    #define OP_DBL_Digs (DECIMAL_DIG)
  #else
    #define OP_DBL_Digs (DBL_DIG + 3)
  #endif
#endif


typedef std::vector<float> REGULATORS_T;

/*****************************************************************************
 * SENSITIVITY ANALYSIS AND  CALIBRATION                                     *
 *****************************************************************************/
// Modify normal values of parameters for S.A. or calibration purposes
//#define CALIBRATION
// Used to print new values of parameters during S.A. or calibration
//#define PRINT_PARAMETER_VALUES

/*****************************************************************************
 * PROFILING FLAGS - FOR TIMING INFORMATION                                  *
 *****************************************************************************/
// Profile major sections of go():
#define PROFILE_MAJOR_STEPS  // Turned this off for SA
// Profile each cell type:
#define PROFILE_CELL_FUNC  // Turned this off for SA
// Profile each section of chemical diffusion:
//#define PROFILE_CHEM_DIFF_DETAILED
// Profile each section of cell seeding:
//#define PROFILE_CELL_PROLIF_DETAILED
// Profile time spent on chemical diffusion by each thread:
//#define PROFILE_THREAD_LEVEL_CHEM_DIFF
// Profile each section of ECM updates:
//#define PROFILE_ECM_UPDATE  // Turned this off for SA
// Profile each section of patch updates
//#define PROFILE_PATCH_UPDATE  // Turned this off for SA
// Profile ECM Function
//#define PROFILE_ECM
// Profile each major steps of diffusion on GPU
#define PROFILE_GPU

/*****************************************************************************
 * PARALLEL OPTION FLAGS                                                     *
 *****************************************************************************/
// Parallel by cell types	(Obsolete)
//#define PAR_BY_CELLTYPE

/*****************************************************************************
 * OPTIMIZATION FLAGS                                                        *
 *****************************************************************************/
// TODO(Nuttiiya)
//#define VECTORIZE
// TODO(Nuttiiya)
//#define OPT_CELL_SEEDING
// TODO(Nuttiiya)
//#define ECM_UNROLL_LOOP
// OBSOLETE
/* With OPT_ECM defined, we assume that we access ECM::HAlife in this order:
 * 1. Decrement life by calling decrement(int n) on HAlife in ECMFunction()
 * 2. Determine number of HA by calling HAlife.size() in fragmentHA()
 * 3. Remove dead HAs in HAlife in updateECM()
 * NOTE: Steps 1 & 2 can occur in either order but must preceed step 3. */
#define OPT_ECM
// Minimize memory required
#define MIN_MEM
// Optimize ECM request with request pool
#define OPT_FECM_REQUEST
// Optimize cell seeding. Should only be set if world is sparsed.
//#define SPARSE_WORLD_OPT
// Optimize pinned memory footprint by alllocate only 1 buffer for input and 1 for output
//#define OPT_PINNED_MEM
// PACK 0, Kernel Spectrum stored in pinned memory
#define PAGELOCK_DIFFUSE

// Optimize data struture for chemical concentrations
//#define OPT_CHEM
// Chemical structure packing options
#ifdef OPT_CHEM
//#define CHEM_PACK8
//#define CHEM_PACK4
#endif



/*****************************************************************************
 * OPTION MACROS                                                             *
 *****************************************************************************/

// Store kernel spectrum on GPU if using Tesla M40
#define M40
// Raw volume output time stride
#define RAW_TIME_STRIDE 10
// ECM initial condition factor
#define RAW_ECM_FACTOR  0.0002
// ECM gradient factor
#define RAW_dECM_FACTOR 0.05  //2.0

/*****************************************************************************
 * OPTION FLAGS                                                              *
 *****************************************************************************/
// Disable randomness in model
//#define DISABLE_RAND
// Enable multiple seeds output
//#define MULT_SEEDS
// Output cell volume as raw file
//#define WRITE_RAW_CELL
// Output ECM volume map as raw file
//#define WRITE_RAW_ECM
// Enable Biomarker Output
//#define BIOMARKER_OUTPUT
// Enable Paraview rendering
//#define PARAVIEW_RENDERING
// Enable In-Situ Visualization
#define VISUALIZATION
// Enable overlap visualization and computation
//#define OVERLAP_VIS
// Enable timing functions for GL rendering
#define TIME_GL
// Run in 3D
#define MODEL_3D
// Initialize Vocal Fold Morphology
#define MODEL_VOCALFOLD
// Human size
#define HUMAN_VF
// Rat size
//#define RAT_VF
// Use PEVOC scale parameters 	(Obsolete?)
//#define PEVOC_SCALE
// Use PDE based chemical diffusion
#define PDE_DIFFUSE
// Use PDE-based Convolution-based GPU diffusion
#define GPU_DIFFUSE
// TODO(Nuttiiya)
//#define COLLECT_CELL_INS_DEL_STATS
// Don't sprout cells in damage zone (this may results in smaller number of cells sprouted)
#define NO_CELLS_IN_DAMZONE

/*****************************************************************************
 * DEBUGGING FLAGS                                                           *
 *****************************************************************************/
// Test diffusion kernel computation
//#define TEST_KERNEL
// Computing diffusion kernel window coverage
//#define COMPUTE_COVERAGE
// Only run diffusions on GPU0
//#define GPU0_ONLY
// Only run diffusions on GPU1
//#define GPU1_ONLY
// Assign first GPU to use
#define GPU1 0
// Assign second GPU to use
#define GPU2 1


/*****************************************************************************
 * VOCAL FOLDS DIMENSION                                                     *
 *****************************************************************************/
#ifdef HUMAN_VF
#define LPx 1.65	//1.6	// vocal fold width (mm)	// This currently gets switched
#define LPy 20.85	//24.9 	// vocal fold length (mm)	//	x and y
#define LPz 15.09	//17.4	// vocal fold height (mm)

#define epitheliumthickness 0.05 /* Epithelium thickness (in mm)
                                 * (Hirano, Minoru. "Phonosurgery: basic and clinical
                                 * investigations." Otologia (Fukuoka) 21.suppl 1 (1975): 239-260.*/
#define capillaryradius 0.075    //0.045
                                /* Capillary radius (in mm)
                                 * (Sato, Kiminori, Minoru Hirano, and Tadashi Nakashima. "Electron microscopic
                                 * and immunohistochemical investigation of Reinke's edema." Annals of Otology,
                                 * Rhinology & Laryngology 108.11 (1999): 1068-1072.) */
#define capillaryXdistance 0.405    //0.4	// mm
#define capillaryYdistance 0.405    //4	// mm

#elif defined(RAT_VF) // HUMAN_VF

#define LPx  0.2  // vocal fold width  (mm)  // Rat's vocal fold dimension (From Aman)
#define LPy  1.0  // vocal fold length (mm)  // Rat's vocal fold dimension (From Aman)
#define LPz  1.0  // vocal fold height (mm)  // Rat's vocal fold dimension (From Aman)

#define epitheliumthickness 0.010 //0.00625
				/* Epithelium thickness (in mm)
                                 * (Hirano, Minoru. "Phonosurgery: basic and clinical
                                 * investigations." Otologia (Fukuoka) 21.suppl 1 (1975): 239-260.*/
#define capillaryradius 0.0035 //0.075    //0.045
                                /* Capillary radius (in mm)
                                 * (Sato, Kiminori, Minoru Hirano, and Tadashi Nakashima. "Electron microscopic
                                 * and immunohistochemical investigation of Reinke's edema." Annals of Otology,
                                 * Rhinology & Laryngology 108.11 (1999): 1068-1072.) */
#define capillaryXdistance 0.015 // making this 2p  //0.01289 //0.405    //0.4   // mm
#define capillaryYdistance 0.015 // making this 2p  //0.01289 //0.405    //4 // mm


#else   // HUMAN_VF

#define LPx  0.825 // mini vocal fold width  (mm) (for visualization testing purposes)
#define LPy 10.500 // mini vocal fold length (mm) (for visualization testing purposes)
#define LPz  5.400 // mini vocal fold height (mm) (for visualization testing purposes)

#define epitheliumthickness 0.015 //0.00625
				/* Epithelium thickness (in mm)
                                 * (Hirano, Minoru. "Phonosurgery: basic and clinical
                                 * investigations." Otologia (Fukuoka) 21.suppl 1 (1975): 239-260.*/
#define capillaryradius 0.045 //0.075    //0.045
                                /* Capillary radius (in mm)
                                 * (Sato, Kiminori, Minoru Hirano, and Tadashi Nakashima. "Electron microscopic
                                 * and immunohistochemical investigation of Reinke's edema." Annals of Otology,
                                 * Rhinology & Laryngology 108.11 (1999): 1068-1072.) */
#define capillaryXdistance 0.2 //0.405    //0.4   // mm
#define capillaryYdistance 0.2 //0.405    //4 // mm


#endif  //HUMAN_VF

#define fractionSLP 0.13 /* Thickness of superficial lamina propria (according to collagen organization)
                          * as fraction of total LP (unitless)
                          * [Ref] Prades, Jean-Michel, et al. "Lamina propria of the human vocal fold: histomorphometric
                          * study of collagen fibers." Surgical and radiologic anatomy 32.4 (2010): 377-382.
                          * [Ref] Kaiser, M. L., et al. "Laryngeal epithelial thickness: a comparison between optical coherence
                          * tomography and histology." Clinical Otolaryngology 34.5 (2009): 460-466.
                          * [Ref] Hirano, Minoru. "Phonosurgery: basic and clinical investigations." Otologia (Fukuoka)
                          * 21.suppl 1 (1975): 239-260. */
#define fractionILP 0.51 /* Thickness of intermediate lamina propria (according to collagen organization)
                          * as fraction of total LP (unitless)*/
#define fractionDLP 0.36 /* Thickness of deep lamina propria (according to collagen organization)
                          * as fraction of total LP (unitless) */
#define fractionProteinLB 0.75		// Hahn, Mariah S., et al. "Midmembranous vocal fold
#define fractionProteinUB 0.85	   // lamina propria proteoglycans across selected species."
																   // Annals of Otology, Rhinology & Laryngology 114.6 (2005): 451-462.


/* Non-uniform cellularity in vocal fold lamina propria with depth 
* VF divided into 5 sections, each 20% of LP thickness (1 superficial, 5 deep)
* [Ref] Catten, Michael, et al. "Analysis of cellular location and concentration in vocal fold
* lamina propria." Otolaryngology--Head and Neck Surgery 118.5 (1998): 663-667. */
#define fibroblastOne 0.212    // Fraction of non-epithelial LP fibroblast in section 1 
#define fibroblastTwo 0.192    // Fraction of non-epithelial LP fibroblast in section 2
#define fibroblastThree 0.169    // Fraction of non-epithelial LP fibroblast in section 3 
#define fibroblastFour 0.177    // Fraction of non-epithelial LP fibroblast in section 4 
#define fibroblastFive 0.250    // Fraction of non-epithelial LP fibroblast in section 5 
#define afibroblastOne 0.347    // Fraction of non-epithelial LP activated fibroblast (myofibroblast) in section 1 
#define afibroblastTwo 0.227    // Fraction of non-epithelial LP activated fibroblast (myofibroblast) in section 2
#define afibroblastThree 0.232    // Fraction of non-epithelial LP activated fibroblast (myofibroblast) in section 3
#define afibroblastFour 0.132    // Fraction of non-epithelial LP activated fibroblast (myofibroblast) in section 4
#define afibroblastFive 0.062    // Fraction of non-epithelial LP activated fibroblast (myofibroblast) in section 5
#define macrophageOne 0.612    // Fraction of non-epithelial LP macrophage in section 1
#define macrophageTwo 0.182    // Fraction of non-epithelial LP macrophage in section 2
#define macrophageThree 0.097    // Fraction of non-epithelial LP macrophage in section 3
#define macrophageFour 0.024    // Fraction of non-epithelial LP macrophage in section 4
#define macrophageFive 0.085    // Fraction of non-epithelial LP macrophage in section 5


/*****************************************************************************
 * ECM MACROS				                                     *
 *****************************************************************************/

#define COL_UNIT                6       // 10^6
#define ELA_UNIT                3       // 10^3
#define HYA_UNIT                3       // 10^3

#define AV_NORM_COL_DLP          400.0f //990.0f	//50.0f                  // 9.9e8 = 990e6 monomers/patch
#define AV_NORM_COL_SLP          380.0f //450.0f                  // 3.8e8 = 380e6 monomers/patch
#define AV_NORM_COL_ILP          450.0f                  // 4.5e8 = 450e6 monomers/patch
#define AV_NORM_ELA_DLP         1900.0f                  // 1.9e6 = 1900e3 monomers/patch
#define AV_NORM_ELA_SLP         1300.0f                  // 1.3e6 = 1300e3 monomers/patch
#define AV_NORM_ELA_ILP          510.0f                  // 5.1e5 =  510e3 monomers/patch
#define AV_NORM_HYA_DLP         3600.0f                  // 3.6e6 = 3600e3 polymers/patch
#define AV_NORM_HYA_SLP         3100.0f                  // 3.1e6 = 3100e3 polymers/patch
#define AV_NORM_HYA_ILP         3600.0f                  // 3.6e6 = 3600e3 polymers/patch

#ifdef VISUALIZATION
#define FIB_COL_PROD_RATE         17.2f					 // 1.72e7 = 17.2e6			//394.0f          // 3.94e5 = 394.0e3
#define FIB_ELA_PROD_RATE         2000.0f          // 3.45e4 = 34.5e3
#define FIB_HYA_PROD_RATE         400.0f          // 6.58e4 = 65.8e3
#else
#define FIB_COL_PROD_RATE         17.2f					 // 1.72e7 = 17.2e6			//394.0f          // 3.94e5 = 394.0e3
#define FIB_ELA_PROD_RATE         34.5f //2000.0f          // 3.45e4 = 34.5e3
#define FIB_HYA_PROD_RATE         65.8f //400.0f          // 6.58e4 = 65.8e3
#endif

#define MAX_MONO_COL           AV_NORM_COL_DLP
#define MAX_MONO_ELA           AV_NORM_ELA_DLP
#define MAX_MONO_HYA					 AV_NORM_HYA_DLP


#define CONV_RATE_COL           AV_NORM_COL_SLP //AV_NORM_COL_ILP   //  1300e3 conversion from tropocollagen to fiber rate
#define CONV_RATE_ELA           AV_NORM_ELA_SLP   //   380e6 conversion from tropoelastin to fiber rate
#define CONV_RATE_HYA						AV_NORM_HYA_SLP		//   3.1e6 conversion from molecules to relative unit

#ifdef HUMAN_VF
#define MAX_COL 5.0f
#define MAX_ELA 2.0f    //4.0f
#define MAX_HYA 1.0f
#elif defined(RAT_VF)
#define MAX_COL 5.0f		//2.0f
#define MAX_ELA 2.0f    //4.0f
#define MAX_HYA 1.0f
#else
#define MAX_COL 2.0f
#define MAX_ELA 2.0f    //4.0f
#define MAX_HYA 1.0f
#endif


/*****************************************************************************
 * MACRO DEFINITIONS                                                         *
 *****************************************************************************/
// Number of ticks (30 min) of model simulation
#define NUM_TICKS	240
// Number of threads for parallelization
#define NUM_THREAD 32
// Maximum number of threads for parallelization
#define MAX_NUM_THREADS 32
// Number of threads allocated to preparing GPU diffusion data
#define PACKTH	2
// Default thread identification number
#define DEFAULT_TID	0
// Constant for ArrayChain. Array size: ~  32k   (32768)
#define DEFAULT_DATA_SMALL   1<<15
// Constant for ArrayChain. Array size: ~ 262k  (262144)
#define DEFAULT_DATA_MEDIUM  1<<18
// Constant for ArrayChain. Array size: ~ 524k  (524288)
#define DEFAULT_DATA_LARGE   1<<19
// Constant for ArrayChain. Array size: ~1048k (1048576)
#define DEFAULT_DATA_XLARGE  1<<20
// Constant for ArrayChain. Array size: ~2097k (2097152)
#define DEFAULT_DATA_XXLARGE 1<<21

//#define PRINT_KERNEL
//#define PRINT_SECRETION


/* NOTE: Please ignore these options (V_a, V_b) since they are for the
 *       currently incorrect diffuseChem() */
//#define V_a
#define V_b



// TODO(Nuttiiya)
#ifdef VECTORIZE
typedef float v4sf __attribute__ ((vector_size(sizeof(float)*4)));
union f4vector
{
	v4sf v;
	float f[4];
};
#endif

// TODO(Nuttiiya)
typedef int sizeType;

/*
 * Description:	Check if a field is modified in this tick (i.e. read and write
 *              entry are not the same). Used for optimizations in updates
 *
 * Return: True if modified
 *
 * Parameters: arr  -- Pointer to attribute to check
 */
template <typename T>
bool isModified(T* arr){
	return arr[read_t] != arr[write_t];
}

#endif /* COMMON_H_ */
