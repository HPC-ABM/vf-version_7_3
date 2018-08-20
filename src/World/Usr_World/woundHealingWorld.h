/*
 * woundHealingWorld.h
 *
 * File Contents: Contains declarations for the WHWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *
 * Created on Jun 17, 2013, 6:59 PM
 *****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 *****************************************************************************/  // TODO(Kim): Update the file comment once we figure out the copyright issues

#ifndef WHWORLD_H
#define	WHWORLD_H

#include "../World.h"
#include "../../FieldVariable/Usr_FieldVariables/WHChemical.h"
#include "../../Agent/Usr_Agents/Platelet.h"
#include "../../Agent/Usr_Agents/Fibroblast.h"
#include "../../Agent/Usr_Agents/Neutrophil.h"
#include "../../Agent/Usr_Agents/Macrophage.h"
#include "../../ECM/ECM.h"
#include "../../ArrayChain/ArrayChain.h"
#include "../../Diffusion/Diffuser.h"

//test

#include "../../common.h"


#include <stdlib.h>
#include <vector>
#include <new>

#ifdef VISUALIZATION
#ifdef OVERLAP_VIS
#ifdef MODEL_3D
#include "../../Visualization/3D/Visualization.h"
#endif  // MODEL_3D
#endif	// OVERLAP_VIS
#endif  // VISUALIZATION

#ifdef GPU_DIFFUSE
// Include CUDA runtime and CUFFT
#include <cuda_runtime.h>
#include <cufft.h>

// Helper functions for CUDA
#include <helper_functions.h>
#include <helper_cuda.h>

#include "../../Diffusion/convolutionFFT_common.h"
#endif  // GPU_DIFFUSE



class Platelet;
class Fibroblast;
class Neutrophil;
class Macrophage; 
class Collagen;
class Elastin;
class Hyaluronan; 
class ECM;

using namespace std;

/*
 * WHWORLD (WOUND HEALING WORLD) CLASS DESCRIPTION:
 * WHWorld is a derived class of the parent class World.
 * The WHWorld class manages the model world.
 * It is used to initialize cells, ECM, patches, and chemicals; to destroy
 * agent ArrayChains; to execute each timestep of the model; to sprout agents;
 * to count patches; and to output data.
 */
class WHWorld: public World {
 public:
    /*
     * Description:	WHWorld constructor. 
     *
     * Return: void
     *
     * Parameters: width    -- Width (x dimension) of the world in millimeters
     *             length   -- Length (y dimension) of the world in millimeters 
     *             height   -- Height (z dimension) of the world in millimeters
     *             plength  -- Length of each patch (grid point) in millimeters
     */
    WHWorld(double width = 5,
    		double length = 7.2,
    		double height = 0.015,
    		double plength = 0.04
    		);

    /*
     * Description:	WHWorld destructor. 
     *
     * Return: void
     *
     * Parameters: void
     *
     * NOTE: These functions are called implicitly:
     *         plats.~ArrayChain()
     *         neus.~ArrayChain()
     *         macs.~ArrayChain()
     *         fibs.~ArrayChain()
     */
    ~WHWorld();
  
    /*
     * Description:	Destructor function for the Platelet ArrayChain 
     *
     * Return: void
     *
     * Parameters: &agent  -- Reference to platelet that will be destroyed
     */
    void destroyPlat(Platelet* &agent);

    /*
     * Description:	Destructor function for the Neutrophil ArrayChain 
     *
     * Return: void
     *
     * Parameters: &agent  -- Reference to neutrophil that will be destroyed
     */
    void destroyNeu(Neutrophil* &agent);

    /*
     * Description:	Destructor function for the Macrophage ArrayChain 
     *
     * Return: void
     *
     * Parameters: &agent  -- Reference to macrophage that will be destroyed
     */
    void destroyMac(Macrophage* &agent);

    /*
     * Description:	Destructor function for the Fibroblast ArrayChain 
     *
     * Return: void
     *
     * Parameters: agent  -- Pointer to fibroblast that will be destroyed
     */
    void destroyFib(Fibroblast* &agent);

    /*
     * Description:	Assign a patch type to each patch in the world in row major
     *              index manner
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializePatches();


    /*
     * Description:	Initializes chemical concentrations in each patch and 
     *              initializes the total concentration of each chemical
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeChemCPU(); 

#ifdef GPU_DIFFUSE
    /*
     * Description:	Initializes chemical concentrations in each patch and 
     *              initializes the total concentration of each chemical
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeChemGPU(); 
#endif

    /*
     * Description:	Initializes chemical concentrations in each patch and 
     *              initializes the total concentration of each chemical
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeChem(); 

    /*
     * Description:	Initializes all fibroblasts to their correct patches
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeFibroblasts();

    /*
     * Description:	Initializes all macrophages to their correct patches
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeMacrophages();

    /*
     * Description:	Initializes all neutrophils to their correct patches
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeNeutrophils();

    /*
     * Description:	Initializes collagen, elastin and hyaluronan to their 
     *              correct patches
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeECM();

    /*
     * Description: Initializes collagen, elastin polymers and hyaluronan per LP section
     *
     * Return: void
     *
     * Parameters:  x_begin         -- beginning coordinate of x for this LP section
     *              x_end           -- end coordinate of x for this LP section
     *              APc             -- Average collagen concentration for this LP section
     *              APe             -- Average elastin concentration for this LP section
     *              APh             -- Average hyaluronan concentration for this LP section
     *              fractionLP      -- fraction of this section related to the whole LP
     */
    void initPolymers(
            int x_begin,
            int x_end,
            float APc,
            float APe,
            float APh,
            float fractionLP);
    void addCol(int in, float npolymers);
    void addEla(int in, float npolymers);
    void addHya(int in, float nh, int multiplier = 1);

    /*
     * Description:	Check whether a given patch qualifies as a center for 
     *              sprouting original hyaluronan.
     *
     * Return: true if the patch can sprout original hyaluronan
     *
     * Parameters: a       -- x-coordinate of target patch
     *             b       -- y-coordinate of target patch
     *             c       -- z-coordinate of target patch
     *             radius  -- specifies how large the zone containing only 
     *                        tissue around a sprouting center must be
     *             gap     -- specifies how far apart the tissue-only zones
     *                        of sprouting centers must be
     */
    bool initHARadius(
    		int a,
    		int b,
    		int c,
    		double radius,
    		double gap
    		);

    int countDamage();

    int getInitialDam();

    /*
     * Description:	Initializes all damaged patches. Sprouts platelets and 
     *              fragments ECM on the damaged patches.
     *
     * Return: void
     *
     * Parameters: void
     */
    void initializeDamage();

    /*
     * Description:	Each call to go() simulates 30 minutes or 'real-world' time
     *              of the biological model. 
     *
     * Return: 0 on success
     *
     * Parameters: void
     */
    int go();             

    /*
     * Description:	Entry function for sprouting cells
     *              Selects the appropriate sprouting function to apply
     *
     * Return: void
     *
     * Parameters: num          -- Number of cells to sprout
     *             patchType    -- Type of patches of sprout on
     *             agentType    -- Type of agent to sprout
     *
     *             Physical boundaries of the sprouting area/volume:
     *             xmin         -- Left
     *             xmax         -- Right
     *             ymin         -- Top
     *             ymax         -- Bottom
     *             zmin         -- Near
     *             zmax         -- Far
     *
     *             bloodOrTiss  -- Pass in true if should be sprouted in blood
     * 							               Pass in false if should be sprouted in 
     *                             tissue (default)
     */
    void sproutAgent(
    		int     num,
    		int     patchType,
    		agent_t agentType,
    		int     xmin, int   xmax,
    		int     ymin, int   ymax,
    		int     zmin, int   zmax,
    		int bloodORtiss = 0
    		);

    /*
     * Description:	Function for sprouting cells in a given area/volume
     *
     * Return: void
     *
     * Parameters: num          -- Number of cells to sprout
     *             patchType    -- Type of patches of sprout on
     *             agentType    -- Type of agent to sprout
     *
     *             Physical boundaries of the sprouting area/volume:
     *             xmin         -- Left
     *             xmax         -- Right
     *             ymin         -- Top
     *             ymax         -- Bottom
     *             zmin         --
     *             zmax         --
     *
     *             bloodOrTiss  -- Pass in true if should be sprouted in blood
     *                             Pass in false if should be sprouted in 
     *                             tissue (default)
     */
    void sproutAgentInArea(
    		int     num,
    		int     patchType,
    		agent_t agentType,
    		int     xmin, int   xmax,
    		int     ymin, int   ymax,
    		int     zmin, int   zmax,
    		int bloodORtiss = 0
    		);

    /*
     * Description:	Function for sprouting cells in the whole world
     *
     * Return: void
     *
     * Parameters: num          -- Number of cells to sprout
     *             patchType    -- Type of patches of sprout on
     *             agentType    -- Type of agent to sprout
     *             bloodOrTiss  -- Pass in true if should be sprouted in blood
     *                             Pass in false if should be sprouted in tissue (default)
     */
    void sproutAgentInWorld(
    		int num,
    		int patchType,
    		int agentType,
    		bool bloodORtiss=0
    		);

    /*
     * Description:	Counts the number of patches of a given type
     *
     * Return: int  -- Number of patches of given patch type
     *
     * Parameters: int  -- Enumic value for the patch type to count 
     */
    int countPatchType(int);

    /****************************************************************
     * HELPER SUBROUTINES                                           *
     ****************************************************************/

    bool searchIntArray(int *arr, int size, int v);

    void sproutCellsInArea(int num, int patchType, int cellType,
                int xmin, int xmax, int ymin, int ymax, int zmin, int zmax, bool bloodORtiss);

    /*
     * Description:	Instantiate a new cellType object and update
     * 			its coordinate to the same as tempPatchPtr.
     *
     * Return: none
     *
     * Parameters: cellType     -- type of cell to create
     * 		   tempPatchPtr -- pointer to patch to put cell on
     */
    void createCell(int cellType, bool bloodORtiss, Patch* tempPatchPtr);

    /*
     * Description:	Converts a length in millimeters to the number of patches
     *
     * Return: The number of patches
     *
     * Parameters: mm  -- Length in millimeters
     */
    int mmToPatch(double mm);

    /*
     * Description:	Set seed for random number generator 
     *
     * Return: void
     *
     * Parameters: seed    -- Desired seed for random number generator
     *
     */
    static void setSeed(unsigned s);

    /*
     * Description: Determines the number of ticks elapsed
     *
     * Return: The number of ticks elapsed
     *
     * Parameters: void
     */

    static double getTick();

    /*
     * Description:	Converts hours and days into ticks
     *
     * Return: The number of ticks
     *
     * Parameters: hour  -- Number of hours
     *             day   -- Number of days
     */
    static int reportTick(int hour = 0, int day = 0);

    /*
     * Description:	Determines the number of hours elapsed
     *
     * Return: The number of hours elapsed
     *
     * Parameters: void
     */
    static double reportHour();

    /*
     * Description:	Determines the number of days elapsed
     *
     * Return: The number of days elapsed
     *
     * Parameters: void
     */
    static double reportDay();

    /*
     * Description:	Determines the number of neighbors with given patch type
     *
     * Return: The number of neighbors with given patch type
     *
     * Parameters: ix         -- x coordinate of current patch
     *             iy         -- y coordinate of current patch
     *             iz         -- z coordinate of current patch
     *             patchType  -- patch type to count
     */
    int countNeighborPatchType(int ix, int iy, int iz, int patchType);

    /*
     * Description:	Reads user input from config file (default: config.txt) to initialize
     *              chemicals, wound, cells.
     *
     * Return: Returns 0 if function proceeded to completion, for testing.
     *         Can be removed later.
     *
     * Parameters: void
     */
    int userInput();

    /*
     * Description:	Initialize chemical mass on each patch.
     * 				Initialize chemical gradients on each patch.
     * 				Initialize world total chemical mass.
     *
     * Return: void
     *
     * Parameters: void
     */
    void initPatchChem();

#ifdef GPU_DIFFUSE
    /****************************************************************
     * HELPER SUBROUTINES - GPU                                     *
     ****************************************************************/

    /*
     * Description:	Calculate and initialize constants used in convolution
     * 				based diffusion on GPU. This includes:
     * 					- Lambda		= D*dt/dx^2
     * 					- Half lives
     * 					- Decay rates
     * 					- Patch baseline cytokine concentrations
     *
     * Return: void
     *
     * Parameters:	lambda		--	pointer to array of 'numChem' floats
     * 				gamma		--	pointer to array of 'numChem' floats (decay)
     * 				numChem		--	number of base chemical types
     * 				kernel_cctx --  pointer to kernel chem context to initialize dt and dx2
     */
    void initKernelConstants(float* lambda, float *gamma, int numChem, c_ctx *kernel_cctx);

#ifdef OPT_CHEM
    /*
     * Description:	Allocate input and output PINNED buffers on host for GPU
     * 				data preparation. Float arrays of (nx*ny*nz) elements are
     * 				allocated for:
     * 					- h_diffusion_ibuffs
     * 					- h_diffusion_obuffs
     *
     * Return: void
     *
     * Parameters:	void
     */
    void allocateHostDiffusionBuffers();
#else	// OPT_CHEM
    /*
     * Description:	Allocate 'typesOfChem' array of (nx*ny*nz) floats for
     * 					- chemAllocation
     * 				Allocate 'typesOfBaseChem' array of (nx*ny*nz) floats for
     * 					- h_diffusion_results
     *
     * Return: void
     *
     * Parameters:	void
     */
    void allocateChemArrays();

#endif	// OPT_CHEM



#endif	// GPU_DIFFUSE

    /****************************************************************
     * OUTPUT SUBROUTINES & VISUALIZATION                           *
     ****************************************************************/

    /*
     * Description:	Outputs cell counts and cytokine levels from the current
     *              tick to the file "Output/Output_Biomarkers.csv".
     *              Used for testing.
     *
     * Return: void
     *
     * Parameters: void
     */
    void outputWorld_csv();

    /*
     * Description:	Outputs all patch assignments (patch type, agent type,
     *              ECM type) to files in output directory.
     *
     * Return: void
     *
     * Parameters: void
     */
    void patchassign_csv();

    /*
     * Description: Update aggregated statistics for output function
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateStats();

    /*
     * Description: Update aggregated cell statistics for visualization
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateCellStats();


    /*
     * Description: Get corner coordinates of retangular prism contain the wound area
     *
     * Return: void
     *
     * Parameters:	w<dim>b -- Begining coordinate of wound in <dim>-dimension
     * 			w<dim>e -- End coordinate of wound in <dim>-dimension
     */
    void getWndPos(int &wxb, int &wxe, int &wyb, int &wye, int &wzb, int &wze);

//#ifdef VISUALIZATION
    /*
     * Description: Increase number of ECM proteins in ecmMap
     *
     * Return: void
     *
     * Parameters:	index		-- World index of patch to increase proteins
     * 							ecmType	-- Type of ECM to be added
     * 							count		-- Amount to be added
     */
    void incECM(int index, ecm_i ecmType, float count);

    /*
     * Description: Decrease number of ECM proteins in ecmMap
     *
     * Return: void
     *
     * Parameters:	index		-- World index of patch to decrease proteins
     * 							ecmType	-- Type of ECM to be removed
     * 							count		-- Amount to be removed
     */
    void decECM(int index, ecm_i ecmType, float count);

    /*
     * Description: Set number of ECM proteins in ecmMap
     *
     * Return: void
     *
     * Parameters:	index		-- World index of patch to det number of proteins
     * 							ecmType	-- Type of ECM to be set
     * 							count		-- Amount to be set to
     */
    void setECM(int index, ecm_i ecmType, float count);
    void resetECM(int index, ecm_i ecmType);
    float getECM(int index, ecm_i ecmType);

    void resetECMmap();


    /****************************************************************
     * MAPS                                                 				*
     ****************************************************************/
    // TODO: Consolidate the maps
    // Map of matured ECM (collagen, elastin, HA)
    float *ecmMap[m_ecmtotal];
    // Map of changes in matured ECM (collagen, elastin, HA)
    float *ecmGradMap[m_ecmtotal];
    // Map of amplified changes and scaled down initial values in matured ECM (collagen, elastin, HA)
    float *ecmPreProcMap[m_ecmtotal];

#ifdef AVEP
    /****************************************************************
     * AVEP DIRTY MAPS                                         		  *
     ****************************************************************/
#ifdef AVEP_INC
    int _iround;

    void resetIncRound();
#endif
    int mapW_AVEP;
    int mapH_AVEP;
    int mapD_AVEP;
    int _mread;
    int _mwrite;

    bool *ecmMapAVEP[2][m_ecmtotal];	// read and write map buffers

    int getMapIndexAVEP(int x, int y, int z);

    bool isDirtyAVEP(int mapx, int mapy, int mapz, ecm_i ecmType);

    void setDirtyAVEP(int x, int y, int z, ecm_i ecmType);

    void resetDirtyAVEP(int mapx, int mapy, int mapz, ecm_i ecmType);

    void resetMapAVEP();

    void printMapAVEP(ecm_i ecmType);

    void updateMapIndexAVEP();

#endif

//#endif

    /****************************************************************
     * STATIC VARIABLES                                             *
     ****************************************************************/
           
    // Keeps track of the current tick
    static double clock;
    // Used to generate random numbers
    static unsigned seed;
    // Width of the basement membrane center (in patches)
    static int bmcx, bmcy;
    // Width of epithelium (in millimeters)
    static float epithickness;
#ifdef MODEL_VOCALFOLD
    // Boundaries of vocal fold lamina propria layers (SLP, ILP, DLP)
    static int SLPxmin, SLPxmax;
    static int ILPxmin, ILPxmax;
    static int DLPxmin, DLPxmax;
    static float VFvolumefraction;
#endif
    // Resonant voice impact(bad)/vibratory(good) stress
    static int RVIS, RVVS;
    // Spontaneous speech impact(bad)/vibratory(good) stress
    static int SSIS, SSVS;
    // Whether there is high TNF damage (which results in ECM fragmentation)
    static bool highTNFdamage;
    // The number of patches per millimeter in the world
    static float patchpermm; 
    // The number if initial tissue patches
    static int initialTissue; 

    /* CALIBRATION Variables */
    // The threshold for TNF damage
    static float thresholdTNFdamage;
    // The threshold for MMP8 damage
    static float thresholdMMP8damage;
    // The number of hours between agent sprouting sessions
    static float sproutingFrequency[6];
    // Constants related to the number of agents to sprout
    static float sproutingAmount[14];
    // The decay rates of the cytokines
    static float cytokineDecay[8];
    // The half lifes of the cytokines in minutes
    static float halfLifes_static[8];


    /****************************************************************
     * CONSTANT VARIABLES                                           *
     ****************************************************************/


    // Maximum number of thousands of collagen fibers allowed per patch
    float maxFCol;
    // Maximum number of thousands of collagen fiber changes per patch
    float maxChangeFCol;
    /* Keeps track of the type of treatment that is being administered
     * (voice rest, spontaneous speech, resonant voice) */
    int treatmentOption;
    // The length of each patch
    double patchlength;
    // x,y coordinates of the wound center
    int woundX, woundY;
    /* wound[0]: wound depth (x-radius) (in millimeters) 
     * wound[1]: wound radius in y-dimension
     * wound[2]: wound radius in z-dimension
     * wound[3]: wound severity */
    float wound[4];
    // Radius of the capillaries
    float capRadius;
    // x-coordinates of each capillary center
    vector<float> capX;
    // y-coordinates of each capillary center
    vector<float> capY;
    // Pointer to instanct of diffusion manager object
    Diffuser *diffuserPtr;
    // Pointer to instance of type to manage chemicals in the world
    WHChemical* WHWorldChem;

#ifdef GPU_DIFFUSE		// Diffusion using Third buffer
    /*typedef struct CCTX		// convolution context
    {
	int KH;
	int KW;
	int KX;
	int KY;
	int DH;
	int DW;
	int FFTH;
	int FFTW;
    } c_ctx;
*/
    typedef struct CCTX c_ctx;
#ifdef OPT_CHEM
    float**   h_diffusion_ibuffs;
    float**   h_diffusion_obuffs;
#else	// OPT_CHEM
    float**   h_diffusion_results;
#endif
    fComplex** d_kernel_spectrum;
    fComplex** h_dKernel_spectrum;
    c_ctx*     chem_cctx;
#endif

    // The number of different chemicals there are in the world
    int typesOfChem;
    // The number of different case chemicals there are in the world
    int typesOfBaseChem;
    // Initial amount of damage
    int initialDam;
    // Current amount of damage
    int totaldamage;
    // Initial amount of each chemical in the world
    vector<float> baselineChem;
#ifdef GPU_DIFFUSE
    // Baseline chem per patch
    float *pbaseline;
#endif
    // Pointer to the array of patches
    Patch* worldPatch;
    // Pointer to array of ECM
    ECM* worldECM;
    // Initial amount of each cell type (agent) in the world
    vector<int> initialCells;
    // ArrayChain to manage all platelet data
    ArrayChain<Platelet*> plats;
    // ArrayChain to manage all fibroblast data
    ArrayChain<Fibroblast*> fibs;
    // ArrayChain to manage all macrophage data
    ArrayChain<Macrophage*> macs;
    // ArrayChain to manage all neutrophil data
    ArrayChain<Neutrophil*> neus;
    // Vector of pointers to local lists of platelet pointers to add to global list
    vector<Platelet*>* localNewPlats[MAX_NUM_THREADS];
    // Vector of pointers to local lists of fibroblast pointers to add to global list
    vector<Fibroblast*>* localNewFibs[MAX_NUM_THREADS];
    // Vector of pointers to local lists of macrophage pointers to add to global list
    vector<Macrophage*>* localNewMacs[MAX_NUM_THREADS];
    // Vector of pointers to local lists of neutrophil pointers to add to global list
    vector<Neutrophil*>* localNewNeus[MAX_NUM_THREADS];
    // Vector of patches which can be centers for sprouting original hyaluronan
    vector<int> initHAcenters;
    // Seeds used to generate random numbers for each thread
    unsigned seeds[MAX_NUM_THREADS];
    // Array of diffusion coefficients, gets allocated in userInput()
    float *D;
    // Array of cytokine half-life (seconds), gets allocated in userInput()
    int *HalfLifes;

    /****************************************************************
     * AGGREGRATED STATS                                            *
     ****************************************************************/
    float totalOC,   totalNC,   totalFC;
    float totalOE,   totalNE,   totalFE;
    float totalHA,   totalFHA;
    int totalCell[p_celltotal];
      
    /****************************************************************
     * VISUALIZATION VARIABLES                                      *
     ****************************************************************/
    int wnd_xb, wnd_xe;
    int wnd_yb, wnd_ye;
    int wnd_zb, wnd_ze;

    /****************************************************************
     * DEBUG VARIABLES                                               *
     ****************************************************************/
    float chemSecreted[4][8];
    int chemSecretedCoord[4][8][3];
    float maxPatchChem[8][3]; // pChem, dChem, tChem
    int maxPatchChemCoord[8][3];
    float maxOldPatchChem[8]; // old pChem
    int maxOldPatchChemCoord[8][3];
    float minOldPatchChem[8]; // old pChem
    int minOldPatchChemCoord[8][3];
    void printChemInfo();

/*
    int deadneus;
    int dead_afibs;
    int deactfibs;
    int actfibs;
    int newfibs;
    int deadfibs;

  int maxHAsize;
  int HAfrags;
  int newOcolls;
*/
#ifdef PROFILE_ECM
  long ECMrepairTime;
  long HAlifeTime   ;
  long ECMdangerTime;
  long ECMscarTime  ;
#endif

 private:
  /****************************************************************
   * MAJOR SECTION SUBROUTINES - begin                            *
   ****************************************************************/

    /*
     * Description: (Stage 0)	Sprout agents on various patches.
     *
     * Return: void
     *
     * Parameters: hours  -- Current hour in model execution
     */
    void seedCells(float hours);

    /*
     * Description:	(Stage 1)	Entry point function for diffusing chemicals
     * 							Selects the appropriate diffusion function to apply
     *
     * Return: void
     *
     * Parameters: coeff  -- Diffusion coefficient
     *             dt     -- Time step
     *             NOTE: Examples of diffusion coefficients can be found at
     * 					 	   http://www.math.ubc.ca/~ais/website/status/diffuse.html
     */
#ifdef MODEL_3D
    void diffuseCytokines(float dt = 0.02);     // 3D stability condition, dt < dx^2/6*D min = 0.020833
                                                //assuming dx = dy = dz = 0.015 mm, D = 0.0018 mm^2/min :
#else
    void diffuseCytokines(float dt = 0.04);//0.03);     // 2D stability condition, dt < dx^2/4*D min = 0.03125
                                                // assuming dx = dy = 0.015 mm, D = 0.0018 mm^2/min :
#endif

    /*
     * Description:	Helper function for cell execution. Execute platelet cell 
     *              function for all living platelets.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executePlats();

    /*
     * Description:	Helper function for cell execution. Execute neutrophil cell
     *              function for all living neutrophils.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeNeus();

    /*
     * Description:	Helper function for cell execution. Execute macrophage cell
     *              function for all living macrophages.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeMacs();

    /*
     * Description:	Helper function for cell execution. Execute all alive fibroblasts.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeFibs();

    /*
     * Description:	(Stage 2)	Execute cell functions for all living cells.
     *
     * Return: void
     *
     * Parameters: void
     */
    void executeCells();

    /*
     * Description:	(Stage 3)	Execute ECM functions.
     *
     * Return: void
     *
     * Parameters: void
     */
    void executeECMs();

    /*
     * Description:	(Stage 3)	Fragment ECM proteins if necessary.
     *
     * Return: void
     *
     * Parameters: void
     */
    void requestECMfragments();

    /*
     * Description:	Helper function for diffuseCytokines(). Update chemical 
     * 				      concentration of all chemicals at each patch over a given 
     *              time step using Netlogo (U. Wilensky, Northwestern 
     *              University, Evanston Illinois) diffusion procedure.
     *
     * Return: void
     *
     * Parameters: void
     */
    void NetlogoDiffuse();

    /*
     * Description:	Helper function for diffuseCytokines(). Update chemical 
     *              concentration of all chemicals at each patch over a given
     * 				      time step using partial differential diffusion equation.
     *
     * Return: void
     *
     * Parameters: coeff  -- Diffusion coefficient (mm^2/min)
     *             dt     -- Time step ( min) 
     *             NOTE: Examples of diffusion coefficients can be found at
     * 					 	   http://www.math.ubc.ca/~ais/website/status/diffuse.html
     */
    void diffuseChem(float dt);

    
#ifdef GPU_DIFFUSE
    /*
     * Description:	Helper function for diffuseCytokines().
     * 			Convolution-based chemical diffusion executed on GPU.
     *
     * Return: void
     *
     * Parameters: void
     * Note:	All parameters are assumed to have been intialized in cctx_t (convolution
     * 		context) via initializeChemGPU()
     */
    void diffuseChemGPU();
#endif

    /*
     * Description:	(Stage 4a)	Update chemicals to reflect next tick's states
     * 			Update in the following manner:
     * 				p<chem> = d<chem> + t<chem>*(1-gamma)
     * 			where gamma is a cytokine specific constant derived from the cytokine's
     * 			halflife.
     * 				gamma = 1 - 2^(-1/halflife)
     * Return: void
     *
     * Parameters: void
     */
    void updateChem();

    
    /*
     * Description:	Helper function for ECM updates. Execute updates for ALL ECM managers.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeAllECMUpdates();

    /*
     * Description:	Helper function for ECM updates. Execute request resets for ALL ECM managers.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeAllECMResetRequests();

    /*
     * Description:	(Stage 4c)	Update ECM managers to reflect next tick's states
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateECMManagers();

    /*
     * Description:	(Stage 4d)	Update patches to reflect next tick's states
     *
     * Return: void
     *
     * Parameters: void
     */
    void updatePatches();

    /*
     * Description:	Helper function for cell updates. Update all platelets.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline updatePlats();

    /*
     * Description:	Helper function for cell updates. Update all neutrophils.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline updateNeus();

    /*
     * Description:	Helper function for cell updates. Update all macrophages.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline updateMacs();

    /*
     * Description:	Helper function for cell updates. Update all fibroblasts.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline updateFibs();

    /*
     * Description:	(Stage 4b)	Update cells to reflect next tick's states
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateCells();

    /*
     * Description:	Update cells to reflect next tick's states. 
     *              This is called instead of updateCells() during setup.
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateCellsInitial();
};
#endif	/* WHWORLD_H */
