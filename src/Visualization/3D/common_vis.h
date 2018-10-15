//
//  common.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 1/3/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_common_h
#define OpenGL_SimpleWH_common_h

#include <stdio.h>
#include <stdlib.h>


//#ifdef __APPLE__
//	#include <GLUT/glut.h>    // Header File For The GLUT Library
//      #include <OPENGL/gl.h>  // Header File For The OpenGL32 Library
//      #include <OPENGL/glu.h> // Header File For The GLu32 Library
//#else
#define GLEW_STATIC
#include "./glew/GL/glew.h"
#include <GL/glut.h>
//#include <GL/glu.h>
//#include <GL/gl.h>
//#endif

#include <unistd.h>     // needed to sleep
#include <string.h>
#include <vector>
#include <array>

//#include "../../Utilities/error_utils.h"

//#define GLM_MESSAGES
#define GLM_FORCE_CXX98
// TODO: Change directory
#include </fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v5/glm/glm/glm.hpp>
#include </fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v5/glm/glm/gtc/matrix_transform.hpp>
#include </fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v5/glm/glm/gtx/compatibility.hpp>
#include </fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v5/glm/glm/gtc/type_ptr.hpp>

#include <math.h>

#include "../../common.h"
#include "draw_util.h"

//// Externs for 3D ECM/Chem Visualization
//extern __device__ surface<void, cudaSurfaceType3D> srfChem0;
//extern __device__ texture<VolumeType, 3, cudaReadModeElementType> texChem0;


/* ASCII code for the escape key. */
#define ESCAPE 27

#define SCREENW 1440
#define SCREENH 900//700

#define WINW    SCREENW//1200
#define WINH    SCREENH//700



// TODO: read from command line input
//#define GRIDW   1160//1660//1440//1200 //740
//#define GRIDH   1660//1160//900//700 //480
#define STEP    1

typedef enum {
  tnf,
  tgf,
  fgf,
  mmp8,
  il1beta,
  il6,
  il8,
  il10,
  numchem
} ChemType;

typedef enum {
  remove_black,
  remove_white,
  add_alpha_channel
} ImLdOption;

typedef enum {
  tx_scar,      // scar
  tx_00h,       //                      HA
  tx_0e0,       //            elastin
  tx_0eh,       //            elastin + HA
  tx_c00,       // collagen
  tx_c0h,       // collagen           + HA
  tx_ce0,       // collagen + elastin
  tx_ceh        // collagen + elastin + HA
} TexP;         // Texture on patch

typedef std::array<glm::vec3,4> QuadV;


#define CHEM_STEP 0.02//3
#define CHEM_REDUCTION_FACTOR 0.98f
#define CHEM_CURRENT_PATCH_FACTOR 0.8f

#define NUM_PATCH_IN_BLOCK 8

/* Half one type of agents, another half is another type */
#define NPOINTS 2000//182000//1000

//GLfloat x[NPOINTS];
//GLfloat y[NPOINTS];
//GLfloat z[NPOINTS];


/* ECM */
#define MAX_NCOLL       15
#define MAX_COLL_SCAR   100

/* VBO Macros */
#define TOTAL_TNF_VBO     0
#define TOTAL_TGF_VBO     1
#define TOTAL_FGF_VBO     2
#define TOTAL_MMP8_VBO    3
#define TOTAL_IL1beta_VBO 4
#define TOTAL_IL6_VBO     5
#define TOTAL_IL8_VBO     6
#define TOTAL_IL10_VBO    7
#define NUM_VBOs      TOTAL_IL10_VBO+1

/* Wound grid */
#define WG_BORDER_X       10   // patches
#define WG_BORDER_Y       10   // patches

/* Heat Map */
#define HML_RESOLUTION    100
#define N_HMLV            5

/* Grid Surface */
#define CHEM_AMP_FACTOR   10//5000
#define CHEM_GRID_MAX     20.0

/* Option Table*/
#define CHEM_OPTION_COLS 1 //2
#define CHEM_OPTION_ROWS 8
#define CELL_OPTION_COLS 3
#define CELL_OPTION_ROWS 1

/* Plot specs */
#define MAX_PLOTS	8
#define VERT_PLOT
#define N_INTV_X    4
#define N_INTV_Y    4

#ifdef VERT_PLOT
#define PLOT_SCREEN_FRACTION  5
#define YOFFSET     20//100

#define MAX_PLOT_POINTS   194//280
#define MIN_X       0
#define MAX_X       (MAX_PLOT_POINTS-1)
#define MIN_Y       0
#define MAX_Y       500
#define PLOT_SHRINK_FACTOR 5000
#define PLOT_AREA_FACTOR  2.5f   // 1/PLOT_AREA_FACTOR
#define PLOTw       (WINW/PLOT_SCREEN_FRACTION)
#define PLOTh       (WINH/(2*PLOT_AREA_FACTOR))
#define PLOT0x      (0 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT0y      ((WINH*3)/(2*PLOT_AREA_FACTOR)) + YOFFSET
#define PLOT1x      (1 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT1y      ((WINH*3)/(2*PLOT_AREA_FACTOR)) + YOFFSET
#define PLOT2x      (0 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT2y      (WINH/(1*PLOT_AREA_FACTOR)) + YOFFSET
#define PLOT3x      (1 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT3y      (WINH/(1*PLOT_AREA_FACTOR)) + YOFFSET
#define PLOT4x      (0 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT4y      (WINH/(2*PLOT_AREA_FACTOR)) + YOFFSET
#define PLOT5x      (1 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT5y      (WINH/(2*PLOT_AREA_FACTOR)) + YOFFSET
#define PLOT6x      (0 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT6y       0 + YOFFSET
#define PLOT7x      (1 * WINW/PLOT_SCREEN_FRACTION)
#define PLOT7y       0 + YOFFSET

#else   // VERT_PLOT

#define MAX_PLOT_POINTS   500//200
#define MIN_X       0
#define MAX_X       (MAX_PLOT_POINTS-1)
#define MIN_Y       0
#define MAX_Y       500
#define PLOT_SHRINK_FACTOR 5000
#define PLOT_AREA_FACTOR  2.5f   // 1/PLOT_AREA_FACTOR
#define CHM_PLOTw   (WINW/4)
#define CHM_PLOTh   (WINH/(2*PLOT_AREA_FACTOR))
#define TNFx        (0 * WINW/4)
#define TNFy        (WINH/(2*PLOT_AREA_FACTOR))
#define TGFx        (1 * WINW/4)
#define TGFy        (WINH/(2*PLOT_AREA_FACTOR))
#define FGFx        (2 * WINW/4)
#define FGFy        (WINH/(2*PLOT_AREA_FACTOR))
#define MMP8x       (3 * WINW/4)
#define MMP8y       (WINH/(2*PLOT_AREA_FACTOR))
#define IL1betax    (0 * WINW/4)
#define IL1betay     0
#define IL6x        (1 * WINW/4)
#define IL6y         0
#define IL8x        (2 * WINW/4)
#define IL8y         0
#define IL10x       (3 * WINW/4)
#define IL10y        0

#endif  // VERT_PLOT

#define ARC4RANDOM_MAX      0x100000000
/* ========================================================================== */
#define PERSPECTIVE_DEPTH   6000.0f
#define PLANE_DEPTH         -1500.0f
//#define COL_PLANE_H         10//100           // Collagen plane vertical distance from tissue plane

//#define GL_X(x)       (x - (GRIDW/2))
//#define GL_Y(y)       (-(y - (GRIDH/2)))
#define GL_Z(z)		(PLANE_DEPTH + z)

/* ========================================================================== */

/* Epithelium bound in y direction */
#define EPI_BOUND   30
#define EPI_SIN_AMP 2
#define PATCHW      0.015f

#define WOUNDW      5.0f      // 3mm
#define WOUNDH      2.0f      // 1mm

/* Direction macros */
#define UP_LEFT     0
#define UP          1
#define UP_RIGHT    2
#define LEFT        3
#define DOWN_RIGHT  (UP_LEFT+4)
#define DOWN        (UP+4)
#define DOWN_LEFT   (UP_RIGHT+4)
#define RIGHT       (LEFT+4)
#define STAY        8

#define USE_TEXTURE
#define SHOW_CHEM
#define SHOW_CHEM_SURF

//#define USE_QUAD_GRID
//#define USE_SPRITE

#define MINCOLOR		0.1f							//The lowest any RGB value can go (RGB values range from 0.0 - 1.0)
#define MAXCOLOR		0.6f							//The highest any RGB value can go 
#define COLORSPEED		0.5f							//A multiplier for how fast the colors change.

#define MINCOLOR2		0.33f							//The lowest any RGB value can go (RGB values range from 0.0 - 1.0)
#define MAXCOLOR2		0.66f							//The highest any RGB value can go 
#define COLORSPEED2		0.5f							//A multiplier for how fast the colors change.

extern int GridH;
extern int GridW;
extern int GridD;

int mmToPatch(float mm);

int GL_X(int x);
int GL_Y(int y);

#endif
