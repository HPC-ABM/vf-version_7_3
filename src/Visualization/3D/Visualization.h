//
//  Visualization.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 12/28/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_Visualization_h
#define OpenGL_SimpleWH_Visualization_h

#include "common_vis.h"

#include "UI.h"
#include "texture_util.h"
//#include "Terrain.h"
//#include "SubPlane.h"
#include "PlotManager.h"
#include "HeatMapManager.h"
//#include "GridSurfaceManager.h"
#include "OptionTable.h"
#include "WorldVisualizer.h"
#include "VolumeManager.h"

#include "../../World/Usr_World/woundHealingWorld.h"


namespace visualization {
  extern float totalChem[MAX_PLOTS];
  extern float totalCell[MAX_PLOTS];

  extern float xcam;
  extern float ycam;
  extern float zcam;

  extern WHWorld* world_ptr;
  extern int iter;
  
  extern int argc;
  extern char** argv;

  extern WorldVisualizer     *EnvVisualizer_ptr;
//  extern Terrain             *ScarTerrain_ptr;
//  extern SubPlane            *ECMplane_ptr;
  extern HeatMapManager      *HeatMapManager_ptr;
  extern PlotManager         *ChemPlotManager_ptr;
  extern PlotManager         *CellPlotManager_ptr;
//  extern GridSurfaceManager  *ChemSurfManager_ptr;
  extern OptionTable         *ChemOptionTable_ptr;
  extern OptionTable         *CellOptionTable_ptr;
  extern VolumeManager       *ECMvolumeManager_Ptr;
  
//  extern GLuint  tissue_texture[3];	/* Storage for 3 textures. */
//  extern GLuint  glass_texture[3];
//  extern GLuint  texture_atlas;		/* Storage for 1 texture. */
//  extern GLuint  scar_texture;
  
#ifdef USE_SPRITE
  extern GLuint  texture_sprite[3];	/* Storage for 3 textures. */
#endif
  
  
  // Text labels for table:
  // Text labels for table:
  extern char *optionHeader[CHEM_OPTION_COLS];
  extern char *optionLabel[CHEM_OPTION_ROWS];
  extern char *optionHeaderCell[CELL_OPTION_COLS];
  extern char *optionLabelCell[CELL_OPTION_ROWS];

  extern float surfColors[CHEM_OPTION_ROWS * CHEM_OPTION_COLS][4]; 
  extern float cellColors[CELL_OPTION_COLS][4];
  
  
  
  /*
   * 1[x, y, r, g, b]-------2[x, y, r, g, b]
   *        |                       |
   *        |                       |
   *        |                       |
   * 0[x, y, r, g, b]-------3[x, y, r, g, b]
   */
  extern GLint *WG_vert;
  extern GLint WG_xb;
  extern GLint WG_yb;
  extern GLint WG_zb;
  extern GLint WG_w;
  extern GLint WG_h;
  extern GLint WG_d;
  extern GLint WG_numvert;
  
  void init(int argc, char** argv, WHWorld *world_ptr);
  void start();

#ifdef OVERLAP_VIS
  extern bool visDone;

  void notifyVisDone();
  void notifyComDone();
  bool waitForVis();
  bool waitForCompute(); 
#endif

  extern bool simulationDone;
  bool isFinished();
  bool isPaused();

  void getTotalChem();
  void getTotalCell();
  void initWndGrid();
  void InitGL(int Width, int Height);
  
  void Draw2D();
  void Draw3D();
  
  void DrawGLScene();
  void ReSizeGLScene(int Width, int Height);

}


#endif
