//
//  PlotManager.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/29/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_PlotManager_h
#define OpenGL_SimpleWH_PlotManager_h

#include "common_vis.h"
#include "draw_util.h"
#include "UI.h"

/*----------------------------------------------------------------------------------------
 *	2D Plot Stuff
 */

struct Plot
{
  int   x;							/* top left x coord of the plot */
  int   y;							/* top left y coord of the plot */
  int   w;							/* the width of the plot */
  int   h;							/* the height of the plot */
  int   VBO_id;
  float *tchemp;                                                /* total chem pointer */
  float vminx;
  float vmaxx;                                                  /* current max value of x */
  float vminy;
  float vmaxy;
  float aminx;
  float amaxx;
  float aminy;
  float amaxy;
  int	state;                                                  /* the state, 1 if show, 0 otherwise */
  std::string title;                                                  /* the text title label of the plot */
  std::string xlabel;                                                 /* the text x-axis label of the plot */
  std::string ylabel;                                                 /* the text y-axis label of the plot */
  GLfloat *va;                                                  /* pointer to vertex array containing data to plot*/
};
typedef struct Plot Plot;

class PlotManager
{
public:
  PlotManager(float *total_a,
			int			init_state,
			std::string title_a	[MAX_PLOTS],
			std::string xlabel_a[MAX_PLOTS],
			std::string ylabel_a[MAX_PLOTS]
			);
  virtual ~PlotManager();
  
  void Render();
  void DrawPlot(Plot *p, int iter, bool isVisible);
  
  void updateXaxisVertsVals(GLfloat xlen,
                            GLfloat offset,
                            int ploth,
                            int l_margin,
                            int b_margin,
                            int axmin, int axmax);
  void updateYaxisVertsVals(int l_margin, int b_margin, int ploth, int plotshift, int aymin, int aymax);
  void DrawChemViewPorts(int iter, bool isVisible);
  
private:
  float *total_a;       // [numChem]    -- Passed in to constructor
  GLfloat **vertices;     // [numChem][MAX_PLOT_POINTS*2]   -- Allocate in constructor
  
  GLint   xaxis_values[N_INTV_X+1];
  GLfloat xaxis_vert_x[N_INTV_X+1];
  GLfloat xaxis_vert_y[N_INTV_X+1];
  
  GLint   yaxis_values[N_INTV_Y];
  GLfloat yaxis_vert_x[N_INTV_Y];
  GLfloat yaxis_vert_y[N_INTV_Y];
  
  GLuint vbo_ids[NUM_VBOs];
  
  Plot plots[MAX_PLOTS];
  
};


#endif
