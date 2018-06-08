//
//  HeatMap.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/29/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_HeatMap_h
#define OpenGL_SimpleWH_HeatMap_h

#include "common_vis.h"
#include "draw_util.h"

/*----------------------------------------------------------------------------------------
 *	Heat Map Stuff
 */
struct HeatMapLegend
{
  int   x;							/* top left x coord of the plot */
  int   y;							/* top left y coord of the plot */
  int   w;							/* the width of the plot */
  int   h;							/* the height of the plot */
  float max;                                                    /* Max value on the heatmap */
  float min;                                                    /* Min value on the heatmap */
  int   init;                                                   /* is initialized */
  int   buff;                                                   /* is buffered */
  GLfloat *data_vert;
  int	state;                                                  /* the state, 1 if show, 0 otherwise */
  char* title;                                                  /* the text title label of the plot */
};
typedef struct HeatMapLegend HeatMapLegend;

class HeatMapManager
{
public:
  HeatMapManager(int nx, int ny, int nz);
  virtual ~HeatMapManager();
  
  void initHMLV(HeatMapLegend *h);
  void initHML(HeatMapLegend *h);
  void DrawLegend();
  void GenerateRegionPointBuffers(int lpx, int lpy, int lpz);
  void Render(float* chemLevels, int gridh, int gridw, int gridd, float xoffset, float yoffset);
  
private:
  int R1xb, R1xe, R1yb, R1ye, R1zb, R1ze;
  int R2xb, R2xe, R2yb, R2ye, R2zb, R2ze;
  int incx1, incy1, incz1;
  int incx2, incy2, incz2;
  int incx3, incy3, incz3;
  bool isInR1(int, int , int);
  bool isInR2(int, int , int);
  
  void DrawLegend(HeatMapLegend *h);
  void bufferColors();
  
  GLfloat HML_vert[HML_RESOLUTION * 2 * (2+3)];   // [x, y, r, g, b]-----[x, y, r, g, b]
  GLfloat HML_values[N_HMLV];
  GLfloat HML_vert_xp[N_HMLV];
  GLfloat HML_vert_yp[N_HMLV];
  
  unsigned int vbo_hml;
  
  HeatMapLegend HML;
  
  // ID's for the VBO's
  GLuint R1vbID;
  GLuint R2vbID;
  GLuint R3vbID;

  GLuint R1cbID;
  GLuint R2cbID;
  GLuint R3cbID;

  std::vector<glm::vec3> R1points;
  std::vector<glm::vec3> R2points;
  std::vector<glm::vec3> R3points;

  std::vector<QuadV> R1vertices;	// rectangular prism
  std::vector<QuadV> R2vertices;	// rectangular prism
  std::vector<QuadV> R3vertices;	// rectangular prism

  std::vector<glm::vec4> R1colors;
  std::vector<glm::vec4> R2colors;
  std::vector<glm::vec4> R3colors;

};

#endif
