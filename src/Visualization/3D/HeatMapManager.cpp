//
//  HeatMap.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/29/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#include "HeatMapManager.h"

#define ADAPTIVE_SAMPLNG
#define SHOW_SAMPLES

HeatMapManager::HeatMapManager(int nx, int ny, int nz){
#ifdef HUMAN_VF
        float cmax = 10.0f;//0.001f;//0.0005f;//1.0f;//0.5f;//0.035f;//0.3f;
#elif defined(RAT_VF)
        float cmax = 0.0000001f;//0.01f;//
#else
        float cmax = 1.0f;
#endif

#ifdef VERT_PLOT
  HML = {WINW - 85, 455, 40, 200, cmax, 0.0, 0, 0, HML_vert, 1, "Chemical Concentration"};
#else
  HML = {50, 50, 40, 200, cmax, 0.0, 0, 0, HML_vert, 1, "Chemical Concentration"};
#endif
  initHML(&HML);
  initHMLV(&HML);
  glGenBuffers(1, &vbo_hml);
  glBindBuffer(GL_ARRAY_BUFFER, vbo_hml);
  glBufferData(GL_ARRAY_BUFFER, HML_RESOLUTION*2*(2+3)*sizeof(float), HML_vert, GL_STATIC_DRAW);  // [x, y, r, g, b]-----[x, y, r, g, b]
  HML.buff = 1;

//#ifdef ADAPTIVE_SAMPLING
  this->GenerateRegionPointBuffers(nx, ny, nz);
//#endif

}

HeatMapManager::~HeatMapManager()
{
  // TODO: CleanUp
}

void HeatMapManager::initHMLV(HeatMapLegend *h){
  float incv = (h->max - h->min)/(N_HMLV-1);
  float ht   = (float)(h->h);
  float incy = ht/((float) (N_HMLV-1));
  printf("h: %f, incy: %f\n", ht, incy);
  for (int i = 0; i < N_HMLV; i++) {
    HML_values[i] = i * incv + h->min;
    HML_vert_xp[i] = h->x - 30;
    HML_vert_yp[i] = (N_HMLV - 1 - i) * incy + h->y + 4;
  }
}

void HeatMapManager::initHML(HeatMapLegend *h){
  int x  = h->x;
  int y  = h->y;
  int wd = h->w;
  int ht = h->h;
  int vert_ind = 0;
  
  float A = h->min;
  float B = h->max;
  float hmv = B;
  float hmy = y;
  float vinc = (B-A)/HML_RESOLUTION;
  float yinc = ht/HML_RESOLUTION;
  
  float *data_vert = h->data_vert;
  
  h->init = 1;
  
  for (int i = 0; i < HML_RESOLUTION; i++) {
    // Update line i
    GLfloat ratio = 2 * (hmv-A) / (B - A);
    GLfloat vB = 0 > 1 - ratio? 0 : 1 - ratio;
    GLfloat vR = 0 > ratio - 1? 0 : ratio - 1;
    GLfloat vG = 1 - vB - vR;
    
    // x1
    data_vert[vert_ind++] = x;
    // y1
    data_vert[vert_ind++] = hmy;
    // r1
    data_vert[vert_ind++] = vR;
    // g1
    data_vert[vert_ind++] = vG;
    // b1
    data_vert[vert_ind++] = vB;
    
    // x2
    data_vert[vert_ind++] = x + wd;
    // y2
    data_vert[vert_ind++] = hmy;
    // r2
    data_vert[vert_ind++] = vR;
    // g2
    data_vert[vert_ind++] = vG;
    // b2
    data_vert[vert_ind++] = vB;
    
    hmv -= vinc;
    hmy += yinc;
    
    //    printf("[%f, %f]----[%f, %f]\t\t[%f, %f, %f]\n", data_vert[i*10 + 0], data_vert[i*10 + 1],
    //           data_vert[i*10 + 5], data_vert[i*10 + 6], data_vert[i*10 + 2], data_vert[i*10 + 3], data_vert[i*10 + 4]);
  }
}

//#ifdef ADAPTIVE_SAMPLING

bool HeatMapManager::isInR1(int x, int y, int z)
{
  return (R1xb < x && x < R1xe && R1yb < y && y < R1ye && R1zb < z && z < R1ze);
}

bool HeatMapManager::isInR2(int x, int y, int z)
{
  return (R2xb < x && x < R2xe && R2yb < y && y < R2ye && R2zb < z && z < R2ze);
}

// TODO: Make regions with respect to wound
void HeatMapManager::GenerateRegionPointBuffers(int lpx, int lpy, int lpz)
{
  // Region 1
  this->R1xb = (lpx/3) * 2;
  this->R1xe =  lpx;
  this->R1yb = (lpy/3) * 1;
  this->R1ye = (lpy/3) * 2;
  this->R1zb = (lpz/3) * 1;
  this->R1ze = (lpz/3) * 2;

  int const_step =10;

//  incx1 = const_step; incy1 = const_step; incz1 = const_step;
  incx1 = 2; incy1 = 2; incz1 = 2;

  for (int z = R1zb; z < R1ze; z += incz1) {
    for (int y = R1yb; y < R1ye; y += incy1) {
      for (int x = R1xb; x < R1xe; x += incx1) {
        this->R1points.push_back(glm::vec3(x, y, z));
//        this->R1vertices.push_back(getQuadVertices(GL_X(x), GL_Y(y), GL_Z(z)));
      }
    }
  }

  // Region 2
  this->R2xb = (lpx/3) * 1;
  this->R2xe =  lpx;
  this->R2yb = (lpy/6) * 1;
  this->R2ye = (lpy/6) * 5;
  this->R2zb = (lpz/6) * 1;
  this->R2ze = (lpz/6) * 5;

//  incx2 = const_step; incy2 = const_step; incz2 = const_step;
  incx2 = 6; incy2 = 6; incz2 = 6;

  for (int z = R2zb; z < R2ze; z += incz2) {
    for (int y = R2yb; y < R2ye; y += incy2) {
      for (int x = R2xb; x < R2xe; x += incx2) {
        if (!isInR1(x, y, z))
          this->R2points.push_back(glm::vec3(x, y, z));
      }
    }
  }

  // Region 3
//  incx3 = const_step; incy3 = const_step; incz3 = const_step;
  incx3 = 10; incy3 = 10; incz3 = 10;
  for (int z = 0; z < lpz; z += incz3) {
    for (int y = 0; y < lpy; y += incy3) {
      for (int x = 0; x < lpx; x += incx3) {
        if (!isInR2(x, y, z))
          this->R3points.push_back(glm::vec3(x, y, z));
      }
    }
  }

}

void HeatMapManager::bufferColors()
{

}

//#endif

void HeatMapManager::Render(float *chemLevels, int gridh, int gridw, int gridd, float xoffset, float yoffset)
{
        glColor4f(1.0f, 1.0f, 1.0f, 1.0f);//0.35f);
        glDisable(GL_DEPTH_TEST);         // Turn Depth Testing Off
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glEnable(GL_BLEND);

        float hmv;
        GLfloat A = this->HML.min;
        GLfloat B = this->HML.max;
//#ifdef HUMAN_VF
//        float B = 10.0f;//0.001f;//0.0005f;//1.0f;//0.5f;//0.035f;//0.3f;
//#else
//        float B = 0.1f;
//#endif
//        float A = 0.0f;

/*#ifdef ADAPTIVE_SAMPLING
#ifdef SHOW_SAMPLES
*/
/*
  glEnable(GL_DEPTH_TEST); // DEBUG
  glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
  glPointSize(3.0);
  glBegin(GL_POINTS);
   // Region 1
  for(std::vector<glm::vec3>::iterator it = R1points.begin(); it != R1points.end(); ++it) {
    int x = (*it)[0];
    int y = (*it)[1];
    int z = (*it)[2];

    glVertex3i(GL_X(x), GL_Y(y), PLANE_DEPTH + z);
  }
  
//  glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
  // Region 2
  for(std::vector<glm::vec3>::iterator it = R2points.begin(); it != R2points.end(); ++it) {
    int x = (*it)[0];
    int y = (*it)[1];
    int z = (*it)[2];

    glVertex3i(GL_X(x), GL_Y(y), PLANE_DEPTH + z);
  }
  
//  glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
  // Region 3
  for(std::vector<glm::vec3>::iterator it = R3points.begin(); it != R3points.end(); ++it) {
    int x = (*it)[0];
    int y = (*it)[1];
    int z = (*it)[2];

    glVertex3i(GL_X(x), GL_Y(y), PLANE_DEPTH + z);
  }
  glEnd();
*/
/*
// RENDER
  glPointSize(10.0f);
  glBegin(GL_POINTS);
  // Region 3
  for(std::vector<glm::vec3>::iterator it = R3points.begin(); it != R3points.end(); ++it) {
    int x = (*it)[0];
    int y = (*it)[1];
    int z = (*it)[2];
    hmv = chemLevels[(z * gridw * gridh) + (y * gridw) + x];
    GLfloat ratio = 2 * (hmv-A) / (B - A);
    GLfloat vB = 0 > 1 - ratio? 0 : 1 - ratio;
    GLfloat vR = 0 > ratio - 1? 0 : ratio - 1;
    GLfloat vG = 1 - vB - vR;

    glColor4f(vR, vG, vB-0.2, (1-0.99*vB));
//    glVertex3i(GL_X(x), GL_Y(y), GL_Z(z));
    drawSolidCube(GL_X(x), GL_Y(y), GL_Z(z), incx3);
  }
  glEnd(); 
 
  glPointSize(0.01f);
  glBegin(GL_POINTS);
  // Region 2
  for(std::vector<glm::vec3>::iterator it = R2points.begin(); it != R2points.end(); ++it) {
    int x = (*it)[0];
    int y = (*it)[1];
    int z = (*it)[2];
    hmv = chemLevels[(z * gridw * gridh) + (y * gridw) + x];
    GLfloat ratio = 2 * (hmv-A) / (B - A);
    GLfloat vB = 0 > 1 - ratio? 0 : 1 - ratio;
    GLfloat vR = 0 > ratio - 1? 0 : ratio - 1;
    GLfloat vG = 1 - vB - vR;

    glColor4f(vR, vG, vB-0.2, (1-0.995*vB));
//  glColor3f(1.0, 0.0, 0.0);
//    glVertex3i(GL_X(x), GL_Y(y), GL_Z(z));
    drawSolidCube(GL_X(x), GL_Y(y), GL_Z(z), incx2);
  }
 
  // Region 1
  for(std::vector<glm::vec3>::iterator it = R1points.begin(); it != R1points.end(); ++it) {
    int x = (*it)[0];
    int y = (*it)[1];
    int z = (*it)[2];
    hmv = chemLevels[(z * gridw * gridh) + (y * gridw) + x];
    GLfloat ratio = 2 * (hmv-A) / (B - A);
    GLfloat vB = 0 > 1 - ratio? 0 : 1 - ratio;
    GLfloat vR = 0 > ratio - 1? 0 : ratio - 1;
    GLfloat vG = 1 - vB - vR;

    glColor4f(vR, vG, vB-0.2, (1-0.998*vB));
//  glColor3f(0.0, 1.0, 0.0);
//    glVertex3i(GL_X(x), GL_Y(y), GL_Z(z));
    drawSolidCube(GL_X(x), GL_Y(y), GL_Z(z), incx1);
  }  
*/
//  glEnd();





/*#endif  // SHOW_SAMPLES
#else	// ADAPTIVE_SAMPLING
*/
	// Note: stride = 
	// 1:   ~2263    ms/iter to render
	// 5:   ~300-400 ms/iter
	// 10:  ~60-120  ms/iter

        int stride = 2;
//	glPointSize(0.01f);
	glBegin(GL_POINTS);
        for (int z = 0; z < gridd; z += stride) {
                for (int y = 0; y < gridh; y += stride) {
                        for (int x = 0; x < gridw; x += stride) {
               			hmv = chemLevels[(z * gridw * gridh) + (y * gridw) + x];
                               // if (hmv > 0.00002) {
                                        GLfloat ratio = 2 * (hmv-A) / (B - A);
                                        GLfloat vB = 0 > 1 - ratio? 0 : 1 - ratio;
                                        GLfloat vR = 0 > ratio - 1? 0 : ratio - 1;
                                        GLfloat vG = 1 - vB - vR;

                                        glColor4f(vR, vG, vB-0.2, (1-(0.9976*vB + 0.83*vG)));
//                                        glVertex3i(GL_X(x), GL_Y(y), GL_Z(z));
                                        drawSolidCube(GL_X(x) + xoffset,
                                                        GL_Y(y) + yoffset,
                                                        PLANE_DEPTH+1+z,
                                                        stride);

                               // }
                        }
                }
        }
//	glEnd();

/*        for (int z = 0; z < gridd; z += stride) {
                for (int y = 0; y < gridh; y += stride) {
                        for (int x = 0; x < gridw; x += stride) {
               			hmv = chemLevels[(z * gridw * gridh) + (y * gridw) + x];
                                {
                                        GLfloat ratio = 2 * (hmv-A) / (B - A);
                                        GLfloat vB = 0 > 1 - ratio? 0 : 1 - ratio;
                                        GLfloat vR = 0 > ratio - 1? 0 : ratio - 1;
                                        GLfloat vG = 1 - vB - vR;

                                        glColor4f(vR, vG, vB-0.2, (1-0.99*vB));
                                        glVertex3i(GL_X(x), GL_Y(y), GL_Z(z));
                                        //drawSolidCube(GL_X(x),
                                        //                GL_Y(y),
                                        //                PLANE_DEPTH+1+z,
                                        //                stride);

                                }
                        }
                }
        }
*/
/*
	glEnable(GL_DEPTH_TEST);



  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_DEPTH_TEST);
  // TODO: Fix max/min order
  GLfloat A = this->HML.max, B = this->HML.min;
//  GLfloat A = 0.0f;
//  GLfloat B = 0.5f; 

  GLfloat chemLevel;
  //  glPointSize( 0.1f );
  glPointSize( 0.01f );
  glBegin(GL_POINTS);
  for (int k = 0; k < gridd; k++)
  for (int i = 0; i < gridh; i++) {
    for (int j = 0; j < gridw; j++) {
      chemLevel = chemLevels[(k * gridw * gridh) + (i * gridw) + j];
      //     glColor3f(chemLevel * 0.5, chemLevel * 0.5, chemLevel);
      if (chemLevel > 0.0002) {//0.2) {
        GLfloat ratio = 2 * (chemLevel-A) / (B - A);
        GLfloat chemB = 0 > 1 - ratio? 0 : 1 - ratio;
        GLfloat chemR = 0 > ratio - 1? 0 : ratio - 1;
        GLfloat chemG = 1 - chemB - chemR;
        
        //              GLfloat chemR = (chemLevel - A) * (br - ar)/(B - A) + ar;
        //              GLfloat chemG = 1 - chemR;
        //              GLfloat chemB = (chemLevel - A) * (bb - ab)/(B - A) + ab;
        glColor4f(chemR, chemG, chemB*0.7, 0.1);
        
        //      glColor4f(0.5, 0.0f, 0.7, chemLevel);
        //      glColor4f(0.3, 0.1f, 0.9, chemLevel);
        //      glVertex3i(GL_X(j), GL_Y(i), 0);
        glVertex3i(GL_X(j), GL_Y(i), PLANE_DEPTH+k);
      }
      
    }
  }
  glEnd();
  
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST); // Re-enable writing to depth buffer

#endif	// ADAPTIVE_SAMPLING
*/
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST); // Re-enable writing to depth buffer

}

// Entry point function to draw heatmap legend
void HeatMapManager::DrawLegend()
{
  DrawLegend(&HML);
}
/*----------------------------------------------------------------------------------------
 *	\brief	This function draws the specified heatmap legend.
 *	\param	h	-	a pointer to the heatmap to draw.
 */
void HeatMapManager::DrawLegend(HeatMapLegend *h)
{
  int tfontx;
  int tfonty;
  int minfontx;
  int minfonty;
  
  if(h)
  {
    if (h->init && h->buff) {
      
      glLineWidth(2);
      
      glBindBuffer(GL_ARRAY_BUFFER, vbo_hml);
      glEnableClientState(GL_VERTEX_ARRAY);
      glEnableClientState(GL_COLOR_ARRAY);
      
      glVertexPointer(2, GL_FLOAT, 5*sizeof(float), 0);
      glColorPointer(3, GL_FLOAT, 5*sizeof(float), (void*) (2*sizeof(float)));
      
      glDrawArrays(GL_LINES, 0, HML_RESOLUTION*2);
      
      glDisableClientState(GL_COLOR_ARRAY);
      glDisableClientState(GL_VERTEX_ARRAY);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
      
      /*
       *  Calculate the x and y coords for the text string in order to center it.
       */
      //      xfontx = h->x + (p->w - glutBitmapLength(GLUT_BITMAP_HELVETICA_10, (unsigned char *) p->xlabel)) / 2 ;
      //      xfonty = h->y + p->h - 5;
      
      tfontx = h->x + (h->w - glutBitmapLength(GLUT_BITMAP_HELVETICA_10, (unsigned char *) h->title)) / 2 ;
      tfonty = h->y + h->h - (h->h + 10);
      
      glColor3f(1,1,1);
      Font(GLUT_BITMAP_HELVETICA_10,h->title, tfontx,tfonty);
      //      Font(GLUT_BITMAP_HELVETICA_10,p->xlabel,xfontx,xfonty);
      
      char v_str[10];
      for (int i = 0; i < N_HMLV; i++) {
        sprintf(v_str, "%.3f", HML_values[i]);
        Font(GLUT_BITMAP_HELVETICA_10,
             v_str,
             HML_vert_xp[i],
             HML_vert_yp[i]);
      }
      
    }
  }
}
