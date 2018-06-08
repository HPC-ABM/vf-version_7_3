//
//  PlotManager.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/29/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#include "PlotManager.h"

PlotManager::PlotManager(float *total_a,
		int	        init_state,
		std::string title_a	[MAX_PLOTS],
		std::string xlabel_a[MAX_PLOTS],
		std::string ylabel_a[MAX_PLOTS]
		)
{

  this->total_a = total_a;

  // Allocation
  this->vertices = (float **) malloc(MAX_PLOTS * sizeof(float *));
  if (!this->vertices) {
    fprintf(stderr, "PlotManager: vertices allocation failed\n");
    exit(-1);
  }

  for (int i = 0; i < MAX_PLOTS; i++) {
    this->vertices[i] = (float *) malloc(MAX_PLOT_POINTS*2 * sizeof(float));
    if (!this->vertices[i]) {
      fprintf(stderr, "PlotManager: vertices[%d] allocation failed\n", i);
      exit(-1);
    }
  }
  
  // Initialization
  int x_a[MAX_PLOTS] = {
		  PLOT0x,
		  PLOT1x,
		  PLOT2x,
		  PLOT3x,
		  PLOT4x,
		  PLOT5x,
		  PLOT6x,
		  PLOT7x
  };

  int y_a[MAX_PLOTS] = {
		  PLOT0y,
		  PLOT1y,
		  PLOT2y,
		  PLOT3y,
		  PLOT4y,
		  PLOT5y,
		  PLOT6y,
		  PLOT7y
  };


  for(int i = 0; i < MAX_PLOTS; i++)
  {
	  plots[i]      = {x_a[i], y_a[i], PLOTw, PLOTh,
			  	  	  i, &total_a[i],
			  	  	  MIN_X, MAX_X, MIN_Y, MAX_Y,
			  	  	  MIN_X, MAX_X, MIN_Y, MAX_Y,
	                  init_state, title_a[i], xlabel_a[i], ylabel_a[i],
	                  vertices[i]};
  }
  
  plots[0].vmaxy = 0;


//  plots[0]      = {TNFx,      (int) TNFy,     CHM_PLOTw, (int) CHM_PLOTh, TOTAL_TNF_VBO, &total_a[tnf],
//                  MIN_X, MAX_X, MIN_Y, 0/*MAX_Y*/, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "TNF",       "Ticks", "Chem Level", vertices[tnf]};
//  plots[1]      = {TGFx,      (int) TGFy,     CHM_PLOTw, (int) CHM_PLOTh, TOTAL_TGF_VBO, &total_a[tgf],
//                  MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "TGF",       "Ticks", "Chem Level", vertices[tgf]};
//  plots[2]      = {FGFx,      (int) FGFy,     CHM_PLOTw, (int) CHM_PLOTh, TOTAL_FGF_VBO, &total_a[fgf],
//                  MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "FGF",       "Ticks", "Chem Level", vertices[fgf]};
//  plots[3]     = {MMP8x,     (int) MMP8y,    CHM_PLOTw, (int) CHM_PLOTh, TOTAL_MMP8_VBO, &total_a[mmp8],
//                  MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "MMP8",      "Ticks", "Chem Level", vertices[mmp8]};
//  plots[4]  = {IL1betax,  (int) IL1betay, CHM_PLOTw, (int) CHM_PLOTh, TOTAL_IL1beta_VBO, &total_a[il1beta],
//                  MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "IL1-beta",  "Ticks", "Chem Level", vertices[il1beta]};
//  plots[5]      = {IL6x,      (int) IL6y,     CHM_PLOTw, (int) CHM_PLOTh, TOTAL_IL6_VBO, &total_a[il6],
//                  MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "IL6",       "Ticks", "Chem Level", vertices[il6]};
//  plots[6]      = {IL8x,      (int) IL8y,     CHM_PLOTw, (int) CHM_PLOTh, TOTAL_IL8_VBO, &total_a[il8],
//                  MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                  1, "IL8",       "Ticks", "Chem Level", vertices[il8]};
//  plots[7]     = {IL10x,     (int) IL10y,    CHM_PLOTw, (int) CHM_PLOTh, TOTAL_IL10_VBO, &total_a[il10],
//                        MIN_X, MAX_X, MIN_Y, MAX_Y, MIN_X, MAX_X, MIN_Y, MAX_Y,
//                        1, "IL10",      "Ticks", "Chem Level", vertices[il10]};
  
  glGenBuffers(NUM_VBOs, vbo_ids);
  for (int i = 0; i < MAX_PLOTS; i++) {
    glBindBuffer(GL_ARRAY_BUFFER, vbo_ids[i]);
    glBufferData(GL_ARRAY_BUFFER, MAX_PLOT_POINTS*2*sizeof(float), NULL, GL_STREAM_DRAW);
  }
  
}



PlotManager::~PlotManager()
{
  // Free vertices
  for (int ic = 0; ic < MAX_PLOTS; ic++) {
    if (this->vertices[ic]) {
      free(this->vertices[ic]);
    }
  }
  
  if (this->vertices) free(this->vertices);
}


void initXverts(GLfloat *va)
{
  int tick = 0;
  for (int i = 0; i < MAX_PLOT_POINTS*2; i += 2) {
    va[i] = tick++;
  }
}

void PlotManager::updateXaxisVertsVals(GLfloat xlen,
                          GLfloat offset,
                          int ploth,
                          int l_margin,
                          int b_margin,
                          int axmin, int axmax)
{
  
  GLfloat xspace = xlen/N_INTV_X;
  // Update x-axis tick labels
  int increment = (axmax - axmin)/(N_INTV_X/*-1*/);
  for (int i = 0; i <= N_INTV_X; i++) {
    xaxis_values[i] = i * increment + axmin;
    xaxis_vert_x[i] = i * xspace    + offset + l_margin;
    xaxis_vert_y[i] = ploth - b_margin;
  }
  
}

void PlotManager::updateYaxisVertsVals(int l_margin, int b_margin, int ploth, int plotshift, int aymin, int aymax)
{
  // Update y-axis tick labels
  int increment = (aymax - aymin)/(N_INTV_Y-1);
  for (int i = 0; i < N_INTV_Y; i++) {
    int currtick = i * (MAX_PLOT_POINTS/4);
    yaxis_values[i] = i * increment + aymin;
    yaxis_vert_x[i] = l_margin;
    yaxis_vert_y[i] = (N_INTV_Y - 1 - i) * (ploth/4) + b_margin;
  }
}

/*----------------------------------------------------------------------------------------
 *	\brief	This function draws the specified 2D plot.
 *	\param	p	-	a pointer to the plot to draw.
 */
void PlotManager::DrawPlot(Plot *p, int iter, bool isVisible)
{
  int xfontx;
  int xfonty;
  int tfontx;
  int tfonty;
  
  int xworld = p->x;
  int yworld = p->y;
  int wworld = p->w;
  int hworld = p->h;
  
  int xlocal = 0;
  int ylocal = 0;//WINH;
  int wlocal = WINW;
  int hlocal = WINH;
  
  float axmin = p->aminx;
  float axmax = p->amaxx;
  float aymin = p->aminy;
  float aymax = p->aminy;
  
  const int l_border = 80;
  const int b_border = 140;
  const int l_margin = 10;
  const int b_margin = 20;
  const int ticksize = 10;

  GLfloat *data_vert = p->va;
  
  float chemLevel = *(p->tchemp);
  
  if(p)
  {

  	if (isVisible)
  	{

  		/*
  		 *	Draw an outline around the plot with width 3
  		 */
  		glLineWidth(1);

  		/*
  		 *	The colours for the outline     */
  		glColor3f(0.9f,0.9f,0.9f);

  		glBegin(GL_LINE_STRIP);
  		glVertex2i( xlocal+wlocal , ylocal      );
  		glVertex2i( xlocal        , ylocal      );
  		glVertex2i( xlocal        , ylocal+hlocal );
  		glEnd();


  		glBegin(GL_LINE_STRIP);
  		glVertex2i( xlocal        , ylocal+hlocal );
  		glVertex2i( xlocal+wlocal , ylocal+hlocal );
  		glVertex2i( xlocal+wlocal , ylocal      );
  		glEnd();

  		glLineWidth(1);

  		/*
  		 * Draw axis
  		 */

  		glColor3f(1,1,1);
  		glBegin(GL_LINE_STRIP);
  		glVertex2i( xlocal+l_border        , ylocal      + b_border);
  		glVertex2i( xlocal+l_border        , ylocal+hlocal - b_border );
  		glVertex2i( xlocal+wlocal -l_border/2 , ylocal+hlocal - b_border );
  		glEnd();


  		/*
  		 *	Calculate the x and y coords for the text string in order to center it.
  		 */
  		xfontx = xlocal + (wlocal - glutBitmapLength(GLUT_BITMAP_HELVETICA_10, (unsigned char *) (p->xlabel.c_str()))) / 2 ;
  		xfonty = ylocal + hlocal - 30;

  		tfontx = xlocal + (wlocal - glutBitmapLength(GLUT_BITMAP_HELVETICA_10, (unsigned char *) (p->title.c_str()))) / 2 ;
  		tfonty = ylocal + hlocal - (hlocal - 60);

  		glColor3f(1,1,1);
  		Font(GLUT_BITMAP_HELVETICA_10,(char *) (p->title).c_str(),  tfontx, tfonty);
  		Font(GLUT_BITMAP_HELVETICA_10,(char *) (p->xlabel).c_str(), xfontx, xfonty);

  	}


    /*
     * Plot data
     */
    
    if (data_vert) {
      
      if (iter == 0) {
        
        initXverts(data_vert);
        updateXaxisVertsVals(wlocal - l_border*2,   // xlen
                             -10,//0,                     // offset
                             hlocal,                // ploth
                             l_border,              // l_margin
                             b_margin*4,              // b_margin
                             axmin,
                             axmax);
        updateYaxisVertsVals(l_margin,
                             b_border,
                             hlocal,                // ploth
                             hlocal,                // plotshift
                             aymin,
                             aymax);
        
      }  else {
        if (p->vminx != p->aminx) {
          p->aminx = p->vminx;
        }
        if (p->vmaxx > axmax) {
          p->amaxx = (p->vmaxx);// * 2;
        }

	    if (!paused) {
          updateXaxisVertsVals(wlocal - l_border*2,   // xlen
                               -10,//0,                     // offset
                               hlocal,                // ploth
                               l_border,              // l_margin
                               b_margin*4,              // b_margin
                               p->aminx,
                               p->amaxx);
        }

        if (p->vminy < aymin) {
          
          if (p->vminy < 0)
            p->aminy = (p->vminy) * 2;
          else
            p->aminy = (p->vminy) - b_margin;
          
        } else if (p->vmaxy > p->amaxy) {
          p->amaxy = (p->vmaxy) + MAX_Y * 20;//* 10;
        }
        
        updateYaxisVertsVals(l_margin,
                             b_border,
                             hlocal,                // ploth
                             hlocal,                // plotshift
                             p->aminy,
                             p->amaxy);
      }

      if (isVisible)
      {
      	// Render x-axis tick values
      	char tick_str[5];
      	for (int i = 0; i <= N_INTV_X; i++) {\
      		sprintf(tick_str, "%d", xaxis_values[i]);
      	Font(GLUT_BITMAP_HELVETICA_10,
      			tick_str,
      			xaxis_vert_x[i],
      			xaxis_vert_y[i]);
      	}
      	// Render y-axis
      	char y_str[15];
      	for (int i = 0; i < N_INTV_Y; i++) {
      		sprintf(y_str, "%d", yaxis_values[i]);
      		Font(GLUT_BITMAP_HELVETICA_10,
      				y_str,
      				yaxis_vert_x[i],
      				yaxis_vert_y[i]);
      	}
      }
      
      // Update y component for this tick
      int vind = iter%MAX_PLOT_POINTS;
      int xind = vind * 2;
      int yind = xind + 1;
      GLfloat xspace = (wlocal - l_border*2)/MAX_PLOT_POINTS;

      if (paused) {
        vind = (iter-1)%MAX_PLOT_POINTS;
        xind = vind * 2;
        yind = xind + 1;
      } 

     // UPDATE Y VALUE
      if (!paused) {
        int ran = rand()%100;
        if (ran < 30){
          chemLevel = (100.0f-(ran/2.0f))/100.0f * chemLevel;
        } else if (ran > 85){
          chemLevel = ((100.0f-ran)/100.0f + 1) * chemLevel;
        }
        data_vert[yind] = chemLevel;
      
        // Update y -- max and min
        if (p->vminy > chemLevel) {
          p->vminy = chemLevel;
        }
        if (p->vmaxy < chemLevel) {
          p->vmaxy = chemLevel;
        }
      
        if (iter >= MAX_PLOT_POINTS) {
          // UPDATE X VALUE
          data_vert[xind] = iter;
          p->vmaxx = iter;
          p->vminx++;
        }
      }
      
      if (isVisible)
      {
      	glMatrixMode(GL_PROJECTION);
      	glPushMatrix();


      	glViewport(
      			xworld + l_margin + ticksize,
      			yworld + b_margin + ticksize,
      			wworld - l_margin * 2 - ticksize,
      			hworld - b_margin * 2 - ticksize
      	);
      	glLoadIdentity();
      	gluOrtho2D(p->aminx,              // left
      			p->amaxx,              // right
      			p->aminy,              // bottom
      			p->amaxy               // top
      	);

      	glBindBuffer(GL_ARRAY_BUFFER, vbo_ids[p->VBO_id]);
      	glBufferSubData(GL_ARRAY_BUFFER, xind * sizeof(GLfloat), 2 * sizeof(GLfloat), &data_vert[xind]);
      	glEnableClientState(GL_VERTEX_ARRAY);
      	glVertexPointer(2, GL_FLOAT, 0, 0);
      	if (iter >= MAX_PLOT_POINTS) {
      		glDrawArrays(GL_LINE_STRIP, vind+1, MAX_PLOT_POINTS - (vind+1));
      		glDrawArrays(GL_LINE_STRIP, 0, vind+1);
      		//        printf("Plotted: %d + %d = %d points\n", MAX_PLOT_POINTS - (vind+1), vind+1, vind+1 + MAX_PLOT_POINTS - (vind+1));
      	} else {
      		glDrawArrays(GL_LINE_STRIP, 0, iter);
      	}
      	glDisableClientState(GL_VERTEX_ARRAY);
      	glBindBuffer(GL_ARRAY_BUFFER, 0);

      	glMatrixMode(GL_PROJECTION);
      	glPopMatrix();
      }
      
    } else {
    	if (isVisible)
    	{
    		int tbifontx = xlocal - 60 + (wlocal - glutBitmapLength(GLUT_BITMAP_HELVETICA_18, (unsigned char *) "-- TBI --")) / 2 ;
    		int tbifonty = ylocal + hlocal/2;
    		Font(GLUT_BITMAP_HELVETICA_18, "-- TBI --", tbifontx, tbifonty);
    	}
    }
    


  }
}

void PlotManager::DrawChemViewPorts(int iter, bool isVisible)
{
  // Update total chem
  // updateTotalChem();
  // Moved to main
  
	for (int i = 0; i < MAX_PLOTS; i++)
	{
		  glViewport(plots[i].x,
				  plots[i].y,
				  plots[i].w,
				  plots[i].h);
		  DrawPlot(&plots[i], iter, isVisible);
	}
  
}


