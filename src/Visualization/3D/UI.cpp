//
//  UI.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/28/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#include "UI.h"

// actual vector representing the camera's direction
float lx=0.0f, ly=1.0f, lz=-1.0f;
// XZ position of the camera
float camX=0.0f, camY=0.0f, camZ=0.0f;
float pcamX=0.0f, pcamY=0.0f, pcamZ=0.0f;

GLfloat cellSize = 0.001f;//0.25f;//2.0f;
GLfloat pcellSize = 0.001f;//0.25f;//2.0f;       // previous cell size

bool paused = false;

// Chem global variables
bool showChem_HM[8]   = {true, false, false, false, false, false, false, false};
bool showChem_Surf[8] = {false, false, false, false, false, false, false, false};
bool showCells      = false;
bool showNeus       = true;
bool showMacs       = true;
bool showFibs       = true;
bool showECMs       = false;
bool showChemCharts = true;
bool showCellCharts = false;

bool showChemOp     = false;

bool lowerECMplane  = false;
bool raiseECMplane  = false;

bool zoomedWound    = false;

#ifdef USE_MOUSE_ONLY

/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y)
{
  /* avoid thrashing this procedure */
  usleep(100);
  
  /* If escape is pressed, kill everything. */
  if (key == ESCAPE)
  {
    /* shut down our window */
    glutDestroyWindow(window);
    
    /* exit the program...normal termination. */
    exit(0);
  }
  
}

#else   // USE_MOUSE_ONLY


/* The function called whenever a key is pressed. */
void keyPressed(unsigned char key, int x, int y)
{
  /* avoid thrashing this procedure */
  usleep(100);
  
  /* If escape is pressed, kill everything. */
  switch(key)
  {
    case ESCAPE:
      /* shut down our window */
      glutDestroyWindow(window);
    
      /* exit the program...normal termination. */
      exit(0);
      break;
    case '+':
      ZoomInButtonCallback();
      break;
    case '_':
      ZoomOutButtonCallback();
      break;
    case 'w':
      ZoomWoundButtonCallBack();
      break;
    case 'u':
      raiseECMplane = (raiseECMplane == true)? false : true;
      break;
    case 'p':
      paused = (paused == true)? false : true;
      break;
    case 'c':
      if (showChemCharts)	{showChemCharts = false; showCellCharts = true;}
      else 					{showChemCharts = true;  showCellCharts = false;}
      break;
  }
  
}

void SpecialInput(int key, int x, int y)
{
  switch(key)
  {
    case GLUT_KEY_UP:
      MoveUpButtonCallBack();
      break;
    case GLUT_KEY_DOWN:
      MoveDownButtonCallBack();
      break;
    case GLUT_KEY_LEFT:
      MoveLeftButtonCallBack();
      break;
    case GLUT_KEY_RIGHT:
      MoveRightButtonCallBack();
      break;
  }
}

void ChemOpButtonCallBack() {
  showChemOp = !showChemOp;
}

#endif  // USE_MOUSE_ONLY


/*----------------------------------------------------------------------------------------
 *	This is an example callback function. Notice that it's type is the same
 *	an the ButtonCallback type. We can assign a pointer to this function which
 *	we can store and later call.
 */
void ZoomInButtonCallback()
{
  GLfloat step = 50.0f;
  if (camZ <= -150) step = 10.0f;
  camZ -= step;
  //camZ -= 50.0f;
  if (camZ < -200) cellSize += 0.2;//0.25;
  
  printf("Zooming in: %f!\n", camZ);
}

void ZoomOutButtonCallback()
{
  printf("Zooming out!\n");
  GLfloat step = 50.0f;
  if (camZ <= -150) step = 10.0f;
  camZ += step;
  if (camZ <= -200) cellSize -= 0.2;//0.25;
}

void MoveLeftButtonCallBack()
{
  GLfloat step = 50.0f;
  if (camZ < -250) step = 10.0f;
  else if (camZ < -500) step = 5.0f;
  else if (camZ < -750) step = 2.0f;
  camX -= step;//50.0f;
  printf("Moving left %f!\n", step);
}

void MoveRightButtonCallBack()
{
  GLfloat step = 50.0f;
  if (camZ < -250) step = 10.0f;
  else if (camZ < -500) step = 5.0f;
  else if (camZ < -750) step = 2.0f;
  camX += step;//50.0f;
  printf("Moving right %f!\n", step);
}

void MoveUpButtonCallBack()
{
  GLfloat step = 50.0f;
  if (camZ < -250) step = 10.0f;
  else if (camZ < -500) step = 5.0f;
  else if (camZ < -750) step = 2.0f;
  camY += step;//50.0f;
  printf("Moving up %f!\n", step);
}

void MoveDownButtonCallBack()
{
  
  GLfloat step = 50.0f;
  if (camZ < -250) step = 10.0f;
  else if (camZ < -500) step = 5.0f;
  else if (camZ < -750) step = 2.0f;
  camY -= step;//50.0f;
  printf("Moving down %f!\n", step);
}

void setAllChemFalse()
{
	for (int ic = 0; ic < 8; ic++)
		showChem_HM[ic] = false;
}

void TNF_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[TNF] = showChem_HM[TNF]? false : true;
  printf("Toggle TNF HM: %d\n", showChem_HM[TNF]);
}

void TNF_SurfButtonCallBack()
{
  showChem_Surf[TNF] = showChem_Surf[TNF]? false : true;
  printf("Toggle TNF Surf: %d\n", showChem_Surf[TNF]);
}

void TGF_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[TGF] = showChem_HM[TGF]? false : true;
  printf("Toggle TGF HM: %d\n", showChem_HM[TGF]);
}

void TGF_SurfButtonCallBack()
{
  showChem_Surf[TGF] = showChem_Surf[TGF]? false : true;
  printf("Toggle TGF Surf: %d\n", showChem_Surf[TGF]);
}

void FGF_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[FGF] = showChem_HM[FGF]? false : true;
  printf("Toggle FGF HM: %d\n", showChem_HM[FGF]);
}

void FGF_SurfButtonCallBack()
{
  showChem_Surf[FGF] = showChem_Surf[FGF]? false : true;
  printf("Toggle FGF Surf: %d\n", showChem_Surf[FGF]);
}

void MMP8_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[MMP8] = showChem_HM[MMP8]? false : true;
  printf("Toggle MMP8 HM: %d\n", showChem_HM[MMP8]);
}

void MMP8_SurfButtonCallBack()
{
  showChem_Surf[MMP8] = showChem_Surf[MMP8]? false : true;
  printf("Toggle MMP8 Surf: %d\n", showChem_Surf[MMP8]);
}
void IL1beta_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[IL1beta] = showChem_HM[IL1beta]? false : true;
  printf("Toggle IL1beta HM: %d\n", showChem_HM[IL1beta]);
}

void IL1beta_SurfButtonCallBack()
{
  showChem_Surf[IL1beta] = showChem_Surf[IL1beta]? false : true;
  printf("Toggle IL1beta Surf: %d\n", showChem_Surf[IL1beta]);
}

void IL6_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[IL6] = showChem_HM[IL6]? false : true;
  printf("Toggle IL6 HM: %d\n", showChem_HM[IL6]);
}

void IL6_SurfButtonCallBack()
{
  showChem_Surf[IL6] = showChem_Surf[IL6]? false : true;
  printf("Toggle IL6 Surf: %d\n", showChem_Surf[IL6]);
}
void IL8_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[IL8] = showChem_HM[IL8]? false : true;
  printf("Toggle IL8 HM: %d\n", showChem_HM[IL8]);
}

void IL8_SurfButtonCallBack()
{
  showChem_Surf[IL8] = showChem_Surf[IL8]? false : true;
  printf("Toggle IL8 Surf: %d\n", showChem_Surf[IL8]);
}

void IL10_HMButtonCallBack()
{
	setAllChemFalse();
  showChem_HM[IL10] = showChem_HM[IL10]? false : true;
  printf("Toggle IL10 HM: %d\n", showChem_HM[IL10]);
}

void IL10_SurfButtonCallBack()
{
  showChem_Surf[IL10] = showChem_Surf[IL10]? false : true;
  printf("Toggle IL10 Surf: %d\n", showChem_Surf[IL10]);
}
void CellsButtonCallBack()
{
  showCells = showCells? false : true;
  printf("Toggle Show Cells: %d\n", showCells);
}

void NeuButtonCallBack()
{
  showNeus = showNeus? false : true;
  printf("Toggle Show Cells: %d\n", showNeus);
}

void MacButtonCallBack()
{
  showMacs = showMacs? false : true;
  printf("Toggle Show Cells: %d\n", showMacs);
}

void FibButtonCallBack()
{
  showFibs = showFibs? false : true;
  printf("Toggle Show Cells: %d\n", showFibs);
}


void ECMButtonCallBack()
{
  showECMs = showECMs? false : true;
  printf("Toggle Show ECMs: %d\n", showECMs);
}

void ZoomWoundButtonCallBack()
{
  zoomedWound = zoomedWound? false : true;
  if (zoomedWound) {
    pcamZ = camZ;
    pcamY = camY;
    pcamX = camX;
    pcellSize = cellSize;
    
#ifdef HUMAN_VF
    camZ = -850.0f;
    camY =  800.0f;
    camX =  860.0f;
#elif defined(RAT_VF)
    camZ = -950.0f;
    camY =  -50.0f;
    camX =  760.0f;
#else
    camZ = -650.0f;
    camY = 400.0f;
    camX = 660.0f;
#endif
    cellSize = 2.0f;
    
  } else {
    camX = pcamX;
    camY = pcamY;
    camZ = pcamZ;
    cellSize = pcellSize;
  }
  printf("Toggle Zoomed Wound: %d\n", zoomedWound);
}

int window = -1;


/*----------------------------------------------------------------------------------------
 *	This is the button visible in the viewport. This is a shorthand way of
 *	initialising the structure's data members. Notice that the last data
 *	member is a pointer to the above function.
 */



#ifdef USE_MOUSE_ONLY
Button ShowECMButton  = {40,  280, 55,25, 0,0, "ECMs", ECMButtonCallBack };
Button ShowCellsButton  = {WINW - 130,  225, 55,25, 0,0, "Cells", CellsButtonCallBack };

Button ZoomInButton         = {WINW - 110,  60, 100,25, 0,0, "Zoom In", ZoomInButtonCallback };
Button ZoomOutButton    = {WINW - 110,  100, 100,25, 0,0, "Zoom Out", ZoomOutButtonCallback };
Button MoveUpButton     = {WINW - 70,   130, 25,25, 0,0, "^", MoveUpButtonCallBack };
Button MoveLeftButton   = {WINW - 110,  160, 25,25, 0,0, "<", MoveLeftButtonCallBack };
Button MoveRightButton  = {WINW - 30,   160, 25,25, 0,0, ">", MoveRightButtonCallBack };
Button MoveDownButton   = {WINW - 70,   190, 25,25, 0,0, "v", MoveDownButtonCallBack };
Button ZoomWoundButton        = {WINW - 60,   225, 55,25, 0,0, "Wound", ZoomWoundButtonCallBack};

Button ShowTNF_HMButton    = {WINW - 130,  260, 55,25, 0,0, "TNF HM", TNF_HMButtonCallBack };
Button ShowTNF_SurfButton  = {WINW - 60,   260, 55,25, 0,0, "TNF Surf", TNF_SurfButtonCallBack };
Button ShowTGF_HMButton    = {WINW - 130,  295, 55,25, 0,0, "TGF HM", TGF_HMButtonCallBack };
Button ShowTGF_SurfButton  = {WINW - 60,   295, 55,25, 0,0, "TGF Surf", TGF_SurfButtonCallBack };
Button ShowFGF_HMButton    = {WINW - 130,  330, 55,25, 0,0, "FGF HM", TNF_HMButtonCallBack };
Button ShowFGF_SurfButton  = {WINW - 60,   330, 55,25, 0,0, "FGF Surf", TNF_SurfButtonCallBack };
Button ShowMMP8_HMButton    = {WINW - 130,  365, 55,25, 0,0, "MMP8 HM", TGF_HMButtonCallBack };
Button ShowMMP8_SurfButton  = {WINW - 60,   365, 55,25, 0,0, "MMP8 Surf", TGF_SurfButtonCallBack };
Button ShowIL1_HMButton    = {WINW - 130,  400, 55,25, 0,0, "IL1b HM", TNF_HMButtonCallBack };
Button ShowIL1_SurfButton  = {WINW - 60,   400, 55,25, 0,0, "IL1b Surf", TNF_SurfButtonCallBack };
Button ShowIL6_HMButton    = {WINW - 130,  435, 55,25, 0,0, "IL6 HM", TGF_HMButtonCallBack };
Button ShowIL6_SurfButton  = {WINW - 60,   435, 55,25, 0,0, "IL6 Surf", TGF_SurfButtonCallBack };
Button ShowIL8_HMButton    = {WINW - 130,  470, 55,25, 0,0, "IL8 HM", TNF_HMButtonCallBack };
Button ShowIL8_SurfButton  = {WINW - 60,   470, 55,25, 0,0, "IL8 Surf", TNF_SurfButtonCallBack };
Button ShowIL10_HMButton    = {WINW - 130,  505, 55,25, 0,0, "IL10 HM", TGF_HMButtonCallBack };
Button ShowIL10_SurfButton  = {WINW - 60,   505, 55,25, 0,0, "IL10 Surf", TGF_SurfButtonCallBack };

#else // USE_MOUSE_ONLY
Button ShowCellsButton  = {WINW - 100,  60, 45,25, 0,0, "Cells", CellsButtonCallBack };
Button ShowECMButton    = {WINW - 50,  60, 45,25, 0,0, "ECMs", ECMButtonCallBack };

Button ShowChemOpButton  = {WINW - 100,  90, 95,25, 0,0, "Chem Options", ChemOpButtonCallBack };

Button optionBoxArr[CHEM_OPTION_COLS * CHEM_OPTION_ROWS] = {{-1, -1, -1, -1, 0, 0, "", TNF_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", TNF_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", TGF_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", TGF_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", FGF_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", FGF_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", MMP8_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", MMP8_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", IL1beta_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", IL1beta_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", IL6_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", IL6_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", IL8_HMButtonCallBack},
//  {-1, -1, -1, -1, 0, 0, "", IL8_SurfButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", IL10_HMButtonCallBack}};
//  {-1, -1, -1, -1, 0, 0, "", IL10_SurfButtonCallBack}};

Button optionBoxArrCell[CELL_OPTION_COLS * CELL_OPTION_ROWS] = {{-1, -1, -1, -1, 0, 0, "", NeuButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", MacButtonCallBack},
  {-1, -1, -1, -1, 0, 0, "", FibButtonCallBack}};

#endif // USE_MOUSE_ONLY

/*----------------------------------------------------------------------------------------
 *	\brief	This function draws the specified button.
 *	\param	b	-	a pointer to the button to draw.
 */
void ButtonDraw(Button *b)
{
  int fontx;
  int fonty;
  
  if(b)
  {
    /*
     *	We will indicate that the mouse cursor is over the button by changing its
     *	colour.
     */
/*    if (b->highlighted)
      glColor3f(0.7f,0.7f,0.8f);
    else
      glColor3f(0.6f,0.6f,0.6f);
*/    
    /*
     *	draw background for the button.
     */
/*    glBegin(GL_QUADS);
    glVertex2i( b->x     , b->y      );
    glVertex2i( b->x     , b->y+b->h );
    glVertex2i( b->x+b->w, b->y+b->h );
    glVertex2i( b->x+b->w, b->y      );
    glEnd();
*/    
    /*
     *	Draw an outline around the button with width 3
     */
    glLineWidth(1);
//    glLineWidth(3);
    
    /*
     *	The colours for the outline are reversed when the button.
     */
    if (b->state)
      glColor3f(0.4f,0.4f,0.4f);
    else
      glColor3f(0.8f,0.8f,0.8f);
    
    glBegin(GL_LINE_STRIP);
    glVertex2i( b->x+b->w, b->y      );
    glVertex2i( b->x     , b->y      );
    glVertex2i( b->x     , b->y+b->h );
    glEnd();
    
    if (b->state)
      glColor3f(0.8f,0.8f,0.8f);
    else
      glColor3f(0.4f,0.4f,0.4f);
    
    glBegin(GL_LINE_STRIP);
    glVertex2i( b->x     , b->y+b->h );
    glVertex2i( b->x+b->w, b->y+b->h );
    glVertex2i( b->x+b->w, b->y      );
    glEnd();
    
    glLineWidth(1);
    
    
    /*
     *	Calculate the x and y coords for the text string in order to center it.
     */
    fontx = b->x + (b->w - glutBitmapLength(GLUT_BITMAP_HELVETICA_10, (unsigned char *) b->label)) / 2 ;
    fonty = b->y + (b->h+10)/2;
    
    /*
     *	if the button is pressed, make it look as though the string has been pushed
     *	down. It's just a visual thing to help with the overall look....
     */
    if (b->state) {
      fontx+=2;
      fonty+=2;
    }
    
    /*
     *	If the cursor is currently over the button we offset the text string and draw a shadow
     */
    if(b->highlighted)
    {
      glColor3f(0,0,0);
      Font(GLUT_BITMAP_HELVETICA_10,b->label,fontx,fonty);
      fontx--;
      fonty--;
    }
    
    glColor3f(1,1,1);
    Font(GLUT_BITMAP_HELVETICA_10,b->label,fontx,fonty);

    // DEBUG
//    Font(GLUT_BITMAP_HELVETICA_10, "HELLO", 100, 100);
//    Font(GLUT_BITMAP_HELVETICA_10, "-HELLO-", fontx, fonty);

  }
}

Mouse TheMouse = {0,0,0,0,0};





/*----------------------------------------------------------------------------------------
 *	\brief	This function is used to see if a mouse click or event is within a button
 *			client area.
 *	\param	b	-	a pointer to the button to test
 *	\param	x	-	the x coord to test
 *	\param	y	-	the y-coord to test
 */
int ButtonClickTest(Button* b,int x,int y)
{
  if( b)
  {
    /*
     *	If clicked within button area, then return true
     */
    if( x > b->x      &&
       x < b->x+b->w &&
       y > b->y      &&
       y < b->y+b->h ) {
      return 1;
    }
  }
  
  /*
   *	otherwise false.
   */
  return 0;
}

/*----------------------------------------------------------------------------------------
 *	\brief	This function draws the specified button.
 *	\param	b	-	a pointer to the button to check.
 *	\param	x	-	the x location of the mouse cursor.
 *	\param	y	-	the y location of the mouse cursor.
 */
void ButtonRelease(Button *b,int x,int y)
{
  if(b)
  {
    /*
     *	If the mouse button was pressed within the button area
     *	as well as being released on the button.....
     */
    if( ButtonClickTest(b,TheMouse.xpress,TheMouse.ypress) &&
       ButtonClickTest(b,x,y) )
    {
      /*
       *	Then if a callback function has been set, call it.
       */
      if (b->callbackFunction) {
        b->callbackFunction();
      }
    }
    
    /*
     *	Set state back to zero.
     */
    b->state = 0;
  }
}

/*----------------------------------------------------------------------------------------
 *	\brief	This function draws the specified button.
 *	\param	b	-	a pointer to the button to check.
 *	\param	x	-	the x location of the mouse cursor.
 *	\param	y	-	the y location of the mouse cursor.
 */
void ButtonPress(Button *b,int x,int y)
{
  if(b)
  {
    /*
     *	if the mouse click was within the buttons client area,
     *	set the state to true.
     */
    if( ButtonClickTest(b,x,y) )
    {
      b->state = 1;
    }
  }
}

/*----------------------------------------------------------------------------------------
 *	\brief	This function is called whenever a mouse button is pressed or released
 *	\param	button	-	GLUT_LEFT_BUTTON, GLUT_RIGHT_BUTTON, or GLUT_MIDDLE_BUTTON
 *	\param	state	-	GLUT_UP or GLUT_DOWN depending on whether the mouse was released
 *						or pressed respectivly.
 *	\param	x		-	the x-coord of the mouse cursor.
 *	\param	y		-	the y-coord of the mouse cursor.
 */
void MouseButton(int button,int state,int x, int y)
{
  /*
   *	update the mouse position
   */
  TheMouse.x = x;
  TheMouse.y = y;
 
  printf("Mouse clicked (%d, %d)\n", x, y);
 
  /*
   *	has the button been pressed or released?
   */
  if (state == GLUT_DOWN)
  {
    /*
     *	This holds the location of the first mouse click
     */
    if ( !(TheMouse.lmb || TheMouse.mmb || TheMouse.rmb) ) {
      TheMouse.xpress = x;
      TheMouse.ypress = y;
    }
   
  char strx[15];
  char stry[15];
  sprintf(strx, "%d", x);
  sprintf(stry, "%d", y);
Font3D(GLUT_BITMAP_HELVETICA_18, strx, 500, 500, PLANE_DEPTH + 100);
Font3D(GLUT_BITMAP_HELVETICA_18, stry, 500+20, 500, PLANE_DEPTH + 100);

 
    /*
     *	Which button was pressed?
     */
    switch(button)
    {
      case GLUT_LEFT_BUTTON:
        TheMouse.lmb = 1;
        // TODO: make more efficient

        ButtonPress(&ShowCellsButton, x, y);
        ButtonPress(&ShowECMButton, x, y);
        
#ifdef USE_MOUSE_ONLY
        ButtonPress(&ZoomInButton,x,y);
        ButtonPress(&ZoomOutButton, x, y);
        ButtonPress(&MoveUpButton, x, y);
        ButtonPress(&MoveLeftButton, x, y);
        ButtonPress(&MoveRightButton, x, y);
        ButtonPress(&MoveDownButton, x, y);
        ButtonPress(&ZoomWoundButton, x, y);
        
        ButtonPress(&ShowTNF_HMButton, x, y);
        ButtonPress(&ShowTGF_HMButton, x, y);
        ButtonPress(&ShowTNF_SurfButton, x, y);
        ButtonPress(&ShowTGF_SurfButton, x, y);
#else
        ButtonPress(&ShowChemOpButton, x, y);
        
        if (showChemOp) {
          for (int i = 0; i < CHEM_OPTION_COLS * CHEM_OPTION_ROWS; i++) {
            ButtonPress(&optionBoxArr[i], x, y);
          }
        }
        
        for (int i = 0; i < CELL_OPTION_COLS * CELL_OPTION_ROWS; i++) {
          ButtonPress(&optionBoxArrCell[i], x, y);
        }

#endif
        
        break;
      case GLUT_MIDDLE_BUTTON:
        TheMouse.mmb = 1;
        break;
      case GLUT_RIGHT_BUTTON:
        printf("Right clicked!\n");
//        ScreenPress(x, y);
        TheMouse.rmb = 1;
        break;
    }
  }
  else
  {
    /*
     *	Which button was released?
     */
    switch(button)
    {
      case GLUT_LEFT_BUTTON:
        TheMouse.lmb = 0;

        ButtonRelease(&ShowCellsButton, x, y);
        ButtonRelease(&ShowECMButton, x, y);
#ifdef USE_MOUSE_ONLY
        ButtonRelease(&ZoomInButton,x,y);
        ButtonRelease(&ZoomOutButton, x, y);
        ButtonRelease(&MoveUpButton, x, y);
        ButtonRelease(&MoveLeftButton, x, y);
        ButtonRelease(&MoveRightButton, x, y);
        ButtonRelease(&MoveDownButton, x, y);
        ButtonRelease(&ZoomWoundButton, x, y);
        
        ButtonRelease(&ShowTNF_HMButton, x, y);
        ButtonRelease(&ShowTGF_HMButton, x, y);
        ButtonRelease(&ShowTNF_SurfButton, x, y);
        ButtonRelease(&ShowTGF_SurfButton, x, y);
#else   // USE_MOUSE_ONLY
        
        ButtonRelease(&ShowChemOpButton, x, y);
        
        if (showChemOp) {
          for (int i = 0; i < CHEM_OPTION_COLS * CHEM_OPTION_ROWS; i++) {
            ButtonRelease(&optionBoxArr[i], x, y);
          }
        }

        for (int i = 0; i < CELL_OPTION_COLS * CELL_OPTION_ROWS; i++) {
          ButtonRelease(&optionBoxArrCell[i], x, y);
        }

#endif  // USE_MOUSE_ONLY
        
        break;
      case GLUT_MIDDLE_BUTTON:
        TheMouse.mmb = 0;
        break;
      case GLUT_RIGHT_BUTTON:
        TheMouse.rmb = 0;
        break;
    }
  }
  
  /*
   *	Force a redraw of the screen. If we later want interactions with the mouse
   *	and the 3D scene, we will need to redraw the changes.
   */
  glutPostRedisplay();
}


#ifdef USE_MOUSE_ONLY


void ScreenPress(int x, int y)
{
  glm::vec3 worldCoord = GetOGLPos(x,y);
  if (worldCoord[0] > 0 && worldCoord[0] < GRIDW && worldCoord[1] > 0 && worldCoord[1] < GRIDH) {
    printf("\tClick INSIDE World plane!!!\n");
    printf("\t\tx: %d, y: %d\n", x, y);
    printf("\t\tword.x: %f, world.y: %f, world.z: %f\n", worldCoord[0], worldCoord[1], worldCoord[2]);
  } else {
    printf("  Click OUTSIDE World plane!!!\n");
    printf("\t\tx: %d, y: %d\n", x, y);
    printf("\t\tword.x: %f, world.y: %f, world.z: %f\n", worldCoord[0], worldCoord[1], worldCoord[2]);
  }
}


/*----------------------------------------------------------------------------------------
 *	\brief	This function draws the specified button.
 *	\param	b	-	a pointer to the button to check.
 *	\param	x	-	the x location of the mouse cursor.
 *	\param	y	-	the y location of the mouse cursor.
 */
void ButtonPassive(Button *b,int x,int y)
{
  if(b)
  {
    /*
     *	if the mouse moved over the control
     */
    if( ButtonClickTest(b,x,y) )
    {
      /*
       *	If the cursor has just arrived over the control, set the highlighted flag
       *	and force a redraw. The screen will not be redrawn again until the mouse
       *	is no longer over this control
       */
      if( b->highlighted == 0 ) {
        b->highlighted = 1;
        glutPostRedisplay();
      }
    }
    else
      
    /*
     *	If the cursor is no longer over the control, then if the control
     *	is highlighted (ie, the mouse has JUST moved off the control) then
     *	we set the highlighting back to false, and force a redraw.
     */
      if( b->highlighted == 1 )
      {
        b->highlighted = 0;
        glutPostRedisplay();
      }
  }
}



/*----------------------------------------------------------------------------------------
 *	\brief	This function is called whenever the mouse cursor is moved AND A BUTTON IS HELD.
 *	\param	x	-	the new x-coord of the mouse cursor.
 *	\param	y	-	the new y-coord of the mouse cursor.
 */
void MouseMotion(int x, int y)
{
  /*
   *	Calculate how much the mouse actually moved
   */
  int dx = x - TheMouse.x;
  int dy = y - TheMouse.y;
  
  /*
   *	update the mouse position
   */
  TheMouse.x = x;
  TheMouse.y = y;
  
  
  /*
   *	Check ZoomInButton to see if we should highlight it cos the mouse is over it
   */
  ButtonPassive(&ZoomInButton,x,y);
  ButtonPassive(&ZoomOutButton, x, y);
  ButtonPassive(&MoveUpButton, x, y);
  ButtonPassive(&MoveLeftButton, x, y);
  ButtonPassive(&MoveRightButton, x, y);
  ButtonPassive(&MoveDownButton, x, y);
  
  /*
   *	Force a redraw of the screen
   */
  glutPostRedisplay();
}

/*----------------------------------------------------------------------------------------
 *	\brief	This function is called whenever the mouse cursor is moved AND NO BUTTONS ARE HELD.
 *	\param	x	-	the new x-coord of the mouse cursor.
 *	\param	y	-	the new y-coord of the mouse cursor.
 */
void MousePassiveMotion(int x, int y)
{
  /*
   *	Calculate how much the mouse actually moved
   */
  int dx = x - TheMouse.x;
  int dy = y - TheMouse.y;
  
  /*
   *	update the mouse position
   */
  TheMouse.x = x;
  TheMouse.y = y;
  
  /*
   *	Check ZoomInButton to see if we should highlight it cos the mouse is over it
   */
  ButtonPassive(&ZoomInButton,x,y);
  ButtonPassive(&ZoomOutButton, x, y);
  ButtonPassive(&MoveUpButton, x, y);
  ButtonPassive(&MoveLeftButton, x, y);
  ButtonPassive(&MoveRightButton, x, y);
  ButtonPassive(&MoveDownButton, x, y);
  
  /*
   *	Note that I'm not using a glutPostRedisplay() call here. The passive motion function
   *	is called at a very high frequency. We really don't want much processing to occur here.
   *	Redrawing the screen every time the mouse moves is a bit excessive. Later on we
   *	will look at a way to solve this problem and force a redraw only when needed.
   */
}


/*----------------------------------------------------------------------------------------
 *	Screen to world coordinates conversion
 */


float linearizeDepth(float nonLinearDepth, float nearPlane, float farPlane)
{
  return (2.0*nearPlane)/(farPlane+nearPlane-(nonLinearDepth*(farPlane-nearPlane)));
}

glm::vec3 GetOGLPos(int x, int y)
{
  //  GLfloat zz;
  //  glReadPixels( 400, 400, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &zz );
  //  printf("\t (1) Winx: 400.0, winY:400.0, Depth component: %f\n", zz);
  
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLfloat winX, winY, winZ;
  GLdouble posX, posY, posZ;
  
  glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
  glGetDoublev( GL_PROJECTION_MATRIX, projection );
  glGetIntegerv( GL_VIEWPORT, viewport );
  
  //  glReadPixels( 400, 400, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &zz );
  //  printf("\t (2)Winx: 400.0, winY:400.0, Depth component: %f\n", zz);
  
  winX = (float)x;
  winY = (float)viewport[3] - (float)y;
  glReadBuffer(GL_BACK);
  //  glClearDepth( 0.44 );
  //  glClear (GL_DEPTH_BUFFER_BIT);
  unsigned char color_rgb[3];
  glReadPixels( x, int(winY), 1, 1, GL_RGB, GL_UNSIGNED_BYTE, color_rgb );
  printf("\tColor_RGB: %d, %d, %d\n", color_rgb[0], color_rgb[1], color_rgb[2]);
  glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
  printf("\t\t\t Depth before linearize: %f\n", winZ);
  winZ = linearizeDepth(winZ, 0.1f, 6000);
  printf("\tWindow height: %d\n", viewport[3]);
  printf("\t Winx: %f, winY: %f, Depth component: %f\n", winX, winY, winZ);
  //  glLoadIdentity();
  //  glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
  //  printf("\tWindow height: %d\n", viewport[3]);
  //  printf("\tAfter loadIden() Winx: %f, winY: %f, Depth component: %f\n", winX, winY, winZ);
  
  //  glReadPixels( 400, 400, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &zz );
  //  printf("\t (3)Winx: 400.0, winY:400.0, Depth component: %f\n", zz);
  
  gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
  //  gluUnProject( winX, winY, 0.5, modelview, projection, viewport, &posX, &posY, &posZ);
  
  printf("[");
  for(int i = 0; i < 16; i += 4) {
    printf("%f, %f, %f, %f\n", modelview[i], modelview[i+1], modelview[i+2], modelview[i+3]);
  }
  printf("]\n");
  
  printf("[");
  for(int i = 0; i < 16; i += 4) {
    printf("%f, %f, %f, %f\n", projection[i], projection[i+1], projection[i+2], projection[i+3]);
  }
  printf("]\n");
  
  printf("["); printf("%d, %d, %d, %d", viewport[0], viewport[1], viewport[2], viewport[3]);
  printf("]\n");
  
  return glm::vec3(posX, posY, posZ);
}

#endif      // USE_MOUSE_ONLY
