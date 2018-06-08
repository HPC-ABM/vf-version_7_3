//
//  UI.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/28/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_UI_h
#define OpenGL_SimpleWH_UI_h

#include "common_vis.h"
#include "draw_util.h"


/*----------------------------------------------------------------------------------------
 *	Global Variables
 */

/*
 *	A structure to represent the mouse information
 */
struct Mouse
{
  int x;		/*	the x coordinate of the mouse cursor	*/
  int y;		/*	the y coordinate of the mouse cursor	*/
  int lmb;	/*	is the left button pressed?		*/
  int mmb;	/*	is the middle button pressed?	*/
  int rmb;	/*	is the right button pressed?	*/
  
  /*
   *	These two variables are a bit odd. Basically I have added these to help replicate
   *	the way that most user interface systems work. When a button press occurs, if no
   *	other button is held down then the co-ordinates of where that click occured are stored.
   *	If other buttons are pressed when another button is pressed it will not update these
   *	values.
   *
   *	This allows us to "Set the Focus" to a specific portion of the screen. For example,
   *	in maya, clicking the Alt+LMB in a view allows you to move the mouse about and alter
   *	just that view. Essentually that viewport takes control of the mouse, therefore it is
   *	useful to know where the first click occured....
   */
  int xpress; /*	stores the x-coord of when the first button press occurred	*/
  int ypress; /*	stores the y-coord of when the first button press occurred	*/
};

/*
 *	rename the structure from "struct Mouse" to just "Mouse"
 */
typedef struct Mouse Mouse;




/*----------------------------------------------------------------------------------------
 *	Button Stuff
 */

/*
 *	We will define a function pointer type. ButtonCallback is a pointer to a function that
 *	looks a bit like this :
 *
 *	void func() {
 *	}
 */
typedef void (*ButtonCallback)();

/*
 *	This is a simple structure that holds a button.
 */
struct Button
{
  int   x;							/* top left x coord of the button */
  int   y;							/* top left y coord of the button */
  int   w;							/* the width of the button */
  int   h;							/* the height of the button */
  int	  state;						/* the state, 1 if pressed, 0 otherwise */
  int	  highlighted;					/* is the mouse cursor over the control? */
  char* label;						/* the text label of the button */
  ButtonCallback callbackFunction;	/* A pointer to a function to call if the button is pressed */
};
typedef struct Button Button;

void ZoomInButtonCallback();
void ZoomOutButtonCallback();
void ZoomWoundButtonCallBack();

void MoveUpButtonCallBack();
void MoveDownButtonCallBack();
void MoveLeftButtonCallBack();
void MoveRightButtonCallBack();

void TNF_HMButtonCallBack();

// actual vector representing the camera's direction
extern float lx, ly, lz;
// XZ position of the camera
extern float camX, camY, camZ;
extern float pcamX, pcamY, pcamZ;

extern GLfloat cellSize;
extern GLfloat pcellSize;       // previous cell size

extern bool paused;

// Chem global variables
//extern bool showTNF_HM   ;
//extern bool showTNF_Surf ;
//extern bool showTGF_HM   ;
//extern bool showTGF_Surf ;
extern bool showChem_HM  [8];
extern bool showChem_Surf[8];
extern bool showCells       ;
extern bool showNeus        ;
extern bool showMacs        ;
extern bool showFibs        ;
extern bool showECMs        ;
extern bool showChemOp      ;

extern bool showChemCharts  ;
extern bool showCellCharts  ;

extern bool lowerECMplane ;
extern bool raiseECMplane ;

extern bool zoomedWound ;


/* The number of our GLUT window */
extern int window;

/*----------------------------------------------------------------------------------------
 *	This is the button visible in the viewport. This is a shorthand way of
 *	initialising the structure's data members. Notice that the last data
 *	member is a pointer to the above function.
 */

extern Button ShowCellsButton;
extern Button ShowECMButton;


#ifdef USE_MOUSE_ONLY

extern Button ZoomInButton;
extern Button ZoomOutButton;
extern Button MoveUpButton;
extern Button MoveLeftButton;
extern Button MoveRightButton;
extern Button MoveDownButton;
extern Button ZoomWoundButton;

extern Button ShowTNF_HMButton;
extern Button ShowTNF_SurfButton;
extern Button ShowTGF_HMButton;
extern Button ShowTGF_SurfButton;
extern Button ShowFGF_HMButton;
extern Button ShowFGF_SurfButton;
extern Button ShowMMP8_HMButton;
extern Button ShowMMP8_SurfButton;
extern Button ShowIL1_HMButton;
extern Button ShowIL1_SurfButton;
extern Button ShowIL6_HMButton;
extern Button ShowIL6_SurfButton;
extern Button ShowIL8_HMButton;
extern Button ShowIL8_SurfButton;
extern Button ShowIL10_HMButton;
extern Button ShowIL10_SurfButton;

#else   // USE_MOUSE_ONLY

extern Button ShowChemOpButton;
extern Button optionBoxArr[CHEM_OPTION_COLS * CHEM_OPTION_ROWS];
extern Button optionBoxArrCell[CELL_OPTION_COLS * CELL_OPTION_ROWS];

#endif  // USE_MOUSE_ONLY

void ButtonDraw(Button *b);

/*
 *	Create a global mouse structure to hold the mouse information.
 */
extern Mouse TheMouse;

void MouseButton(int button,int state,int x, int y);
int ButtonClickTest(Button* b,int x,int y);
void ButtonRelease(Button *b,int x,int y);
void ButtonPress(Button *b,int x,int y);

void keyPressed(unsigned char key, int x, int y);


#ifdef USE_MOUSE_ONLY

glm::vec3 GetOGLPos(int x, int y);
void ButtonPassive(Button *b,int x,int y);

void ScreenPress(int x, int y);

void MouseMotion(int x, int y);
void MousePassiveMotion(int x, int y);

#else     // USE_MOUSE_ONLY

void SpecialInput(int key, int x, int y);

#endif    // USE_MOUSE_ONLY




#endif
