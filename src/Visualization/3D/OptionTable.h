//
//  OptionTable.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 10/19/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_OptionTable_h
#define OpenGL_SimpleWH_OptionTable_h

#include "common_vis.h"
#include "UI.h"

#define OPTION_BOX_W  10
#define OPTION_BOX_H  10

///*
// *	We will define a function pointer type. ButtonCallback is a pointer to a function that
// *	looks a bit like this :
// *
// *	void func() {
// *	}
// */
//typedef void (*OptionCallback)();
//
///*
// *	This is a simple structure that holds a button.
// */
//struct OptionBox
//{
//  int   x;							/* top left x coord of the button */
//  int   y;							/* top left y coord of the button */
//  int   w;							/* the width of the button */
//  int   h;							/* the height of the button */
//  int	state;						/* the state, 1 if pressed, 0 otherwise */
//  int	highlighted;					/* is the mouse cursor over the control? */
//  OptionCallback callbackFunction;                      /* A pointer to a function to call if the button is pressed */
//};
//typedef struct OptionBox OptionBox;

class OptionTable{
  
public:
  OptionTable(GLfloat ulX,
              GLfloat ulY,
              GLfloat tableW,
              GLfloat tableH,
              int nRowsOption,
              int nColsOption,
              char **textHeaders,
              char **textLabels,
//              OptionBox *boxArray);
              Button *boxArray,
	      GLfloat colors[][4]);
  
  virtual ~OptionTable();
  
  void Render();
  
private:

  GLfloat ulX, ulY;   // upper left corner coordinate
  GLfloat tableW, tableH;
  
  int nRowsOption;
  int nColsOption;
  
  char **textHeaders;
  char **textLabels;
//  OptionBox *boxArray;
  Button *boxArray;
  GLfloat (*colors)[4];  

  void InitBoxLocations();
  void DrawHeaders();
  void DrawLabels();
  void DrawOptionBoxes();
  
  int nRows;
  int nCols;
  
};


#endif
