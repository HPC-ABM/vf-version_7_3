//
//  draw_util.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/24/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_draw_util_h
#define OpenGL_SimpleWH_draw_util_h

#include "common_vis.h"

void Font(void *font,char *text,int x,int y);

void Font3D(void *font, char *text, int x, int y, int z);

void drawSquare(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width);

void drawSquare(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb,
                GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width);

void drawSquare2D(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb,
                GLfloat cx, GLfloat cy, GLfloat width);

void drawTransparentRect(GLfloat cx, GLfloat cy, GLfloat width, GLfloat height);

void drawTexturedRect(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb,
                      GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat height);

void drawRectangle(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat height);

//void drawRectangleRepeat(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb, GLfloat txw,
//                   GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat height);

void drawSolidCube(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width);

void drawCube(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width);

void drawTransparentBox(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat depth, GLfloat height);

#endif
