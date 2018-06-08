//
//  OptionTable.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 10/19/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#include "OptionTable.h"


OptionTable::OptionTable(GLfloat ulX,
                         GLfloat ulY,
                         GLfloat tableW,
                         GLfloat tableH,
                         int nRowsOption,
                         int nColsOption,
                         char **textHeaders,
                         char **textLabels,
                         Button *boxArray,
                         GLfloat colors[][4])
{
  this->ulX         = ulX;
  this->ulY         = ulY;
  this->tableW      = tableW;
  this->tableH      = tableH;
  this->nRowsOption = nRowsOption;
  this->nColsOption = nColsOption;
  this->textHeaders = textHeaders;
  this->textLabels  = textLabels;
  this->boxArray    = boxArray;
  
  this->nRows       = 1 + nRowsOption;
  this->nCols       = 1 + nColsOption;

  this->colors      = colors;
  
  InitBoxLocations();
}

OptionTable::~OptionTable()
{
  // TODO: Cleanup
}

//float *getSurfColorPtr()
//{
//  return this->surfColor;
//}


void OptionTable::Render()
{
  // Draw outline
  glColor3f(0.5f, 0.5f, 0.5f);
  drawTransparentRect(ulX + (tableW/2),
                      ulY + (tableH/2),
                      tableW,
                      tableH);

  glColor3f(1.0f, 1.0f, 1.0f);
  // Render Headers
  DrawHeaders();
  // Render Labels
  DrawLabels();
  // Render option boxes
  DrawOptionBoxes();
}

void OptionTable::InitBoxLocations()
{
  if (boxArray) {
    GLfloat posHorizantalInc  = ((float) tableW)/((float) nCols);
    GLfloat halfSpaceW        = (posHorizantalInc/2) - OPTION_BOX_W;
    GLfloat horizontal_offset = posHorizantalInc + ulX + halfSpaceW;
    
    GLfloat posVerticalInc    = ((float) tableH)/((float) nRows);
    GLfloat halfSpaceH        = (posVerticalInc/2) - OPTION_BOX_H;
    GLfloat charH             = 0.0;
    GLfloat vertical_offset   = ulY + posVerticalInc + charH + halfSpaceH;
    for (int i = 0; i < nRowsOption; i++) {
      for (int j = 0; j < nColsOption; j++) {
        int index = i * nColsOption + j;
        boxArray[index].x = horizontal_offset + posHorizantalInc * j;
        boxArray[index].y = vertical_offset   + posVerticalInc * i;
        boxArray[index].w = OPTION_BOX_W;
        boxArray[index].h = OPTION_BOX_H;
      }
    }
  }
}

void OptionTable::DrawHeaders()
{
  if (textHeaders) {
    GLfloat posInc = ((float) tableW)/((float) nCols);
    GLfloat horizontal_offset = posInc + ulX;
    
    for (int i = 0; i < nColsOption; i++) {
      Font(GLUT_BITMAP_HELVETICA_10, textHeaders[i], horizontal_offset + posInc * i, ulY + 20);
    }
  }
}

void OptionTable::DrawLabels()
{
  if (textLabels) {
    GLfloat posInc = ((float) tableH)/((float) nRows);
    GLfloat charH = 12.0;
    GLfloat vertical_offset = ulY + posInc + charH;
    
    for (int i = 0; i < nRowsOption; i++) {
      Font(GLUT_BITMAP_HELVETICA_10, textLabels[i], ulX + 10, vertical_offset + posInc * i);
    }
  }
}

void OptionTable::DrawOptionBoxes()
{
  if (boxArray) {
    int nTotal = nColsOption * nRowsOption;
    Button *b;
    int ci = 0;
    for (int i = 0; i < nTotal; i++) {
      ci = i/2;
      b = &boxArray[i];
      
      /*
       *	We will indicate that the mouse cursor is over the button by changing its
       *	colour.
       */
      if (b->highlighted) {
        if (colors) {
          glColor4f(colors[i][0], colors[i][1], colors[i][2], colors[i][3]);
        } else {
          glColor3f(0.7f,0.7f,0.8f);
        }
      } else {
        glColor3f(0.3f,0.3f,0.3f);
      }

      /*
       *	draw background for the button.
       */
      glBegin(GL_QUADS);
      glVertex2i( b->x     , b->y      );
      glVertex2i( b->x     , b->y+b->h );
      glVertex2i( b->x+b->w, b->y+b->h );
      glVertex2i( b->x+b->w, b->y      );
      glEnd();
      
      /*
       *	Draw an outline around the button with width 3
       */
      glLineWidth(0.2);
      
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
      
    }
  }
}




