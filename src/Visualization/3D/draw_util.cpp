//
//  draw_util.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/24/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//


#include "draw_util.h"

/*----------------------------------------------------------------------------------------
 *	\brief	This function draws a text string to the screen using glut bitmap fonts.
 *	\param	font	-	the font to use. it can be one of the following :
 *
 *					GLUT_BITMAP_9_BY_15
 *					GLUT_BITMAP_8_BY_13
 *					GLUT_BITMAP_TIMES_ROMAN_10
 *					GLUT_BITMAP_TIMES_ROMAN_24
 *					GLUT_BITMAP_HELVETICA_10
 *					GLUT_BITMAP_HELVETICA_12
 *					GLUT_BITMAP_HELVETICA_18
 *
 *	\param	text	-	the text string to output
 *	\param	x		-	the x co-ordinate
 *	\param	y		-	the y co-ordinate
 */
void Font(void *font,char *text,int x,int y)
{
  glRasterPos2i(x, y);
  
  while( *text != '\0' )
  {
    glutBitmapCharacter( font, *text );
    ++text;
  }
}

void Font3D(void *font,char *text,int x,int y, int z)
{
  glRasterPos3i(x, y, z);
  
  while( *text != '\0' )
  {
    glutBitmapCharacter( font, *text );
    ++text;
  }
}

void drawTexturedRect(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb,
                      GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat height)
{
  GLfloat halfWidth = width / 2.0f;
  GLfloat halfHeight = height / 2.0f;
  GLfloat posx = GL_X(cx + halfWidth);
  GLfloat posy = GL_Y(cy + halfHeight);
  GLfloat negx = GL_X(cx - halfWidth);
  GLfloat negy = GL_Y(cy - halfHeight);
  
  glBegin(GL_QUADS);
  glTexCoord2f( txl, txt ); glVertex3f( negx, negy, cz );
  glTexCoord2f( txr, txt ); glVertex3f( posx, negy, cz );
  glTexCoord2f( txr, txb ); glVertex3f( posx, posy, cz );
  glTexCoord2f( txl, txb ); glVertex3f( negx, posy, cz );
  glEnd( );
}


void drawTransparentRect(GLfloat cx, GLfloat cy, GLfloat width, GLfloat height)
{
  GLfloat halfWidth = width / 2.0f;
  GLfloat halfHeight = height / 2.0f;
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfHeight;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfHeight;
  
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex2f( negx, negy );
  glVertex2f( posx, negy );
  glVertex2f( posx, posy );
  glVertex2f( negx, posy );
  glEnd( );
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void drawTransparentBox(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat depth, GLfloat height){
  GLfloat halfWidth = width / 2.0f;
  GLfloat halfDepth = depth / 2.0f;
  GLfloat halfHeight = height / 2.0f;
  
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfDepth;
  GLfloat posz = cz + halfHeight;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfDepth;
  GLfloat negz = cz - halfHeight;
  
  /*      bnw __________ bne  */
  /*        /|         /|     */
  /*       / |        / |     */
  /*   fnw/__________fne|     */
  /*     |bsw|_______|__|bse  */
  /*     |  /        |  /     */
  /*     | /         | /      */
  /*     |/__________|/       */
  /*   fsw           fse      */
  
  GLfloat fsw[3] = {negx, negy, posz};
  GLfloat fnw[3] = {negx, posy, posz};
  GLfloat fne[3] = {posx, posy, posz};
  GLfloat fse[3] = {posx, negy, posz};
  
  GLfloat bsw[3] = {negx, negy, negz};
  GLfloat bnw[3] = {negx, posy, negz};
  GLfloat bne[3] = {posx, posy, negz};
  GLfloat bse[3] = {posx, negy, negz};
  
  
  /*
   * Draw the cube. A cube consists of six quads, with four coordinates (glVertex3f)
   * per quad.
   *
   */
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  /* Front Face */
  glVertex3f( fsw[0], fsw[1], fsw[2] );     // fsw
  glVertex3f( fse[0], fse[1], fse[2] );     // fse
  glVertex3f( fne[0], fne[1], fne[2] );     // fne
  glVertex3f( fnw[0], fnw[1], fnw[2] );     // fnw
  
  /* Back Face */
  glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glVertex3f( bne[0], bne[1], bne[2] );    // bne
  glVertex3f( bse[0], bse[1], bse[2] );    // bse
  
  /* Top Face */
  glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glVertex3f( fnw[0], fnw[1], fnw[2] );    // fnw
  glVertex3f( fne[0], fne[1], fne[2] );    // fne
  glVertex3f( bne[0], bne[1], bne[2] );    // bne
  
  /* Bottom Face */
  /* Top Right Of The Texture and Quad */
  glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glVertex3f( bse[0], bse[1], bse[2] );    // bse
  glVertex3f( fse[0], fse[1], fse[2] );    // fse
  glVertex3f( fsw[0], fsw[1], fsw[2] );    // fsw
  
  /* Right face */
  glVertex3f( bse[0], bse[1], bse[2] );     // bse
  glVertex3f( bne[0], bne[1], bne[2] );     // bne
  glVertex3f( fne[0], fne[1], fne[2] );     // fne
  glVertex3f( fse[0], fse[1], fse[2] );     // fse
  
  /* Left Face */
  glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glVertex3f( fsw[0], fsw[1], fsw[2] );    // fsw
  glVertex3f( fnw[0], fnw[1], fnw[2] );    // fnw
  glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glEnd( );
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}


void drawRectangle(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat height)
{
  GLfloat halfWidth = width / 2.0f;
  GLfloat halfHeight = height / 2.0f;
  GLfloat posx = GL_X(cx + halfWidth);
  GLfloat posy = GL_Y(cy + halfHeight);
  GLfloat negx = GL_X(cx - halfWidth);
  GLfloat negy = GL_Y(cy - halfHeight);
  
  glBegin(GL_QUADS);
  glVertex3f( negx, negy, cz );
  glVertex3f( posx, negy, cz );
  glVertex3f( posx, posy, cz );
  glVertex3f( negx, posy, cz );
  glEnd( );
}

//void drawRectangleRepeat(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb, GLfloat txw,
//                    GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width, GLfloat height)
//{
//  GLint nFullSquaresX = width / txw;
//  GLint nFullSquaresY = height / txw;
//  GLint nFullSquares = nFullSquaresX * nFullSquaresY;
//  
//  GLfloat halfWidth = width / 2.0f;
//  GLfloat posx = cx + halfWidth;
//  GLfloat posy = cy + halfWidth;
//  GLfloat negx = cx - halfWidth;
//  GLfloat negy = cy - halfWidth;
//  
//  glBegin(GL_QUADS);
//  glTexCoord2f( txl, txt ); glVertex3f( negx, negy, cz );
//  glTexCoord2f( txr, txt ); glVertex3f( posx, negy, cz );
//  glTexCoord2f( txr, txb ); glVertex3f( posx, posy, cz );
//  glTexCoord2f( txl, txb ); glVertex3f( negx, posy, cz );
//  glEnd( );
//}

void drawSquare2D(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb,
                  GLfloat cx, GLfloat cy, GLfloat width)
{
  GLfloat halfWidth = width / 2.0f;
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfWidth;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfWidth;
  
  
  
  glBegin(GL_QUADS);
  glTexCoord2f( txl, txt ); glVertex2f( negx, negy );
  glTexCoord2f( txr, txt ); glVertex2f( posx, negy );
  glTexCoord2f( txr, txb ); glVertex2f( posx, posy );
  glTexCoord2f( txl, txb ); glVertex2f( negx, posy );
  glEnd( );
}

QuadV getQuadVertices(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width)
{
  QuadV ret;
  GLfloat halfWidth = width / 2.0f;
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfWidth;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfWidth;
  
  ret[0] = glm::vec3( negx, negy, cz );
  ret[1] = glm::vec3( posx, negy, cz );
  ret[2] = glm::vec3( posx, posy, cz );
  ret[3] = glm::vec3( negx, posy, cz );

  return ret;
}

void drawSquare(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width)
{
  GLfloat halfWidth = width / 2.0f;
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfWidth;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfWidth;
  
  glBegin(GL_QUADS);
  glVertex3f( negx, negy, cz );
  glVertex3f( posx, negy, cz );
  glVertex3f( posx, posy, cz );
  glVertex3f( negx, posy, cz );
  glEnd( );
}

void drawSquare(GLfloat txl, GLfloat txr, GLfloat txt, GLfloat txb,
                GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width)
{
  
  GLfloat halfWidth = width / 2.0f;
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfWidth;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfWidth;
  
  
  
  glBegin(GL_QUADS);
  glTexCoord2f( txl, txt ); glVertex3f( negx, negy, cz );
  glTexCoord2f( txr, txt ); glVertex3f( posx, negy, cz );
  glTexCoord2f( txr, txb ); glVertex3f( posx, posy, cz );
  glTexCoord2f( txl, txb ); glVertex3f( negx, posy, cz );
  glEnd( );
}

void drawSolidCube(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width) {
  
  GLfloat halfWidth = width / 2.0f;
   
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfWidth;
  GLfloat posz = cz + halfWidth;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfWidth;
  GLfloat negz = cz - halfWidth;
  
  /*      bnw __________ bne  */
  /*        /|         /|     */
  /*       / |        / |     */
  /*   fnw/__________fne|     */
  /*     |bsw|_______|__|bse  */
  /*     |  /        |  /     */
  /*     | /         | /      */
  /*     |/__________|/       */
  /*   fsw           fse      */
  
  GLfloat fsw[3] = {negx, negy, posz};
  GLfloat fnw[3] = {negx, posy, posz};
  GLfloat fne[3] = {posx, posy, posz};
  GLfloat fse[3] = {posx, negy, posz};
  
  GLfloat bsw[3] = {negx, negy, negz};
  GLfloat bnw[3] = {negx, posy, negz};
  GLfloat bne[3] = {posx, posy, negz};
  GLfloat bse[3] = {posx, negy, negz};
  
  
  /*
   * Draw the cube. A cube consists of six quads, with four coordinates (glVertex3f)
   * per quad.
   *
   */
  glBegin(GL_QUADS);
  /* Front Face */
  glVertex3f( fsw[0], fsw[1], fsw[2] );     // fsw
  glVertex3f( fse[0], fse[1], fse[2] );     // fse
  glVertex3f( fne[0], fne[1], fne[2] );     // fne
  glVertex3f( fnw[0], fnw[1], fnw[2] );     // fnw
  
  /* Back Face */
  glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glVertex3f( bne[0], bne[1], bne[2] );    // bne
  glVertex3f( bse[0], bse[1], bse[2] );    // bse
  
  /* Top Face */
  glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glVertex3f( fnw[0], fnw[1], fnw[2] );    // fnw
  glVertex3f( fne[0], fne[1], fne[2] );    // fne
  glVertex3f( bne[0], bne[1], bne[2] );    // bne
  
  /* Bottom Face */
  glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glVertex3f( bse[0], bse[1], bse[2] );    // bse
  glVertex3f( fse[0], fse[1], fse[2] );    // fse
  glVertex3f( fsw[0], fsw[1], fsw[2] );    // fsw
  
  /* Right face */
  glVertex3f( bse[0], bse[1], bse[2] );     // bse
  glVertex3f( bne[0], bne[1], bne[2] );     // bne
  glVertex3f( fne[0], fne[1], fne[2] );     // fne
  glVertex3f( fse[0], fse[1], fse[2] );     // fse
  
  /* Left Face */
  glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glVertex3f( fsw[0], fsw[1], fsw[2] );    // fsw
  glVertex3f( fnw[0], fnw[1], fnw[2] );    // fnw
  glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glEnd( );

  
}
void drawCube(GLfloat cx, GLfloat cy, GLfloat cz, GLfloat width) {
  
  GLfloat halfWidth = width / 2.0f;
  
  GLfloat txbl[2] = {0.25, 0.75};
  GLfloat txtl[2] = {0.25, 1.00};
  GLfloat txtr[2] = {0.50, 1.00};
  GLfloat txbr[2] = {0.50, 0.75};
  
  GLfloat posx = cx + halfWidth;
  GLfloat posy = cy + halfWidth;
  GLfloat posz = cz + halfWidth;
  GLfloat negx = cx - halfWidth;
  GLfloat negy = cy - halfWidth;
  GLfloat negz = cz - halfWidth;
  
  /*      bnw __________ bne  */
  /*        /|         /|     */
  /*       / |        / |     */
  /*   fnw/__________fne|     */
  /*     |bsw|_______|__|bse  */
  /*     |  /        |  /     */
  /*     | /         | /      */
  /*     |/__________|/       */
  /*   fsw           fse      */
  
  GLfloat fsw[3] = {negx, negy, posz};
  GLfloat fnw[3] = {negx, posy, posz};
  GLfloat fne[3] = {posx, posy, posz};
  GLfloat fse[3] = {posx, negy, posz};
  
  GLfloat bsw[3] = {negx, negy, negz};
  GLfloat bnw[3] = {negx, posy, negz};
  GLfloat bne[3] = {posx, posy, negz};
  GLfloat bse[3] = {posx, negy, negz};
  
  
  /*
   * Draw the cube. A cube consists of six quads, with four coordinates (glVertex3f)
   * per quad.
   *
   */
  glBegin(GL_QUADS);
  /* Front Face */
  glTexCoord2f( txtl[0], txtl[1] ); glVertex3f( fsw[0], fsw[1], fsw[2] );     // fsw
  glTexCoord2f( txtr[0], txtr[1] ); glVertex3f( fse[0], fse[1], fse[2] );     // fse
  glTexCoord2f( txbr[0], txbr[1] ); glVertex3f( fne[0], fne[1], fne[2] );     // fne
  glTexCoord2f( txbl[0], txbl[1] ); glVertex3f( fnw[0], fnw[1], fnw[2] );     // fnw
  
  /* Back Face */
  glTexCoord2f( txbl[0], txbl[1] ); glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glTexCoord2f( txtl[0], txtl[1] ); glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glTexCoord2f( txtr[0], txtr[1] ); glVertex3f( bne[0], bne[1], bne[2] );    // bne
  glTexCoord2f( txbr[0], txbr[1] ); glVertex3f( bse[0], bse[1], bse[2] );    // bse
  
  /* Top Face */
  glTexCoord2f( txtr[0], txtr[1] ); glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glTexCoord2f( txbr[0], txbr[1] ); glVertex3f( fnw[0], fnw[1], fnw[2] );    // fnw
  glTexCoord2f( txbl[0], txbl[1] ); glVertex3f( fne[0], fne[1], fne[2] );    // fne
  glTexCoord2f( txtl[0], txtl[1] ); glVertex3f( bne[0], bne[1], bne[2] );    // bne
  
  /* Bottom Face */
  /* Top Right Of The Texture and Quad */
  glTexCoord2f( txtl[0], txtl[1] ); glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glTexCoord2f( txtr[0], txtr[1] ); glVertex3f( bse[0], bse[1], bse[2] );    // bse
  glTexCoord2f( txbr[0], txbr[1] ); glVertex3f( fse[0], fse[1], fse[2] );    // fse
  glTexCoord2f( txbl[0], txbl[1] ); glVertex3f( fsw[0], fsw[1], fsw[2] );    // fsw
  
  /* Right face */
  glTexCoord2f( txbl[0], txbl[1] ); glVertex3f( bse[0], bse[1], bse[2] );     // bse
  glTexCoord2f( txtl[0], txtl[1] ); glVertex3f( bne[0], bne[1], bne[2] );     // bne
  glTexCoord2f( txtr[0], txtr[1] ); glVertex3f( fne[0], fne[1], fne[2] );     // fne
  glTexCoord2f( txbr[0], txbr[1] ); glVertex3f( fse[0], fse[1], fse[2] );     // fse
  
  /* Left Face */
  glTexCoord2f( txbr[0], txbr[1] ); glVertex3f( bsw[0], bsw[1], bsw[2] );    // bsw
  glTexCoord2f( txbl[0], txbl[1] ); glVertex3f( fsw[0], fsw[1], fsw[2] );    // fsw
  glTexCoord2f( txtl[0], txtl[1] ); glVertex3f( fnw[0], fnw[1], fnw[2] );    // fnw
  glTexCoord2f( txtr[0], txtr[1] ); glVertex3f( bnw[0], bnw[1], bnw[2] );    // bnw
  glEnd( );

  
}
