//
//  common.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 1/3/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//


int GridH = -1;
int GridW = -1;
int GridD = -1;

int GL_X(int x) {
  return (x - (GridW/2));
}

int GL_Y(int y) {
 return (-(y - (GridH/2)));
}


/*int mmToPatch(float mm){
  return mm / PATCHW;  // each patch is 15um x 15 um -> 0.015 mm x 0.015 mm
}*/
