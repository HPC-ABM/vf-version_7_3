//
//  texture_util.h
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/22/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#ifndef OpenGL_SimpleWH_texture_util_h
#define OpenGL_SimpleWH_texture_util_h

#include "common_vis.h"

/* ================================================================================================ */
/* Image type - contains height, width, and data */
struct Image {
  unsigned long sizeX;
  unsigned long sizeY;
  unsigned char *data;
};
typedef struct Image Image;


static unsigned int getint(FILE *fp);
static unsigned int getshort(FILE *fp);

int ImageLoad(char *filename, Image *image);
int ImageLoad32(char *filename, Image *image, ImLdOption option);

GLvoid LoadGLTextures_Sprite(GLuint *texture_sprite);
GLuint LoadPointTextureSprite();
GLvoid LoadGLTextures(char *filename, GLuint *texture, int size);
GLvoid LoadGLTexture(char *filename, GLuint *texture, ImLdOption option);
GLvoid LoadGLTexture24(char *filename, GLuint *texture);

#endif
