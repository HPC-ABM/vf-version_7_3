//
//  texture_util.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 9/22/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#include "texture_util.h"

/*
 * getint and getshort arehelp functions to load the bitmap byte by byte on
 * SPARC platform.
 * I've got them from xv bitmap load routinebecause the original bmp loader didn't work
 * I've tried to change as less code as possible.
 */

static unsigned int getint(FILE *fp)
//FILE *fp;
{
  int c, c1, c2, c3;
  
  // get 4 bytes
  c = getc(fp);
  c1 = getc(fp);
  c2 = getc(fp);
  c3 = getc(fp);
  
  return ((unsigned int) c) +
  (((unsigned int) c1) << 8) +
  (((unsigned int) c2) << 16) +
  (((unsigned int) c3) << 24);
}

static unsigned int getshort(FILE *fp)
//FILE *fp;
{
  int c, c1;
  
  //get 2 bytes
  c = getc(fp);
  c1 = getc(fp);
  
  return ((unsigned int) c) + (((unsigned int) c1) << 8);
}


// quick and dirty bitmap loader...for 24 bit bitmaps with 1 plane only.
// See http://www.dcs.ed.ac.uk/~mxr/gfx/2d/BMP.txt for more info.

int ImageLoad(char *filename, Image *image) {
  FILE *file;
  unsigned long size;                 // size of the image in bytes.
  unsigned long i;                    // standard counter.
  unsigned short int planes;          // number of planes in image (must be 1)
  unsigned short int bpp;             // number of bits per pixel (must be 24)
  char temp;                          // used to convert bgr to rgb color.
  
  // make sure the file is there.
  if ((file = fopen(filename, "rb"))==NULL) {
    printf("File Not Found : %s\n",filename);
    return 0;
  }
  
  // seek through the bmp header, up to the width/height:
  fseek(file, 18, SEEK_CUR);
  
  // No 100% errorchecking anymore!!!
  
  // read the width
  image->sizeX = getint (file);
  printf("Width of %s: %lu\n", filename, image->sizeX);
  
  // read the height
  image->sizeY = getint (file);
  printf("Height of %s: %lu\n", filename, image->sizeY);
  
  // calculate the size (assuming 24 bits or 3 bytes per pixel).
  size = image->sizeX * image->sizeY * 3;
  
  // read the planes
  planes = getshort(file);
  if (planes != 1) {
    printf("Planes from %s is not 1: %u\n", filename, planes);
    return 0;
  }
  
  // read the bpp
  bpp = getshort(file);
  if (bpp != 24) {
    printf("Bpp from %s is not 24: %u\n", filename, bpp);
    return 0;
  }
  
  // seek past the rest of the bitmap header.
  fseek(file, 24, SEEK_CUR);
  
  // read the data.
  image->data = NULL;
  image->data = (unsigned char *) malloc(size);
  if (image->data == NULL) {
    printf("Error allocating memory for color-corrected image data");
    return 0;
  }
  printf("Allocated %lu bytes\n", size);
  
  if ((i = fread(image->data, size, 1, file)) != 1) {
    printf("Error reading image data from %s.\n", filename);
    return 0;
  }
  
  
  for (i=0;i<size;i+=3) { // reverse all of the colors. (bgr -> rgb)
    if (image->sizeX == 64 && i < 90)
      printf("reversing bmp: image[%lu] = %02X, %02X, %02X\n", i, image->data[i], image->data[i+1], image->data[i+2]);
    temp = image->data[i];
    image->data[i] = image->data[i+2];
    image->data[i+2] = temp;
    if (image->sizeX == 64 && image->data[i] != image->data[i+2])// && i < 90)
      printf("\tto bmp: image[%lu] = %02X, %02X, %02X\n", i, image->data[i], image->data[i+1], image->data[i+2]);
  }
  
  // we're done.
  return 1;
}



int ImageLoad32(char *filename, Image *image, ImLdOption option) {
  FILE *file;
  unsigned long orig_size;
  unsigned long size;                 // size of the image in bytes.
  unsigned long i;                    // standard counter.
  unsigned short int planes;          // number of planes in image (must be 1)
  unsigned short int bpp;             // number of bits per pixel (must be 24)
  char temp;                          // used to convert bgr to rgb color.
  
  // make sure the file is there.
  if ((file = fopen(filename, "rb"))==NULL) {
    printf("File Not Found : %s\n",filename);
    return 0;
  }
  
  // seek through the bmp header, up to the width/height:
  fseek(file, 18, SEEK_CUR);
  
  // No 100% errorchecking anymore!!!
  
  // read the width
  image->sizeX = getint (file);
  printf("Width of %s: %lu\n", filename, image->sizeX);
  
  // read the height
  image->sizeY = getint (file);
  printf("Height of %s: %lu\n", filename, image->sizeY);
  
  // calculate the size (assuming 32 bits or 4 bytes per pixel).
  size = image->sizeX * image->sizeY * 4;
  if (option == add_alpha_channel) orig_size = image->sizeX * image->sizeY * 3;
  else orig_size = size;
  
  // read the planes
  planes = getshort(file);
  if (planes != 1) {
    printf("Planes from %s is not 1: %u\n", filename, planes);
    return 0;
  }
  
  // read the bpp
  bpp = getshort(file);
  if (bpp != 32) {
    if (bpp != 24 && option != add_alpha_channel) {
      printf("Bpp from %s is not 32: %u\n", filename, bpp);
      printf("\t\tLoad option: %d\n", option);
      return 0;
    }
  }
  
  // seek past the rest of the bitmap header.
  fseek(file, 24, SEEK_CUR);
  
  // read the data.
  image->data = (unsigned char *) malloc(size);
  if (image->data == NULL) {
    printf("Error allocating memory for color-corrected image data");
    return 0;
  }
  
  
  size_t rc;
  unsigned char *cptr = image->data;
  
  switch (option) {
    case remove_black:
      if ((rc = fread(image->data, orig_size, 1, file)) != 1) {
        printf("Error reading image data from %s.\n", filename);
        return 0;
      }
      
      for (i=0;i<size;i+=4) { // reverse all of the colors. (bgra -> rgba)
        temp = image->data[i];
        image->data[i] = image->data[i+2];
        image->data[i+2] = temp;
        
        if (image->data[i] == 0 && image->data[i+1] == 0 && image->data[i+2] == 0) {
          image->data[i+3] = 0;
        } else {
          image->data[i+3] = 0xFF;
        }
      }
      break;
      
    case remove_white:
      if ((rc = fread(image->data, orig_size, 1, file)) != 1) {
        printf("Error reading image data from %s.\n", filename);
        return 0;
      }
      
      for (i=0;i<size;i+=4) { // reverse all of the colors. (bgra -> rgba)
        temp = image->data[i];
        image->data[i] = image->data[i+2];
        image->data[i+2] = temp;
        
        if (image->data[i] == 0xFF && image->data[i+1] == 0xFF && image->data[i+2] == 0xFF) {
          image->data[i+3] = 0;
        } else {
          image->data[i+3] = 0xFF;
        }
      }
      break;
      
    case add_alpha_channel:
      
      for (i=0;i<size;i+=4) { // reverse all of the colors. (bgra -> rgba)
        if ((rc = fread(cptr, 3, 1, file)) != 1) {
          printf("Error reading image data[%lu] from %s.\n", i, filename);
          return 0;
        }
        temp = image->data[i];
        image->data[i] = image->data[i+2];
        image->data[i+2] = temp;
        image->data[i+3] = 0xFF;
        cptr += 4;
      }
      break;
      
    default:
      break;
  }
  // we're done.
  return 1;
}

#ifdef USE_SPRITE
GLvoid LoadGLTextures_Sprite(GLuint *texture_sprite) {
  // Load Texture
  Image *image1;
  
  // allocate space for texture
  image1 = (Image *) malloc(sizeof(Image));
  if (image1 == NULL) {
    printf("Error allocating space for image");
    exit(0);
  }
  
  if (!ImageLoad("/fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v6/images/white_particle.bmp", image1)) {
    exit(1);
  }
  
  // Create Textures
  glGenTextures(3, &texture_sprite[0]);
  
  // texture 1 (poor quality scaling)
  glBindTexture(GL_TEXTURE_2D, texture_sprite[0]);   // 2d texture (x and y size)
  
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); // cheap scaling when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST); // cheap scaling when image smalled than texture
  
  // 2d texture, level of detail 0 (normal), 3 components (red, green, blue), x size from image, y size from image,
  // border 0 (normal), rgb color data, unsigned byte data, and finally the data itself.
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // texture 2 (linear scaling)
  glBindTexture(GL_TEXTURE_2D, texture_sprite[1]);   // 2d texture (x and y size)
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); // scale linearly when image smalled than texture
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // texture 3 (mipmapped scaling)
  glBindTexture(GL_TEXTURE_2D, texture_sprite[2]);   // 2d texture (x and y size)
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST); // scale linearly + mipmap when image smalled than texture
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // 2d texture, 3 colors, width, height, RGB in that order, byte data, and the data.
  gluBuild2DMipmaps(GL_TEXTURE_2D, 3, image1->sizeX, image1->sizeY, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
};
#endif

GLuint LoadPointTextureSprite() {
  // Load Texture
  Image *image1;

  GLuint txID;

  // allocate space for texture
  image1 = (Image *) malloc(sizeof(Image));
  if (image1 == NULL) {
    printf("Error allocating space for image");
    exit(0);
  }

  if (!ImageLoad32(
        "/fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v6/images/white_particle.bmp",
        image1,
        add_alpha_channel)) {
    exit(1);
  }

  // Create Textures
  glGenTextures(1, &txID);
  glBindTexture( GL_TEXTURE_2D, txID);

  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
  // 2d texture, 3 colors, width, height, RGB in that order, byte data, and the data.
  gluBuild2DMipmaps(GL_TEXTURE_2D, 4, image1->sizeX, image1->sizeY, GL_RGBA, GL_UNSIGNED_BYTE, image1->data);

  return txID;
}


// Load Bitmaps And Convert To Textures
//GLvoid LoadGLTextures(GLvoid) {
GLvoid LoadGLTextures(char *filename, GLuint *texture, int size) {
  if (size != 3) {
    printf("LoadTextures: Texture size error\n");
    return;
  }
  // Load Texture
  Image *image1;
  
  // allocate space for texture
  image1 = (Image *) malloc(sizeof(Image));
  if (image1 == NULL) {
    printf("Error allocating space for image");
    exit(0);
  }
  
  if (!ImageLoad(filename, image1)) {
    exit(1);
  }
  
  // Create Textures
  glGenTextures(3, &texture[0]);
  //  printf("Tissue texture: %d, %d, %d\n", texture[0], texture[1], texture[2]);
  
  //  glActiveTexture(GL_TEXTURE0);
  
  // texture 1 (poor quality scaling)
  glBindTexture(GL_TEXTURE_2D, texture[0]);   // 2d texture (x and y size)
  
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); // cheap scaling when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST); // cheap scaling when image smalled than texture
  
  // 2d texture, level of detail 0 (normal), 3 components (red, green, blue), x size from image, y size from image,
  // border 0 (normal), rgb color data, unsigned byte data, and finally the data itself.
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // texture 2 (linear scaling)
  glBindTexture(GL_TEXTURE_2D, texture[1]);   // 2d texture (x and y size)
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); // scale linearly when image smalled than texture
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // texture 3 (mipmapped scaling)
  glBindTexture(GL_TEXTURE_2D, texture[2]);   // 2d texture (x and y size)
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST); // scale linearly + mipmap when image smalled than texture
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // 2d texture, 3 colors, width, height, RGB in that order, byte data, and the data.
  gluBuild2DMipmaps(GL_TEXTURE_2D, 3, image1->sizeX, image1->sizeY, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
};


// Load Bitmaps And Convert To Textures
//GLvoid LoadGLTextures(GLvoid) {
GLvoid LoadGLTexture24(char *filename, GLuint *texture) {

  // Load Texture
  Image *image1;
  
  // allocate space for texture
  image1 = (Image *) malloc(sizeof(Image));
  if (image1 == NULL) {
    printf("Error allocating space for image");
    exit(0);
  }
  
  if (!ImageLoad(filename, image1)) {
    exit(1);
  }
  
  // Create Textures
  glGenTextures(1, texture);
//  //  printf("Tissue texture: %d, %d, %d\n", texture[0], texture[1], texture[2]);
//  
//  //  glActiveTexture(GL_TEXTURE0);
//  
//  // texture 1 (poor quality scaling)
//  glBindTexture(GL_TEXTURE_2D, texture[0]);   // 2d texture (x and y size)
//  
//  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST); // cheap scaling when image bigger than texture
//  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST); // cheap scaling when image smalled than texture
//  
//  // 2d texture, level of detail 0 (normal), 3 components (red, green, blue), x size from image, y size from image,
//  // border 0 (normal), rgb color data, unsigned byte data, and finally the data itself.
//  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
//  
//  // texture 2 (linear scaling)
//  glBindTexture(GL_TEXTURE_2D, texture[1]);   // 2d texture (x and y size)
//  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
//  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); // scale linearly when image smalled than texture
//  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // texture 3 (mipmapped scaling)
  glBindTexture(GL_TEXTURE_2D, *texture);   // 2d texture (x and y size)
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST); // scale linearly + mipmap when image smalled than texture
  glTexImage2D(GL_TEXTURE_2D, 0, 3, image1->sizeX, image1->sizeY, 0, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
  
  // 2d texture, 3 colors, width, height, RGB in that order, byte data, and the data.
  gluBuild2DMipmaps(GL_TEXTURE_2D, 3, image1->sizeX, image1->sizeY, GL_RGB, GL_UNSIGNED_BYTE, image1->data);
};


GLvoid LoadGLTexture(char *filename, GLuint *texture, ImLdOption option) {
  // Load Texture
  Image *image1;
  
  // allocate space for texture
  image1 = (Image *) malloc(sizeof(Image));
  if (image1 == NULL) {
    printf("Error allocating space for image");
    exit(0);
  }
  
  if (!ImageLoad32(filename, image1, option)) {
    exit(1);
  }
  
  // Create Textures
  glGenTextures(1, texture);
  
  // texture 2 (linear scaling)
  glBindTexture(GL_TEXTURE_2D, *texture);   // 2d texture (x and y size)
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR); // scale linearly when image smalled than texture
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
  glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
  glTexImage2D(GL_TEXTURE_2D, 0, 4, image1->sizeX, image1->sizeY, 0, GL_RGBA, GL_UNSIGNED_BYTE, image1->data);
  
//  // texture 3 (mipmapped scaling)
//  glBindTexture(GL_TEXTURE_2D, texture_scar[2]);   // 2d texture (x and y size)
//  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR); // scale linearly when image bigger than texture
//  glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR_MIPMAP_NEAREST); // scale linearly + mipmap when image smalled than texture
//  glTexImage2D(GL_TEXTURE_2D, 0, 4, image1->sizeX, image1->sizeY, 0, GL_RGBA, GL_UNSIGNED_BYTE, image1->data);
//  
//  // 2d texture, 3 colors, width, height, RGB in that order, byte data, and the data.
//  gluBuild2DMipmaps(GL_TEXTURE_2D, 4, image1->sizeX, image1->sizeY, GL_RGBA, GL_UNSIGNED_BYTE, image1->data);
};

