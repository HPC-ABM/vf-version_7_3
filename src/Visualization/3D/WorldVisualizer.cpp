/*
 * WorldVisualizer.cpp
 * Created by Nuttiiya on 11/30/2016
 *
 * Copyright (c) 2016 Nuttiiya Seekhao. All rights reserved
 *
 **/

#include "WorldVisualizer.h"

WorldVisualizer::WorldVisualizer()
: m_GLVertexBufferCap(0)
, m_GLVertexBufferDam(0)
, m_GLColorBufferCap(0)
, m_GLColorBufferDam(0)
{}

WorldVisualizer::~WorldVisualizer()
{
  DestroyWorldVisualizer();
}
inline void DeleteVertexBuffer( GLuint& vboID )
{
  if ( vboID != 0 )
  {
    glDeleteBuffersARB( 1, &vboID );
    vboID = 0;
  }
}

inline void CreateVertexBuffer( GLuint& vboID )
{
  // Make sure we don't loose the reference to the previous VBO if there is one
  DeleteVertexBuffer( vboID );
  glGenBuffersARB( 1, &vboID );
}


void WorldVisualizer::DestroyWorldVisualizer()
{
  DeleteVertexBuffer(m_GLVertexBufferCap);
  DeleteVertexBuffer(m_GLColorBufferCap);
  DeleteVertexBuffer(m_GLVertexBufferDam);
  DeleteVertexBuffer(m_GLColorBufferDam);
}

void WorldVisualizer::LoadWorldEnv(Patch* worldPatch, int width, int height, int depth, GLuint txID)
{

  if (worldPatch == NULL || width < 0 || height < 0 || depth < 0) {
    printf("WorldVisualizer::LoadWorldEnv(%p, %d, %d, %d) failed\n",
      worldPatch, width, height, depth);
    exit(-1);
  }

  this->pointtxID = txID;

//  this->numVerts = width * height * depth;
//  m_PositionBuffer.resize( numVerts );
//  m_ColorBuffer.resize( numVerts );

  printf("width: %d\theight: %d\tdepth: %d\n", width, height, depth);
  this->patchPtr = worldPatch;
  this->width    = width;
  this->height   = height;
  this->depth    = depth;

  int index = 0;
  for (int z = 0; z < depth; z++) {
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        GLuint index = z*width*height + y*width + x;
        int  patchType = worldPatch[index].getType();
        bool isDamZone = worldPatch[index].isInDamZone();
	bool isDamaged = worldPatch[index].isDamaged();
       
        if (patchType == capillary){
          float X = (float) GL_X(x);
          float Y = (float) GL_Y(y);
          float Z = (float) PLANE_DEPTH + z;
          m_ColorBufferCap.push_back(glm::vec4(1.0f, 0.0f, 0.0f, 0.6f)); //red
          m_PositionBufferCap.push_back(glm::vec3( X, Y, Z ));
	   
        } else if (isDamZone) {
          float X = (float) GL_X(x);
          float Y = (float) GL_Y(y);
          float Z = (float) PLANE_DEPTH + z;
          m_PatchIndexBufferDam.push_back(index);
          m_PositionBufferDam.push_back(glm::vec3( X, Y, Z ));
          if (!isDamaged)
            m_ColorBufferDam.push_back(glm::vec4(1.0f, 1.0f, 1.0f, 1.0f)); // void
          else {
#ifdef HUMAN_VF
            m_ColorBufferDam.push_back(glm::vec4(0.687255f, 0.510784f, 0.510784f, 0.009f)); //dark pink
#elif defined(RAT_VF)
            m_ColorBufferDam.push_back(glm::vec4(1.0f, 1.0f, 1.0f, 0.1f)); //dark pink
#else
            m_ColorBufferDam.push_back(glm::vec4(0.687255f, 0.510784f, 0.510784f, 0.009f)); //dark pink
#endif
          }
        }
      }
    }
  }

  this->numVertsCap = m_PositionBufferCap.size();
  this->numVertsDam = m_PositionBufferDam.size();


  printf("World Environment has been loaded!\n");
//exit(-1);
  this->GenerateVertexBuffer();
}

void WorldVisualizer::BufferDam()
{
  for (int i = 0; i < this->numVertsDam; i++) {
    int  isDamaged = (this->patchPtr)[m_PatchIndexBufferDam[i]].isDamaged();
    if (!isDamaged)
    {
      m_ColorBufferDam.at(i) = glm::vec4(0.0f, 0.0f, 0.0f, 0.0f); // void
    }
    else
    {
#ifdef HUMAN_VF
      m_ColorBufferDam.at(i) = glm::vec4(0.687255f, 0.510784f, 0.510784f, 0.009f); //dark pink
#elif defined(RAT_VF)
      m_ColorBufferDam.at(i) = glm::vec4(1.0f, 1.0f, 1.0f, 0.1f); //dark pink
#else
      m_ColorBufferDam.at(i) = glm::vec4(0.687255f, 0.510784f, 0.510784f, 0.009f); //dark pink
#endif
    }
  }

  glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, sizeof(glm::vec4) * numVertsDam,
                     &(m_ColorBufferDam[0]));
//  glBufferDataARB( GL_ARRAY_BUFFER_ARB, sizeof(glm::vec4) * m_ColorBufferDam.size(),
//                   &(m_ColorBufferDam[0]), GL_STATIC_DRAW_ARB );
}


void WorldVisualizer::GenerateVertexBuffer()
{
  if (this->width < 0 || this->height < 0 || this->depth < 0 || this->patchPtr == NULL)
  {
    fprintf(stderr, "Error: WorldVisualizer::GenerateVertexBuffer(): Environment not loaded\n");
    exit(-1);
  }

  // First generate the buffer object ID's
  CreateVertexBuffer(m_GLVertexBufferCap);
  CreateVertexBuffer(m_GLColorBufferCap);

  CreateVertexBuffer(m_GLVertexBufferDam);
  CreateVertexBuffer(m_GLColorBufferDam);

  // Copy the host data into the vertex buffer objects
  glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_GLVertexBufferCap );
  glBufferDataARB( GL_ARRAY_BUFFER_ARB, sizeof(glm::vec3) * m_PositionBufferCap.size(), &(m_PositionBufferCap[0]), GL_STATIC_DRAW_ARB );

  glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_GLColorBufferCap );
  glBufferDataARB( GL_ARRAY_BUFFER_ARB, sizeof(glm::vec4) * m_ColorBufferCap.size(), &(m_ColorBufferCap[0]), GL_STATIC_DRAW_ARB );


  glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_GLVertexBufferDam );
  glBufferDataARB( GL_ARRAY_BUFFER_ARB, sizeof(glm::vec3) * m_PositionBufferDam.size(), &(m_PositionBufferDam[0]), GL_STATIC_DRAW_ARB );

  glBindBufferARB( GL_ARRAY_BUFFER_ARB, m_GLColorBufferDam );
  glBufferDataARB( GL_ARRAY_BUFFER_ARB, sizeof(glm::vec4) * m_ColorBufferDam.size(), &(m_ColorBufferDam[0]), GL_DYNAMIC_DRAW_ARB );

}

void WorldVisualizer::startSprite()
{
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, this->pointtxID);
    glEnable( GL_BLEND );
    glBlendFunc( GL_SRC_ALPHA, GL_ONE );
    float quadratic[] =  { 0.0001f, 0.005f, 0.001f };
    glPointParameterfvARB( GL_POINT_DISTANCE_ATTENUATION_ARB, quadratic );

    float maxSize = 0.0f;
    glGetFloatv( GL_POINT_SIZE_MAX_ARB, &maxSize );
    glPointSize( maxSize );

    glPointParameterfARB( GL_POINT_FADE_THRESHOLD_SIZE_ARB, 10.0f );

    glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, 1.0f );
    glPointParameterfARB( GL_POINT_SIZE_MAX_ARB, /*maxSize*/100000 );
    glTexEnvf( GL_POINT_SPRITE_ARB, GL_COORD_REPLACE_ARB, GL_TRUE );
    glEnable( GL_POINT_SPRITE_ARB );

}

void endSprite()
{
    glDisable( GL_POINT_SPRITE_ARB );
    glDisable(GL_TEXTURE_2D);
    glDisable( GL_BLEND );
}

void WorldVisualizer::Render()
{

  /********************************
   * Capillary rendering          *
   ********************************/
/*
//    this->startSprite();
//    glColor3f(0.0f, 1.0f, 0.0f);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_DEPTH_TEST);

    glPointSize(0.05f);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_GLVertexBufferCap);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_GLColorBufferCap);
    glColorPointer(4, GL_FLOAT, 0, 0);
    
    
    glDrawArrays(GL_POINTS, 0, numVertsCap);
    
    glDisableClientState(GL_VERTEX_ARRAY); 
    glDisableClientState(GL_COLOR_ARRAY);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

//    endSprite();
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST); // Re-enable writing to depth buffer
*/
  /********************************
   * Damage rendering             *
   ********************************/
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDisable(GL_DEPTH_TEST);

    glPointSize(2.0f);
    glEnableClientState(GL_VERTEX_ARRAY);
    glEnableClientState(GL_COLOR_ARRAY);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_GLVertexBufferDam);
    glVertexPointer(3, GL_FLOAT, 0, 0);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, m_GLColorBufferDam);
    this->BufferDam();
//    glBufferSubDataARB(GL_ARRAY_BUFFER_ARB, 0, numVerts * 4 * sizeof(float), &(m_ColorBuffer[0]));
    glColorPointer(4, GL_FLOAT, 0, 0);
    
    
    glDrawArrays(GL_POINTS, 0, numVertsDam);
    
    glDisableClientState(GL_VERTEX_ARRAY); 
    glDisableClientState(GL_COLOR_ARRAY);
    
    glBindBufferARB(GL_ARRAY_BUFFER_ARB, 0);

//    endSprite();
  glDisable(GL_BLEND);
  glEnable(GL_DEPTH_TEST); // Re-enable writing to depth buffer
  glPointSize(1.0f);
}



