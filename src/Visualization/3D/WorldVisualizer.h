#ifndef WorldVisualizer_h
#define WorldVisualizer_h

#include "common_vis.h"
#include "../../Patch/Patch.h"

#pragma once

#define BUFFER_OFFSET(i) ((char*)NULL + (i))

class WorldVisualizer
{
public:
  WorldVisualizer();
  virtual ~WorldVisualizer();

  void LoadWorldEnv(Patch* worldPatch, int width, int height, int depth, GLuint txID);
  // TODO: Update functions

  void Render();
  void DestroyWorldVisualizer();
private:
  typedef std::vector<glm::vec3>  PositionBuffer;
  typedef std::vector<glm::vec4>  ColorBuffer;
  typedef std::vector<GLuint>     IndexBuffer;

  /* Capillaries */
  PositionBuffer  m_PositionBufferCap;
  ColorBuffer     m_ColorBufferCap;
  IndexBuffer     m_IndexBufferCap;

  // ID's for the VBO's
  GLuint  m_GLVertexBufferCap;
  GLuint  m_GLColorBufferCap;
  GLuint  m_GLIndexBufferCap;

  GLint numVertsCap;
  GLuint pointtxID;


  /* Damages */
  PositionBuffer  m_PositionBufferDam;
  ColorBuffer     m_ColorBufferDam;
  IndexBuffer     m_IndexBufferDam;
  IndexBuffer     m_PatchIndexBufferDam;

  // ID's for the VBO's
  GLuint  m_GLVertexBufferDam;
  GLuint  m_GLColorBufferDam;
  GLuint  m_GLIndexBufferDam;

  GLint numVertsDam;


  GLint width;
  GLint height;
  GLint depth;

 
  Patch *patchPtr;

  void startSprite();
  void GenerateVertexBuffer();
  void BufferDam();
};

#endif	// WorldVisualizer_h


