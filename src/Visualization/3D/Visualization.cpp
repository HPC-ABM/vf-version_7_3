//
//  Visualization.cpp
//  OpenGL_SimpleWH
//
//  Created by NungnunG on 12/28/15.
//  Copyright (c) 2015 NungnunG. All rights reserved.
//

#include "Visualization.h"

using namespace visualization;

WHWorld* visualization::world_ptr = NULL;
int visualization::iter           = 0;

int visualization::argc           = 0;
char** visualization::argv        = NULL;

#ifdef OVERLAP_VIS
bool visualization::visDone = false;
#endif

float visualization::totalChem[MAX_PLOTS] = {0.0};
float visualization::totalCell[MAX_PLOTS] = {0.0};
float visualization::xcam      = 0;
float visualization::ycam      = 0;
float visualization::zcam      = 0;

WorldVisualizer    *visualization::EnvVisualizer_ptr    = NULL;
//Terrain             *visualization::ScarTerrain_ptr      = NULL;
//SubPlane            *visualization::ECMplane_ptr         = NULL;
HeatMapManager      *visualization::HeatMapManager_ptr   = NULL;
PlotManager         *visualization::ChemPlotManager_ptr      = NULL;
PlotManager         *visualization::CellPlotManager_ptr      = NULL;
//GridSurfaceManager  *visualization::ChemSurfManager_ptr  = NULL;
OptionTable         *visualization::ChemOptionTable_ptr  = NULL;
OptionTable         *visualization::CellOptionTable_ptr  = NULL;
//VolumeManager       *visualization::ECMvolumeManager_Ptr = NULL;

//GLuint  visualization::tissue_texture[3] = {0};	/* Storage for 3 textures. */
//GLuint  visualization::glass_texture[3]  = {0};
//GLuint  visualization::texture_atlas     =  0 ;	/* Storage for 1 texture. */
//GLuint  visualization::scar_texture      =  0 ;

#ifdef USE_SPRITE
GLuint  visualization::texture_sprite[3] = {0};	/* Storage for 3 textures. */
#endif

#ifdef HUMAN_VF
float4	particleColor = make_float4(MAXCOLOR, MINCOLOR, MINCOLOR, 0.25f);//holds color for particle
float4	particleColor2 = make_float4(MAXCOLOR2, MINCOLOR2, MINCOLOR2, 0.2f);//holds color for particle
float4	particleColor3 = make_float4(0.3f, 0.1f, 0.9f, 0.15f);//holds color for particle

#elif defined(RAT_VF)
float4	particleColor2 = make_float4(MAXCOLOR, MINCOLOR+0.4, MINCOLOR+0.2, 0.9f);//holds color for particle
float4	particleColor = make_float4(MAXCOLOR2+0.3, MINCOLOR2, MINCOLOR2, 0.9f);//holds color for particle
float4	particleColor3 = make_float4(0.3f, 0.1f, 0.9f, 0.45f);//holds color for particle

#else
float4	particleColor = make_float4(MAXCOLOR, MINCOLOR, MINCOLOR, 0.25f);//holds color for particle
float4	particleColor2 = make_float4(MAXCOLOR2, MINCOLOR2, MINCOLOR2, 0.2f);//holds color for particle
float4	particleColor3 = make_float4(0.3f, 0.1f, 0.9f, 0.15f);//holds color for particle

#endif

// Text labels for table:

char* visualization::optionHeader[CHEM_OPTION_COLS] = {"Heat Map"};//, "3D Surface"};
char* visualization::optionLabel[CHEM_OPTION_ROWS] = {"TNF", "TGF", "FGF", "MMP8",
                                        "IL-1beta", "IL-6", "IL-8", "IL-10"};
char* visualization::optionHeaderCell[CELL_OPTION_COLS] = {"Neutrophil", "Macrophage", "Fibroblast"};
char* visualization::optionLabelCell[CELL_OPTION_ROWS] = {"Cells"};

// Color for table elements

float visualization::surfColors[CHEM_OPTION_ROWS * CHEM_OPTION_COLS][4] = {
                                 {0.9f, 0.9f, 0.95f, 1.0f},
                                 {0.6f, 0.6f, 0.82f, 1.0f},
                                 {0.3f, 0.3f, 0.70f, 1.0f},
                                 {0.6f, 0.3f, 0.70f, 1.0f},
                                 {0.9f, 0.3f, 0.70f, 1.0f},
                                 {0.9f, 0.6f, 0.30f, 1.0f},
                                 {0.3f, 0.6f, 0.30f, 1.0f},
                                 {0.3f, 0.9f, 0.30f, 1.0f}};


float visualization::cellColors[CELL_OPTION_COLS][4]= {
  {particleColor.x,  particleColor.y,  particleColor.z,  particleColor.w},
  {particleColor2.x, particleColor2.y, particleColor2.z, particleColor2.w},
  {particleColor3.x, particleColor3.y, particleColor3.z, particleColor3.w}
};



/*
 * 1[x, y, r, g, b]-------2[x, y, r, g, b]
 *        |                       |
 *        |                       |
 *        |                       |
 * 0[x, y, r, g, b]-------3[x, y, r, g, b]
 */
GLint* visualization::WG_vert = NULL;
GLint visualization::WG_xb = 0;
GLint visualization::WG_yb = 0;
GLint visualization::WG_zb = 0;
GLint visualization::WG_w = 0;
GLint visualization::WG_h = 0;
GLint visualization::WG_d = 0;
GLint visualization::WG_numvert = 0;


void visualization::init(int argc, char** argv, WHWorld *world_ptr)
{
  if (!world_ptr) {
	fprintf(stderr, "visualization::init() NULL world_ptr\n");
	exit(-1);
  }

  xcam = 0.0;
  ycam = 0.0;
  zcam = 0.0;

  visualization::world_ptr = world_ptr;
  visualization::iter      = 0;
  visualization::argc      = argc;
  visualization::argv      = argv;

//  visualization::totalChem = (float *) malloc(world_ptr->typesOfBaseChem * sizeof(float));
//  if (!visualization::totalChem) {
//	fprintf(stderr, "visualization::init() fail to allocate totaChem\n");
//	exit(-1);
//  } else {
//	fprintf(stderr, "Allocated %d entries in visualization::totalChem\n",
//		world_ptr->typesOfBaseChem);
//  }
}



void visualization::initWndGrid(){
  int wnd_xb, wnd_xe, wnd_yb, wnd_ye, wnd_zb, wnd_ze;
  world_ptr->getWndPos(wnd_xb, wnd_xe, wnd_yb, wnd_ye, wnd_zb, wnd_ze);
  
  
//  wnd_xb -= WG_BORDER_X;
//  wnd_xe += WG_BORDER_X;
//  wnd_ye += WG_BORDER_Y;
  
  int gw = wnd_xe - wnd_xb;
  int gh = wnd_ye - wnd_yb;
  int gd = wnd_ze - wnd_zb;
  
  WG_w = gw;
  WG_h = gh;
  WG_d = gd;
  WG_xb = wnd_xb;
  WG_yb = wnd_yb;
  WG_zb = wnd_zb;
  WG_numvert = gw * gh * gd;
}

#ifdef OVERLAP_VIS
void visualization::notifyVisDone()
{
  visualization::visDone = true;
}

void visualization::notifyComDone()
{
  visualization::visDone = false;
}

bool isVisDone()
{
  return visualization::visDone;
}

bool visualization::waitForVis()
{
  bool flag = !visualization::visDone;
  while (flag){
    flag = !isVisDone();
  }
  return flag;
}

bool visualization::waitForCompute()
{
  bool flag = visualization::visDone;
  while (flag){
    flag = isVisDone();
  }
  return flag;
}
#endif	// OVERLAP_VIS


const char *volumeFilename = "/fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v6/src/Visualization/3D/data/Bucky.raw";
cudaExtent volumeSize = make_cudaExtent(32, 32, 32);

cudaMemcpy3DParms copyParams = {0};

// Define the files that are to be save and the reference images for validation
const char *sOriginal[] =
{
    "volume.ppm",
    NULL
};

const char *sReference[] =
{
    "ref_volume.ppm",
    NULL
};


//char *volumeFilename = "/fs/HPC_ABMs/GitProjects/vocalcord-cpuabm-v6/src/Visualization/3D/data/mrt16_angio.raw";
//cudaExtent volumeSize = make_cudaExtent(416, 512, 112);


uint width = 800, height = 800;
//uint width = 512, height = 512;
dim3 blockSize(16, 16);
dim3 gridSize;

float3 viewRotation;
float3 viewTranslation = make_float3(0.0, 0.0, -4.0f);
float invViewMatrix[12];

float density = 0.5f;//0.2f;//0.05f
float brightness = 1.0f;
float transferOffset = 0.0f;
float transferScale = 1.0f;
bool linearFiltering = true;

GLuint pbo = 0;     // OpenGL pixel buffer object
GLuint tex = 0;     // OpenGL texture object
struct cudaGraphicsResource *cuda_pbo_resource; // CUDA Graphics Resource (to transfer PBO)

StopWatchInterface *timer = 0;

// Auto-Verification Code
const int frameCheckNumber = 2;
int fpsCount = 0;        // FPS count for averaging
int fpsLimit = 1;        // FPS limit for sampling
int g_Index = 0;
unsigned int frameCount = 0;

int *pArgc;
char **pArgv;


void initPixelBuffer();

void computeFPS()
{
    frameCount++;
    fpsCount++;

    if (fpsCount == fpsLimit)
    {
        char fps[256];
        float ifps = 1.f / (sdkGetAverageTimerValue(&timer) / 1000.f);
//        sprintf(fps, "Volume Render: %3.1f fps", ifps);

        glutSetWindowTitle(fps);
        fpsCount = 0;

        fpsLimit = (int)MAX(1.f, ifps);
        sdkResetTimer(&timer);
    }
}

// render image using CUDA
void render()
{

    copyInvViewMatrix(invViewMatrix, sizeof(float4)*3);

    // map PBO to get CUDA device pointer
    uint *d_output;
    // map PBO to get CUDA device pointer
    checkCudaErrors(cudaGraphicsMapResources(1, &cuda_pbo_resource, 0));
    size_t num_bytes;
    checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void **)&d_output, &num_bytes,
                                                         cuda_pbo_resource));
    //printf("CUDA mapped PBO: May access %ld bytes\n", num_bytes);

    // clear image
    checkCudaErrors(cudaMemset(d_output, 0, width*height*4));

    // call CUDA kernel, writing results to PBO
//    render_kernel(gridSize, blockSize, d_output, width, height, density, brightness, transferOffset, transferScale);

    render_kernel_dim(
    		gridSize,
    		blockSize,
    		d_output,
    		volumeSize.width,
    		volumeSize.height,
    		volumeSize.depth,
    		width, height,
    		density, brightness, transferOffset, transferScale);


    getLastCudaError("kernel failed");

    checkCudaErrors(cudaGraphicsUnmapResources(1, &cuda_pbo_resource, 0));
}

// display results using OpenGL (called by GLUT)
void display()
{

    sdkStartTimer(&timer);

    // use OpenGL to build view matrix
    GLfloat modelView[16];
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    glRotatef(-viewRotation.x, 1.0, 0.0, 0.0);
    glRotatef(-viewRotation.y, 0.0, 1.0, 0.0);
    glTranslatef(-viewTranslation.x, -viewTranslation.y, -viewTranslation.z);
    glGetFloatv(GL_MODELVIEW_MATRIX, modelView);
    glPopMatrix();

//  	glColor3f(1,1,1);
  //    drawTransparentBox(GL_X(nx/2) +  fibOffset, GL_Y(ny/2) +  fibOffsetY, PLANE_DEPTH + nz/2, nx, ny, nz);
//    drawTransparentBox(GL_X(142/2), GL_Y(142/2), -2 + 28/2, 142, 142, 28);


    invViewMatrix[0] = modelView[0];
    invViewMatrix[1] = modelView[4];
    invViewMatrix[2] = modelView[8];
    invViewMatrix[3] = modelView[12];
    invViewMatrix[4] = modelView[1];
    invViewMatrix[5] = modelView[5];
    invViewMatrix[6] = modelView[9];
    invViewMatrix[7] = modelView[13];
    invViewMatrix[8] = modelView[2];
    invViewMatrix[9] = modelView[6];
    invViewMatrix[10] = modelView[10];
    invViewMatrix[11] = modelView[14];

    render();

    // for white bg
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);           // This Will Clear The Background Color To White
    // display results
    glClear(GL_COLOR_BUFFER_BIT);


    // draw image from PBO
    glDisable(GL_DEPTH_TEST);

    // for white bg
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
#if 0
    // draw using glDrawPixels (slower)
    glRasterPos2i(0, 0);
    glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
    glDrawPixels(width, height, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);
#else
    // draw using texture

    // copy from pbo to texture
    glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

    // draw textured quad
    glEnable(GL_TEXTURE_2D);
    glBegin(GL_QUADS);
    glTexCoord2f(0, 0);
    glVertex2f(0, 0);
    glTexCoord2f(1, 0);
    glVertex2f(1, 0);
    glTexCoord2f(1, 1);
    glVertex2f(1, 1);
    glTexCoord2f(0, 1);
    glVertex2f(0, 1);
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, 0);
#endif




    glutSwapBuffers();
    glutReportErrors();

    sdkStopTimer(&timer);

    computeFPS();

    char winTitle[256];
    sprintf(winTitle, "Collagen: Iteration %d", visualization::iter);

    glutSetWindowTitle(winTitle);

    if (!paused)
    {
    	/*****************************************/
#ifdef OVERLAP_VIS
    	visualization::notifyVisDone();
    	cout << "visualization::start() waiting for compute ..." << endl;
    	visualization::waitForCompute();
    	cout << "visualization::start() done waiting for compute" << endl;
#else	// OVERLAP_VIS
    	// Update new locations
    	visualization::world_ptr->go();
#endif	// OVERLAP_VIS
    	/*****************************************/

    	bufferECMmap(copyParams);
    	visualization::iter++;
    }
}

void idle()
{
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
    switch (key)
    {
        case 27:
            #if defined (__APPLE__) || defined(MACOSX)
                exit(EXIT_SUCCESS);
            #else
                glutDestroyWindow(glutGetWindow());
                return;
            #endif
            break;

        case 'f':
            linearFiltering = !linearFiltering;
            setTextureFilterMode(linearFiltering);
            break;

        case '+':
            density += 0.01f;
            break;

        case '-':
            density -= 0.01f;
            break;

        case ']':
            brightness += 0.1f;
            break;

        case '[':
            brightness -= 0.1f;
            break;

        case ';':
            transferOffset += 0.01f;
            break;

        case '\'':
            transferOffset -= 0.01f;
            break;

        case '.':
            transferScale += 0.01f;
            break;

        case ',':
            transferScale -= 0.01f;
            break;

        default:
            break;
    }

    printf("density = %.2f, brightness = %.2f, transferOffset = %.2f, transferScale = %.2f\n", density, brightness, transferOffset, transferScale);
    glutPostRedisplay();
}

int ox, oy;
int buttonState = 0;

void mouse(int button, int state, int x, int y)
{
    if (state == GLUT_DOWN)
    {
        buttonState  |= 1<<button;
    }
    else if (state == GLUT_UP)
    {
        buttonState = 0;
    }

    ox = x;
    oy = y;
    glutPostRedisplay();
}

void motion(int x, int y)
{
    float dx, dy;
    dx = (float)(x - ox);
    dy = (float)(y - oy);

    if (buttonState == 4)
    {
        // right = zoom
        viewTranslation.z += dy / 100.0f;
    }
    else if (buttonState == 2)
    {
        // middle = translate
        viewTranslation.x += dx / 100.0f;
        viewTranslation.y -= dy / 100.0f;
    }
    else if (buttonState == 1)
    {
        // left = rotate
        viewRotation.x += dy / 5.0f;
        viewRotation.y += dx / 5.0f;
    }

    ox = x;
    oy = y;
    glutPostRedisplay();
}

//int iDivUp(int a, int b)
//{
//    return (a % b != 0) ? (a / b + 1) : (a / b);
//}

void reshape(int w, int h)
{
    width = w;
    height = h;
    initPixelBuffer();

    // calculate new grid size
    gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y));

    glViewport(0, 0, w, h);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0, 1.0, 0.0, 1.0, 0.0, 1.0);
}

void cleanup()
{
    sdkDeleteTimer(&timer);

    freeCudaBuffers();

    if (pbo)
    {
        cudaGraphicsUnregisterResource(cuda_pbo_resource);
        glDeleteBuffersARB(1, &pbo);
        glDeleteTextures(1, &tex);
    }
    // cudaDeviceReset causes the driver to clean up all state. While
    // not mandatory in normal operation, it is good practice.  It is also
    // needed to ensure correct operation when the application is being
    // profiled. Calling cudaDeviceReset causes all profile data to be
    // flushed before the application exits
    cudaDeviceReset();
}

void initGL(int *argc, char **argv)
{
    // initialize GLUT callback functions
    glutInit(argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutInitWindowSize(width, height);
    glutCreateWindow("CUDA volume rendering");

    glewInit();

    if (!glewIsSupported("GL_VERSION_2_0 GL_ARB_pixel_buffer_object"))
    {
        printf("Required OpenGL extensions missing.");
        exit(EXIT_SUCCESS);
    }
}

void initPixelBuffer()
{
    if (pbo)
    {
        // unregister this buffer object from CUDA C
        checkCudaErrors(cudaGraphicsUnregisterResource(cuda_pbo_resource));

        // delete old buffer
        glDeleteBuffersARB(1, &pbo);
        glDeleteTextures(1, &tex);
    }

    // create pixel buffer object for display
    glGenBuffersARB(1, &pbo);
    glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, pbo);
    glBufferDataARB(GL_PIXEL_UNPACK_BUFFER_ARB, width*height*sizeof(GLubyte)*4, 0, GL_STREAM_DRAW_ARB);
    glBindBufferARB(GL_PIXEL_UNPACK_BUFFER_ARB, 0);

    // register this buffer object with CUDA
    checkCudaErrors(cudaGraphicsGLRegisterBuffer(&cuda_pbo_resource, pbo, cudaGraphicsMapFlagsWriteDiscard));

    // create texture for display
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, NULL);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Load raw data from disk
void *loadRawFile(char *filename, size_t size)
{
    FILE *fp = fopen(filename, "rb");

    if (!fp)
    {
        fprintf(stderr, "Error opening file '%s'\n", filename);
        return 0;
    }

    void *data = malloc(size);
    size_t read = fread(data, 1, size, fp);
    fclose(fp);

#if defined(_MSC_VER_)
    printf("Read '%s', %Iu bytes\n", filename, read);
#else
    printf("Read '%s', %zu bytes\n", filename, read);
#endif

    return data;
}

// General initialization call for CUDA Device
int chooseCudaDevice(int argc, const char **argv, bool bUseOpenGL)
{
    int result = 0;

    if (bUseOpenGL)
    {
        result = findCudaGLDevice(argc, argv);
    }
    else
    {
        result = findCudaDevice(argc, argv);
    }

    return result;
}



void visualization::start()
{

//  char *ref_file = NULL;

#if defined(__linux__)
  setenv ("DISPLAY", ":0", 0);
#endif

  volumeSize.width  = world_ptr->nx;
  volumeSize.height = world_ptr->ny;
  volumeSize.depth  = world_ptr->nz;

  // First initialize OpenGL context, so we can properly set the GL for CUDA.
  // This is necessary in order to achieve optimal performance with OpenGL/CUDA interop.
  initGL(&argc, argv);

//  // load volume data
//  char *path = sdkFindFilePath(volumeFilename, argv[0]);
//
//  if (path == 0)
//  {
//      printf("Error finding file '%s'\n", volumeFilename);
//      exit(EXIT_FAILURE);
//  }

  size_t size = volumeSize.width*volumeSize.height*volumeSize.depth*sizeof(VolumeType);
//  void *h_volume = world_ptr->ecmMap[m_col];//loadRawFile(path, size);
  void *h_volume = world_ptr->ecmPreProcMap[m_col];

  initCuda(h_volume, volumeSize, copyParams);
  //free(h_volume);

  sdkCreateTimer(&timer);

  printf("Press '+' and '-' to change density (0.01 increments)\n"
         "      ']' and '[' to change brightness\n"
         "      ';' and ''' to modify transfer function offset\n"
         "      '.' and ',' to modify transfer function scale\n\n");

  // calculate new grid size
  gridSize = dim3(iDivUp(width, blockSize.x), iDivUp(height, blockSize.y));

  // This is the normal rendering path for VolumeRender
  glutDisplayFunc(display);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutReshapeFunc(reshape);
  glutIdleFunc(idle);

  initPixelBuffer();

//#if defined (__APPLE__) || defined(MACOSX)
//  atexit(cleanup);
//#else
//  glutCloseFunc(cleanup);
//#endif

  glutMainLoop();


//  /* Initialize GLUT state - glut will take any command line arguments that pertain to it or
// *  *    X Windows - look at its documentation at http://reality.sgi.com/mjk/spec3/spec3.html */
//  glutInit(&argc, argv);
//
////  printf("executed glutInit()\n");
//
// /* Select type of Display mode:
// *  Double buffer
// *  RGBA color
// *  Alpha components supported
// *  Depth buffered for automatic clipping */
//  glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_ALPHA | GLUT_DEPTH | GLUT_MULTISAMPLE);
//
//  /* get a 640 x 480 window */
//  glutInitWindowSize(WINW, WINH);
//
//  /* the window starts at the upper left corner of the screen */
//  glutInitWindowPosition(0, 0);
//
//  /* Open a window */
//  window = glutCreateWindow("Wound Healing ABM");
//
//  GLenum res = glewInit();
//  if (res != GLEW_OK)
//  {
//    fprintf(stderr, "Error: '%s'\n", glewGetErrorString(res));
//    exit(-1);
//  }
//
//  /* Register the function to do all our OpenGL drawing. */
//  glutDisplayFunc(&DrawGLScene);
//
//  /* Go fullscreen.  This is as soon as possible. */
//  glutFullScreen();
//
//  /* Even if there are no events, redraw our gl scene. */
//  glutIdleFunc(&DrawGLScene);
//
//  /* Register the function called when our window is resized. */
//  glutReshapeFunc(&ReSizeGLScene);
//
//  /* Register the function called when the keyboard is pressed. */
//  glutKeyboardFunc(&keyPressed);
//  glutSpecialFunc(&SpecialInput);
//
//  /* Initialize our window. */
//  InitGL(WINW, WINH);
////      InitGL(640, 480);
//  /* Set scar terrain ptr in myWorld (Need to be called AFTER InitGL) */
////  world_ptr->setScarTerrainPtr(ScarTerrain_ptr);
//  /* Register mouse function */
//  glutMouseFunc(MouseButton);
//#ifdef USE_MOUSE_ONLY
//  glutMotionFunc(MouseMotion);
//  glutPassiveMotionFunc(MousePassiveMotion);
//#endif
//
//
//  /* Start Event Processing Engine */
//  glutMainLoop();
}

/*----------------------------------------------------------------------------------------
 *
 */
/* A general OpenGL initialization function.  Sets all of the initial parameters. */
void visualization::InitGL(int Width, int Height)	        // We call this right after our OpenGL window is created.
{ 
  // Get world dimension
  GridW = world_ptr->nx;
  GridH = world_ptr->ny;
  GridD = world_ptr->nz;
  
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);           // This Will Clear The Background Color To Black
  //glClearColor(1.0f, 1.0f, 1.0f, 1.0f);           // This Will Clear The Background Color To
  glClear ( GL_COLOR_BUFFER_BIT ) ;
  glClearDepth(1.0);                              // Enables Clearing Of The Depth Buffer
  glDepthFunc(GL_LESS);                           // The Type Of Depth Test To Do
  glEnable(GL_DEPTH_TEST);                        // Enables Depth Testing
  glShadeModel(GL_SMOOTH);                        // Enables Smooth Color Shading

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();                               // Reset The Projection Matrix

  // Calculate The Aspect Ratio Of The Window
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,PERSPECTIVE_DEPTH);//100.0f);

  glMatrixMode(GL_MODELVIEW); 

  // load point texture
#ifdef USE_SPRITE
  LoadGLTextures_Sprite(texture_sprite);
#endif

  GLuint txID = LoadPointTextureSprite();
  /* Create and initialize world environment visualizer */
  visualization::EnvVisualizer_ptr = new WorldVisualizer();
  EnvVisualizer_ptr->LoadWorldEnv(world_ptr->worldPatch, GridW, GridH, GridD, txID);

  /* Initialize surface grid */
//  visualization::ChemSurfManager_ptr = new GridSurfaceManager();
  
  /* Initialize plots */
  std::string chemplot_title [MAX_PLOTS] = {"TNF", "TGF", "FGF", "MMP8", "IL1", "IL6", "IL8", "IL10"};
  std::string chemplot_xlabel[MAX_PLOTS], chemplot_ylabel[MAX_PLOTS];
  for(int i = 0; i < MAX_PLOTS; i++)
  {
	  chemplot_xlabel[i] = "Ticks";
	  chemplot_ylabel[i] = "Chem Level";
  }

  std::string cellplot_title [MAX_PLOTS] = {
		  "Neutrophils", "Activated Neutrophils", "Macrophages", "Activated Macrophages",
		  "Fibroblasts", "Activated Fibroblasts", "Platelets",   "Damage"};
  std::string cellplot_xlabel[MAX_PLOTS], cellplot_ylabel[MAX_PLOTS];
  for(int i = 0; i < MAX_PLOTS; i++)
  {
	  cellplot_xlabel[i] = "Ticks";
	  cellplot_ylabel[i] = "Total #Cells";
  }

  visualization::ChemPlotManager_ptr = new PlotManager(
		  	  	  	  	  	  	  	  	  	visualization::totalChem,
		  	  	  	  	  	  	  			1,
		  	  	  	  	  	  	  			chemplot_title,
		  	  	  	  	  	  	  			chemplot_xlabel,
		  	  	  	  	  	  	  			chemplot_ylabel);
  visualization::CellPlotManager_ptr = new PlotManager(
		  	  	  	  	  	  	  	  	  	visualization::totalCell,
		  	  	  	  	  	  	  			1,
		  	  	  	  	  	  	  			cellplot_title,
		  	  	  	  	  	  	  			cellplot_xlabel,
		  	  	  	  	  	  	  			cellplot_ylabel);
  
  /* Initalize chem and cell option table */

#ifdef VERT_PLOT
  visualization::ChemOptionTable_ptr = new OptionTable(WINW - 185.0f,
                                        180.0f,
                                        180.0f,
                                        250.0f,
                                        CHEM_OPTION_ROWS,
                                        CHEM_OPTION_COLS,
                                        visualization::optionHeader,
                                        visualization::optionLabel,
                                        optionBoxArr,
                                        visualization::surfColors);
  
  visualization::CellOptionTable_ptr = new OptionTable(WINW - 255.0f,
                                        123.0f,
                                        250.0f,
                                        50.0f,
                                        CELL_OPTION_ROWS,
                                        CELL_OPTION_COLS,
                                        visualization::optionHeaderCell,
                                        visualization::optionLabelCell,
                                        optionBoxArrCell,
                                        visualization::cellColors);
#else
  ChemOptionTable_ptr = new OptionTable(WINW - 190.0f,
                                        275.0f,
                                        180.0f,
                                        250.0f,
                                        CHEM_OPTION_ROWS,
                                        CHEM_OPTION_COLS,
                                        optionHeader,
                                        optionLabel,
                                        optionBoxArr,
                                        visualization::surfColors);
  
  CellOptionTable_ptr = new OptionTable(WINW - 250.0f,
                                        530.0f,
                                        250.0f,
                                        50.0f,
                                        CELL_OPTION_ROWS,
                                        CELL_OPTION_COLS,
                                        optionHeaderCell,
                                        optionLabelCell,
                                        optionBoxArrCell,
                                        visualization::cellColors);
#endif

  /* Initialize heatmaps */
  visualization::HeatMapManager_ptr = new HeatMapManager(GridW, GridH, GridD);
  
//  /* Initialize ECM volume */
//  visualization::ECMvolumeManager_Ptr = new VolumeManager(GridW, GridH, GridD, NULL);
//  ECMvolumeManager_Ptr->preprocessECMdata();


  
  /* Initialize ECM vertices */
  initWndGrid();
  
  /* Initialize ECM plane */
/*
  printf("Subplane x: %d -> %d\ty: %d -> %d\n", GL_X(visualization::WG_xb),
         GL_X(visualization::WG_xb + visualization::WG_w), GL_Y(visualization::WG_yb),
         GL_Y(visualization::WG_yb + visualization::WG_h));

  visualization::ECMplane_ptr = new SubPlane(
			      (GLfloat) visualization::WG_xb,
                              (GLfloat) visualization::WG_yb,
                              (GLfloat) visualization::WG_w,
                              (GLfloat) visualization::WG_h,
                              1.0f,//4.0f
                              visualization::texture_atlas,
                              visualization::glass_texture[1],
                              visualization::world_ptr->worldECM);
  
  visualization::ScarTerrain_ptr = new Terrain();
  if ( !visualization::ScarTerrain_ptr->LoadHeightmap(visualization::world_ptr->worldECM,
                                                      visualization::WG_xb,
                                                      visualization::WG_yb,
                                                      visualization::WG_w,
                                                      visualization::WG_h ) )
  {
    fprintf(stderr, "Failed to load heightmap for terrain!\n" );
  }
  if ( !visualization::ScarTerrain_ptr->LoadTexture( visualization::scar_texture ) )
  {
    fprintf(stderr, "Failed to load terrain texture for texture id: %d!\n", visualization::scar_texture);
  }
*/
  
}


void visualization::Draw2D()
{

  ShowCellsButton.highlighted = showCells;
 /* 
  if (showECMs) {
    ShowECMButton.highlighted = 1;
#ifdef USE_MOUSE_ONLY
    visualization::ECMplane_ptr->DrawLegend(20.0f, 340.0f, 40.0f);
#else   // USE_MOUSE_ONLY
    
#ifdef VERT_PLOT
    visualization::ECMplane_ptr->DrawLegend(WINW - 130.0f, WINH - 215.0f, 40.0f);
#else   // VERT_PLOT
    //    ECMplane_ptr->DrawLegend(WINW - 195.0f, 90.0f, 40.0f);
    ECMplane_ptr->DrawLegend(10.0f, 320.0f, 40.0f);
#endif  // VERT_PLOT
    
#endif  // USE_MOUSE_ONLY
  } else {
    ShowECMButton.highlighted = 0;
  }
*/
  glColor3f(1,1,1);
#ifdef DEBUG_VIS  
  // Print mouse
  char strx[15];
  char stry[15];
  sprintf(strx, "x: %d", TheMouse.x);
  sprintf(stry, "y: %d", TheMouse.y);
  
  Font(GLUT_BITMAP_HELVETICA_10, strx,30, 850);
  Font(GLUT_BITMAP_HELVETICA_10, stry,30+60, 850);
#endif  

  // Damage stats
  int initialDam = visualization::world_ptr->getInitialDam();
  int currentDam = visualization::world_ptr->countDamage();
  char str_idam[50];
  char str_cdam[50];
  char str_pheal[50];
  sprintf(str_idam,  "Initial Damage (patches): %d", initialDam);
  sprintf(str_cdam,  "Current Damage (patches): %d", currentDam);
  sprintf(str_pheal, "Percent Healed (%%):       %3.1f",
	((float) (initialDam - currentDam) / (float) initialDam) * 100.0f);
  
#ifdef VERT_PLOT
  Font(GLUT_BITMAP_8_BY_13, str_idam,  15, 20);
  Font(GLUT_BITMAP_8_BY_13, str_cdam,  15, 40);
  Font(GLUT_BITMAP_8_BY_13, str_pheal, 15, 60);
#else   // VERT_PLOT
  
#endif  // VERT_PLOT
  
#ifdef USE_MOUSE_ONLY
  ZoomWoundButton.highlighted = zoomedWound;
  ShowTNF_HMButton.highlighted = showChem_HM[TNF];
  ShowTGF_HMButton.highlighted = showChem_HM[TGF];  
//  ShowTNF_HMButton.highlighted = showChem_HM[FGF];
//  ShowTNF_HMButton.highlighted = showChem_HM[MMP8];
//  ShowTNF_HMButton.highlighted = showChem_HM[IL1beta];
//  ShowTGF_HMButton.highlighted = showChem_HM[IL6];
//  ShowTGF_HMButton.highlighted = showChem_HM[IL8];
//  ShowTGF_HMButton.highlighted = showChem_HM[IL10];

  ShowTNF_SurfButton.highlighted = showChem_Surf[TNF];
  ShowTGF_SurfButton.highlighted = showChem_Surf[TGF];
  if (showTNF_HM || showTGF_HM) {
    HeatMapManager_ptr->DrawLegend();
  }
  
#else     // USE_MOUSE_ONLY

  ShowChemOpButton.highlighted = showChemOp;

  bool isHM = false;
  for (int ic = 0; ic < CHEM_OPTION_ROWS; ic++) {
	  optionBoxArr[ic].highlighted	= showChem_HM  [ic];
	  if (showChem_HM[ic]) isHM 		= true;

//    optionBoxArr[ic*2].highlighted	= showChem_HM  [ic];
//    if (showChem_HM[ic]) isHM 		= true;
//
//    optionBoxArr[ic*2 + 1].highlighted	= showChem_Surf[ic];
  }

  optionBoxArrCell[0].highlighted = showNeus;
  optionBoxArrCell[1].highlighted = showMacs;
  optionBoxArrCell[2].highlighted = showFibs; 

  if (isHM) {
    HeatMapManager_ptr->DrawLegend();
  }
  
#endif    // USE_MOUSE_ONLY

  ButtonDraw(&ShowECMButton);
  ButtonDraw(&ShowCellsButton); 

#ifdef USE_MOUSE_ONLY
  ButtonDraw(&ZoomInButton);
  ButtonDraw(&ZoomOutButton);
  ButtonDraw(&MoveUpButton);
  ButtonDraw(&MoveLeftButton);
  ButtonDraw(&MoveRightButton);
  ButtonDraw(&MoveDownButton);
  ButtonDraw(&ZoomWoundButton);
  
  ButtonDraw(&ShowTNF_HMButton);
  ButtonDraw(&ShowTNF_SurfButton);
  ButtonDraw(&ShowTGF_HMButton);
  ButtonDraw(&ShowTGF_SurfButton);
  
  ButtonDraw(&ShowFGF_HMButton);
  ButtonDraw(&ShowFGF_SurfButton);
  ButtonDraw(&ShowMMP8_HMButton);
  ButtonDraw(&ShowMMP8_SurfButton);
  ButtonDraw(&ShowIL1_HMButton);
  ButtonDraw(&ShowIL1_SurfButton);
  ButtonDraw(&ShowIL6_HMButton);
  ButtonDraw(&ShowIL6_SurfButton);
  ButtonDraw(&ShowIL8_HMButton);
  ButtonDraw(&ShowIL8_SurfButton);
  ButtonDraw(&ShowIL10_HMButton);
  ButtonDraw(&ShowIL10_SurfButton);
#else
  ButtonDraw(&ShowChemOpButton);
#endif 
 

  HeatMapManager_ptr->DrawLegend();
  
  char str_tc[15];
  char str_dy[15];
  char str_hr[15];
  char str_mn[15];
  int total_mn = visualization::iter*30;
  int total_hr  = total_mn / 60;
  sprintf(str_tc, "%d", visualization::iter);
  sprintf(str_dy, "%d", total_hr/24);
  sprintf(str_hr, "%d", total_hr%24);
  sprintf(str_mn, "%d", total_mn%60);
  glColor3f(0.7,0.7,0.7);
  Font(GLUT_BITMAP_HELVETICA_12,"elapsed time:       d       h      m", WINW - 190, 40);
  Font(GLUT_BITMAP_HELVETICA_12,"tick:  ",   WINW - 140, 20);
  glColor3f(1,1,1);
  Font(GLUT_BITMAP_HELVETICA_12,str_dy,      WINW - 105,  40);
  Font(GLUT_BITMAP_HELVETICA_12,str_hr,      WINW - 73,  40);
  Font(GLUT_BITMAP_HELVETICA_12,str_mn,      WINW - 42,  40);
  
  Font(GLUT_BITMAP_HELVETICA_12,str_tc,      WINW - 100,  20);
}

void getColor()
{
	static unsigned int state		= 0;	//holds which state we are in
	switch(state) {							// (R, G, B)
	case(0):								//(High, Increasing, Low)
		particleColor.y  += COLORSPEED*1.0f/255.0f;
		if(particleColor.y >= MAXCOLOR)
		{
			state = 1;
			particleColor.y = MAXCOLOR;
		}
		break;
	case(1):								//(Decreasing, High, Low)
		particleColor.x  -= COLORSPEED*1.0f/255.0f;
		if(particleColor.x <= MINCOLOR)
		{
			state = 2;
			particleColor.x = MINCOLOR;
		}
		break;
	case(2):								//(Low, High, Increasing)
		particleColor.z  += COLORSPEED*1.0f/255.0f;
		if(particleColor.z >= MAXCOLOR)
		{
			state = 3;
			particleColor.z = MAXCOLOR;
		}
		break;
	case(3):								//(Low, Decreasing, High)
		particleColor.y  -= COLORSPEED*1.0f/255.0f;
		if(particleColor.y <= MINCOLOR)
		{
			state = 4;
			particleColor.y = MINCOLOR;
		}
		break;
	case(4):								//(Increasing, Low, High)
		particleColor.x  += COLORSPEED*1.0f/255.0f;
		if(particleColor.x >= MAXCOLOR)
		{
			state = 5;
			particleColor.x = MAXCOLOR;
		}
		break;
	case(5):								//(High, Low, Decreasing)
		particleColor.z  -= COLORSPEED*1.0f/255.0f;
		if(particleColor.z <= MINCOLOR)
		{
			state = 0;
			particleColor.z = MINCOLOR;
		}
		break;
    } 

}

void getColor2()
{
	static unsigned int state2		= 0;	//holds which state we are in
	switch(state2) {							// (R, G, B)
	case(0):								//(High, Increasing, Low)
		particleColor2.y  += COLORSPEED2*1.0f/255.0f;
		if(particleColor2.y >= MAXCOLOR2)
		{
			state2 = 1;
			particleColor2.y = MAXCOLOR2;
		}
		break;
	case(1):								//(Decreasing, High, Low)
		particleColor2.x  -= COLORSPEED2*1.0f/255.0f;
		if(particleColor2.x <= MINCOLOR2)
		{
			state2 = 2;
			particleColor2.x = MINCOLOR2;
		}
		break;
	case(2):								//(Low, High, Increasing)
		particleColor2.z  += COLORSPEED2*1.0f/255.0f;
		if(particleColor2.z >= MAXCOLOR2)
		{
			state2 = 3;
			particleColor2.z = MAXCOLOR2;
		}
		break;
	case(3):								//(Low, Decreasing, High)
		particleColor2.y  -= COLORSPEED2*1.0f/255.0f;
		if(particleColor2.y <= MINCOLOR2)
		{
			state2 = 4;
			particleColor2.y = MINCOLOR2;
		}
		break;
	case(4):								//(Increasing, Low, High)
		particleColor2.x  += COLORSPEED2*1.0f/255.0f;
		if(particleColor2.x >= MAXCOLOR2)
		{
			state2 = 5;
			particleColor2.x = MAXCOLOR2;
		}
		break;
	case(5):								//(High, Low, Decreasing)
		particleColor2.z  -= COLORSPEED2*1.0f/255.0f;
		if(particleColor2.z <= MINCOLOR2)
		{
			state2 = 0;
			particleColor2.z = MINCOLOR2;
		}
		break;
    } 

}
void visualization::Draw3D()
{
  // Set the camera
#ifdef HUMAN_VF
  float eyePos[]          ={  50.0f + camX, -50.0f + camY, 1300.0f + camZ};
  float lookAtDir[]       ={-900.0f, -50.0f, PLANE_DEPTH+1};
  float upDir[]           ={0.0f, 1.0f, 0.0f};

  int fibOffset  =     0; //400;
  int neuOffset  =  -400;//400; //0;
  int macOffset  =  -800;//-600; //-400;
  int chemOffset = -1400;//800

  int neuOffsetY  =     0;
  int macOffsetY  =     0;
  int fibOffsetY  =     0;
  int chemOffsetY =     0;

#elif defined(RAT_VF)
  float eyePos[]          ={  -50.0f + camX, -50.0f + camY, -850.0f + camZ};
  float lookAtDir[]       ={ -250.0f, -50.0f, PLANE_DEPTH+1};
  float upDir[]           ={    0.0f, 1.0f, 0.0f};

  int neuOffset  =     0;//   0;
  int macOffset  =     0;
  int fibOffset  =  -200;
  int chemOffset =  -200;//-600;

  int neuOffsetY  =     0;//   0;
  int macOffsetY  =  -200;
  int fibOffsetY  =  -200;
  int chemOffsetY =     0;//-600;

#else
  float eyePos[]          ={ 250.0f + camX, -50.0f + camY,  500.0f + camZ};
  float lookAtDir[]       ={-900.0f, -50.0f, PLANE_DEPTH+1};
  float upDir[]           ={0.0f, 1.0f, 0.0f};

  int fibOffset  =     0; //400;
  int neuOffset  =  -200;//400; //0;
  int macOffset  =  -400;//-600; //-400;
  int chemOffset =  -700;//800

  int neuOffsetY  =     0;
  int macOffsetY  =     0;
  int fibOffsetY  =     0;
  int chemOffsetY =     0;

#endif

  gluLookAt(eyePos[0],       eyePos[1], eyePos[2],
                        lookAtDir[0],   lookAtDir[1],   lookAtDir[2],
                        upDir[0],       upDir[1],       upDir[2]); 
  
  int nx = world_ptr->nx;
  int ny = world_ptr->ny;
  int nz = world_ptr->nz;
  
  /*========= Damage and Capillary ============*/

  visualization::EnvVisualizer_ptr->Render();
  // Draw rectangle perimeterizing world
  glColor3f(0.3f, 0.3f, 0.3f);
  glLineWidth(0.15);
  drawTransparentBox(GL_X(nx/2) + chemOffset, GL_Y(ny/2) + chemOffsetY, PLANE_DEPTH + nz/2, nx, ny, nz);
  drawTransparentBox(GL_X(nx/2) +  neuOffset, GL_Y(ny/2) +  neuOffsetY, PLANE_DEPTH + nz/2, nx, ny, nz);
  drawTransparentBox(GL_X(nx/2) +  macOffset, GL_Y(ny/2) +  macOffsetY, PLANE_DEPTH + nz/2, nx, ny, nz);
  drawTransparentBox(GL_X(nx/2) +  fibOffset, GL_Y(ny/2) +  fibOffsetY, PLANE_DEPTH + nz/2, nx, ny, nz);

  // ========= Render Scar Terrain ========== //
  
  
//  visualization::ScarTerrain_ptr->UpdateHeightVertexData();
//  visualization::ScarTerrain_ptr->Render(); 
  
  // ========= Cytokine Surface Mesh ========== //
 /* 
  heightFromPlane += 0.1f;
  bool isSurf = false;

  glLineWidth(0.5);
  for (int ic = 0; ic < CHEM_OPTION_ROWS; ic++) {
    if (showChem_Surf[ic]) {
      int coli = ic * 2 + 1;
      glColor4f(visualization::surfColors[coli][0], visualization::surfColors[coli][1],
		visualization::surfColors[coli][2], visualization::surfColors[coli][3]);
      visualization::ChemSurfManager_ptr->Render(
	                                    visualization::world_ptr->chemAllocation[ic],
                                            heightFromPlane-2.8f,
                                            GridH, GridW);
      isSurf = true;
    }
  }
  if (isSurf) {
    visualization::ChemSurfManager_ptr->RenderBox(heightFromPlane, GridH, GridW);
  }   
  */
  // ========= Render ECM Plane ========== // 
 /*

  if (showECMs)
  {
    visualization::ECMplane_ptr->Render(heightFromPlane, raiseECMplane);
  }
  */
  
  // ========= Render Cells ========== //
#define SHOW_CELL
#ifdef SHOW_CELL  
//  heightFromPlane -= 0.3f;

  int aX = 0, aY = 0, aZ = 0;
//  glClear(GL_COLOR_BUFFER_BIT);

#ifdef RAT_VF
  glPointSize( 2.0f);
#else
  glPointSize( 1.0f );
#endif
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE);					//Use additive blending


  //getColor();
  glColor4f(particleColor.x, particleColor.y, particleColor.z, particleColor.w);  
    glBegin(GL_POINTS);    
    if (showNeus) {
      //glColor4f(visualization::cellColors[0][0], visualization::cellColors[0][1],
      //          visualization::cellColors[0][2], visualization::cellColors[0][3]);
      int neusSize = world_ptr->neus.size();
      for (int i=0 ; i < neusSize; i++){
        Neutrophil *neu = visualization::world_ptr->neus.getDataAt(i);
        if (!neu) continue;
      
        aX = neu->getX();
        aY = neu->getY();
        aZ = neu->getZ();
         		      
	if (neu->isActivated()) {
		glColor4f(particleColor.x+0.2, particleColor.y+0.2, particleColor.z+0.1, particleColor.w+0.2);
	} else {
		glColor4f(particleColor.x, particleColor.y, particleColor.z, particleColor.w);
	}
        glVertex3i(GL_X(aX) + neuOffset, GL_Y(aY) + neuOffsetY, PLANE_DEPTH + aZ);
      
      }
    }
    glEnd();

  //getColor2();
  glColor4f(particleColor2.x, particleColor2.y, particleColor2.z, particleColor2.w);
    glBegin(GL_POINTS);
    if (showMacs) {
      //glColor4f(visualization::cellColors[1][0], visualization::cellColors[1][1],
      //          visualization::cellColors[1][2], visualization::cellColors[1][3]);
      int macsSize = world_ptr->macs.size();
      for (int i=0 ; i < macsSize; i++){
        Macrophage *mac = visualization::world_ptr->macs.getDataAt(i);
        if (!mac) continue;
      
        aX = mac->getX();
        aY = mac->getY();
        aZ = mac->getZ();
         		      
	if (mac->isActivated()) {
		glColor4f(particleColor2.x+0.2, particleColor2.y+0.2, particleColor2.z+0.2, particleColor2.w+0.3);
	} else {
		glColor4f(particleColor2.x, particleColor2.y, particleColor2.z, particleColor2.w);
	}
        glVertex3i(GL_X(aX) + macOffset, GL_Y(aY) + macOffsetY, PLANE_DEPTH + aZ);
      
      }
    }
    glEnd();

    glPointSize( 1.0f);
  glColor4f(particleColor3.x, particleColor3.y, particleColor3.z, particleColor3.w);
    glBegin(GL_POINTS);
    if (showFibs) {
      //glColor4f(visualization::cellColors[2][0], visualization::cellColors[2][1],
      //          visualization::cellColors[2][2], visualization::cellColors[2][3]); 
      int fibsSize = world_ptr->fibs.size();
      for (int i=0 ; i < fibsSize; i++){
        Fibroblast *fib = visualization::world_ptr->fibs.getDataAt(i);
        if (!fib) continue;
      
        aX = fib->getX();
        aY = fib->getY();
        aZ = fib->getZ();
        
	if (fib->isActivated()) {
		glColor4f(particleColor3.x+0.3, particleColor3.y+0.3, particleColor3.z+0.3, particleColor3.w+0.4);
	} else {
		glColor4f(particleColor3.x, particleColor3.y, particleColor3.z, particleColor3.w);
	}
        glVertex3i(GL_X(aX) + fibOffset, GL_Y(aY) + fibOffsetY, PLANE_DEPTH + aZ);
      
      }
    }
    glEnd();


  glDisable(GL_BLEND);
/*
#ifdef USE_SPRITE
  
  glActiveTexture( GL_TEXTURE1 );
  glEnable( GL_TEXTURE_2D );
  glBindTexture(GL_TEXTURE_2D, texture_sprite[0]);
  
  glPointSize( cellSize);//6.0f );
  
#else
  glPointSize( 0.001f);//cellSize );
#endif
 

  if (showNeus || showMacs || showFibs) {
    int aX = 0, aY = 0, aZ = 0;
    
    glBegin(GL_POINTS);
    
    if (1) {//showNeus) {
      glColor4f(visualization::cellColors[0][0], visualization::cellColors[0][1],
                visualization::cellColors[0][2], visualization::cellColors[0][3]);
      int neusSize = world_ptr->neus.size();
      for (int i=0 ; i < neusSize; i++){
        Neutrophil *neu = visualization::world_ptr->neus.getDataAt(i);
        if (!neu) continue;
      
        aX = neu->getX();
        aY = neu->getY();
        aZ = neu->getZ();
         		      
        glVertex3i(GL_X(aX) + neuOffset, GL_Y(aY), PLANE_DEPTH + aZ);
      
      }
    }
    
    if (1) {//showMacs) {
      glColor4f(visualization::cellColors[1][0], visualization::cellColors[1][1],
                visualization::cellColors[1][2], visualization::cellColors[1][3]);
      int macsSize = world_ptr->macs.size();
      for (int i=0 ; i < macsSize; i++){
        Macrophage *mac = visualization::world_ptr->macs.getDataAt(i);
        if (!mac) continue;
      
        aX = mac->getX();
        aY = mac->getY();
        aZ = mac->getZ();
         		      
        glVertex3i(GL_X(aX) + macOffset, GL_Y(aY), PLANE_DEPTH + aZ);
      
      }
    }

    
    if (showFibs) {
      glColor4f(visualization::cellColors[2][0], visualization::cellColors[2][1],
                visualization::cellColors[2][2], visualization::cellColors[2][3]); 
      int fibsSize = world_ptr->fibs.size();
      for (int i=0 ; i < fibsSize; i++){
        Fibroblast *fib = visualization::world_ptr->fibs.getDataAt(i);
        if (!fib) continue;
      
        aX = fib->getX();
        aY = fib->getY();
        aZ = fib->getZ();
         		      
        glVertex3i(GL_X(aX) + fibOffset, GL_Y(aY), PLANE_DEPTH + aZ);
      
      }
    }

    
    glEnd();
    
  }
 
  
#ifdef USE_SPRITE
  
  glActiveTexture( GL_TEXTURE0 );
  glEnable( GL_TEXTURE_2D );
  
#endif
*/ 
#endif // SHOW_CELL
  // ========= Cytokine Heat Maps ========== //
  visualization::getTotalChem();
  visualization::getTotalCell();
 
#ifdef SHOW_CHEM
 
  for (int ic = 0; ic < CHEM_OPTION_ROWS; ic++) {
  if (showChem_HM[ic])
      visualization::HeatMapManager_ptr->Render(
	visualization::world_ptr->WHWorldChem->pChem[ic],
	GridH,
	GridW,
    GridD,
	chemOffset,
	chemOffsetY);
 
  }

#endif

}

/* ================================================================================================ */

void visualization::getTotalChem()
{
	for (int ic = 0; ic < 8; ic++)
		totalChem[ic] = visualization::world_ptr->WHWorldChem->total[ic];
}

void visualization::getTotalCell()
{
	for (int i = 0; i < 8; i++)
		totalCell[i] = (float) visualization::world_ptr->totalCell[i];
}


/* The function called when our window is resized (which shouldn't happen, because we're fullscreen) */
void visualization::ReSizeGLScene(int Width, int Height)
{
  if (Height==0)				// Prevent A Divide By Zero If The Window Is Too Small
    Height=1;
  
  glViewport(0, 0, Width, Height);		// Reset The Current Viewport And Perspective Transformation
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  
  //  gluOrtho2D(0, Width, Height, 0);
  gluPerspective(45.0f,(GLfloat)Width/(GLfloat)Height,0.1f,PERSPECTIVE_DEPTH);
  glMatrixMode(GL_MODELVIEW);
}


#ifdef TIME_GL
#define TIME_BUFFER_SIZE    24
long gl_time[TIME_BUFFER_SIZE];
long go_time[TIME_BUFFER_SIZE];
long total_time[TIME_BUFFER_SIZE];
#endif

/* The main drawing function. */
//void visualization::DrawGLVolScene()
//{
//
//}

void visualization::DrawGLScene()
{

// DEBUG
cout << "DrawGLScene: begin" << endl; 
#ifdef TIME_GL
  struct timeval start, end, go_start, go_end;
  gettimeofday(&start, NULL);
#endif
 
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);	// Clear The Screen And The Depth Buffer
  //  glViewport(0,0/*900/3*/,1440,900/* - 900/3*/);   // full screen resolution: 900x1440
  
  double w = glutGet( GLUT_WINDOW_WIDTH );
  double h = glutGet( GLUT_WINDOW_HEIGHT );
  glViewport(0,0,w,h);   // full screen resolution: 900x1440
  /*
   *	Set perspective viewing transformation
   */
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  //  gluPerspective(45,(WINH==0)?(1):((float)WINW/WINH),0.1f, 6000);//PERSPECTIVE_DEPTH);
  gluPerspective(45,(w==0)?(1):((float)w/h),0.1f, PERSPECTIVE_DEPTH);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  
  
  Draw3D();
  
  
  /*
   *	Disable depth test and lighting for 2D elements
   */
//  glDisable(GL_DEPTH_TEST);
  //  glDisable(GL_LIGHTING);
  
  /*
   *	Set the orthographic viewing transformation
   */

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glOrtho(0,w,h,0,-1,1);
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glViewport(0,0,w,h);

  /*
   *	Draw the 2D overlay
   */
  Draw2D();

  if (showChemOp) {
    visualization::ChemOptionTable_ptr->Render();
  }

  if (showCells) {
    visualization::CellOptionTable_ptr->Render();
  }
  

  /*
   *	Draw chem viewports
   */

  if (showChemCharts) {
    visualization::ChemPlotManager_ptr->DrawChemViewPorts(visualization::iter, showChemCharts);
  } else if (showCellCharts){
	visualization::CellPlotManager_ptr->DrawChemViewPorts(visualization::iter, showCellCharts);
  }
  glColor3f(1,1,1);
  
  
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
  
  glEnable(GL_DEPTH_TEST);
 

// DEBUG
cout << "DrawGLScene: Before swapbuffer()" << endl; 
  // swap the buffers to display, since double buffering is used.
  glutSwapBuffers();


cout << "DrawGLScene: After swapbuffer()" << endl; 



#ifdef TIME_GL
  gettimeofday(&go_start, NULL);
#endif 
  if (!paused)
  {
/*****************************************/
#ifdef OVERLAP_VIS 
    visualization::notifyVisDone();
    cout << "visualization::start() waiting for compute ..." << endl;
    visualization::waitForCompute();
    cout << "visualization::start() done waiting for compute" << endl;
#else	// OVERLAP_VIS
    // Update new locations
    visualization::world_ptr->go();
#endif	// OVERLAP_VIS
/*****************************************/
    visualization::iter++;
  }
#ifdef TIME_GL
  gettimeofday(&go_end, NULL);
#endif  




#ifdef TIME_GL
  gettimeofday(&end, NULL);

  int ind = (iter - 1) % TIME_BUFFER_SIZE;
  total_time[ind] = (end.tv_sec * 1000 + end.tv_usec / 1000)
		  - (start.tv_sec * 1000 + start.tv_usec / 1000);

  go_time[ind] = (go_end.tv_sec * 1000 + go_end.tv_usec / 1000)
		  - (go_start.tv_sec * 1000 + go_start.tv_usec / 1000);

  gl_time[ind] = total_time[ind] - go_time[ind];
  long sum_total = 0, sum_go = 0, sum_gl = 0;
  for (int i = 0; i < TIME_BUFFER_SIZE; i++) {
    sum_total += total_time[i];
    sum_go    += go_time[i];
    sum_gl    += gl_time[i];
  }

  printf("Iter %d times and average of the last %d iterations: \n", iter - 1, TIME_BUFFER_SIZE);
  printf("\tcompute time: %ld\taverage: %ld\n", go_time[ind], sum_go/TIME_BUFFER_SIZE);
  printf("\trender time:  %ld\taverage: %ld\n", gl_time[ind], sum_gl/TIME_BUFFER_SIZE);
  printf("\tTOTAL time:   %ld\taverage: %ld\n", total_time[ind], sum_total/TIME_BUFFER_SIZE);
#endif

}




/* ================================================================================================ */

