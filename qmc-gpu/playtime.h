#include <cg/cgGL.h>

void getError(char *msg);
void idle();
void reshape(int w, int h);
void keyboard(unsigned char key, int x, int y);
void display();
void init();
void initializeCg();
void cgErrorCallback();
void perspective(bool tranformed);
void specTest();
void PrintArray(Array2D<GLfloat> matrix);
void testTiming();
void matrixMultiplyTimings(int d);

CGcontext g_cgContext;
CGprofile g_cgProfile;

static int window_w, window_h;