

#ifndef _henon_sim_h_ 
#define _henon_sim_h_ 1


#include "blochlib.h"
#include <windows.h>
#include <commctrl.h>
#include <winuser.h>

extern double a;
extern double b;
extern int xdim;
extern int ydim;
extern bool running;
extern bool calcattr;
extern bool calclyp;
extern bool calcbasin;
extern int which;

void func(Vector<double> &X);
void Jacobi(Vector<double> &X);
int simmain(HWND master,HWND prog, RECT &inR);
bool DrawPoint(HWND cwnd, RECT &curSC, double cX, double cY);
bool DrawBoundingBox(HWND hwnd);
bool ReDrawAllPoints(HWND hwnd, int which);

#endif