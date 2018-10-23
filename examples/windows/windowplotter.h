

#ifndef _winHENON_H_ 
#define _winHENON_H_ 1

#include<windows.h>
#include "resource.h"

extern int xdim;
extern int ydim;

extern double a;
extern double b;

extern bool running;

extern bool calcattr;
extern bool calclyp;
extern bool calcbasin;

extern int which; //0=basin, 1=lyp, 2=attr

extern char *pr_ClassWinName;;		//the basic registry class name
extern  HINSTANCE pr_hInst;							//The main window instance...
extern HDC myHDC;
extern RECT winRECT;	//the main window rectangle size
extern HWND hwndPB; 	//the progress bar sub window


extern int simmain(HWND master,HWND prog, RECT &inR);
extern bool DrawPoint(HDC draw, RECT &curRC);
extern bool ReDrawAllPoints(HWND hwnd, int which);


/** Destroy Handler **/
LRESULT CALLBACK KillProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);

// ABout dialog
int CALLBACK AboutProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);

// sets the X and Y divisions
BOOL CALLBACK EditXYdims(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);

// sets the A and B factors
BOOL CALLBACK EditABvals(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);

/** Handles Menu stuff **/
LRESULT CALLBACK MenuProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);

/***  Main Message Handler ***/
LRESULT CALLBACK MainProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam);

//redraws (or draws) the scroll bar on the bottom of the screen
bool ReDrawScrollBar(HWND hwnd);

//draws a white rectangle that will cover everything in the screen
bool DrawWhiteBG(HWND hwnd);
#endif


