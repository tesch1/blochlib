

#include<windows.h>
#include "resource.h"
#include "windowplotter.h"
#include <commctrl.h>
#include <stdio.h>


 /* Global variable declarations */
int xdim=10;
int ydim=10;

double a=1.28;
double b=0.3;

bool calcattr=false;
bool calclyp=false;
 bool calcbasin=true;

int which=0;
bool running=false;

static char *pr_ClassWinName = "Windows Plotter";		//the basic registry class name
static HINSTANCE pr_hInst=NULL;							//The main window instance...
static HDC myHDC=NULL;	//the main drawing window
static RECT winRECT={0,0,10,10}; //main window rect size
static HWND hwndPB; 	//the progress bar sub window


/** Destroy Handler **/
LRESULT CALLBACK KillProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	switch(Message)		//Mesage being sent
	{
	case WM_CLOSE:		//somebody hit the [X] on the top right
		DestroyWindow(hwnd);	//kill the window
		break;
	case WM_DESTROY:	//everything has finished kill everything
		PostQuitMessage(0);	//LAST function to call...does the house cleaning and exits totaly
		break;
	default:
		return DefWindowProc(hwnd, Message, wParam, lParam); //does NOTHING...i.e. Default
	}
	return 0;
}

/** Handles the ABout Dialog Box **/
// NOTE: the "WM_COMMAND" is passed from the main window to this one from 'MenuProc'
// DO _NOT_ call "DefWindProc" for Dialogs...it is auto...
int CALLBACK AboutProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	switch(Message)
	{
		case WM_INITDIALOG:		//initiallize (do nothin here)
			return 1;
		case WM_COMMAND:		//a  command was sent from the dialog
			switch(LOWORD(wParam))	//get the commad type
			{
				case IDOK:			//They hit OK 
					EndDialog(hwnd,IDOK);	//close the dialog box
					return IDOK;			//return the ok value 

				case IDCANCEL:		//they hit cancel
					EndDialog(hwnd, IDCANCEL);	//close the dialog boc
					return IDCANCEL;	//return the cancel value
			}
			break;
	}
	return 0;
}

/** Handled the 'X,Y' dim menu entering */
BOOL CALLBACK EditXYdims(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	const int bufsize=256;
	char inputX[bufsize];
	char inputY[bufsize];
	int err=1; //error flag
	int sig=1; //digned value flag
	switch(Message)
	{
		case WM_INITDIALOG: //open up the edit dlg box
			//this sets the SIZE limit to any entered text in the box
			SendDlgItemMessage(hwnd, IDC_XDIM, EM_SETLIMITTEXT, (WPARAM)bufsize-1, (LPARAM)0);			
			SetDlgItemInt(hwnd, IDC_XDIM, xdim, sig);			
			SendDlgItemMessage(hwnd, IDC_YDIM, EM_SETLIMITTEXT, (WPARAM)bufsize-1, (LPARAM)0);
			SetDlgItemInt(hwnd, IDC_YDIM, ydim, sig);			
			return true;
		case WM_COMMAND: //get the 'OK' or 'Cancel' hits
			switch(LOWORD(wParam))
			{
				case IDOK:
					xdim=GetDlgItemInt(hwnd,IDC_XDIM, &err, sig);
					ydim=GetDlgItemInt(hwnd,IDC_YDIM, &err, sig);
					if(xdim<=0 || ydim<=0 || !err)
					{
						MessageBox(hwnd,"X and Y dims must be >=1", "ERROR",  MB_ICONEXCLAMATION | MB_OK );
					}else{
						EndDialog(hwnd,IDOK);
					}
					return true;
				case IDCANCEL:
					if(xdim<=0 || ydim<=0)
					{
						MessageBox(hwnd,"X and Y dims must be >=1", "ERROR",  MB_ICONEXCLAMATION | MB_OK );
					}else{
						EndDialog(hwnd,IDCANCEL);
					}
					return true;
			}
			break;
	}
	return false;
}

/** Handled the 'A,B' dim menu entering */
BOOL CALLBACK EditABvals(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	const int bufsize=256;
	char inputX[bufsize];
	char inputY[bufsize];
	sprintf(inputX,"%f", a);
	sprintf(inputY,"%f", b);
	
	int err=1; //error flag
	int sig=1; //digned value flag
	switch(Message)
	{
		case WM_INITDIALOG: //open up the edit dlg box
			//this sets the SIZE limit to any entered text in the box
			SendDlgItemMessage(hwnd, IDC_AVAL, EM_SETLIMITTEXT, (WPARAM)bufsize-1, (LPARAM)0);			
			SetDlgItemText(hwnd, IDC_AVAL, inputX);			
			SendDlgItemMessage(hwnd, IDC_BVAL, EM_SETLIMITTEXT, (WPARAM)bufsize-1, (LPARAM)0);
			SetDlgItemText(hwnd, IDC_BVAL, inputY);			
			return true;
		case WM_COMMAND: //get the 'OK' or 'Cancel' hits
			switch(LOWORD(wParam))
			{
				case IDOK:
					GetDlgItemText(hwnd,IDC_AVAL, inputY, bufsize); a=atof(inputY);
					GetDlgItemText(hwnd,IDC_BVAL, inputY, bufsize); b=atof(inputY);
					if(!err)
					{
						MessageBox(hwnd,"X and Y dims must be >=1", "ERROR",  MB_ICONEXCLAMATION | MB_OK );
					}else{
						EndDialog(hwnd,IDOK);
					}
					return true;
				case IDCANCEL:
					EndDialog(hwnd,IDCANCEL);
					return true;
			}
			break;
	}
	return false;
}

/** Handled the 'A,B' dim menu entering */
BOOL CALLBACK SimStopper(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	int err=1; //error flag
	int sig=1; //digned value flag
	switch(Message)
	{
		case WM_INITDIALOG: //open up the edit dlg box
			//this sets the SIZE limit to any entered text in the box			
			return true;
		case WM_COMMAND: //get the 'OK' or 'Cancel' hits
			switch(LOWORD(wParam))
			{
				case IDSTOP:
					EndDialog(hwnd,IDSTOP);
					return true;
			}
			break;
	}
	return false;
}

/** Handles Menu stuff **/
LRESULT CALLBACK MenuProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	int ret=0;
	static PAINTSTRUCT ps;
	bool running=false;
	static	HMENU hmenu;
	switch( LOWORD(wParam) )	//handle things from out menu...
	{
        case CM_FILE_EXIT:	//Check to see if the 'exit' selection was used
           PostMessage(hwnd,WM_CLOSE,0,0);
        break;
		case CM_XYDIMS: //set sdim and ydim
			ret=DialogBox(pr_hInst, MAKEINTRESOURCE(IDD_DIMSENTER), hwnd, EditXYdims);
			if(!ret){
				MessageBox(hwnd, "Failed to open Dialog", "ERROR!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
			}
		break;
		case CM_ABSET: //set a and b
			ret=DialogBox(pr_hInst, MAKEINTRESOURCE(IDD_ABSETTER), hwnd, EditABvals);
			if(!ret){
				MessageBox(hwnd, "Failed to open Dialog", "ERROR!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
			}
		break;

		case CM_RUN: //The Runner.....
		
			const int maxct=1000;
			SendMessage(hwndPB,PBM_SETRANGE , 0, MAKELPARAM(0, maxct)); //the maxct for the prog
   			SendMessage(hwndPB,PBM_SETSTEP, (WPARAM) 1, 0); //sets the step size for the prog
   			//erase any info already there
   			
   			
   			DrawWhiteBG(hwnd);
   			running=true;
   			simmain(hwnd, hwndPB,  winRECT); //calculate the henon bits
			
		break;
		case CM_LYP:
			calclyp=true;
			calcattr=false;
			calcbasin=false;
			which=1;
			hmenu=GetMenu(hwnd);
			CheckMenuItem(hmenu, CM_LYP, MF_BYCOMMAND |  MF_CHECKED); 
			CheckMenuItem(hmenu, CM_ATTR, MF_BYCOMMAND |  MF_UNCHECKED); 
			CheckMenuItem(hmenu, CM_BASIN, MF_BYCOMMAND |  MF_UNCHECKED); 
   			
		break;
		case CM_ATTR: 
			calclyp=false;
			calcattr=true;
			calcbasin=false;
			which=2;
			hmenu=GetMenu(hwnd);
			CheckMenuItem(hmenu, CM_LYP, MF_BYCOMMAND |  MF_UNCHECKED); 
			CheckMenuItem(hmenu, CM_ATTR, MF_BYCOMMAND |  MF_CHECKED); 
			CheckMenuItem(hmenu, CM_BASIN, MF_BYCOMMAND |  MF_UNCHECKED); 
   			
		break;
		case CM_BASIN: 
			calclyp=false;
			calcattr=false;
			calcbasin=true;
			which=0;
			hmenu=GetMenu(hwnd);
			CheckMenuItem(hmenu, CM_LYP, MF_BYCOMMAND |  MF_UNCHECKED); 
			CheckMenuItem(hmenu, CM_ATTR, MF_BYCOMMAND |  MF_UNCHECKED); 
			CheckMenuItem(hmenu, CM_BASIN, MF_BYCOMMAND |  MF_CHECKED); 
   			
		break;
		case CM_ABOUT:	//May About Menu item was hit..send it to the About Dialoag Handlers
			ret=DialogBox(pr_hInst, MAKEINTRESOURCE(IDD_ABOUT), hwnd, AboutProc);
			//Use the 'ret' value above to do things AFTER they hit something like 'OK' or 'CANCEL'
			//in the dialog box 
			if(!ret){
				MessageBox(hwnd, "Failed to open Dialog", "ERROR!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
			}
		break;
			
    }
	return 0;
}

//redraws (or draws) the scroll bar on the bottom of the screen
bool ReDrawScrollBar(HWND hwnd)
{
		GetClientRect(hwnd, &winRECT); 
		int cyVScroll = GetSystemMetrics(SM_CYVSCROLL); 

		hwndPB=CreateWindowEx(0,PROGRESS_CLASS, (LPSTR) NULL, 
		 WS_CHILD | WS_VISIBLE, winRECT.left, 
		 winRECT.bottom - cyVScroll, 
		 winRECT.right, cyVScroll, 
		 hwnd, (HMENU) 0, pr_hInst, NULL); 
		if(!hwndPB){
			MessageBox(hwnd, "Failed to open Prog Bar", "ERROR!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
			return false;
		}
		return true;
}

//draws a white rectangle that will cover everything in the screen
bool DrawWhiteBG(HWND hwnd)
{
	static PAINTSTRUCT ps;
	GetClientRect(hwnd, &winRECT);
	myHDC=GetDC(hwnd);
	BeginPaint(hwnd, &ps);
	Rectangle(myHDC,winRECT.left, winRECT.top, winRECT.right, winRECT.bottom);
 	EndPaint(hwnd, &ps);
 	ReleaseDC(hwnd,myHDC);
 	return true;
}

/***  Main Message Handler ***/
LRESULT CALLBACK MainProc(HWND hwnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
	static char ModName[MAX_PATH];

	switch(Message)
	{
	case WM_CREATE:
	//draw the scroll bar from the start
		InitCommonControls();
		 DrawWhiteBG(hwnd);
		ReDrawScrollBar(hwnd);
		break;
	case WM_SIZE:
	
		if(WM_LBUTTONUP)
		{
			GetClientRect(hwnd, &winRECT); 
			 DrawWhiteBG(hwnd);
			 ReDrawScrollBar(hwnd);
			ReDrawAllPoints(hwnd, which);
		}
		break;
/*	case WM_LBUTTONDOWN:	//left mouse down...
		GetModuleFileName(pr_hInst,ModName,MAX_PATH);		//a function that gets the program instance name (from the main app instatnce)
		MessageBox(hwnd,ModName, "The Program Be...",  MB_ICONINFORMATION | MB_OK );
		break;
*/
	case WM_COMMAND:		//the MEnu was used
		MenuProc(hwnd,Message, wParam, lParam);		//send it to our menu contoler
		break;
	
	default:
		return KillProc(hwnd,Message,wParam,lParam);
	}
	return 0;
}



/*
	!!THE MAIN!!
		hInstance-->the handled in mem to this program
		hPrevInstance-->the last instance of this program (OLD kept for compat issues...Win16 only...for Win32 always NULL)
		lpCmdLine-->program arguments (does NOT include the program name)
		nShowCmd-->value to pass to "ShowWindow()"

*/
int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,LPSTR lpCmdLine, int nShowCmd)
{
	WNDCLASSEX WndClass;	//our main window to create
	HWND hwnd;				//
	MSG msg;			//our message to pass

	pr_hInst=hInstance;		//set our application window to the global var
	WndClass.cbSize=sizeof(WNDCLASSEX);	//the mem size of the window
	WndClass.style=NULL;				//Class Styles _NOT_ Window Style (typically NULL)
	WndClass.lpfnWndProc=MainProc;		//the ptr to the main message handler
	WndClass.cbClsExtra=0;				//Any extra mem to allocate for this window
	WndClass.cbWndExtra=0;				//Any extra mem to allocate for this window
	WndClass.hInstance=pr_hInst;		//handel to the application instance (i.e. the hInstance)
	WndClass.hIcon=LoadIcon(pr_hInst,"windowplotter.ico");	//icon resorce (for Alt+Tab requests)
	WndClass.hCursor=LoadCursor(NULL, IDC_ARROW);	//icon for the cursor to show while over the window
	WndClass.hbrBackground=(HBRUSH)(COLOR_WINDOW-1);	//the backjound window color
	WndClass.lpszMenuName=MAKEINTRESOURCE(ID_MAIN_MENU);							//a Menu resource (if any)
	WndClass.lpszClassName=pr_ClassWinName;				//the class name _NOT_ the window name
	WndClass.hIconSm=LoadIcon(NULL,IDI_APPLICATION);	//the 16x16 icon in the top left corner

	//not we need to register our window...
	if(!RegisterClassEx(&WndClass)){
		MessageBox(0,"Failed to Register the Window!", "Error!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
		return 0;		//exit
	}

	//now we must create a HWND to draw (this does not draw it yet)
	hwnd=CreateWindowEx(
		WS_EX_CLIENTEDGE,		//the appearence of the the window Extended WINDOW STYLE (this is basic one)
		pr_ClassWinName,		//Must match the one we registerd above
		"Henon Map Exporler",	//display title

		WS_OVERLAPPEDWINDOW | WS_VSCROLL,		//a windows style parameter
		CW_USEDEFAULT,				//the initial X placement (here default, could be an int)
		CW_USEDEFAULT,				//the initial Y placement (here default, could be an int)
		700,500,					//Width, Height
		NULL,						//The parent handel (this is the top so NO handel)
		NULL,						//The menu handle (there is not one here)
		pr_hInst,					//application instance (i.e. this)
		NULL);						//ptr to the creation window (this again is the top, so NO handel)

	//check to see if out creation failed
	if(hwnd==NULL){
		MessageBox(0,"Failed to Create the Window!", "Error!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
		return 0;		//exit
	}
	
	ShowWindow(hwnd, nShowCmd);		//draw the window!
	UpdateWindow(hwnd);				//update the window

	/* *** THE MIGHTY EVENT LOOP CATCHER ***/
	while(GetMessage(&msg, NULL, 0,0))	//grab any events in the area....
	//(the message var, another Window to watch, a min filter, a max filter)
	{
		TranslateMessage(&msg);		//translate it to its proper form (i.e. like a character if a key was pushed)
		DispatchMessage(&msg);		//send the message along its way....
	}
	return msg.wParam;
}

