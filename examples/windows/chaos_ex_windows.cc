

#ifndef _henon_cc_ 
#define _henon_cc_ 1

#include "chaos_ex_windows.h"

double Mxdiv=2.5, Mydiv=2.5;
Vector<coord<double, 2> > basindata;
Vector<coord<double, 2> > lypdata;
Vector<coord<double, 2> > attrdata;

ofstream inf("inf");
	
void func(coord<double, 2> &X)
{
	static double tx;
	tx=X.x();
	X.x()=a-tx*tx+b*X.y();
	X.y()=tx;
	inf<<X<<endl;
}

rmatrix Jc(2,2,0);
void Jacobi(coord<double, 2> &pt)
{
	Jc(0,0)=-2*pt.x();
	Jc(0,1)=b;
	Jc(1,0)=1;
	Jc(1,1)=0;
}

bool DrawPoint(HWND cwnd, RECT &winR, double cX, double cY)
{
	//color fo the fill and the outline
	static HPEN hpen, hpenOld;
    static HBRUSH hbrush, hbrushOld;

    hpen = CreatePen(PS_SOLID, 0, RGB(100, 100, 150));
    hbrush = CreateSolidBrush( RGB(100, 100, 150));

    // Select the new pen and brush, and then draw.
   static PAINTSTRUCT ps;
	HDC inhdc=GetDC(cwnd);
	 SelectObject(inhdc, hpen);
     SelectObject(inhdc, hbrush);
	
	ps.hdc=inhdc;
	ps.rcPaint=winR;
	static double rw=1, rh=1;
	int wid=(winR.right-winR.left);
	int hei=(winR.bottom-winR.top);
	int Mw=2*Mxdiv;
	int Mh=2*Mydiv;
	cX*=(double(wid)/double(Mw))*0.95;
	cY*=(double(hei)/double(Mh))*0.93;
	cX+=double(wid)/2.*0.9+0.06*hei;
	cY+=double(hei)/2.*0.9+0.025*wid;
	//char inn[256];
	//sprintf(inn, "%f %f", cX, cY);
	//MessageBox(cwnd,inn, "ERROR!", MB_ICONEXCLAMATION | MB_OK | MB_SYSTEMMODAL);
						
	BeginPaint(cwnd, &ps);
	Rectangle(inhdc, cX-rw, cY-rh, cX+rw, cY+rh);
	EndPaint(cwnd, &ps);
	//ReleaseDC(cwnd, inhdc); 

	return true;
}

bool DrawBoundingBox(HWND hwnd)
{
	HDC hdc=GetDC(hwnd);
	static RECT rs;
	int cyVScroll = GetSystemMetrics(SM_CYVSCROLL); 
	GetClientRect(hwnd, &rs);
	MoveToEx(hdc,rs.left+5, rs.top+5, (LPPOINT) NULL); 
	LineTo(hdc, rs.left+5, rs.bottom-cyVScroll-5); 
	LineTo(hdc,rs.right-5, rs.bottom-cyVScroll-5); 
	LineTo(hdc,rs.right-5, rs.top+5); 
	LineTo(hdc,rs.left+5, rs.top+5);  
	ReleaseDC(hwnd, hdc); 
	return true;
}

bool ReDrawAllPoints(HWND hwnd, int which)
{
	static RECT rr;
	GetClientRect(hwnd, &rr);
	if(which==0)
	{
		for(int i=0;i<basindata.size();++i)
		{
			 DrawPoint(hwnd, rr, basindata[i].x(),  basindata[i].y());
		}
	}
	else if(which==1)
	{
		for(int i=0;i<lypdata.size();++i)
		{
			 DrawPoint(hwnd, rr, lypdata[i].x(),  lypdata[i].y());
		}
	}
	else if(which==2)
	{
		for(int i=0;i<attrdata.size();++i)
		{
			 DrawPoint(hwnd, rr,attrdata[i].x(),  attrdata[i].y());
		}
	}
}

	

int simmain(HWND master,HWND prog, RECT &inR)
{
	DrawBoundingBox(master);
	if(calcbasin) basindata.resize(xdim*ydim);
	if(calclyp) lypdata.resize(xdim*ydim);
	if(calcattr) attrdata.resize(xdim*ydim);
	coord<double, 2> IC(2,0);
	IC[0]=0.78;
	IC[1]=0.94;
	
	rmatrix Jt(2,2,0);
	Jt.identity();
	
	coord<double, 2> lyp(2,0);
	double xdiv=Mxdiv, ydiv=Mydiv;
	int i=0, j=0,kk=0;;
	double xs=-xdiv, ys=-ydiv;
	double cap=100;
	int ProgStepSize=(ydim*xdim)/1000;
	int ct=0, basinct=0, lypct=0, attrct=0;
	while(j<ydim)
	{
		i=0;
		while(i<xdim)
		{
			kk=0;
			ct++;
			//update progress bar
			if(ct==ProgStepSize)
			{
				SendMessage(prog, PBM_STEPIT, 0, 0); 
				ct=0;
			}
			xs=-xdiv+xdiv*2.0*double(j)/double(ydim);
			IC[0]=xs;
			ys=-ydiv+ydiv*2.0*double(i)/double(xdim);
			IC[1]=ys;
			while(kk<10)
			{
			//	if(!running) return 2;
				func(IC);
				if(calcbasin && norm(IC)>cap)
				{
					basindata[basinct]=coord<double, 2>(xs,ys);
					++basinct;
					DrawPoint(master, inR, xs,ys);
					break; 
				}				
				kk++;
			}
			if(calcattr && norm(IC)<cap || calclyp)
			{
				kk=0;
				if(calclyp) Jt.identity();
				lyp=0;
				while(kk<50)
				{
					//if(!running) return 2;
					if(norm(IC)>cap) break;
					func(IC);
					if(calclyp)
					{
						Jacobi(IC);
						Jt=Jt*Jc;
						Jt=GramSchmidt(Jt);
						lyp[0]+=log(norm(Jt.col(0)));
						Jt.putCol(0, Jt.col(0)/norm(Jt.col(0)));
						lyp[1]+=log(norm(Jt.col(1)));
						Jt.putCol(1, Jt.col(1)/norm(Jt.col(1)));
					}
					if(kk==49 && calcattr)
					{
						attrdata[attrct]=coord<double, 2>(IC);
						//inf<<attrdata[attrct-1]<<" "<<basindata.size()<<endl; 
						DrawPoint(master, inR, IC.x(),IC.y());
						++attrct;
					}
					++kk;
				}
				if(calclyp)
				{
					lypdata[lypct]=(lyp/(50));
					DrawPoint(master, inR, lyp.x(),lyp.y());
					++lypct;
				}
			}
			i++;
		}		
		j++;
	}
	if(calclyp)	lypdata.resizeAndPreserve(lypct);
	if(calcattr) attrdata.resizeAndPreserve(attrct);
	if(calcbasin) basindata.resizeAndPreserve(basinct);
	return 1;
}

#endif

