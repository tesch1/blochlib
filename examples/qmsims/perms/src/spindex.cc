

#ifndef _spin_tensors_index__cc_
#define _spin_tensors_index__cc_ 1

#include "spindex.h"

//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/*

This is a set of classes that generate and maintain massive lists of
various combinations of spin tensors to certain 'orders'
i.e.
1st order--> Ie, Ix, Iy, Iz
2nd Order--> Ix^2, Iy^2, Iz^2, Ix*Iy.... (al powers of 2)
3rd order-->all powers of three
etc

there are 2 brands...

1) the Cartiesian tensors (Ix, Iy, Iz representation)
2) the Spherical Tensors (T00, T10, T11, T1m1..)





*/

bool SpinIdxCartData::operator==(const SpinIdxCartData &rhs)
{
	if(&rhs==this)	return true;
	if(rhs.axis_.size() != axis_.size())	return false;
	for(int i=0;i<axis_.size();++i) if(axis_(i) != rhs.axis_(i)) return false;
	if(rhs.spins_.size() != spins_.size())	return false;
	for(int i=0;i<spins_.size();++i) if(spins_(i) != rhs.spins_(i)) return false;
	if(MulOps!=rhs.MulOps) 	return false;
	return true;
}

bool SpinIdxCartData::operator!=(const SpinIdxCartData &rhs)
{
	if(&rhs==this)	return false;
	if(rhs.axis_.size() != axis_.size())	return true;
	for(int i=0;i<axis_.size();++i) if(axis_(i) != rhs.axis_(i)) return true;
	if(rhs.spins_.size() != spins_.size())	return true;
	for(int i=0;i<spins_.size();++i) if(spins_(i) != rhs.spins_(i)) return true;
	if(MulOps!=rhs.MulOps) 	return true;
	return false;
}

//uses the length of the 'axis' to determin this inequlity
bool operator<(const SpinIdxCartData &lhs,const SpinIdxCartData &rhs)
{
	if(&rhs==&lhs)	return false;
	if(lhs.axis_.size()<rhs.axis_.size())	return true;
	return false;
}

//uses the length of the 'axis' to determin this inequlity
bool operator>(const SpinIdxCartData &lhs,const SpinIdxCartData &rhs)
{
	if(&rhs==&lhs)	return false;
	if(lhs.axis_.size()>rhs.axis_.size())	return true;
	return false;
}


void SpinIdxCartData::operator=(const SpinIdxCartData &rhs)
{
	if(&rhs==this)	return;
	axis_=rhs.axis_;
	spins_=rhs.spins_;
	MulOps=rhs.MulOps;
}

/******* Set parameters Vectors ***/
void SpinIdxCartData::setSpin(int i)
{	spins_.resize(1);		spins_=i; }

void SpinIdxCartData::setSpin(const Vector<int> &in)
{	spins_=in;	}

void SpinIdxCartData::setAxes(const Vector<int> &in)
{	axis_=in;	}

void SpinIdxCartData::setOptions(int ops)
{	MulOps=ops;	}

/********* Nameing generator ***/
std::string SpinIdxCartData::name() const
{
	std::string nam="", tms="", tms2="";

	if(MulOps & SpinIdxOps::Commutators && axis_.size()>1){
	//the commutator "[ spinten name, spinten name ]"
	// for the first spin is different then the rest...
		tms="";
		switch(axis_[0])
		{
			case Tensors::x:	tms+="_x";	break;
			case Tensors::y:	tms+="_y";	break;
			case Tensors::z:	tms+="_z";	break;
			default: break;
		}
		if(spins_[0]!=-1){
			tms="I"+tms+"^"+itost(spins_[0]);
		}else{
			tms="F"+tms;
		}
		nam+=tms;
	//the second spin
		for(int i=1;i<axis_.size();++i){
			tms="";
			switch(axis_[i])
			{
				case Tensors::x:	tms2="_x";	break;
				case Tensors::y:	tms2="_y";	break;
				case Tensors::z:	tms2="_z";	break;
				default: break;
			}
			if(spins_[i]!=-1){
				tms="I"+tms+"^"+itost(spins_[i]);
			}else{
				tms="F"+tms;
			}
			nam="["+tms+","+nam+"]";
		}
	}else{
		for(int i=0;i<axis_.size();++i){
			tms="";
			switch(axis_[i])
			{
				case Tensors::x:	tms+="_x";	break;
				case Tensors::y:	tms+="_y";	break;
				case Tensors::z:	tms+="_z";	break;
				default: break;
			}
			if(spins_[i]!=-1){
				tms="I"+tms+"^"+itost(spins_[i]);
			}else{
				tms="F"+tms;
			}
			if(MulOps & SpinIdxOps::Add && i!=axis_.size()-1) tms+="+";
			nam+=tms;
		}
	}


	return nam;
}

/********* Tensor generator ***/
matrix SpinIdxCartData::tensor(SpinSys& in)
{
	matrix mat=(in.Fe());
	if(MulOps & SpinIdxOps::Add) mat=in.F0();
	matrix spinI;

	for(int i=0;i<axis_.size();++i){
		//spinI=in.F0();
		if(spins_[i]!=-1){
			switch(axis_[i])
			{
				case Tensors::x:	spinI = in.Ix(spins_[i]);	break;
				case Tensors::y:	spinI = in.Iy(spins_[i]);	break;
				case Tensors::z:	spinI = in.Iz(spins_[i]);	break;
				default: break;
			}
		}else{
			switch(axis_[i])
			{
				case Tensors::x:	spinI = in.Fx();	break;
				case Tensors::y:	spinI = in.Fy();	break;
				case Tensors::z:	spinI = in.Fz();	break;
				default: break;
			}
		}
		if(MulOps & SpinIdxOps::Add){ mat+=spinI;}
		else if(MulOps & SpinIdxOps::Commutators && i==0){	mat=spinI;	}
		else if(MulOps & SpinIdxOps::Commutators){	mat=spinI*mat-mat*spinI;	}
		else {	mat*=spinI;	}
	}
	return mat;
}

matrix SpinIdxCartData::operator()(SpinSys &in)
{	return tensor(in);	}

//I/O
int SpinIdxCartData::binarySize()
{
	return 3*sizeof(char)+
			4*sizeof(int)+ //<spins ele+axis ele+options>
	       spins_.size()*sizeof(int)+ //<spins list>
	       axis_.size()*sizeof(int); //axis list
}

//format is...
// "SCD<binarysize><numele><spin1>...<spin2><axis1>...<axis2>"
void SpinIdxCartData::write(fstream &oo) //binary write
{
	static const char *sym="SCD";
	int siz=binarySize();
	oo.write(&sym[0], sizeof(char));
	oo.write(&sym[1], sizeof(char));
	oo.write(&sym[2], sizeof(char));
	oo.write((char *)&siz, sizeof(int));

	oo.write((char *)&MulOps, sizeof(int));

	siz=spins_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<spins_.size();++i){
		siz=spins_[i];
		oo.write((char *)&siz, sizeof(int));
	}

	siz=axis_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<axis_.size();++i){
		siz=axis_[i];
		oo.write((char *)&siz, sizeof(int));
	}
}

void SpinIdxCartData::read(fstream &oo) //binary read
{
	static char sym[3];
	int siz;
	int pos=oo.tellg();
	oo.read(&sym[0], sizeof(char));
	oo.read(&sym[1], sizeof(char));
	oo.read(&sym[2], sizeof(char));
	if(sym[0]=='S' && sym[1]=='C' && sym[2]=='D')
	{
		oo.read((char *)&siz, sizeof(int)); //num bytes
		oo.read((char *)&MulOps, sizeof(int)); //options

		oo.read((char *)&siz, sizeof(int)); //spins size
		spins_.resize(siz);
		for(int i=0;i<spins_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			spins_[i]=siz;
		}

		oo.read((char *)&siz, sizeof(int)); //spins size
		axis_.resize(siz);
		for(int i=0;i<axis_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			axis_[i]=siz;
		}

	}else{
		oo.seekg(pos); //go back the 3 steps we just did
	}
}


std::ostream &operator<<(std::ostream &oo,const SpinIdxCartData &out)
{
	oo<<" Spin Index Cartesian Tensor: "<<out.name()<<std::endl;
	return oo;
}

std::fstream &operator<<(std::fstream &oo, SpinIdxCartData &out)
{
	out.write(oo);
	return oo;
}

std::fstream &operator>>(std::fstream &oo, SpinIdxCartData &out)
{
	out.read(oo);
	return oo;
}

/***************************************************************/
/******** SPHERICAL DATA ELEMENT **************/
/***************************************************************/

bool SpinIdxSphData::operator==(const SpinIdxSphData &rhs)
{
	if(&rhs==this)	return true;
	if(MulOps!=rhs.MulOps) 	return false;
	if(rhs.axis_.size() != axis_.size())	return false;
	for(int i=0;i<axis_.size();++i) if(axis_(i) != rhs.axis_(i)) return false;
	if(rhs.spins1_.size() != spins1_.size())	return false;
	for(int i=0;i<spins1_.size();++i) if(spins1_(i) != rhs.spins1_(i)) return false;
	if(rhs.spins2_.size() != spins2_.size())	return false;
	for(int i=0;i<spins2_.size();++i) if(spins2_(i) != rhs.spins2_(i)) return false;
	return true;
}

bool SpinIdxSphData::operator!=(const SpinIdxSphData &rhs)
{
	if(&rhs==this)	return false;
	if(MulOps!=rhs.MulOps) 	return true;
	if(rhs.axis_.size() != axis_.size())	return true;
	for(int i=0;i<axis_.size();++i) if(axis_(i) != rhs.axis_(i)) return true;
	if(rhs.spins1_.size() != spins1_.size())	return true;
	for(int i=0;i<spins1_.size();++i) if(spins1_(i) != rhs.spins1_(i)) return true;
	if(rhs.spins2_.size() != spins2_.size())	return true;
	for(int i=0;i<spins2_.size();++i) if(spins2_(i) != rhs.spins2_(i)) return true;
	return false;
}

bool operator<(const SpinIdxSphData &lhs,const SpinIdxSphData &rhs)
{
	if(&rhs==&lhs)	return false;
	if(lhs.axis_.size()<rhs.axis_.size())	return true;
	return false;
}

bool operator>(const SpinIdxSphData &lhs,const SpinIdxSphData &rhs)
{
	if(&rhs==&lhs)	return false;
	if(lhs.axis_.size()>rhs.axis_.size())	return true;
	return false;
}

void SpinIdxSphData::operator=(const SpinIdxSphData &rhs)
{
	if(&rhs==this)	return;
	axis_=rhs.axis_;
	spins1_=rhs.spins1_;
	spins2_=rhs.spins2_;
	MulOps=rhs.MulOps;
}

/******* Set parameters Vectors ***/
void SpinIdxSphData::setSpin(int i)
{
	spins1_=i;
	spins2_=i;
}

/******* Set parameters Vectors ***/
void SpinIdxSphData::setSpin(int i, int j)
{
	spins1_=i;
	spins2_=j;
}

void SpinIdxSphData::setSpin(const Vector<int> &in, const Vector<int> &sp2)
{	spins1_=in;	spins2_=sp2;	}

void SpinIdxSphData::setAxes(const Vector<int> &in)
{		axis_=in;	}

void SpinIdxSphData::setOptions(int ops)
{		MulOps=ops;	}

/********* Nameing generator ***/
std::string SpinIdxSphData::name() const
{
	std::string nam="", tms="";
	RunTimeAssert(spins1_.size()==spins2_.size());

	if(MulOps & SpinIdxOps::Commutators && axis_.size()>1){
		if(spins1_.size()>1) nam+="(";
		tms="";
		switch(axis_[0])
		{
				case Tensors::T1m1:	tms+="T_{1,-1}";	break;
				case Tensors::T10:	tms+="T_{1,0}";	break;
				case Tensors::T11:	tms+="T_{1,1}";	break;
				case Tensors::T2m2:	tms+="T_{2,-2}";	break;
				case Tensors::T2m1:	tms+="T_{2,-1}";	break;
				case Tensors::T20:	tms+="T_{2,0}";	break;
				case Tensors::T21:	tms+="T_{2,1}";	break;
				case Tensors::T22:	tms+="T_{2,2}";	break;
		}
		if(spins1_[0]!=-1 && axis_[0] & Tensors::T1){
			tms=tms+"^"+itost(spins1_[0])+" ";
		}
		if(spins1_[0]!=-1 && axis_[0] & Tensors::T2){
			tms=tms+"^{"+itost(spins1_[0])+","+itost(spins2_[0])+"} ";
		}
		nam+=tms;
		for(int i=1; i<axis_.size();++i){
			tms="";
			switch(axis_[i])
			{
				case Tensors::T1m1:	tms+="T_{1,-1}";	break;
				case Tensors::T10:	tms+="T_{1,0}";	break;
				case Tensors::T11:	tms+="T_{1,1}";	break;
				case Tensors::T2m2:	tms+="T_{2,-2}";	break;
				case Tensors::T2m1:	tms+="T_{2,-1}";	break;
				case Tensors::T20:	tms+="T_{2,0}";	break;
				case Tensors::T21:	tms+="T_{2,1}";	break;
				case Tensors::T22:	tms+="T_{2,2}";	break;
			}
			if(spins1_[i]!=-1 && axis_[i] & Tensors::T1){
				tms=tms+"^"+itost(spins1_[i])+" ";
			}
			if(spins1_[i]!=-1 && axis_[i] & Tensors::T2){
				tms=tms+"^{"+itost(spins1_[i])+","+itost(spins2_[i])+"} ";
			}
			nam="["+tms+","+nam+"]";
		}
	}else{
		for(int i=0;i<axis_.size();++i){
			tms="";
			switch(axis_[i])
			{
				case Tensors::T1m1:	tms+="T_{1,-1}";	break;
				case Tensors::T10:	tms+="T_{1,0}";	break;
				case Tensors::T11:	tms+="T_{1,1}";	break;
				case Tensors::T2m2:	tms+="T_{2,-2}";	break;
				case Tensors::T2m1:	tms+="T_{2,-1}";	break;
				case Tensors::T20:	tms+="T_{2,0}";	break;
				case Tensors::T21:	tms+="T_{2,1}";	break;
				case Tensors::T22:	tms+="T_{2,2}";	break;
			}
			if(spins1_[i]!=-1 && axis_[i] & Tensors::T1){
				nam+=tms+"^"+itost(spins1_[i])+" ";
			}
			if(spins1_[i]!=-1 && axis_[i] & Tensors::T2){
				nam+=tms+"^{"+itost(spins1_[i])+","+itost(spins2_[i])+"} ";
			}

			if(MulOps & SpinIdxOps::Add && i != axis_.size()-1) nam+="+";
		}
	}
	return nam;
}

/********* Tensor generator ***/
matrix SpinIdxSphData::tensor(SpinSys& in)
{
	matrix mat=(in.Fe());
	if(MulOps & SpinIdxOps::Add) mat=in.F0();

	matrix spinI=(in.F0());
//	RunTimeAssert(spins1_.size()==axis_.size());
	RunTimeAssert(spins2_.size()==spins1_.size());

	for(int i=0;i<axis_.size();++i){
		spinI=in.F0();
		if(spins1_[i]!=-1){
			if(axis_[i] & Tensors::T1){
				switch(axis_[i])
				{
					case Tensors::T1m1:	spinI += T1(in, spins1_[i], -1);	break;
					case Tensors::T10:	spinI += T1(in, spins1_[i], 0);	break;
					case Tensors::T11:	spinI += T1(in, spins1_[i], 1);	break;
				}
			}else{
				switch(axis_[i])
				{
					case Tensors::T2m2:	spinI += T2(in, spins1_[i],spins2_[i], -2);	break;
					case Tensors::T2m1: 	spinI += T2(in, spins1_[i],spins2_[i], -1);	break;
					case Tensors::T20:	spinI += T2(in, spins1_[i],spins2_[i], 0);	break;
					case Tensors::T21:	spinI += T2(in, spins1_[i],spins2_[i], 1);	break;
					case Tensors::T22:	spinI += T2(in, spins1_[i],spins2_[i], 2);	break;
				}
			}
		}else{
			if(axis_[i] & Tensors::T1){

				switch(axis_[i])
				{
					case Tensors::T1m1:
						for(int k=0;k<in.size();++k)		spinI+=T1(in, k, -1);
						break;

					case Tensors::T10:
						for(int k=0;k<in.size();++k)		spinI+=T1(in, k, 0);
						break;

					case Tensors::T11:
						for(int k=0;k<in.size();++k)		spinI+=T1(in, k, 1);
						break;
				}
			}else{
				switch(axis_[i])
				{
					case Tensors::T2m2:
						for(int j=0;j<in.size();++j)
							for(int k=0;k<in.size();++k)		spinI+=T2(in, k,j, -2);

						break;

					case Tensors::T2m1:
						for(int j=0;j<in.size();++j)
							for(int k=0;k<in.size();++k)		spinI+=T2(in, k,j, -1);

						break;
					case Tensors::T20:
						for(int j=0;j<in.size();++j)
							for(int k=0;k<in.size();++k)		spinI+=T2(in, k,j, 0);

						break;
					case Tensors::T21:
						for(int j=0;j<in.size();++j)
							for(int k=0;k<in.size();++k)		spinI+=T2(in, k,j, 1);

						break;
					case Tensors::T22:
						for(int j=0;j<in.size();++j)
							for(int k=0;k<in.size();++k)		spinI+=T2(in, k,j, 2);

						break;
				}
			}
		}
		if(MulOps & SpinIdxOps::Add){ mat+=spinI;	}
		else if(MulOps & SpinIdxOps::Commutators && i==0){	mat=spinI;	}
		else if(MulOps & SpinIdxOps::Commutators){ mat=spinI*mat-mat*spinI;	}
		else{ mat*=spinI;	}
	}

	return mat;
}

matrix SpinIdxSphData::operator()(SpinSys &in)
{	return tensor(in);	}


//I/O
int SpinIdxSphData::binarySize()
{
	return 3*sizeof(char)+ //"SSC"
		5*sizeof(int)+ //<spins ele+axis ele+numbytes+options>
	       spins1_.size()*sizeof(int)+ //<spins list>
	       spins2_.size()*sizeof(int)+ //<spins list>
	       axis_.size()*sizeof(int); //axis list
}

//format is...
// "SCD<binarysize><numele><numsp1><spinarray1><numaxis><axis1>...<axis2>"
void SpinIdxSphData::write(fstream &oo) //binary write
{
	static char *sym="SSD";
	int siz=binarySize();
	oo.write(&sym[0], sizeof(char));
	oo.write(&sym[1], sizeof(char));
	oo.write(&sym[2], sizeof(char));
	oo.write((char *)&siz, sizeof(int));

	oo.write((char *)&MulOps, sizeof(int));

	siz=spins1_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<spins1_.size();++i){
		siz=spins1_[i];
		oo.write((char *)&siz, sizeof(int));
	}

	siz=spins2_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<spins2_.size();++i){
		siz=spins2_[i];
		oo.write((char *)&siz, sizeof(int));
	}

	siz=axis_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<axis_.size();++i){
		siz=axis_[i];
		oo.write((char *)&siz, sizeof(int));
	}
}

void SpinIdxSphData::read(fstream &oo) //binary read
{
	static char sym[3];
	int siz;
	int pos=oo.tellg();

	oo.read(&sym[0], sizeof(char));
	oo.read(&sym[1], sizeof(char));
	oo.read(&sym[2], sizeof(char));
	if(sym[0]=='S' && sym[1]=='S' && sym[2]=='D')
	{
		oo.read((char *)&siz, sizeof(int)); //num bytes
		oo.read((char *)&MulOps, sizeof(int));

		oo.read((char *)&siz, sizeof(int)); //spins size
		spins1_.resize(siz);
		for(int i=0;i<spins1_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			spins1_[i]=siz;
		}

		oo.read((char *)&siz, sizeof(int)); //spins size
		spins2_.resize(siz);
		for(int i=0;i<spins2_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			spins2_[i]=siz;
		}

		oo.read((char *)&siz, sizeof(int)); //spins size
		axis_.resize(siz);
		for(int i=0;i<axis_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			axis_[i]=siz;
		}

	}else{
		oo.seekg(pos); //go back the 3 steps we just did
	}
}


std::ostream &operator<<(std::ostream &oo,const SpinIdxSphData &out)
{
	oo<<" Spin Index Spherical Tensor: "<<out.name()<<std::endl;
	return oo;
}


std::fstream &operator<<(std::fstream &oo,SpinIdxSphData &out)
{
	out.write(oo);
	return oo;
}

std::fstream &operator>>(std::fstream &oo,SpinIdxSphData &out)
{
	out.read(oo);
	return oo;
}

/***************************************************************/
/******** SPHERICAL DATA ELEMENT **************/
/***************************************************************/

bool SolidSysData::operator==(const SolidSysData &rhs)
{
	if(&rhs==this)	return true;
	if(MulOps!=rhs.MulOps) 	return false;
	if(rhs.Htypes_.size() != Htypes_.size())	return false;
	for(int i=0;i<Htypes_.size();++i) if(Htypes_(i) != rhs.Htypes_(i)) return false;
	if(rhs.indexs_.size() != indexs_.size())	return false;
	for(int i=0;i<indexs_.size();++i) if(indexs_(i) != rhs.indexs_(i)) return false;
	return true;
}

bool SolidSysData::operator!=(const SolidSysData &rhs)
{
	if(&rhs==this)	return false;
	if(MulOps!=rhs.MulOps) 	return true;
	if(rhs.Htypes_.size() != Htypes_.size())	return true;
	for(int i=0;i<Htypes_.size();++i) if(Htypes_(i) != rhs.Htypes_(i)) return true;
	if(rhs.indexs_.size() != indexs_.size())	return true;
	for(int i=0;i<indexs_.size();++i) if(indexs_(i) != rhs.indexs_(i)) return true;
	return false;
}

bool operator<(const SolidSysData &lhs,const SolidSysData &rhs)
{
	if(&rhs==&lhs)	return false;
	if(lhs.Htypes_.size()<rhs.Htypes_.size())	return true;
	return false;
}

bool operator>(const SolidSysData &lhs,const SolidSysData &rhs)
{
	if(&rhs==&lhs)	return false;
	if(lhs.Htypes_.size()>rhs.Htypes_.size())	return true;
	return false;
}

void SolidSysData::operator=(const SolidSysData &rhs)
{
	if(&rhs==this)	return;
	Htypes_=rhs.Htypes_;
	indexs_=rhs.indexs_;
	MulOps=rhs.MulOps;
}

/******* Set parameters Vectors ***/
void SolidSysData::setHtypes(int i)
{
	indexs_=i;
}

void SolidSysData::setHtypes(const Vector<int> &in)
{	indexs_=in;		}

void SolidSysData::setIndex(const Vector<int> &in)
{		Htypes_=in;	}

void SolidSysData::setOptions(int ops)
{		MulOps=ops;	}

/********* Nameing generator ***/
std::string SolidSysData::name(SolidSys &sys) const
{
	std::string nam="", tms="";
	RunTimeAssert(indexs_.size()==Htypes_.size());

	if(MulOps & SpinIdxOps::Commutators && Htypes_.size()>1){
		tms="";
		switch(Htypes_[0])
		{
			case HamilTypes::Csa:	tms+=sys.csa[indexs_[0]].name(); 	break;
			case HamilTypes::Dip:	tms+=sys.dip[indexs_[0]].name();	break;
			case HamilTypes::J:		tms+=sys.jcop[indexs_[0]].name();	break;
			case HamilTypes::Q:		tms+=sys.qua[indexs_[0]].name();	break;
		}
		nam+=tms;
		for(int i=1; i<Htypes_.size();++i){
			tms="";
			switch(Htypes_[i])
			{
				case HamilTypes::Csa:	tms+=sys.csa[indexs_[i]].name(); 	break;
				case HamilTypes::Dip:	tms+=sys.dip[indexs_[i]].name();	break;
				case HamilTypes::J:		tms+=sys.jcop[indexs_[i]].name();	break;
				case HamilTypes::Q:		tms+=sys.qua[indexs_[i]].name();	break;
			}
			nam="["+tms+","+nam+"]";
		}
	}else{
		for(int i=0;i<Htypes_.size();++i){
			tms="";
			switch(Htypes_[i])
			{
				case HamilTypes::Csa:	tms+=sys.csa[indexs_[i]].name(); 	break;
				case HamilTypes::Dip:	tms+=sys.dip[indexs_[i]].name();	break;
				case HamilTypes::J:		tms+=sys.jcop[indexs_[i]].name();	break;
				case HamilTypes::Q:		tms+=sys.qua[indexs_[i]].name();	break;
			}
			if(MulOps & SpinIdxOps::Add && i != Htypes_.size()-1) tms+="+";
			nam+=tms;
		}
	}
	return nam;
}

/********* Tensor generator ***/
matrix SolidSysData::tensor(SolidSys& sys)
{
	matrix mat=(sys.Fe());
	if(MulOps & SpinIdxOps::Add) mat=sys.F0();

	matrix spinI=(sys.F0());
	RunTimeAssert(Htypes_.size()==indexs_.size());

	for(int i=0;i<Htypes_.size();++i){
		switch(Htypes_[i])
		{
			case HamilTypes::Csa:	spinI=sys.csa[indexs_[i]].H(sys, sys.theRotations); 	break;
			case HamilTypes::Dip:	spinI=sys.dip[indexs_[i]].H(sys, sys.theRotations);	break;
			case HamilTypes::J:		spinI=sys.jcop[indexs_[i]].H(sys, sys.theRotations);	break;
			case HamilTypes::Q:		spinI=sys.qua[indexs_[i]].H(sys, sys.theRotations);	break;
		}
		if(MulOps & SpinIdxOps::Add){ mat+=spinI;	}
		else if(MulOps & SpinIdxOps::Commutators && i==0){	mat=spinI;	}
		else if(MulOps & SpinIdxOps::Commutators){ mat=spinI*mat-mat*spinI;	}
		else{ mat*=spinI;	}
	}
	return mat;
}

matrix SolidSysData::operator()(SolidSys &in)
{	return tensor(in);	}


//I/O
int SolidSysData::binarySize()
{
	return 3*sizeof(char)+ //"SOS"
		4*sizeof(int)+ //<axis ele+numbytes+options>
	       indexs_.size()*sizeof(int)+ //<spins list>
	       Htypes_.size()*sizeof(int); //axis list
}

//format is...
// "SOS<binarysize><numele><numsp1><spinarray1><numaxis><axis1>...<axis2>"
void SolidSysData::write(fstream &oo) //binary write
{
	static char *sym="SOS";
	int siz=binarySize();
	oo.write(&sym[0], sizeof(char));
	oo.write(&sym[1], sizeof(char));
	oo.write(&sym[2], sizeof(char));
	oo.write((char *)&siz, sizeof(int));

	oo.write((char *)&MulOps, sizeof(int));

	siz=indexs_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<indexs_.size();++i){
		siz=indexs_[i];
		oo.write((char *)&siz, sizeof(int));
	}

	siz=Htypes_.size();
	oo.write((char *)&siz, sizeof(int));
	for(int i=0;i<Htypes_.size();++i){
		siz=Htypes_[i];
		oo.write((char *)&siz, sizeof(int));
	}
}

void SolidSysData::read(fstream &oo) //binary read
{
	static char sym[3];
	int siz;
	int pos=oo.tellg();

	oo.read(&sym[0], sizeof(char));
	oo.read(&sym[1], sizeof(char));
	oo.read(&sym[2], sizeof(char));
	if(sym[0]=='S' && sym[1]=='O' && sym[2]=='S')
	{
		oo.read((char *)&siz, sizeof(int)); //num bytes
		oo.read((char *)&MulOps, sizeof(int));

		oo.read((char *)&siz, sizeof(int)); //spins size
		indexs_.resize(siz);
		for(int i=0;i<indexs_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			indexs_[i]=siz;
		}

		oo.read((char *)&siz, sizeof(int)); //spins size
		Htypes_.resize(siz);
		for(int i=0;i<Htypes_.size();++i){
			oo.read((char *)&siz, sizeof(int));
			Htypes_[i]=siz;
		}

	}else{
		oo.seekg(pos); //go back the 3 steps we just did
	}
}


std::ostream &operator<<(std::ostream &oo,const SolidSysData &out)
{
	oo<<" Spin Index Solid System: of size"<<out.size()<<std::endl;
	return oo;
}


std::fstream &operator<<(std::fstream &oo,SolidSysData &out)
{
	out.write(oo);
	return oo;
}

std::fstream &operator>>(std::fstream &oo,SolidSysData &out)
{
	out.read(oo);
	return oo;
}

/***************************************************************/
/******** The Permutation Tensor Generator **************/
/***************************************************************/
const int TensorGen::cartarr_[3]={Tensors::z,Tensors::y,Tensors::x};
const int TensorGen::totalspharr_[8]={Tensors::T1m1,Tensors::T10,Tensors::T11,
							Tensors::T2m2,Tensors::T2m1,Tensors::T20,
							Tensors::T21,Tensors::T22};
const int TensorGen::doublespharr_[5]={Tensors::T2m2,Tensors::T2m1,Tensors::T20,
							Tensors::T21,Tensors::T22};
const int TensorGen::singlespharr_[3]={Tensors::T1m1,Tensors::T10,Tensors::T11};

//generates the initial permutation lists
int TensorGen::permutationLength()
{
	int mulfact=0;
	int tmpf=0;
	if(options_ & TensorGenOps::I){
		if(spinperms1_.size()==0 || spinperms1_[0]==-1){
			mulfact=3;tmpf=mulfact;
		}else{
			mulfact=3*spinperms1_.size();
			tmpf=mulfact;
		}
		while(mulfact<order_){	 mulfact+=tmpf;	}
	}else if(options_ & TensorGenOps::T1){
		if(spinperms1_.size()==0 || spinperms1_[0]==-1){
			mulfact=3;tmpf=mulfact;
		}else{
			mulfact=3*spinperms1_.size();
			tmpf=mulfact;
		}
		while(mulfact<order_){	 mulfact+=tmpf;	}
	}else if(options_ & TensorGenOps::T2){
		if(spinperms1_.size()==0 || spinperms1_[0]==-1){
			mulfact=5;tmpf=mulfact;
		}else{
			mulfact=5*spinperms1_.size();
			tmpf=mulfact;
		}
		while(mulfact<order_){	 mulfact+=tmpf;	}
	}else if(options_ & TensorGenOps::T){
		if(spinperms1_.size()==0 || spinperms1_[0]==-1){
			mulfact=8;tmpf=mulfact;
		}else{
			mulfact=8*spinperms1_.size();
			tmpf=mulfact;
		}
		while(mulfact<order_){	 mulfact+=tmpf;	}
	}
	return mulfact;
}

//the length of the perm list for the Hamiltonian type
// the 'spin' list allows me to specifiy only interactions
//coming from those spins...if the first element is
// '-1' then ALL the interactions are used
int TensorGen::permutationLength(SolidSys &sys)
{
	int mulfact=0, i,j;
	if(options_ & TensorGenOps::Hamiltonians){
		if(spinperms1_.size()==0 || spinperms1_[0]==-1){
			mulfact=sys.csa.size()+sys.dip.size()+sys.jcop.size()+sys.qua.size();
		}else{
			for(i=0;i<spinperms1_.size();++i){
				for(j=0;j<sys.csa.size();++j){
					if(sys.csa[j].on()==spinperms1_[i]) mulfact++;
				}
				for(j=0;j<sys.dip.size();++j){
					if(sys.dip[j].on1()==spinperms1_[i] || sys.dip[j].on2()==spinperms1_[i] ) mulfact++;
				}
				for(j=0;j<sys.jcop.size();++j){
					if(sys.jcop[j].on1()==spinperms1_[i] || sys.jcop[j].on2()==spinperms1_[i] ) mulfact++;
				}
				for(j=0;j<sys.qua.size();++j){
					if(sys.qua[j].on()==spinperms1_[i]) mulfact++;
				}
			}
		}
	}else{
		std::cerr<<std::endl<<"Error: permutationLength(SolidSys &)"<<std::endl;
		std::cerr<<" tensor type is Not 'Hamiltonian'...wrong function"<<std::endl;
		std::cerr<<" returning 'permutationLength()'"<<std::endl;
		return permutationLength();
	}
	int tmpf=mulfact;
	while(mulfact<order_) mulfact+=tmpf;
	return mulfact;
}

//generates the initial permutation lists
void TensorGen::genPermutationList(int mulfact, int *theperms, int *spinx, int *spinx2)
{
	if(options_ & TensorGenOps::I){
		int ct=0;
		while(ct<mulfact){
			if(spinperms1_.size()==0 || spinperms1_[0]==-1){
				theperms[ct]=cartarr_[ct%3];
				spinx[ct]=-1;
				++ct;
			}else{
				for(int j=0;j<3;++j){
					for(int i=0;i<spinperms1_.size();++i){
						spinx[ct]=spinperms1_[i];
						theperms[ct]=cartarr_[j];
						++ct;
					}
				}
			}
		}
	}else if(options_ & TensorGenOps::T1){
		int ct=0;
		while(ct<mulfact){
			if(spinperms1_.size()==0 || spinperms1_[0]==-1){
				theperms[ct]=singlespharr_[ct%3];
				spinx[ct]=-1;
				++ct;
			}else{
				for(int j=0;j<3;++j){
					for(int i=0;i<spinperms1_.size();++i){
						spinx[ct]=spinperms1_[i];
						theperms[ct]=singlespharr_[j];
						++ct;
					}
				}
			}
		}
	}else if(options_ & TensorGenOps::T2){
		int ct=0;
		while(ct<mulfact){
			if(spinperms1_.size()==0 || spinperms1_[0]==-1){
				theperms[ct]=doublespharr_[ct%5];
				spinx[ct]=-1;
				spinx2[ct]=-1;
				++ct;
			}else{
				for(int j=0;j<5;++j){
					for(int i=0;i<spinperms1_.size();++i){
						spinx[ct]=spinperms1_[i];
						spinx2[ct]=spinperms2_[i];
						theperms[ct]=doublespharr_[j];
						++ct;
					}
				}
			}
		}
	}else if(options_ & TensorGenOps::T){
		int ct=0;
		while(ct<mulfact){
			if(spinperms1_.size()==0 || spinperms1_[0]==-1){
				theperms[ct]=totalspharr_[ct%8];
				spinx[ct]=-1;
				spinx2[ct]=-1;
				++ct;
			}else{
				for(int j=0;j<8;++j){
					for(int i=0;i<spinperms1_.size();++i){
						spinx[ct]=spinperms1_[i];
						spinx2[ct]=spinperms2_[i];
						theperms[ct]=totalspharr_[j];
						++ct;
					}
				}
			}
		}
	}
}

//generates the initial permutation For the hamiltonian types
void TensorGen::genPermutationList(SolidSys &sys, int mulfact, int *htypes, int *indexs)
{
	int ct=0, j,i;
	while(ct<mulfact){
		if(options_ & TensorGenOps::Hamiltonians){
			if(spinperms1_.size()==0 || spinperms1_[0]==-1){
				for(j=0;j<sys.csa.size();++j){
					htypes[ct]=HamilTypes::Csa;
					indexs[ct]=j;
					++ct;
				}
				for(j=0;j<sys.dip.size();++j){
					htypes[ct]=HamilTypes::Dip;
					indexs[ct]=j;
					++ct;
				}
				for(j=0;j<sys.jcop.size();++j){
					htypes[ct]=HamilTypes::J;
					indexs[ct]=j;
					++ct;
				}
				for(j=0;j<sys.qua.size();++j){
					htypes[ct]=HamilTypes::Q;
					indexs[ct]=j;
					++ct;
				}
			}else{
				for(i=0;i<spinperms1_.size();++i){
					for(j=0;j<sys.csa.size();++j){
						if(sys.csa[j].on()==spinperms1_[i]){
							htypes[ct]=HamilTypes::Csa;
							indexs[ct]=j;
							++ct;
						}
					}
					for(j=0;j<sys.dip.size();++j){
						if(sys.dip[j].on1()==spinperms1_[i] || sys.dip[j].on2()==spinperms1_[i] ){
							htypes[ct]=HamilTypes::Dip;
							indexs[ct]=j;
							++ct;
						}
					}
					for(j=0;j<sys.jcop.size();++j){
						if(sys.jcop[j].on1()==spinperms1_[i] || sys.jcop[j].on2()==spinperms1_[i] ){
							htypes[ct]=HamilTypes::J;
							indexs[ct]=j;
							++ct;
						}
					}
					for(j=0;j<sys.qua.size();++j){
						if(sys.qua[j].on()==spinperms1_[i]){
							htypes[ct]=HamilTypes::Q;
							indexs[ct]=j;
							++ct;
						}
					}
				}
			}
		}else{
			std::cerr<<std::endl<<"Error: genPermutationList(SolidSys &, int)"<<std::endl;
			std::cerr<<" tensor type in Not 'Hamiltonian'...wrong function"<<std::endl;
			std::cerr<<" can do nothing.."<<std::endl;
			return;
		}
	}
}


int TensorGen::MulOption()
{
	if(options_ & TensorGenOps::Commutators){	return TensorGenOps::Commutators;	}
	else if(options_ & TensorGenOps::Add){	return TensorGenOps::Add;	}
	else{	return TensorGenOps::Multiply;	}
}

//generates the next subset of length 'k' from
// the aster list 'master' of length n
void TensorGen::nextKSubset(int n, int k, int *subset,bool &more)
{
	static int m2, m;
	if(n<0 || k<0) return;
	if( !more ){
		m2 = 0;
		m = k;
	}else{
		if( m2 < n-m )  m = 0;
		++m;
		m2 = subset[k-m];
	}
	for(int i=0;i<m;++i) subset[k+i-m] = m2 + i+1;
	more = subset[0] != n-k+1;
}


void TensorGen::generateVectors()
{
	if(spinperms1_.size() != spinperms2_.size()){
		std::cerr<<std::endl<<"Error: TensorGen::generate() "<<std::endl;
		std::cerr<<" both spins lists must be the same length... "<<std::endl;
		exit(0);
	}
	//figure out the multiply options
	int mulops=MulOption();
	int permlength=permutationLength();
	int *theperms=new int[permlength];
	int *spinidx=new int[permlength];
	int *spinidx2=new int[permlength];
	genPermutationList(permlength,theperms, spinidx, spinidx2);


	int maxsize=100000;
	if(options_ & TensorGenOps::I) cart_.resize(maxsize);
	else sphere_.resize(maxsize);

//the order loop loop from 1..order_
	int curorder=1;
	do{
		bool more=false;
		int *subset=new int[curorder];
	//this is the k-subset loop;
		do{
			nextKSubset(permlength, curorder, subset, more);
			int *perms=new int[curorder];
			for(int i=0;i<curorder;++i) perms[i]=subset[i]-1;
			std::sort(perms, perms+curorder);

			do{ //this is the permutation on the k-subset loop

			//fill up our tensor vector
				Vector<int> tmax(curorder);
				Vector<int> spi(curorder);
				for(int i=0;i<curorder;++i){
					tmax[i]=theperms[perms[i]];
					spi[i]=spinidx[perms[i]];
				}
			//add the tensor to the correct list (based the users choise)
				if(options_ & TensorGenOps::I){
					SpinIdxCartData tmpc;
					if(options_ & TensorGenOps::SingleSpin){
						tmpc=SpinIdxCartData(spinperms1_[0], tmax, mulops);
					}else if(options_ & TensorGenOps::TotalSpace){ // the total spins space
						tmpc=SpinIdxCartData(-1, tmax, mulops);
					}else if(options_ & TensorGenOps::ChooseSpins){ // the user chossen set of spins
						tmpc=SpinIdxCartData(spi, tmax, mulops);
					}
					if(size_<maxsize){
						cart_[size_]=(tmpc); ++size_;
					}else{
						maxsize*=2;
						cart_.resizeAndPreserve(maxsize);
						cart_[size_]=(tmpc); ++size_;
					}
				}else{
					SpinIdxSphData tmpc;
					if(options_ & TensorGenOps::SingleSpin){
						tmpc=SpinIdxSphData(spinperms1_[0], tmax, mulops);
					}else if(options_ & TensorGenOps::TotalSpace){ // the total spins space
						tmpc=SpinIdxSphData(-1, tmax, mulops);
					}else if(options_ & TensorGenOps::ChooseSpins){ // the user chossen set of spins
						Vector<int> spi2(curorder);
						for(int i=0;i<curorder;++i) spi2[i]=spinidx2[perms[i]];
						tmpc=SpinIdxSphData(spi, spi2,tmax, mulops);
					}
					if(size_<maxsize){
						sphere_[size_]=(tmpc); ++size_;
					}else{
						maxsize*=2;
						sphere_.resizeAndPreserve(maxsize);
						sphere_[size_]=(tmpc); ++size_;
					}
				}
			}while(next_permutation(perms, perms+curorder)); //end the permutation loop
			delete [] perms;
		}while(more); //end subset loop

		delete [] subset;
		++curorder;
	}while(curorder<=order_); //end order loop
	delete [] theperms;
	delete [] spinidx;
	if(options_ & TensorGenOps::I)	 cart_.resizeAndPreserve(size_);
	else sphere_.resizeAndPreserve(size_);

}

void TensorGen::generate()
{
	size_=0;
	cart_.resize(0);
	sphere_.resize(0);

	if(order_<=0) return;
//Single Spin Tensor Generations...
	generateVectors();
	//TensorGen::sort();
}


//generates the indexing list for the hamiltonian types
void TensorGen::generateVectors(SolidSys &sys)
{
	if(spinperms1_.size() != spinperms2_.size()){
		std::cerr<<std::endl<<"Error: TensorGen::generate() "<<std::endl;
		std::cerr<<" both spins lists must be the same length... "<<std::endl;
		exit(0);
	}
	//figure out the multiply options
	int mulops=MulOption();
	int permlength=permutationLength(sys);
	int *thepermH=new int[permlength]; //Hytpes
	int *thepermI=new int[permlength]; //indexes
	genPermutationList(sys, permlength, thepermH,thepermI);

	//std::sort(theperms, theperms+permlength);
	std::cout<<std::endl<<"TensorGen Info: There will be ~"<<sum_factorial(order_)<<" terms to calculate "<<std::endl;

	int maxsize=100000;
	solidsys_.resize(maxsize);

//the order loop loop from 1..order_
	int curorder=1;
	do{
		bool more=false;
		int *subset=new int[curorder];
	//this is the k-subset loop;
		do{
			nextKSubset(permlength, curorder, subset, more);
			int *perms=new int[curorder];
			for(int i=0;i<curorder;++i) perms[i]=subset[i]-1;
			std::sort(perms, perms+curorder);

			do{ //this is the permutation on the k-subset loop

			//fill up our tensor vector
				Vector<int> Hp(curorder),Ip(curorder);
				for(int i=0;i<curorder;++i){
					Hp[i]=thepermH[perms[i]];
					Ip[i]=thepermI[perms[i]];
				}
			//add the tensor to the correct list (based the users choise)
				SolidSysData tmpc;
				tmpc=SolidSysData(Hp, Ip, mulops);
				if(size_<maxsize){
					solidsys_[size_]=(tmpc); ++size_;
				}else{
					maxsize*=2;
					solidsys_.resizeAndPreserve(maxsize);
					solidsys_[size_]=(tmpc); ++size_;
				}

			}while(next_permutation(perms, perms+curorder)); //end the permutation loop
			delete [] perms;
		}while(more); //end subset loop

		delete [] subset;
		++curorder;
	}while(curorder<=order_); //end order loop
	delete [] thepermH;
	delete [] thepermI;
	solidsys_.resizeAndPreserve(size_);

}

void TensorGen::generate(SolidSys &sys)
{
	size_=0;
	cart_.resize(0);
	sphere_.resize(0);
	solidsys_.resize(0);
	//cout<<"Hamils: "<<TensorGenOps::Hamiltonians<<endl;
	//cout<<"options & Hamils: "<<int(options_ & TensorGenOps::Hamiltonians)<<endl;
	if(order_<=0) return;
//Single Spin Tensor Generations...
	if(options_ & TensorGenOps::Hamiltonians){
		generateVectors(sys);
	}else{
		generateVectors();
	}
}


//NO Nedd for 'sort anymore, the generating alogorithm sorts automatically
void TensorGen::sort()
{
	/*std::sort(cart_.data(), cart_.data()+cart_.size());
	std::sort(sphere_.data(), sphere_.data()+sphere_.size());
	*/
}

//size of the subelemenets
int TensorGen::size(int i)const
{
	if(options_ & TensorGenOps::I) return cart_[i].size();
	else if(options_ & TensorGenOps::Hamiltonians) return solidsys_[i].size();
	else return sphere_[i].size();
}

//set the parameters
void TensorGen::setSpins(int i)
{
	spinperms1_.resize(1);
	spinperms1_=i;
//	generate();
}

void TensorGen::setSpins(int i, int j)
{
//	if(TensorGenOps::TwoSpins){
//		spinperms1_.resize(2);
//		spinperms2_.resize(2);
//	}
	spinperms1_=i;
	spinperms2_=j;
//	generate();
}

void TensorGen::setSpins(Vector<int> i)
{
	spinperms1_=i;
	spinperms2_=i;
//	generate();
}

void TensorGen::setSpins(Vector<int> i, Vector<int> j)
{
	if(i.size() != j.size()){
		std::cerr<<std::endl<<"Error: TensorGen::setSpins(Vector<int> i, Vector<int> j)"<<std::endl;
		std::cerr<<" input vectors must be the same size...."<<std::endl;
		exit(0);
	}
	spinperms1_=i;
	spinperms2_=j;
//	generate();
}


//get the string name of the operator
std::string TensorGen::name(int i)
{
	RunTimeAssert(i<size_);
	if(options_ & TensorGenOps::I) return cart_[i].name();
	else return sphere_[i].name();
}

//get the string name of the operator
std::string TensorGen::name(SolidSys &sys, int i)
{
	RunTimeAssert(i<size_);
	if(options_ & TensorGenOps::Hamiltonians) return solidsys_[i].name(sys);
	else return name(i);
}
//get the operator from the list
matrix TensorGen::tensor(SpinSys &in, int i)
{
	RunTimeAssert(i<size_);
	if(options_ & TensorGenOps::I) return cart_[i].tensor(in);
	else	return sphere_[i].tensor(in);
}

matrix TensorGen::operator()(SpinSys &in,int i)
{	return tensor(in,i);	}


//get the operator from the list
matrix TensorGen::tensor(SolidSys &in, int i)
{
	RunTimeAssert(i<size_);
	if(options_ & TensorGenOps::Hamiltonians) return solidsys_[i].tensor(in);
	else	return tensor(*in.spinSys(), i);
}

matrix TensorGen::operator()(SolidSys &in,int i)
{	return tensor(in,i);	}


//i/o
int TensorGen::binarySize()
{
	int tmp=sizeof(int); //cart or sphere & spinoptions
	if(options_ & TensorGenOps::I){
		tmp+=cart_.size()*sizeof(int); //num elements
		for(int i=0;i<cart_.size();++i) tmp+=cart_[i].binarySize();
	}else if(options_ & TensorGenOps::Hamiltonians){
		tmp+=solidsys_.size()*sizeof(int); //num elements
		for(int i=0;i<solidsys_.size();++i) tmp+=solidsys_[i].binarySize();
	}else{
		tmp+=sphere_.size()*sizeof(int); //num elements
		for(int i=0;i<sphere_.size();++i) tmp+=sphere_[i].binarySize();
	}
	return tmp;
}

//binary i/o...format
//BEGIN TensorGen
//<cart or sph &spin options><numele><the indivdual data elements>
//END TensorGen
void TensorGen::write(std::fstream &oo)
{
	oo.seekg(0, ios::end); //go to end of file
	oo<<"\nBEGIN TensorGen\n";
	oo.write((char *)&options_, sizeof(int));
	int siz;
	if(options_ & TensorGenOps::I){
		siz=cart_.size();
		oo.write((char *)&siz, sizeof(int));
		for(int i=0;i<siz;++i){
			oo<<cart_[i];
		}
	}else if(options_ & TensorGenOps::Hamiltonians){
		siz=solidsys_.size();
		oo.write((char *)&siz, sizeof(int));
		for(int i=0;i<siz;++i){
			oo<<solidsys_[i];
		}
	}else{
		siz=sphere_.size();
		oo.write((char *)&siz, sizeof(int));
		for(int i=0;i<siz;++i){
			oo<<sphere_[i];
		}
	}
	oo<<"\nEND TensorGen\n";
}

void TensorGen::read(std::fstream &oo)
{
	char liner[100];
	while(!oo.eof()){
		oo.getline(liner, 100, '\n');
		if(std::string(liner)=="BEGIN TensorGen")	break;
	}

	oo.read((char *)&options_, sizeof(int));
	if(options_ & TensorGenOps::I){
		oo.read((char *)&size_, sizeof(int));
		cart_.resize(size_);
		for(int i=0;i<size_;++i){
			oo>>cart_[i];
		}
	 }else if(options_ & TensorGenOps::Hamiltonians){
		oo.read((char *)&size_, sizeof(int));
		solidsys_.resize(size_);
		for(int i=0;i<size_;++i){
			oo>>solidsys_[i];
		}
	}else{
		oo.read((char *)&size_, sizeof(int));
		sphere_.resize(size_);
		for(int i=0;i<size_;++i){
			oo>>sphere_[i];
		}
	}
	while(!oo.eof()){
		oo.getline(liner, 100, '\n');
		if(std::string(liner)=="END TensorGen")	break;
	}
}
//binary out operator
std::fstream &operator<<(std::fstream &oo, TensorGen &out)
{
	out.write(oo);
	return oo;
}

//binary read operator
std::fstream &operator>>(std::fstream &oo, TensorGen &out)
{
	out.read(oo);
	return oo;
}



#endif


