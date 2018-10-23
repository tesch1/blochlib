

#include "blochlib.h"
//the required 2 namespaces
using namespace BlochLib;
using namespace std;

/**
a little file that goes into the
invokeing of functions from another proc
**/



template<class T>
struct RemoveReference
{
  typedef T type_t;
};

template<class T>
struct RemoveReference<T&>
{
  typedef T type_t;
};


//class for void functions (no inputs)
template<class Return>
class FunctionHolder0
{
	Return (*func_t) () ;

	public:

		FunctionHolder0(Return (*func)())
		{
			func_t=func;
		}
		Return run()
		{
			return (*func_t)();
		}
};

//class for one input
template<class Return, class T1>
class FunctionHolder1
{
	Return (*func_t) (T1) ;

	public:
		FunctionHolder1(Return (*func)(T1))
		{
			func_t=func;
		}

		Return run(RemoveReference<T1> &thing)
		{
			return (*func_t)(thing);
		}
};

//class for two inputs
template<class Return, class T1, class T2>
class FunctionHolder2
{
	Return (*func_t) (T1, T2) ;

	public:
		FunctionHolder2(Return (*func)(T1, T2))
		{
			func_t=func;
		}

		Return run(RemoveReference<T1> &thing, RemoveReference<T2> &thing2)
		{
			return (*func_t)(thing, thing2);
		}
};


//receive function for Void funcitons
template<class Return>
Return recieve(int FromProc, int Tag, Return (*func_t)())
{
	return (*func_t)();
}

//receive function for single input funcitons
template<class Return, class T1>
Return recieve(int FromProc, int Tag, Return (*func_t)(T1))
{
	T1 thing;
	MPIworld.get(thing, FromProc, Tag);
	return (*func_t)(thing);
}

//receive function for double input funcitons
template<class Return, class T1, class T2>
Return recieve(int FromProc, int Tag, Return (*func_t)(T1, T2))
{
	T1 thing;
	MPIworld.get(thing, FromProc, Tag);
	T2 thing2;
	MPIworld.get(thing2, FromProc, Tag+1);
	return (*func_t)(thing, thing2);
}


//send function for Void funcitons
void send(int To, int Tag)
{	}

//receive function for single input funcitons
template<class T1>
void send(int To, int Tag, T1 &thing)
{
	MPIworld.put(thing, To, Tag);
}

//receive function for single input funcitons
template<class T1, class T2>
void send(int FromProc, int Tag, T1 &thing, T2 &thing2)
{
	MPIworld.put(thing, To, Tag);
	MPIworld.put(thing2, To, Tag);
}


//send function for Void funcitons
void sendAll(int Tag)
{	}

//receive function for single input funcitons
template<class T1>
void sendAll(int Tag, T1 &thing)
{
	for(int i=0;i<MPIworld.size();++i){
		if(i!=MPIworld.rank()) MPIworld.send(thing, i, Tag);
	}
}

//receive function for double input funcitons
template<class T1, class T2>
void sendAll(int Tag, T1 &thing, T2 &thing2)
{
	for(int i=0;i<MPIworld.size();++i){
		if(i!=MPIworld.rank()) {
			MPIworld.send(thing, i, Tag);
			MPIworld.send(thing2, i, Tag+1);
		}
	}
}

void moo(int kk){
	cout<<endl<<"I was called on: "<<MPIworld.rank()<<" with value: "<<kk<<endl;
	sleep(MPIworld.rank()-1);
}

int main(int argc,char* argv[])
{

	MPIworld.start(argc, argv);
	std::cout<<MPIworld.name()<<"::"<<MPIworld.rank()<<"/"<<MPIworld.size()<<std::endl<<endl;

//ALL procs must register the func before we can use it
	int WorkTag=2, EndTag=1, done =-1, cur=0;
	if(MPIworld.rank()==0){
		Vector<int> V=Range(1,10);
		int ct=0, rr=-1;

		for(int qq=1;qq<MPIworld.size();++qq){
			MPIworld.put(ct, qq, WorkTag); ++ct;
			if(ct>V.size()) break;
		}
		int get;
		while(ct<V.size()){
			get=MPIworld.getAny(rr);
			MPIworld.put(ct,get , WorkTag);
			++ct;
		}
		//get the last returns
		//for(int qq=1;qq<MPIworld.size();++qq)
		//	MPIworld.get(rr, MPI_ANY_SOURCE, MPI_ANY_TAG);

		//put the final kills
		for(int qq=1;qq<MPIworld.size();++qq)
			MPIworld.put(done, qq, EndTag);


	}else{
		while(1){
			MPIworld.get(cur,0, MPI_ANY_TAG);
			if(cur==-1) break;
			moo(cur);
			MPIworld.put(cur,0, WorkTag);

		}
	}

	MPIworld.end();
}
