

#ifndef _listblochparsGrad_h_
#define _listblochparsGrad_h_

#include "bloch/blochParams.h"
#include "bloch/listblochpars.h"
#include "bloch/gradientgrid.h"

BEGIN_BL_NAMESPACE


/* maintains a list of bloch parameters with a Gradient! */

template<class Engine_t>
class ListBlochParamsGrad : public ListBlochParams {
	private:

		GradientGrid<Engine_t> *grads_;
		const void IterErr();
		static double off_;
		bool apply_;
	public:

		ListBlochParamsGrad(): ListBlochParams(),
			grads_(NULL)
		{}

		ListBlochParamsGrad(int nsp):
			ListBlochParams(nsp),
			grads_(NULL)
		{}

		ListBlochParamsGrad(int nsp, string inspin):
			ListBlochParams(nsp, inspin),
			grads_(NULL)
		{}

		ListBlochParamsGrad(int nsp, string inspin, GradientGrid<Engine_t> &gr):
			ListBlochParams(nsp, inspin)
		{
			grads_=&gr;
		}

		ListBlochParamsGrad(const BlochParams &one):
			ListBlochParams(one),
			grads_(NULL)
		{}

		ListBlochParamsGrad(int nsp, const BlochParams &dup):
			ListBlochParams(nsp,dup),
			grads_(NULL)
		{}

		ListBlochParamsGrad( const ListBlochParams &dup):
			ListBlochParams(dup),
			grads_(NULL)
		{}

		ListBlochParamsGrad( const ListBlochParams &dup,GradientGrid<Engine_t> &gr):
			ListBlochParams(dup),
			grads_(&gr)
		{}

		~ListBlochParamsGrad(){ grads_=NULL;	}

		void SetGrads(GradientGrid<Engine_t> &in);

		double &offset();
		inline void off(){	apply_=false;	}
		inline void on(){		apply_=true;	}
		inline void GradOff(){	apply_=false;	}
		inline void GradOn(){		apply_=true;	}

		void operator++();
		void reset();
		operator bool(){ return (ListBlochParams::operator bool())&&(*grads_); }

		ListBlochParamsGrad &operator=(const ListBlochParamsGrad &rhs);

//		void push_back(const BlochParams &rhs);

/*		friend std::ostream& operator<<(std::ostream &oo,ListBlochParams &out);
		void print(std::ostream &oo);
		void print(){ print(std::cout);	}	//text print
		bool write(std::fstream &oo);	//binary write
		bool read(std::fstream &in);	//binart read
*/
};

template<class Engine_t>
double ListBlochParamsGrad<Engine_t>::off_=0;

END_BL_NAMESPACE


#include "listblochparsGrad_meth.h"

#endif

