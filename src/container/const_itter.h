

/* const_itter class and methods for array */

#include<MasterHeader.h>

#ifndef _const_itter_h_
#define _const_itter_h_

BEGIN_BL_NAMESPACE

struct endTag { }; //know about the 'end' of the data lists


template<class T, int N>
class constArrayIter {
	private:
		constArrayIter() { }		//kill the copy constructor from abuse

	private:
		NVector<int,N> strides_, extent_, base_;
		T* stack_[N];
		T* last_[N];
		int stride_;
		int maxRank_;

	protected:
		NVector<int,N> pos_;
		T * data_;
		T * first_;

	public:
		constArrayIter(const Array<T,N>& array){
			// Making internal copies of these avoids keeping
			// a pointer to the array and doing indirection.
			strides_ = array.stride();
			extent_ = array.extent();
			first_ = const_cast<T*>(array.dataFirst());
			data_ = first_;
			base_=array.base();

			maxRank_ = order_(0);
			stride_ = strides_(maxRank_);

			for (int i=0; i < N; ++i){
				stack_[i] = data_;
				last_[i] = data_ + array.extent(i)* strides_(i);
			}
			pos_=array.base()
		}

		constArrayIter(const Array<T,N>& array, endTag){
			// The endTag type is provided by the end() method
			// in Array<T,N>, and indicates that an end iterator
			// is to be constructed.

			// Use 0 pointer to mark end of array.
			// This also handles the case of empty arrays, which
			// have their data pointer set to 0.
			data_ = 0;
		}

		T operator*() const{
			if(data_ == 0){
				std::cerr<<std::endl<<"Error: constArrayIter->()"<<std::endl;
				std::cerr<<"Attempted to dereference invalid iterator "<<std::endl;
				std::cerr<< "(empty array or past end of array)"<<std::endl;
				    BLEXCEPTION(__FILE__,__LINE__)
			}
			return *data_;
		}

		const T* operator->() const{
			if(data_ == 0){
				std::cerr<<std::endl<<"Error: constArrayIter->()"<<std::endl;
				std::cerr<<"Attempted to dereference invalid iterator "<<std::endl;
				std::cerr<< "(empty array or past end of array)"<<std::endl;
				    BLEXCEPTION(__FILE__,__LINE__)
			}
			return data_;
		}

		constArrayIter<T,N>& operator++();

		// This operator returns void, which breaks the STL forward
		// iterator requirements.  Unfortunately many people have
		// gotten into the habit of writing iter++ when they really
		// mean ++iter.  iter++ implemented the proper way requires
		// returning a copy of the original state of the iterator,
		// before increment.  This would be very inefficient, since
		// the iterator contains a lot of data.  Hence the void
		// return: it will be efficient even if you write iter++.
		// Maybe this is a bad idea, let me know if this causes
		// you problems.
		void operator++(int)  { ++(*this); }

		const NVector<int,N>& position() const
		{
			if(data_ != 0){
				std::cerr<<std::endl<<"Error: constArrayIter::position()"<<std::endl;
				std::cerr<<"Attempted on an invalid iterator "<<std::endl;
				std::cerr<< "(empty array or past end of array)"<<std::endl;
				    BLEXCEPTION(__FILE__,__LINE__)
			}
			return pos_;
		}

		bool operator==(const constArrayIter<T,N>& x) const{
			return data_ == x.data_;
		}

		bool operator!=(const constArrayIter<T,N>& x) const{
			return data_ != x.data_;
		}
};


//the 'non' const itterator class (basically const itterator, with out the const flags)
template<class T, int N>
class arrayIterator : public constArrayIter<T,N> {
  public:
    arrayIterator(Array<T,N>& x): constArrayIter<T,N>(x) { }

    arrayIterator(Array<T,N>& x, endTag y): constArrayIter<T,N>(x,y){ }

    arrayIterator<T,N>& operator++(){
        constArrayIter<T,N>::operator++();
        return *this;
    }

    T& operator*(){      return *data_;    }

    T* operator->(){        return data_;  }
};


//the beef of itterator the '++'....

template<class T, int N>
constArrayIter<T,N>& constArrayIter<T,N>::operator++()
{
    if(data_ == 0){
		std::cerr<<std::endl<<"Error: constArrayIter++"<<std::endl;
		std::cerr<<"Attempted to iterate past the end of an array."<<std::endl;
		    BLEXCEPTION(__FILE__,__LINE__)
	}

    data_ += stride_;

    if (data_ != last_[0]){
        // We hit this case almost all the time.
        ++pos_[maxRank_];
        return *this;
    }

    // We've hit the end of a row/column/whatever.  Need to
    // increment one of the loops over another dimension.

    int j = 1;
    for (; j < N; ++j){
        int r = order_(j);
        data_ = stack_[j];
        data_ += strides_[r];
        ++pos_(r);

        if (data_ != last_[j])
            break;
    }

    // All done?
    if (j == N){
        // Setting data_ to 0 indicates the end of the array has
        // been reached, and will match the end iterator.
        data_ = 0;
        return *this;
    }

    stack_[j] = data_;

    // Now reset all the last pointers
    for (--j; j >= 0; --j) {
        int r2 = j;
        stack_[j] = data_;
        last_[j] = data_ + extent_[r2] * strides_[r2];
        pos_[r2] = base_[r2];
    }

    return *this;
}


END_BL_NAMESPACE

#endif
