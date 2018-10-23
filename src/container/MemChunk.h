/*
 * Copyright (c)2000-2001 Bo Blanton::UC Berkeley Dept of Chemistry
 * author:: Bo Blanton
 * contact:: magneto@dirac.cchem.berkeley.edu
 * last modified:: 06-25-01
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */



/*
	MemChunk.h-->a chunk of memmory allocator, destructor, copier

templated to handle many a type...

based on blitz++' (http:://oonumerics.com/blitz/)
*/


#ifndef _MemChunk_h_
#define _MemChunk_h_ 1


#include <iostream>
#include <new>
#include "utils/blassert.h"

BEGIN_BL_NAMESPACE


//#define ToDeBug 1

enum preMemoryPolicy {
  duplicateData,
  deleteDataWhenDone,
  neverDeleteData
};



template<class In_type>
class MemChunkRef;

template<class In_type>
class MemChunk {
	friend class MemChunkRef<In_type>;

	public:
		typedef In_type theType;

	private:   // Data members
	    theType *  data_;
    	theType *  dataAddress_;
    	int     references_;
		int  length_;

    	//here we essentially kill out any 'copy' type constructors
    	//...they mess with the memory alloc so they do nothing here
		MemChunk(const MemChunk<theType>&){ }
		void operator=(const MemChunk<theType>&){ }
		void operator=( MemChunk<theType>){ }

	protected:
		MemChunk(){
			length_ = 0;
			data_ = 0;
			dataAddress_ = 0;
			references_ = 0;
		}

		explicit MemChunk(int items){
			length_ = items;
			allocate((length_%2==0)?length_:(length_+1));

#ifdef ToDeBug
std::cout << "MemChunk: allocated " << setw(8) << length_
     << " at " << ((void *)dataAddress_) << std::endl;

if(dataAddress_==NULL){
	BLEXCEPTION("Assertion failed:")
}
#endif

			references_ = 0;
		}

		MemChunk(int length, theType* data){
			length_ = (length_%2==0)?length_:(length_+1);//length;
			data_ = data;
			dataAddress_ = 0;
			references_ = 0;
		}

		~MemChunk(){
			if (dataAddress_){
#ifdef ToDeBug
std::cout << "MemChunk:     freed " << setw(8) << length_
	 << " at " << ((void *)dataAddress_) << std::endl;
#endif
				deallocate();
			}
		}

		void addReference(){
			++references_;

#ifdef ToDeBug
std::cout << "MemChunk:    reffed " << setw(8) << length_
     << " at " << ((void *)dataAddress_) << " (r="
     << (int)references_ << ")" << std::endl;
#endif

		}
		theType* end()const{
			return &data_[length-2];//&data_[length_-1];

		}

		theType* begin()const{
			return &data_[0];
		}


		theType*  data(){	return data_;	}

		const theType*  data() const{	return data_;	}

		int length() const{		return length_-1;		}

		void removeReference(){
			--references_;

#ifdef ToDeBug
std::cout << "MemChunk: dereffed  " << setw(8) << length_
	 << " at " << ((void *)dataAddress_) << " (r=" << (int)references_
	 << ")" << std::endl;
#endif

		}

		int references() const{	return references_;	}
		int refs() const{ return references_;		}
		inline void allocate(int length){
			data_ =  new theType[length];
    		dataAddress_ = data_;
		}
		void deallocate(){		delete [] dataAddress_;		}

};

//a class to hold 'dangling' bits of memory

template<class In_type>
class UnownedMemChunk : public MemChunk<In_type> {
public:
    UnownedMemChunk(int length, In_type*  data): MemChunk<In_type>(length,data)   {}

    virtual ~UnownedMemChunk()  {}
};

//a class to hold 'Null' bits of memory (i.e. no memory)

template<class In_type>
class NullMemChunk : public MemChunk<In_type> {
public:
    NullMemChunk(){
        // This ensures that the delete operator will not be invoked
        // on an instance of NullMemChunk in removeReference().
        this->addReference();
    }

    virtual ~NullMemChunk(){}
};


template<class In_type>
class MemChunkRef {
	public:
		typedef In_type theType;
	private:
	    MemChunk<theType>* block_;
    	static NullMemChunk<theType> nullBlock_;
    	void operator=(const MemChunkRef<theType>&){ }  //keep assignments out of this class
    	void operator=(MemChunkRef<theType>&){ }  //keep assignments out of this class


	public:

		MemChunkRef(){
			block_ = &nullBlock_;
			block_->addReference();
			data_ = 0;
		}

		MemChunkRef(MemChunkRef<theType>& ref){
			block_ = ref.block_;
			block_->addReference();
			data_ = block_->data();
		}


		MemChunkRef(const MemChunkRef<theType>& ref){
			block_ = ref.block_;
			block_->addReference();
			data_ = block_->data();
		}

		MemChunkRef(MemChunkRef<theType>& ref, int offset){
			block_ = ref.block_;
			block_->addReference();
			data_ = block_->data() + offset;
		}

		MemChunkRef(int length, theType* data,	preMemoryPolicy deletionPolicy=duplicateData){
			// Create a memory block using already allocated memory.

			// Note: if the deletionPolicy is duplicateData, this must
			// be handled by the leaf class.  In MemChunkRef,
			// this is treated as neverDeleteData; the leaf class (e.g. Array)
			// must duplicate the data.

			if ((deletionPolicy == neverDeleteData)	  || (deletionPolicy == duplicateData))
				block_ = new UnownedMemChunk<theType>(length, data);
			else if (deletionPolicy == deleteDataWhenDone)
				block_ = new MemChunk<theType>(length, data);
			block_->addReference();

#ifdef ToDeBug
std::cout << "MemChunkRef: created MemChunk at "
	 << ((void*)block_) << std::endl;
#endif

			data_ = data;
		}

		explicit MemChunkRef(int items){
			block_ = new MemChunk<theType>(items);
			block_->addReference();
			data_ = block_->data();

#ifdef ToDeBug
std::cout << "MemChunkRef: created MemChunk at "
	 << ((void*)block_) << std::endl;
#endif

		}

		void createChunk(int length, theType* data,	preMemoryPolicy deletionPolicy=duplicateData)
		{
			//remove older olds refs
			blockRemoveReference();
			// Create a memory block using already allocated memory.

			// Note: if the deletionPolicy is duplicateData, this must
			// be handled by the leaf class.  In MemChunkRef,
			// this is treated as neverDeleteData; the leaf class (e.g. Array)
			// must duplicate the data.

			if ((deletionPolicy == neverDeleteData)	  || (deletionPolicy == duplicateData))
				block_ = new UnownedMemChunk<theType>(length, data);
			else if (deletionPolicy == deleteDataWhenDone)
				block_ = new MemChunk<theType>(length, data);

			block_->addReference();

#ifdef ToDeBug
std::cout << "MemChunkRef: created MemChunk at "
	 << ((void*)block_) << std::endl;
#endif
			data_ = data;
		}


		void blockRemoveReference()
		{
			block_->removeReference();
			if ((block_->refs() == 0) && (block_ != &nullBlock_)){
#ifdef ToDeBug
std::cout << "MemChunk: no more refs, delete MemChunk object at "
	 << ((void*)block_) << std::endl;
#endif

				delete block_;
			}
		}

	   ~MemChunkRef(){		blockRemoveReference();		}

		int numReferences() 	const{	return block_->references();	}
		int numRefs() 			const{	return block_->references();	}
		int References() 		const{	return block_->references();	}
		int refs() 				const{	return block_->references();	}

		theType* end()const{
			return block_->end();
		}

		theType* begin()const{
			return block_->begin();
		}
		//void push_back(const theType &in){
		//	block_->push_back(in);
		//}

	protected:

		//our data chunk
		theType *  data_;

		void changeToNullBlock(){		//change to a Null block
			blockRemoveReference();
			block_ = &nullBlock_;
			block_->addReference();
			data_ = 0;
		}

		void changeBlock(MemChunkRef<theType>& ref, int offset){
			blockRemoveReference();
			block_ = ref.block_;
			block_->addReference();
			data_ = block_->data() + offset;
		}

		void newBlock(int items){
			blockRemoveReference();
			block_ = new MemChunk<theType>(items);
			block_->addReference();
			data_ = block_->data();

#ifdef ToDeBug
		std::cout << "MemChunkRef: created MemChunk at "
			 << ((void*)block_) << std::endl;
#endif
		}

};


///some 'outside' methods that rely on all the above classes at being defined

// Call up an instance of the 'Null' (that static member inside 'MemChunkRef'..such that here
// we have the static member declared..yippie....
template<class In_type>  NullMemChunk<In_type> MemChunkRef<In_type>::nullBlock_;


END_BL_NAMESPACE



#endif



