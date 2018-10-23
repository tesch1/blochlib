

#ifndef _full_binary_tree_h_
#define _full_binary_tree_h_ 1

#include <iostream>

//a simple full binary tree class...for 2 values 'On' and 'Off'
//where we simply wish to specifiy a 'depth' and it fills out the
//pointer list to that 'depth' recursively

//this is meant for setting up a nice iterator to go through
//large lists of permutations where starting at the Root we travers down
//each possible branch of the tree, changing an EXTERNAL class parameter
//based on the current value in this tree

//so for the example of this pC7 testing, since we can possibly extend our sequence to
// upwards of 100 different pC7 cycles, and we can choose from
// an arbitrary combo of pC7 and pC7bar...this tree allows
// for an easy set up of the 2^(99) combos we can have at a depth of 100 (alot really
// and a bit to many to have a normal machine test that many seqeunce)



/****** EXAMPLE USAGE ****
	FullBinaryTree<1, 0> moo(5);
	FullBinaryTree<1, 0>::pathiterator paths(moo);
	while(paths){
		FullBinaryTree<1, 0>::nodeiterator nodes(paths);
		while(nodes){
			cout<<nodes.value()<<" ";
			++nodes;
		}
		cout<<endl;
		++paths;
	}

//this should give you all
//the unique path ways in the tree (32 of them)
1 1 1 1 1 1
1 1 1 1 1 0
1 1 1 1 0 0
1 1 1 1 0 1
1 1 1 0 0 0
1 1 1 0 0 1
1 1 1 0 1 1
1 1 1 0 1 0
1 1 0 0 0 0
1 1 0 0 0 1
1 1 0 0 1 1
1 1 0 0 1 0
1 1 0 1 1 1
1 1 0 1 1 0
1 1 0 1 0 0
1 1 0 1 0 1
1 0 0 0 0 0
1 0 0 0 0 1
1 0 0 0 1 1
1 0 0 0 1 0
1 0 0 1 1 1
1 0 0 1 1 0
1 0 0 1 0 0
1 0 0 1 0 1
1 0 1 1 1 1
1 0 1 1 1 0
1 0 1 1 0 0
1 0 1 1 0 1
1 0 1 0 0 0
1 0 1 0 0 1
1 0 1 0 1 1
1 0 1 0 1 0

**/




template <int On, int Off>
class PathIterator;

template <int On, int Off>
class NodeIterator;

template <int On=1, int Off=0>
class FullBinaryTree
{
public:
    typedef FullBinaryTree<On, Off> node;

	typedef PathIterator<On, Off> iterator;
	typedef NodeIterator<On, Off> nodeiterator;

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::m_data
//  Description:    data stored in current tree node
// ----------------------------------------------------------------
    int m_data;

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::m_parent
//  Description:    pointer to the parent of the node.
// ----------------------------------------------------------------
    node *m_parent;

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::m_left
//  Description:    pointer to the left child of the node.
// ----------------------------------------------------------------
    node *m_left;

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::m_right
//  Description:    pointer to the right child of the node.
// ----------------------------------------------------------------
    node *m_right;

//-----------------------------------------------------------------
//	Name:			FullBinaryTree::l_r_flag
//	Description:	USED FOR ITERATORS...if 0 i choose the left tree
//                  if 1 i choose the right tree
//-----------------------------------------------------------------
	bool l_r_flag;

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::FullBinaryTree
//  Description:    Constructor. Creates a full tree node.
//                  uses the constructor below to assign
//                  'this' as the parent of the rest
//  Arguments:      - depth: max number of levels
//                  - start: The starting value
//  Return Value:   None.
// ----------------------------------------------------------------
	FullBinaryTree( int depth, int start=On )
		: m_data( start ), l_r_flag(true)
	{
		if(depth>=1)
		{
			m_parent = 0;
			m_left =new FullBinaryTree(this, depth-1, start) ;
			m_right = new FullBinaryTree(this, depth-1, (start==On)?Off:On );
		//	cout<<depth<<" "<<m_left<<" "<<m_right<<" "<<m_parent<<" "<<m_data<<endl;

		}else{
			m_parent= 0;
			m_left= 0;
			m_right=0;
		}
	};

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::FullBinaryTree
//  Description:    Constructor. for the recursion creationof the full tree.
//  Arguments:      - data: data to initialise node with
//  Return Value:   None.
// ----------------------------------------------------------------
	FullBinaryTree(FullBinaryTree *p_tree, int depth, int start=On )
		: m_data( start ), l_r_flag(true)
	{
		if(depth>=1)
		{
			m_parent = p_tree;
			m_left =new FullBinaryTree(this, depth-1, start) ;
			m_right = new FullBinaryTree(this, depth-1, (start==On)?Off:On );
		//	cout<<depth<<" "<<m_left<<" "<<m_right<<" "<<m_parent<<" "<<m_data<<endl;

		}else{
			m_parent= p_tree;
			m_left= 0;
			m_right=0;
		}
	};

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::~FullBinaryTree
//  Description:    Destructor. Deletes all child nodes.
//  Arguments:      None.
//  Return Value:   None.
// ----------------------------------------------------------------
	~FullBinaryTree()
	{
		// delete both children. Since the destructor is recursive,
		// each child node we delete automatically deletes all
		// of it's children. Neat, huh?
		if( m_left != 0 )	delete m_left;
		if( m_right != 0 )	delete m_right;
	}


// ----------------------------------------------------------------
//  Name:           FullBinaryTree::isLeft
//  Description:    determines if node is a left subtree. Note that
//                  a result of false does NOT mean that it is a
//                  right child; it may be a root instead.
//  Arguments:      None.
//  Return Value:   True or False.
// ----------------------------------------------------------------
    inline bool isLeft()
    {
        if( isRoot() )
            return false;
        if( m_parent->m_left == this )
            return true;
        return false;
    }


// ----------------------------------------------------------------
//  Name:           FullBinaryTree::isRight
//  Description:    determines if node is a right subtree. Note that
//                  a result of false does NOT mean that it is a
//                  left child; it may be a root instead.
//  Arguments:      None.
//  Return Value:   True or False.
// ----------------------------------------------------------------
    inline bool isRight()
    {
        if( isRoot() )
            return false;
        if( m_parent->m_right == this )
            return true;
        return false;
    }


// ----------------------------------------------------------------
//  Name:           FullBinaryTree::isRoot
//  Description:    determines if node is a root node.
//  Arguments:      None.
//  Return Value:   True or False.
// ----------------------------------------------------------------
    inline bool isRoot()
    {
        if( m_parent == 0 )
            return true;
        return false;
    }

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::hasChild
//  Description:    determines if node has a child
//  Arguments:      None.
//  Return Value:   True or False.
// ----------------------------------------------------------------
    inline bool hasChild()
    {
        if( m_left == 0 && m_right==0)	return false;
        return true;
    }

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::count
//  Description:    recursively counts all children.
//  Arguments:      None.
//  Return Value:   the number of children.
// ----------------------------------------------------------------
    int count()
    {
        int temp = 1;
        if( m_left != 0 )
            temp += m_left->count();
        if( m_right != 0 )
            temp += m_right->count();
        return temp;
    }


// ----------------------------------------------------------------
//  Name:           FullBinaryTree::depth
//  Description:    recursively determines the lowest relative
//                  depth of all its children.
//  Arguments:      None.
//  Return Value:   the lowest depth of the current node.
// ----------------------------------------------------------------
    int depth()
    {
        int left = -1;
        int right = -1;
        if( m_left != 0 )
            left = m_left->depth();
        if( m_right != 0 )
            right = m_right->depth();
        if( left > right )
            return left + 1;
        return right + 1;
    }

// ----------------------------------------------------------------
//  Name:           FullBinaryTree::resetLRflags
//  Description:    recursively resets all the l_r_flags to true
//  Arguments:      None.
//  Return Value:   None
// ----------------------------------------------------------------
	void resetLRflags()
	{
		l_r_flag=true;
		if( m_left != 0 )		m_left->resetLRflags();
		if( m_right != 0 )		 m_right->resetLRflags();
	}
};


/********************************************************/
/* FullBinaryTreeIterator ****/
/********************************************************/
// THIS HAS 2 Distinct Iterator Types..

// a) the 'PathIterator' -> advances through all the possible paths in the FullTree
// b) the 'NodeIterator' -> advances through the nodes of a Path..

/**********************************************************************/
/****  Path Iterator...iterates down paths ***/
/**********************************************************************/

//the path iterator simply sets the proper 'l_r_flag' inside
// the FullTree
template<int On, int Off>
class PathIterator{
	friend class NodeIterator<On, Off>;
	private:
		FullBinaryTree<On,Off> *m_tree;
		FullBinaryTree<On,Off> *m_cur;
		bool notend;

	public:

		typedef NodeIterator<On, Off> iterator;

		PathIterator(FullBinaryTree<On, Off> &tree):
			m_tree(&tree)
		{
			m_tree->resetLRflags();
			if(m_tree->depth()>1) notend=true;
			m_cur=m_tree;
		}

		PathIterator(FullBinaryTree<On, Off> *tree):
			m_tree(tree)
		{
			m_tree->resetLRflags();
			if(m_tree->depth()>1) notend=true;
			m_cur=m_tree;
		}

		~PathIterator()
		{
			m_tree =0;
			m_cur=0;
		}


	//NOTE:: this algorithm is lacking in finesse...
	//it is brut force and ugly...but it works
		void advancePath()
		{
			//need to get to the last point in the last path iteration
			iterator on(m_tree);
			while(on){	++on;	}
			m_cur=on.current(); //the last point on the chain

			//now we go backwards until we fing a node that
			//has not has l_r_flag switched, and switch it
			m_cur=m_cur->m_parent;
			while(m_cur!=NULL){
				if(!m_cur->l_r_flag){
					m_cur=m_cur->m_parent;
				}else {
					m_cur->l_r_flag=false;
					break;
				}
			}
			//if we attempted to go back to the root parent...
			//we must leave
			if(m_cur==NULL){ notend=false;	}
		}

		void operator++(){	advancePath();	}
		void operator++(int){	advancePath();	}

		operator bool(){	return notend;	}

		void reset(){	m_tree->resetLRflags();	 }

		iterator begin(){	return iterator(m_tree);	}
};


/**********************************************************************/
/****  Node Iterator...iterates down a specific path ***/
/**********************************************************************/

template<int On, int Off>
class NodeIterator{
	private:
		FullBinaryTree<On,Off> *m_tree;
		FullBinaryTree<On,Off> *m_cur;
		bool notend;
		int ct;

	public:

		typedef NodeIterator iterator;

		NodeIterator(FullBinaryTree<On, Off> *tree):
			m_tree(tree), notend(false), ct(0)
		{
			if(m_tree !=NULL){ m_cur=m_tree; notend=true;	}
		}

		NodeIterator(PathIterator<On, Off> &tree):
			m_tree(tree.m_tree), notend(false), ct(0)
		{
			if(m_tree !=NULL){ m_cur=m_tree; notend=true;	}
		}

		~NodeIterator()
		{
			m_tree =0;
			m_cur = 0;
		}


		void advancePath(){
			if(m_cur->hasChild()){
				if(m_cur->l_r_flag){
					m_cur=m_cur->m_left;
				}else{
					m_cur=m_cur->m_right;
				}
				++ct;
			}else{
				notend=false;
			}
		}

		int value(){	return m_cur->m_data;	}
		FullBinaryTree<On,Off> *current()const {	return m_cur;	}

		void operator++(){	advancePath();	}
		void operator++(int){	advancePath();	}
		int curpos()const {	return ct;	}

		void reset(){	m_cur=m_tree; notend=(m_cur!=NULL)?true:false;	ct=0; }

		operator bool(){	return notend;	}
};


/*

template <class BidirectionalIterator>
bool next_permutation(BidirectionalIterator first,
                      BidirectionalIterator last) {
    if (first == last) return false;
    BidirectionalIterator i = first;
    ++i;
    if (i == last) return false;
    i = last;
    --i;

    for(;;) {
        BidirectionalIterator ii = i--;
        if (*i < *ii) {
            BidirectionalIterator j = last;
            while (!(*i < *--j));
            iter_swap(i, j);
            reverse(ii, last);
            return true;
        }
        if (i == first) {
            reverse(first, last);
            return false;
        }
    }
}


*/





#endif
