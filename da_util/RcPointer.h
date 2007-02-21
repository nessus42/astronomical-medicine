// -*- Mode: C++; fill-column: 73; fill-prefix: "//# " -*-
#ifndef RcPointer_H
#define RcPointer_H

#define RcPointerRcsId_H \
"$Id: RcPointer.h,v 5.4 1997/05/08 16:42:39 nessus Exp $"

// Description : An implementation of reference-counted pointers.
// Author :      Douglas Alan <nessus@mit.edu>.  Plagiarized from an article
//               published in the March-April issue of *C++ Report*, "Memory
//               Management and Smart Pointers".
// Copyright   : Massachusetts Institute of Technology, 1993

#include <EacUtil.h>
#include <da_usual.h>

template <class T> class RcPointer;
template <class T> class ConstRcPointer;
template <class T> class RciPointer;
template <class T> class ConstRciPointer;
class RciBasePointer;
class ConstRciBasePointer;

//#############################################################################
//#####                                                                   #####
//#####            RcPointer<T>: template handle class                    #####
//#####                                                                   #####
//#############################################################################

//# An 'RcPointer' is a reference-counted pointer template class.  It is used
//# like this:

//#      RcPointer<Foo> pf = new Foo(bar);
//#      RcPointer<Foo> pf2;
//#      pf2 = pf;
//#      pf2->operation();

//# In the above code there are no restrictions on the data type Foo -- it can
//# be any class or any built-in data type.  'pf' will act almost exactly like
//# a Foo*, except that a reference count to *pf is maintained.  When 'pf' is
//# copied (or a copy of 'pf' is copied), the reference count is incremented.
//# When the destructor of 'pf' (or a copy of 'pf') is called, the reference
//# count is decremented.  When the count goes to 0, *pf will be deleted.

//# There are a couple of caveats: (1) For each object that needs a reference
//# count, a small structure to maintain the reference count is allocated on
//# the heap.  (2) RcPointer does not understand inheritance; i.e. If Bar is
//# derived from Foo, you cannot assign an RcPointer<Bar> to a variable
//# declared as an RcPointer<Foo>.

//# At some point in the near future, I will define another template class that
//# fixes both of these problems, but it will have its own caveat: the
//# parameter to the template class will need to be a class, and it will need
//# to define three methods, incrementRefCount(), decrementRefCount(), and
//# refCount(), along with an appropriate private data member to keep track of
//# the reference count.  I know of no way to fix caveat (2) without making
//# this requirement.  If you know of a way, please tell me.


//-----------------------------------------------------------------------------
// RcPointerStore: local template struct
//-----------------------------------------------------------------------------

// This struct 'RcPointerStore' is for use only by the implementation of
// RcPointer:

template <class T>
struct RcPointerStore {
private:
  friend class RcPointer<T>;
  friend class ConstRcPointer<T>;

  unsigned rc;
  T*       p;

  RcPointerStore(T* pt): rc(1), p(pt) { assert(pt); }
  ~RcPointerStore() { delete p; }
};


//*****************************************************************************
//*****                                                                   *****
//*****              RcPointer<T>: template handle class                  *****
//*****                                                                   *****
//*****************************************************************************

template <class T>
class RcPointer {
  friend class ConstRcPointer<T>;

  //* Instance variables:
  RcPointerStore<T>*		_store;

  //* Private methods:
  void 				free() { if (_store && --_store->rc == 0)
					   delete _store;
				       }
public:
  
  //* Default constructor:
  RcPointer(): _store(0) {}

  //* Copy constructor:
  RcPointer(const RcPointer<T>& src): _store(src._store)
    { if (_store) _store->rc++; }

  //* Conversion constructors:
  RcPointer(T* p)
     : _store(p ? daCheckHeap(new RcPointerStore<T>(p)) : 0) {}

  //* Destructor:
  ~RcPointer() { free(); }

  //* Nonvirtual methods:
  const RcPointer<T>&		operator=(const RcPointer<T>& src)
                                  { if (src._store) src._store->rc++;
				    free();
				    _store = src._store;
				    return *this; }

  T* 				operator->() const
                                  { assert(_store); return _store->p; }

  T&				operator*() const
                                  { assert(_store); return *_store->p; }

  operator 			const void*() const
                                  { return _store ? _store->p : 0; }
};


//*****************************************************************************
//*****                                                                   *****
//*****      ConstRcPointer<T>: template const handle class               *****
//*****                                                                   *****
//*****************************************************************************

template <class T>
class ConstRcPointer {
  friend class RcPointer<T>;

  //* Instance variables:
  RcPointerStore<T>*		_store;

  //* Private methods:
  void 				free() { if (_store && --_store->rc == 0)
					   delete _store; }

public:
  
  //* Default constructor:
  ConstRcPointer(): _store(0) {}

  //* Copy constructor:
  ConstRcPointer(const ConstRcPointer<T>& src): _store(src._store)
    { if (_store) _store->rc++; }

  //* Conversion constructors:
  ConstRcPointer(const RcPointer<T>& src): _store(src._store)
    { if (_store) _store->rc++; }

  ConstRcPointer(const T* p)
    : _store(p ? daCheckHeap(new RcPointerStore<T>((T*) p)): 0) {}

  //* Destructor:
  ~ConstRcPointer() { free(); }

  //* Nonvirtual methods:
  const ConstRcPointer<T>&      operator=(const ConstRcPointer<T>& src)
                                  { if (src._store) src._store->rc++;
				    free();
				    _store = src._store;
				    return *this; }

  const ConstRcPointer<T>&      operator=(const RcPointer<T>& src)
                                   { if (src._store) src._store->rc++;
				     free();
				     _store = src._store;
				     return *this; }

  const ConstRcPointer<T>&     operator=(const T* p)
                                   { free();
				     if (p) {
				       _store = new RcPointerStore<T>((T*) p);
				       daCheckHeap(_store);
				     } else _store = 0;
				     return *this;
				   }

  T* 				operator->() const
                                  { assert(_store); return _store->p; }

  T&				operator*() const
                                  { assert(_store); return *_store->p; }

  operator 			const void*() const
                                  { return _store ? _store->p : 0; }
};


//*****************************************************************************
//*****                                                                   *****
//*****                RciObject: base class                              *****
//*****                                                                   *****
//*****************************************************************************

class RciObject {
  friend class RciBasePointer;
  friend class ConstRciBasePointer;

  //* Instance variables:
  int				_refCount;
  
  //* Private methods:
  void				incrementRefCount() const
                                   { ++((RciObject*) this)->_refCount; }
  int				decrementRefCount() const
                                   { return --((RciObject*) this)->_refCount; }

protected:

  //* Deactive delete:
  void				operator delete(void*, size_t)
                                   { rciObjectDeletionError(); }

public:

  //* Constructors, etc:
  RciObject(): _refCount(0) {}
  virtual ~RciObject() {}

  //* Nonvirtual methods:
  int				refCount() const { return _refCount; }

  
  // We define operator new() here to be the same as ::operator new().  We do
  // this to avoid a warning message for having defined delete() above without
  // having also defined new():
  void*				operator new(size_t size)
                                   { return ::operator new(size); }
  void*				operator new(size_t size, void* p)
                                   { return ::operator new(size, p); }
};


//*****************************************************************************
//*****                                                                   *****
//*****                RciBasePointer: handle class                       *****
//*****                                                                   *****
//*****************************************************************************

class RciBasePointer {
  friend class ConstRciBasePointer;
  
  //* Instance variables:
  RciObject*			_p;

  //* Private methods:
  void 				free(){ if (_p && _p->decrementRefCount() == 0)
					  ::delete _p; }

public:
  
  //* Default constructor:
  RciBasePointer(): _p(0) {}

  //* Copy constructor:
  RciBasePointer(const RciBasePointer& src): _p(src._p)
    { if (_p) _p->incrementRefCount(); }

  //* Converting constructors:
  RciBasePointer(RciObject* p) : _p(p)
    { if (p) p->incrementRefCount(); }

  //* Destructor:
  ~RciBasePointer() { free(); }

  //* Nonvirtual methods:

  const RciBasePointer&		operator=(const RciBasePointer& src)
                                  { if (src._p) src._p->incrementRefCount();
				    free();
				    _p = src._p;
				    return *this; }

  const RciBasePointer&		operator=(RciObject* src)
                                  { if (src) src->incrementRefCount();
				    free();
				    _p = src;
				    return *this; }

  RciObject* 			operator->() const { return _p; }
  RciObject&			operator*()  const { return *_p; }
  operator			RciObject*() const { return _p; }
  RciObject*			rawPointer() const { return _p; }
};


//*****************************************************************************
//*****                                                                   *****
//*****          RciPointer<T>: template handle class                     *****
//*****                                                                   *****
//*****************************************************************************

template <class T>
class RciPointer {
  friend class ConstRciPointer<T>;

  //* Instance variables:
  RciBasePointer			_p;

public:

  //* Constructors, etc:
  RciPointer() {}	                	        // Default ctor
  RciPointer(const RciPointer<T>& src): _p(src._p) {}   // Copy ctor
  RciPointer(T* p): _p((RciObject*)(void*) p) {}       // Converting ctor

  //* Methods:
  const RciPointer<T>&		operator=(const RciPointer<T>& src)
                                   { _p = src._p; return *this;}
  const RciPointer<T>&		operator=(T* src) 
                                   { _p = (RciObject*)(void*) src;
				     return *this; }
  T* 				operator->() const
                                   { return (T*) _p.operator->(); }
  T&				operator*() const { return (T&) *_p; }
  T*				rawPointer() const
                                   { return (T*) _p.rawPointer(); }
  operator			T*() const { return rawPointer(); }
};


//*****************************************************************************
//*****                                                                   *****
//*****      ConstRciPointer<T>: template const handle class              *****
//*****                                                                   *****
//*****************************************************************************

template <class T>
class ConstRciPointer {

  //* Instance variables:
  RciBasePointer	        	_p;

public:
  
  //* Constructors, etc:
  ConstRciPointer() {};
  ConstRciPointer(const ConstRciPointer<T>& src): _p(src._p) {}
  ConstRciPointer(const RciPointer<T>& src): _p(src._p) {}
  ConstRciPointer(const T* p): _p((RciObject*)(void*) p) {}

  //* Methods:
  const ConstRciPointer<T>&	operator=(const ConstRciPointer<T>& src)
                                   { _p = src._p; return *this; }
  const ConstRciPointer<T>&	operator=(const RciPointer<T>& src)
                                   { _p = src._p; return *this; }
  const ConstRciPointer<T>&     operator=(const T* src)
                                   { _p = (RciObject*)(void*)src;
				     return *this; }
  const T* 		        operator->() const
                                   { return (T*) _p.operator->(); }
  const T&		        operator*() const { return (T&) *_p; }
  const T*			rawPointer() const
                                   { return (T*) _p.rawPointer(); }
  operator		        const T*() const { return rawPointer(); }
};

#endif // RcPointer_h
