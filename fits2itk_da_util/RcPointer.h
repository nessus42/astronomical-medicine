// -*- Mode: C++; fill-column: 79 -*-
#ifndef RcPointer_H
#define RcPointer_H
//=============================================================================
// Description: An implementation of reference-counted pointers
// Module:      RcPointer.h
// Language:    C++
// Author :     Douglas Alan <doug AT alum.mit.edu>.  Plagiarized from an
//              article published in the March-April issue of *C++ Report*,
//              "Memory Management and Smart Pointers".
//
// Copyright (c) 1993 Douglas Alan
//
// This is free software available under the terms of the "The MIT License".
// See LICENSE.txt for for details
//=============================================================================

template <class T> class RcPointer;
template <class T> class ConstRcPointer;
template <class T> class RcMallocPointer;
template <class T> class ConstRcMallocPointer;
template <class T> class RciPointer;
template <class T> class ConstRciPointer;


class RciBasePointer;

void rciObjectDeletionError();

//#############################################################################
//#####                                                                   #####
//#####            RcPointer<T>                                           #####
//#####            ConstRcPointer<T>                                      #####
//#####            RciPointer<T>                                          #####
//#####            ConstRciPointer<T>                                     #####
//#####            RciObject                                              #####
//#####                                                                   #####
//#####            RcMallocPointer<T>                                     #####
//#####            ConstMallocRcPointer<T>                                #####
//#####                                                                   #####
//#####               Classes for doing automated memory management       #####
//#####               via reference counting                              #####
//#####                                                                   #####
//#####                                                                   #####
//#############################################################################

//! @class RcPointer<T>
//! 
//! @brief A reference-counted pointer template class.

//! It is used //! like this:

//!      RcPointer<Foo> pf = new Foo(bar);
//!      RcPointer<Foo> pf2;
//!      pf2 = pf;
//!      pf2->operation();

//! In the above code there are no restrictions on the data type Foo -- it can
//! be any class or any built-in data type.  'pf' will act almost exactly like
//! a Foo*, except that a reference count to *pf is maintained.  When 'pf' is
//! copied (or a copy of 'pf' is copied), the reference count is incremented.
//! When the destructor of 'pf' (or a copy of 'pf') is called, the reference
//! count is decremented.  When the count goes to 0, *pf will be deleted.

//! There are a couple of caveats: (1) For each object that needs a reference
//! count, a small structure to maintain the reference count is allocated on
//! the heap.  (2) RcPointer does not understand inheritance; i.e. If Bar is
//! derived from Foo, you cannot assign an RcPointer<Bar> to a variable
//! declared as an RcPointer<Foo>.

//! We have also defined another reference counted pointer class,
//! RciPointer<T>, that fixes both of these problems, but it has its own
//! caveat: For an object to use an RciPointer to manage its lifetime, its
//! class must be derived from RciBase.  I know of no way to fix caveat #2
//! above without making this requirement, or a similar one.  (If you know of a
//! better way around caveat #2, please tell me.)

//! Note: There are const versions of the classes too, and versions that work
//! on malloc'ed data, rather than new'ed data.

//! TODO: RcMallocPointer<T> and ConstMallocRcPointer<T> should be eliminated
//! via the clever use of template specialization.  As it is, they are (with
//! the exception of the one line of code that does a free() vs a delete)
//! completely cut-and-paste versions of RcPointer<T> and ConstRcPointer<T>.

//----------

//! @class ConstRcPointer<T>
//! 
//! @brief See the documentation for RcPointer<T>.

//----------

//! @class RcMallocPointer<T>
//! 
//! @brief See the documentation for RcPointer<T>.

//----------

//! @class ConstRcMallocPointer<T>
//! 
//! @brief See the documentation for RcPointer<T>.

//----------

//! @class RciPointer<T>
//! 
//! @brief See the documentation for RcPointer<T>.

//----------

//! @class ConstRciPointer<T>
//! 
//! @brief See the documentation for RcPointer<T>.

//----------


//-----------------------------------------------------------------------------
//-----                                                                   -----
//-----          RcPointerStore: local template struct                    -----
//-----                                                                   -----
//-----------------------------------------------------------------------------

//! Only to be used by the implementations of RcPointer and ConstRcPointer.

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
     : _store(p ? new RcPointerStore<T>(p) : 0) {}

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

  T*			        rawPointer() const
                                  { return _store ? _store->p : 0; }

  				operator bool() const
                                  { return _store ? _store->p : false; }


  // TODO: Remove this commented out method in this class and in classes below.
  // operator void*() was implemented previously (I think) so that we could use
  // an RcPointer in a boolean context, but that C++ has bools, we can now just
  // use operator bool() instead.

// operator 			void*() const
//                                { return _store ? _store->p : 0; }

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
    : _store(p ? new RcPointerStore<T>((T*) p): 0) {}

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

  const ConstRcPointer<T>&      operator=(const T* p)
                                   { free();
				     if (p) {
				       _store = new RcPointerStore<T>((T*) p);
				     } else _store = 0;
				     return *this;
				   }

  const T*		 	operator->() const
                                  { assert(_store); return _store->p; }

  const T&			operator*() const
                                  { assert(_store); return *_store->p; }

  const T*		        rawPointer() const
                                  { return _store ? _store->p : 0; }

  				operator bool() const
                                  { return _store ? _store->p : false; }

//   operator 			const void*() const
//                                   { return _store ? _store->p : 0; }

};


//-----------------------------------------------------------------------------
//-----                                                                   -----
//-----          RcMallocPointerStore: local template struct              -----
//-----                                                                   -----
//-----------------------------------------------------------------------------

//! Only to be used by the implementations of RcMallocPointer and
//! ConstRcMallocPointer.

template <class T>
struct RcMallocPointerStore {
private:
  friend class RcMallocPointer<T>;
  friend class ConstRcMallocPointer<T>;

  unsigned rc;
  T*       p;

  RcMallocPointerStore(T* pt): rc(1), p(pt) { assert(pt); }
  ~RcMallocPointerStore() { free(p); }
};


//*****************************************************************************
//*****                                                                   *****
//*****              RcMallocPointer<T>: template handle class            *****
//*****                                                                   *****
//*****************************************************************************

template <class T>
class RcMallocPointer {
  friend class ConstRcMallocPointer<T>;

  //* Instance variables:
  RcMallocPointerStore<T>*		_store;

  //* Private methods:
  void 				free() { if (_store && --_store->rc == 0)
					   delete _store;
				       }
public:
  
  //* Default constructor:
  RcMallocPointer(): _store(0) {}

  //* Copy constructor:
  RcMallocPointer(const RcMallocPointer<T>& src): _store(src._store)
    { if (_store) _store->rc++; }

  //* Conversion constructors:
  RcMallocPointer(T* p)
     : _store(p ? new RcMallocPointerStore<T>(p) : 0) {}

  //* Destructor:
  ~RcMallocPointer() { free(); }

  //* Nonvirtual methods:
  const RcMallocPointer<T>&    operator=(const RcMallocPointer<T>& src)
                                  { if (src._store) src._store->rc++;
				    free();
				    _store = src._store;
				    return *this; }

  T* 				operator->() const
                                  { assert(_store); return _store->p; }

  T&				operator*() const
                                  { assert(_store); return *_store->p; }

  T*			        rawPointer() const
                                  { return _store ? _store->p : 0; }

  				operator bool() const
                                  { return _store ? _store->p : false; }

//   operator 			const void*() const
//                                   { return _store ? _store->p : 0; }
};


//*****************************************************************************
//*****                                                                   *****
//*****      ConstRcMallocPointer<T>: template const handle class         *****
//*****                                                                   *****
//*****************************************************************************

template <class T>
class ConstRcMallocPointer {
  friend class RcMallocPointer<T>;

  //* Instance variables:
  RcMallocPointerStore<T>*	_store;

  //* Private methods:
  void 				free() { if (_store && --_store->rc == 0)
					   delete _store; }

public:
  
  //* Default constructor:
  ConstRcMallocPointer(): _store(0) {}

  //* Copy constructor:
  ConstRcMallocPointer(const ConstRcMallocPointer<T>& src): _store(src._store)
    { if (_store) _store->rc++; }

  //* Conversion constructors:
  ConstRcMallocPointer(const RcMallocPointer<T>& src): _store(src._store)
    { if (_store) _store->rc++; }

  ConstRcMallocPointer(const T* p)
    : _store(p ? new RcMallocPointerStore<T>((T*) p): 0) {}

  //* Destructor:
  ~ConstRcMallocPointer() { free(); }

  //* Nonvirtual methods:
  const ConstRcMallocPointer<T>&  operator=(const ConstRcMallocPointer<T>& src)
                                  { if (src._store) src._store->rc++;
				    free();
				    _store = src._store;
				    return *this; }

  const ConstRcMallocPointer<T>&   operator=(const RcMallocPointer<T>& src)
                                   { if (src._store) src._store->rc++;
				     free();
				     _store = src._store;
				     return *this; }

  const ConstRcMallocPointer<T>&   operator=(const T* p)
                                   { free();
				     if (p) {
				       _store =
					 new RcMallocPointerStore<T>((T*) p);
				     } else _store = 0;
				     return *this;
				   }

  const T*			operator->() const
                                  { assert(_store); return _store->p; }

  const T&			operator*() const
                                  { assert(_store); return *_store->p; }

  const T*		        rawPointer() const
                                  { return _store ? _store->p : 0; }

  				operator bool() const
                                  { return _store ? _store->p : false; }

//   operator 			const void*() const
//                                   { return _store ? _store->p : 0; }

};


//*****************************************************************************
//*****                                                                   *****
//*****                RciObject: base class                              *****
//*****                                                                   *****
//*****************************************************************************

class RciObject {
  friend class RciBasePointer;

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


//-----------------------------------------------------------------------------
//-----                                                                   -----
//-----          RcBasePointer: local handle class                        -----
//-----                                                                   -----
//-----------------------------------------------------------------------------

//! Only to be used by the implementations of RciPointer and ConstRciPointer.

class RciBasePointer {
  
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
  RciPointer(T* p): _p((RciObject*)(void*) p) {}        // Converting ctor

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
