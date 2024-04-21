// arrays.h: various templates for safe and convenient arrays

#if !defined(AFX_ARRAYS_H__D99CA647_F19D_463F_9511_12D9FE4974F2__INCLUDED_)
#define AFX_ARRAYS_H__D99CA647_F19D_463F_9511_12D9FE4974F2__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<stdexcept>
#include<new>

namespace arrays
{
	
	template<class Type>
		class array
	{
	protected:
		Type* pointer;
		int n;
	public:
		array<Type>() {n=0; pointer=0;};
		array<Type> (int n) {this->n=0; makeSize(n);};
		void makeSize (int n);
		const int getSize() const
			{return n;};
		Type& operator[] (int index);
		const Type& operator[] (int index) const;
		virtual ~array<Type>(){nullify();};
	protected:
		void nullify();
	private:
		array<Type> (const array<Type>&);
		array<Type>& operator= (const array<Type>&);
	}; // template array<Type>
	
	
	template<class Type>
		class array2D
	{
	protected:
		array< array<Type> > outerArray;
		//	int nOuter, nInner;
	public:
		array2D<Type>();
		array2D<Type>(int n, int m);
		virtual ~array2D<Type>(){}; // everything is taken care of by array!!!
		array<Type>& operator[] (int n) {return outerArray[n];};
		const array<Type>& operator[] (int n) const {return outerArray[n];};
		void makeSize(int n, int m);
	}; // template array2D<Type>
	
	
	template<class Type>
		class array3D
	{
	protected:
		array< array2D<Type> > outerArray;
		//	int n1, n2, n3;
	public:
		array3D<Type>();
		array3D<Type>(int n1, int n2, int n3);
		virtual ~array3D<Type>(){};
		array2D<Type>& operator[] (int n) {return outerArray[n];};
		const array2D<Type>& operator[] (int n) const {return outerArray[n];};
		void makeSize(int n1, int n2, int n3);		
	}; // template array3D<Type>
	
	
	
	template<class Type>
		class array4D
	{
	protected:
		array< array3D<Type> > outerArray;
		//	int n1, n2, n3, n4;
	public:
		array4D<Type>();
		array4D<Type>(int n1, int n2, int n3, int n4);
		virtual ~array4D<Type>(){};
		array3D<Type>& operator[] (int n){return outerArray[n];};
		const array3D<Type>& operator[] (int n) const {return outerArray[n];};
		void makeSize(int n1, int n2, int n3, int n4);
	}; // template array4D<Type>
	


// useful for arrays of non-trivial classes
// when usage of default-ctor is undesired.
	template<class Pointer>
		class pointerArray : public array<Pointer>
	{
	public:
		pointerArray<Pointer>()
			:array<Pointer>()
		{};
		pointerArray<Pointer> (int n)
			: array<Pointer>()
			{makeSize(n);};

		void makeSize(int n)
		{
			array<Pointer>::makeSize(n);
			for (int i=0; i<n; i++) pointer[i]=0;
		};

		virtual ~pointerArray<Pointer>()
			{for (int i=0; i<n; i++) delete pointer[i];};
	}; // template pointerarray<Pointer> : public array<Pointer>



	// IMPLEMENTATIONS!!!!!!!!!!
	
	
	// Implementations of array
	template<class Type>
		void array<Type>::makeSize(int n)
	{
		if (n<=0)
			throw std::logic_error("Array size must be non-negative!");
		if (this->n)
			nullify();
		this->n = n;
		try
		{pointer = new Type[n];}
		catch (...)
		{
			pointer=0;
			throw;
		};
	};
	
	template<class Type>
		Type& array<Type>::operator[] (int index)
	{
		if (index<0)
			throw std::logic_error("Negative index encountered");
		if (index>=n)
			throw std::logic_error("Out of array bounds exception");
		return pointer[index];
	};
	
	
	template<class Type>
		const Type& array<Type>::operator[] (int index) const
	{
		if (index<0)
			throw std::logic_error("Negative index encountered");
		if (index>=n)
			throw std::logic_error("Out of array bounds exception");
		return pointer[index];
	};
	
	
	template<class Type>
		void array<Type>::nullify()
	{
		delete [] pointer;
		pointer=0;
		n=0;
	};
	// \array
	
	
	
	// Implementations of array2D
	template<class Type>
		array2D<Type>::array2D<Type>()
		:
	outerArray()
	{
		//	nOuter=nInner=0;
	};
	
	template<class Type>
		array2D<Type>::array2D<Type>(int n, int m)
		:
	outerArray(n)
	{
		for (int i=0; i<n; i++)
			outerArray[i].makeSize(m);
		//	nOuter=n;
		//	nInner=m;
	};
	
	template<class Type>
		void array2D<Type>::makeSize(int n, int m)
	{
		outerArray.makeSize(n);
		for (int i=0; i<n; i++)
			outerArray[i].makeSize(m);
		//	nOuter=n; 
		//	nInner=m;
	};
	
	// \array2d
	
	
	
	// Implementations of array3D
	template<class Type>
		array3D<Type>::array3D<Type>()
		:
	outerArray()
	{
		//	n1=n2=n3=0;
	};
	
	template<class Type>
		array3D<Type>::array3D<Type>(int n1, int n2, int n3)
		:
	outerArray(n1)
	{
		for (int i=0; i<n1; i++)
			outerArray[i].makeSize(n2, n3);
		//	this->n1=n1;
		//	this->n2=n2;
		//	this->n3=n3;
	};
	
	template<class Type>
		void array3D<Type>::makeSize(int n1, int n2, int n3)
	{
		outerArray.makeSize(n1);
		for (int i=0; i<n1; i++)
			outerArray[i].makeSize(n2, n3);
		//	this->n1=n1;
		//	this->n2=n2;
		//	this->n3=n3;
	};
	// \array3D
	
	
	// Implementations of array4D
	template<class Type>
		array4D<Type>::array4D<Type>()
		:
	outerArray()
	{
		//	n1=n2=n3=n4=0;
	}; // array4D::default-ctor
	
	template<class Type>
		array4D<Type>::array4D<Type>(int n1, int n2, int n3, int n4)
		:
	outerArray(n1)
	{
		for (int i=0; i<n1; i++)
			outerArray[i].makeSize(n2, n3, n4);
		//	this->n1=n1;
		//	this->n2=n2;
		//	this->n3=n3;
		//	this->n4=n4;
	}; // array4D::ctor(int^4)
	
	template<class Type>
		void array4D<Type>::makeSize(int n1, int n2, int n3, int n4)
	{
		outerArray.makeSize(n1);
		for (int i=0; i<n1; i++)
			outerArray[i].makeSize(n2, n3, n4);
		//	this->n1=n1;
		//	this->n2=n2;
		//	this->n3=n3;
		//	this->n4=n4;
	}; // array4D::makeSize
	// \array4D
	
} // namespace arrays


#endif // !defined(AFX_ARRAYS_H__D99CA647_F19D_463F_9511_12D9FE4974F2__INCLUDED_)
