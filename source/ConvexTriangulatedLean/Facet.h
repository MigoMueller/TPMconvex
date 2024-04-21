// little utility class for shape models
// Author: michael.mueller@dlr.de
// Date:   Nov. 19 2003

// Not supposed to be a base-class for anything, non-virtual dtor!!!

#if !defined(AFX_FACET_H__04F1946E_BD82_4FB5_B604_CCF29D1974AA__INCLUDED_)
#define AFX_FACET_H__04F1946E_BD82_4FB5_B604_CCF29D1974AA__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<stdexcept>

class Facet
{
protected:
	int a [3];
public:
	~Facet(){};
	
	// throws range_error when an index is <1
	Facet(int i1, int i2, int i3)
		throw(std::range_error)
	{
		if (i1<1 || i2<1 || i3<1)
			throw (std::range_error("Index < 1 encountered"));
		a[0]=i1; a[1]=i2; a[2]=i3;
	};

	Facet(int rhs[3])
	{
		a[0]=rhs[0];
		a[1]=rhs[1];
		a[2]=rhs[2];
	};
		
	// index zero-based; range_error may be thrown
	int operator[]  (int index) const
		throw(std::range_error)
	{
		if (index<0 || index >= 3)
			throw (std::range_error("Facet::operator[]: Bad index encountered"));
		return a[index];
	}; 

	int get0() const
		{return a[0];};
	int get1() const
		{return a[1];};
	int get2() const
		{return a[2];};

private:
	Facet(){};
//	Facet(const Facet& rhs){}; // copy-ctor IS needed!
//	Facet operator=(const Facet& rhs){return *this;};
}; // class Facet





#endif // !defined(AFX_FACET_H__04F1946E_BD82_4FB5_B604_CCF29D1974AA__INCLUDED_)
