// ObjFile.h: interface for the ObjFile class.
//
//////////////////////////////////////////////////////////////////////
// Author: michael.mueller@dlr.de
// Date:   2003 Nov 6
//
// Reads in OBJ-Files (Wavefront-format)
// Purpose: Pass asteroid shape model to thermal IR-routine

// Actually, only a small subset of the actual OBJ-syntax is supported:
// 1. 3D-vertices
// 2. Triangular facets (positive numbering only)
// This comprises all shape models from S. Hudson's web-page.

// Not supposed to be a base-class for anything, non-virtual dtor!!!


// Update 2004 Apr 29:
// Replaced the 'home-grown' ConcatList by STL-vector
// --> Increased stability on run-time 
// --> More concise to read
// --> Interfaces remain all the same

// Update 2004 May 11:
// Replaced the interface to getting vertices and facets.
// Are now 'templatized' such fitting any sort of iterator (C-arrays or STL-containers)
// which support dereferencing and ++.
// Convention fits STL-conventions: begin and end are passed to the functions,
// where begin is the first valid element and end is the first INvalid element.
// (In C-arrays, array[N]: begin=array[0], end=array[N])

// Update 2004 May 12:
// Added public sub-class 'proxy', which serves for storing one ObjFile,
// reading it in the first time it is referenced.
// A conversion operator to const ObjFile& is supplied, so proxy can be used like ObjFile itself.
//
// Background: Useful as substitute for static const ObjFile - accelerates startup of program,
// since ObjFile will only be read in when it's actually needed.

// body file: ObjFile.cpp

#if !defined(AFX_OBJFILE_H__A9876350_F06B_4A25_AFF9_E32E0B2782F8__INCLUDED_)
#define AFX_OBJFILE_H__A9876350_F06B_4A25_AFF9_E32E0B2782F8__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include<iostream>
#include<fstream>
#include<stdexcept>
#include<exception>
#include<new> // for exception bad_alloc
#include<string>
#include<deque>
#include<algorithm> // for std::copy used in getVertexArray
#include "vector3.h"
#include "Facet.h"
#include "TriangulatedShape.h" // for TriangulatedShape::FacetTypeBasic


class ObjFile  
{
public:
	// throws invalid_argument if filename is not a valid OBJ-File (check .what for an explanation)
	ObjFile(const char* filename)
		throw (std::invalid_argument, std::bad_alloc, std::range_error);

	inline int howManyVertices() const throw()
		{return nVertices;};
	inline int howManyFacets() const throw()
		{return nFacets;};
	~ObjFile(); 

	class proxy
	{
	public:
		proxy(const std::string& filename) : filename(filename), obj(0) {};
		inline operator const ObjFile& () const; // defined below declaration of ObjFile 
		~proxy() {delete obj;};
	private:
		const std::string filename;
		mutable ObjFile* obj; // mutable: operator (const ObjFile&) can be const!
		proxy();
		proxy(const proxy&);
		const proxy operator=(const proxy&);
	}; // class ObjFile::proxy
protected:
	std::deque<vector3*> vertices;
	std::deque<Facet*>   facets;
	void CleanUp();
	std::ifstream in;
	inline void skipComments();
	int nVertices, nFacets;
public:
	// Two template functions to pass vertices and facets.
	// Both must be cluttered into the class declaration, else VC++6.0 won't compile (BUG!)

	// begin and end must point to vector3*
	// throws logic_error if container is too short or too long
	// as in STL: begin is the first valid member of container, end the first invalid one
	template<class iterator>
	inline void getVertexArray(const iterator begin, const iterator end) const
		throw(std::logic_error)
	{ // check container for correct size, then copy
		int i=0;
		for (iterator currentVertex = begin; currentVertex!=end; ++currentVertex)
			++i;
		if (i!=nVertices)
			throw std::logic_error("ObjFile::getVertexArray: container too short or too long!");
		std::copy(vertices.begin(), vertices.end(), begin);
	}; // ObjFile::getVertexArray
	


	// begin and end must point to TriangulatedShape::FacetTypeBasic*
	// throws logic_error if container is too short or too long, bad_alloc if out-of-memory
	template<class iteratorF>
	void getFacetArray(const iteratorF begin, const iteratorF end, const std::vector<vector3*>& vertexList) const 
		throw(std::bad_alloc, std::logic_error)
	{ // check container for correct size, then create TriangulatedShape::FacetTypeBasic's
		int i=0;
		iteratorF currentFacet = begin;
		for (; currentFacet!=end; ++currentFacet)
			++i;
		if (i!=nFacets)
			throw std::logic_error("ObjFile::getFacetArray: container too short or too long!");
		try
		{ // catch out-of-memory
			currentFacet=begin;
			for (i=0; i<nFacets; i++, currentFacet++)
			{ 
				int i1= (*facets[i])[0];
				int i2= (*facets[i])[1];
				int i3= (*facets[i])[2];
				try 
				{*currentFacet = new TriangulatedShape::FacetTypeBasic (i1, i2, i3, vertexList);}
				catch (...)
				{throw i;};
			};
		}
		catch(int n)
		{
			currentFacet=begin;
			for (i=0; i<n; i++, currentFacet++)
				delete *currentFacet;
			throw std::bad_alloc(); 
		};
	}; // ObjFile::getFacetArray
private: // prohibitive declarations
	ObjFile(){};
	ObjFile operator=(const ObjFile&){return *this;};
	ObjFile(const ObjFile& rhs){};
}; 
// class ObjFile



inline ObjFile::proxy::operator const ObjFile& () const
{
	if (!obj)
		obj=new ObjFile(filename.c_str());
	return *obj;
};




#endif // !defined(AFX_OBJFILE_H__A9876350_F06B_4A25_AFF9_E32E0B2782F8__INCLUDED_)
