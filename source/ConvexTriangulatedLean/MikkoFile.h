// MikkoFile.h: interface for the MikkoFile class.
// reads in an asteroid's shape model in Mikko Kaasalainen's format
// similar to ObjFile!

// Michael.Mueller@dlr.de
// Dec. 03 2003


// Works fine in the present version!
// Memory management, in particular, works ok:
// All memory is released upon destruction,
// the vectors contained in vector3** vertices, however,
// are NOT destructed as soon as they were passed to a shape model.
// (which, in turn, works only once)

// UPDATE:	
// Memory IS released in dtor, also pass***-functions are adapted to
// STL-containers.

// Convention:
// Vertices are copied pointer-wise, and are destroyed in dtor.
// Calling routine MUST NOT delete them!
// The TriangulatedShape::FacetTypeBasic's, however, must be deleted by caller.

// body file: MikkoFile.cpp

// Not supposed to be a base-class for anything, non-virtual dtor!!!

#if !defined(AFX_MIKKOFILE_H__B5B54BF7_003B_40DA_8543_46F8990D674E__INCLUDED_)
#define AFX_MIKKOFILE_H__B5B54BF7_003B_40DA_8543_46F8990D674E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include<iostream>
#include<fstream>
#include<stdexcept>
#include<new>
#include<algorithm>

class vector3;
class Facet;
class FacetTypeBasic;
class TriangulatedShape;


class MikkoFile  
{
protected:
	vector3** vertices;
	Facet**   facets;
	int       nVertices, nFacets;
public: 
	// throws invalid_argument if filename is not a valid MikkoFile
	// bad_alloc if there is not enough memory
	MikkoFile(const char* filename)
		throw(std::invalid_argument, std::bad_alloc);

	inline int howManyVertices() const throw()
		{return nVertices;};
	inline int howManyFacets() const throw()
		{return nFacets;};

	void writeObjFile (const char* const filename) const;

	~MikkoFile() throw();

	// Below you'll find two templates for passing vertices and facets.
	// Unfortunately, VC++ 6.0 won't let me define the template AFTER it's declaration (BUG!), sorry!

	// begin and end must point to vector3*
	// throws logic_error if container is too short or too long
	template<class iterator>
	void getVertexArray(const iterator begin, const iterator end) const
		throw(std::bad_alloc, std::logic_error)
	{ // check container for correct size, then copy
		int i=0;
		for (iterator currentVertex = begin; currentVertex!=end; ++currentVertex)
			++i;
		if (i!=nVertices)
			throw std::logic_error("MikkoFile::getVertexArray: Container too short or too long!");
		std::copy(vertices, &vertices[nVertices], begin);
	}; // MikkoFile::getVertexArray
	
	
	// begin and end must point to TriangulatedShape::FacetTypeBasic*
	// throws logic_error if container is too short or too long,
	// throws bad_alloc when out of memory
	template<class iteratorF>
	inline void getFacetArray(const iteratorF begin, const iteratorF end) const
		throw (std::bad_alloc, std::logic_error)
	{
		//check container for correct size
		int i=0; 
		iteratorF currentFacet = begin;
		for ( ; currentFacet!=end; ++currentFacet)
			++i;
		if (i!=nFacets)
			throw std::logic_error("MikkoFile::getFacetArray: Container too short or too long!");
		vector3* dA=0;
		try
		{ // be prepared for out of memory
			currentFacet=begin;
			for (i=0; i<nFacets; i++)
			{ 
				int i1= (*facets[i])[0];
				int i2= (*facets[i])[1];
				int i3= (*facets[i])[2];
				// vertex numbers (ccw as seen from outside, according to Mikko)
				
				const vector3& v1= *vertices[i1-1];
				const vector3& v2= *vertices[i2-1];
				const vector3& v3= *vertices[i3-1];
				// note: vertices is a zero-based array, Mikko however counts them 1-based!				
				try 
				{ // outbound surface element
					dA=new vector3(v2-v1, v3-v1); 
					*dA/=2;
					*currentFacet = new TriangulatedShape::FacetTypeBasic (i1, i2, i3, dA);
					currentFacet++;
					dA=0;
				}
				catch (...)
				{ 
					delete dA; 
					throw i;
				};
			};
		}
		catch(int n)
		{// out of memory, delete the facets already allocated
			currentFacet=begin;
			for (i=0; i<n; i++, currentFacet++)
				delete *currentFacet;
			throw std::bad_alloc(); 
		};
	}; // MikkoFile::getFacetArray

private:
	MikkoFile();
	MikkoFile(const MikkoFile&);
	MikkoFile& operator=(const MikkoFile&);
}; // class MikkoFile


#endif // !defined(AFX_MIKKOFILE_H__B5B54BF7_003B_40DA_8543_46F8990D674E__INCLUDED_)
