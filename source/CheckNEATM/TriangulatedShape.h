// TriangulatedShape inherits shape
// Still an abstract class

// Author: michael.mueller@dlr.de
// Date:   Nov. 19 2003

// Version 2:
// Date:   2004 May 11

// body file: TriangulatedShape.cpp
// body file: FacetTypeBasic.cpp


// Interface with subclasses:
// Arrays for vertices and facets are allocated and initialized with null.
// Subclasses can fill them in.
// Their deletion is taken care of in here!

// Arrays are 1-based!



// Version 2:
// Arrays are replaced by std::vectors, for which memory is reserved here
// Keep them 1-based!


// Convention: Any dangling pointer is set to null!!!!!!!


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// Another convention:
// The vertex vectors are NOT destructed in here,
// rather they are only taken as pointers.
// The *file-classes (which deliver them, anyway) take care of this.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#if !defined(AFX_TRIANGULATEDSHAPE_H__EAFBC63C_E375_42F0_A478_221C2E1E9787__INCLUDED_)
#define AFX_TRIANGULATEDSHAPE_H__EAFBC63C_E375_42F0_A478_221C2E1E9787__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<new> 
#include<vector>
#include"shape.h"



class TriangulatedShape : public shape
{
public:
	class FacetTypeBasic; // to be declared just below
protected:
	const int nVertices, nFacets;
	// 1-based arrays of *vertices...
	// allocated + destructed in here, so subclasses needn't take care of this
	//vector3**        vertices;
	//FacetTypeBasic** facets;

	std::vector<vector3*>        vertices;
	std::vector<FacetTypeBasic*> facets;

	// intrinsicDiameter must be set (nonzero) in ctor of subclasses!!!
	// use void TriangulatedShape::calculateIntrinsicDiameter()! (sets scaleFactor, too)
	double intrinsicDiameter;
	double scaleFactor; // = shape::diameter/intrinsicDiameter
public:
	TriangulatedShape(int nVertices, int nFacets,
		double pV, 
		const SpinState& spinAxis,
		double H, double G=0.15)
		throw (std::bad_alloc, std::invalid_argument);
	virtual inline void setPV(double newAlbedo) throw(std::invalid_argument);
	virtual ~TriangulatedShape() throw();
protected:
//	void makeFacetsLookOut();
	void barycenter2zero();
	vector3 barycenter();
	void calculateIntrinsicDiameter() throw (std::invalid_argument, std::logic_error);
	void setScaleFactor()
		{scaleFactor=diameter/intrinsicDiameter;};
public:
	class FacetTypeBasic
	{
	protected:
		// one-based, zeroth elements set to 0
		int vertexID[4];
		// outbound surface element
		const vector3* const dA;	
	public:
		// Calculates dA from the vertices saved in vertexlist (assumed to be one-based!).
		// The shape model must be star-shaped w.r.t. the coordinate system's origin!
		FacetTypeBasic(int vertex1, int vertex2, int vertex3,
			const std::vector<vector3*>& vertexlist)
			throw (std::invalid_argument);
		// dA will be deleted in dtor -- so don't do this by yourself!!!
		FacetTypeBasic(int i1, int i2, int i3, const vector3* const dA)
			throw(std::invalid_argument);

		// passes the array of vertices --> vertex numbers mustn't be changed
		const int* const getVertices() const
			{return vertexID;};

		const vector3& getDA() const throw()
			{return *dA;};

		virtual ~FacetTypeBasic(){delete dA;};
	private:
		FacetTypeBasic();
		FacetTypeBasic(const FacetTypeBasic&);
		FacetTypeBasic& operator=(const FacetTypeBasic&);
	}; //class TriangulatedShape::FacetTypeBasic
}; // class TriangulatedShape



inline void TriangulatedShape::setPV(double newAlbedo) throw (std::invalid_argument)
{
	shape::setPV(newAlbedo);
	// relies on intrinsic Diameter being nonzero!!!!
	setScaleFactor();
};





#endif // !defined(AFX_TRIANGULATEDSHAPE_H__EAFBC63C_E375_42F0_A478_221C2E1E9787__INCLUDED_)
