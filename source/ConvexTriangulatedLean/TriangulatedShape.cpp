#include "TriangulatedShape.h"
using namespace std;



TriangulatedShape::TriangulatedShape(int nVertices, int nFacets,
					 double pV, 
					 const SpinState& spinAxis,
					 double H, double G)// =0.15 by default
			throw (bad_alloc, invalid_argument)
	: shape(pV, spinAxis, H, G),
	nVertices(nVertices), nFacets(nFacets)
{
	intrinsicDiameter=0;
	scaleFactor=0;
	// fill vectors up with zeroes
	facets.resize(nFacets+1);
	vertices.resize(nVertices+1);
}; // ctor



// delete facets -- the std::vector's as such take care of themselves 
// note: vertices are NOT destroyed!
TriangulatedShape::~TriangulatedShape() throw()
{
	vector<FacetTypeBasic*>::iterator it, end=facets.end();
	for (it=facets.begin(); it!=end; it++)
		delete *it;
}; // TriangulatedShape::dtor



/*vector3 TriangulatedShape::barycenter()
{
	vector3 dummy(0,0,0);
	for (int i=1; i<= nVertices; i++)
		dummy+= *vertices[i];
	dummy/=nVertices;
	return dummy;
}; // vector3 TriangulatedShape::barycenter()
*/


/*void TriangulatedShape::barycenter2zero()
{
	vector3 center=barycenter();
	for (int i=1; i<=nVertices; i++)
		(*vertices[i])-=center;
	//center=barycenter(); //is now zero!
}; // TriangulatedShape::barycenter2zero()
*/


// calculates first the asteroid's volume
//		(by adding up tetrahedra's volumes as defined by 0 and the three vertices)
// and makes this an effective diameter 
//		(V=4/3*PI*r^3)
void TriangulatedShape::calculateIntrinsicDiameter() throw (invalid_argument, logic_error)
{
	if (nVertices<1 || nFacets<1 || nFacets!=(2*nVertices-4))
		throw invalid_argument("TriangulatedShape::calculateIntrinsicDiameter: problem with vertex/facet numbers!");
	if (intrinsicDiameter!=0)
		throw logic_error("TriangulatedShape::calculateIntrinsicDiameter: was already called, must be called only once!");
	double volume=0;
	const int* facetsVertices=0;
	for (int i=1; i<=nFacets; i++)
	{
		// array of the current facet's vertices
		facetsVertices=facets[i]->getVertices();
		const vector3& v1=*vertices[facetsVertices[1]];
		const vector3& v2=*vertices[facetsVertices[2]];
		const vector3& v3=*vertices[facetsVertices[3]];
		// += (v1xv2).v3
		volume+=fabs( vector3(v1,v2).dot(v3)); 
	};
	if (volume==0) throw invalid_argument("TriangulatedShape::calculateIntrinsicDiameter: volume is zero ?!?");
	//volume/=6; // now volume correct
	//intrinsicDiameter = pow( 6*volume/vector3::PI, 1./3.);
	intrinsicDiameter = pow(volume/vector3::PI, 1./3.); // optimized -> 6/6=1!
	setScaleFactor();
};



TriangulatedShape::FacetTypeBasic::FacetTypeBasic(int vertex1, int vertex2, int vertex3, 
							    const vector<vector3*>& vertexList)
							   throw(invalid_argument)
							   :dA(0)
{
	if(vertex1<1 && vertex2<1 && vertex3<1)
		throw invalid_argument("Vertex number < 1 encountered in FacetTypeBasic");
	vertexID[1]=vertex1; vertexID[2]=vertex2; vertexID[3]=vertex3; 
	vertexID[0]=0; // just in case

	// assuming vertexList is one-based
	const vector3& v1 = * vertexList[vertex1];
	const vector3& v2 = * vertexList[vertex2];
	const vector3& v3 = * vertexList[vertex3];

	// cast const'ness away, then determine value = +- cross product/2
	vector3*& dA = const_cast<vector3*&> (this->dA);
	dA=new vector3(v2-v1, v3-v1);
	*dA/=2;

	// Flip if dA looks inward: (assuming that coordinate system's origin is INSIDE)
	vector3 midpoint (v1);	midpoint+=v2;	midpoint+=v3;	midpoint/=3;
	if (vector3::dot(*dA, midpoint)<0)
		*dA*=-1;
}; // ctor FacetTypeBasic (int^3, vertexlist)



TriangulatedShape::FacetTypeBasic::FacetTypeBasic(int i1, int i2, int i3, 
							   const vector3* const dA)
							   throw(invalid_argument)
							   : dA(dA)
{
	if(i1<1 || i2<1 || i3<1)
		throw invalid_argument("Vertex number < 1 encountered in TriangulatedShape::FacetTypeBasic");
	vertexID[1]=i1; vertexID[2]=i2; vertexID[3]=i3; 
	vertexID[0]=0; // just in case
}; // ctor (int^3, dA)
