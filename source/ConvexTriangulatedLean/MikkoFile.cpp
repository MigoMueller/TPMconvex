// MikkoFile.cpp: implementation of the MikkoFile class.
//
//////////////////////////////////////////////////////////////////////

#include "MikkoFile.h"
#include "vector3.h"
#include "Facet.h"
#include "TriangulatedShape.h"
using namespace std;

MikkoFile::~MikkoFile() throw()
{
	int i;
	for(i=0;i<nVertices;i++)
		delete vertices[i];
	for (i=0; i<nFacets; i++)
		delete facets[i];
	delete [] facets;
	delete [] vertices;
}; // dtor



MikkoFile::MikkoFile(const char* filename)
	throw(invalid_argument, bad_alloc)
{
	vertices=0;
	facets=0;
	nVertices=nFacets=0;
	try
	{
		ifstream in(filename);
		if (!in)
			throw invalid_argument("Cannot open file for reading");
		in>>nVertices>>nFacets;
		if (!in)
			throw invalid_argument("Couldn't read in number of vertices / facets");

		try
		{// guarantee that exception thrown is std::bad_alloc
			vertices = new vector3* [nVertices];
			facets   = new Facet*  [nFacets];
		}
		catch (...) 
		{throw bad_alloc();};
		
		{ // set pointers to 0
			int i;
			for (i=0; i<nVertices; i++)
			{
				vertices[i]=0;
				facets[i]=0;
			}
			for (i=nVertices; i<nFacets; i++)
				facets[i]=0;
		};
		{ // read in vertices
			double x,y,z;
			for (int i=0; i<nVertices; i++)
			{
				in>>x>>y>>z;
				if (!in)
					throw invalid_argument("Not a valid MikkoFile (error in vertex-line)");
				try {vertices[i]=new vector3(x,y,z);} catch(...){throw bad_alloc();};
			}
		}
		{ // read in facets
			int three, i1, i2, i3;
			for (int i=0; i<nFacets; i++)
			{
				in>>three>>i1>>i2>>i3;
				if (!in)
					throw invalid_argument("Not a valid MikkoFile (error in facet-line)");
				if ( i1<1 || i2<1 || i3<1 )
					throw invalid_argument("Not a valid MikkoFile (vertex# < 1 encountered)");
				if (i1>nVertices || i2>nVertices || i3>nVertices)
					throw invalid_argument("Not a valid MikkoFile (vertex# > max encountered)");
				if (three != 3)
					throw invalid_argument("Not a valid MikkoFile (facets must be triangular)");
				try {facets[i]=new Facet(i1, i2, i3);} catch(...){throw bad_alloc();};
			}
		};
		// Done!
	} // try
	catch(...)
	{
		if (vertices)
		{
			for (int i=0; i<nVertices; i++)
				delete vertices[i];
			delete [] vertices;
		}
		if (facets)
		{
			for (int i=0; i<nFacets; i++)
				delete facets[i];
			delete [] facets;
		}
		throw;
	};
}; // ctor(const char* filename)



void MikkoFile::writeObjFile(const char *const filename) const
{
	ofstream out(filename);
	int i;
	for (i=0; i<nVertices; i++)
		out<<"v  "<<vertices[i]->getX()<<"  "<<vertices[i]->getY()<<"  "<<vertices[i]->getZ()<<endl;
	for (i=0; i<nFacets; i++)
		out<<"f  "<<facets[i]->get0()<<"  "<<facets[i]->get1()<<"  "<<facets[i]->get2()<<endl;
}; //MikkoFile::writeObjFile;


