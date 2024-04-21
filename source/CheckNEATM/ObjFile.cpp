#include "ObjFile.h"
using namespace std;


ObjFile::ObjFile(const char* filename)
throw (invalid_argument, bad_alloc, range_error)
{
	try
	{
		vector3* dummyVector=0;
		Facet* dummyFacet=0;
		
		in.open(filename, ios::in);
		if (!in)
			throw invalid_argument("Cannot open file for reading");
		
		char c;
		double d1, d2, d3;
		int i1, i2, i3;
		bool keepLooping=true;
		nVertices=0;
		nFacets=0;

		do
		{
			skipComments();
			in>>c;
			if (in)
			{
				switch(c)
				{
				case 'v':
					in>>d1>>d2>>d3;
					if(in)
					{
						dummyVector = new vector3(d1, d2, d3);
						vertices.push_back(dummyVector);
						nVertices++;					
					}
					else
						throw invalid_argument(
							"Not a valid OBJ-file (Error in vertex-line)");
					break;
				case 'f':
					in>>i1>>i2>>i3;
					if(in)
					{
						dummyFacet = new Facet(i1, i2, i3);
						facets.push_back(dummyFacet);
						nFacets++;
					}
					else
						throw invalid_argument(
							"Not a valid OBJ-File (Error in facet-line)");
					break;
				default:
					if (in.eof())
						keepLooping=false;
					else
						throw invalid_argument(
							"Not a valid OBJ-file (with vertices and triangular faces, only, that is)");
				}// switch(c)
			}// if(in)
			else 
			{
				if(in.eof())
					keepLooping=false;
				else
					throw invalid_argument(
						"Not a valid OBJ-file (with vertices and triangular faces, only, that is)");	
			}// else if(in)
		} while (keepLooping);
	} // try
	catch(bad_alloc& exc)
	{
		cerr<<exc.what()<<endl;
		CleanUp();
		throw; // nothing I can do here
	}
	catch(invalid_argument &exc)
	{ 
		cerr<<exc.what()<<endl;
		CleanUp();
		throw; 
	}
	catch(exception& exc)
	{
		cerr<<"Unhandled exception in ObjFile:"<<endl;
		cerr<<exc.what()<<endl;
		CleanUp();
		throw;
	};
} // ObjFile(const char* filename)



ObjFile::~ObjFile()
{
	CleanUp();
} // ~ObfFile


inline void ObjFile::skipComments()
{
	static char dummystr[556];
	while(in.peek()=='#') 
		in.getline(dummystr, 555);
	return;
}; // ObjFile::skip_comments



void ObjFile::CleanUp()
{
	if (!facets.empty())
	{
		deque<Facet*>::iterator it, end=facets.end();
		for (it=facets.begin(); it!=end; it++)
			delete *it;
	};
	if (!vertices.empty())
	{
		deque<vector3*>::iterator it, end=vertices.end();
		for (it=vertices.begin(); it!=end; it++)
			delete *it;
	};
	facets.clear();
	vertices.clear();
	in.clear();
	in.close();
} // ObjFile::CleanUp



