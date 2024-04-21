#include "ConvexFile.h"
#include<fstream>
#include<string>
#include<sstream>
#include"bfstream.h"
using namespace std;


ConvexFile::ConvexFile (const char filename[])
: dA(), intrinsicDiameter(0) // will be read in afterwards
{
	bifstream in(filename); // binary file
	if (!in)
		throw invalid_argument("Could not open " + (string)filename + " for reading!");
	int nFacets;
	in>>nFacets;
	if (!in || nFacets<1)
		throw invalid_argument("Parse error while reading #facets!");

	// cast constness away
	double& intrinsicDiameter = const_cast<double&> (this->intrinsicDiameter);
	in>>intrinsicDiameter;
	if (!in || intrinsicDiameter <= 0)
		throw invalid_argument("Parse error while reading in intrinsic diameter!");

	vector<vector3>& dA = const_cast<vector<vector3>&> (this->dA);
	dA.reserve(nFacets);

	double x,y,z;
	for (int i=0; i<nFacets; i++)
	{
		in>>x>>y>>z;
		if (!in)
		{
			ostringstream dummy;
			dummy<<"Parse error in facet line "<<i<<"!"<<ends;
			throw invalid_argument(dummy.str());
		};
		dA.push_back(vector3(x,y,z));
	};
	if (!in.eof())
	{
		char c;
		in>>c;
		if (!!in)
			throw invalid_argument("Parse error: found more facets than there should be!");
	};
};



/*void ConvexFile::skipComments(ifstream &in)
{
	while (in.peek() == '#')
	{
		char dummy[555];
		in.getline(dummy, 554);
	};
}; // skipComments
*/
