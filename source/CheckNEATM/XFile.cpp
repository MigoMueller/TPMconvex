// XFile.cpp: implementation of the XFile class.
//
//////////////////////////////////////////////////////////////////////

#include "XFile.h"
//using namespace arrays;
#include<fstream>
#include<strstream>
#include<sstream>
#include<iostream>
using namespace std;


XFile::XFile(const char filename[])
{
	ifstream inFile(filename);
	if (!inFile)
	{
		strstream dummy;
		dummy<<"Couldn't open "<<filename<<" for opening!";
		throw invalid_argument(dummy.str());
	};

	stringstream in;
	//in.rdbuf(inFile.rdbuf()); // wouldn't work -- why?
	in << inFile.rdbuf();
	
	while (in.peek() == '#')
	{
		char dummy[513];
		in.getline(dummy, 512);
		if (!in)
		{
			strstream dummy;
			dummy<<"Error while parsing comments in file "<<filename;
			throw invalid_argument(dummy.str());
		};
	};

	in>>gamma>>blueness>>nSun>>nEarth>>nAzi;
	if (!in)
	{
		strstream dummy;
		dummy<<"Error while reading in header of file "<<filename;
		throw invalid_argument(dummy.str());
	};
	sunStep   = 1./nSun;
	earthStep = 1./nEarth;
	aziStep   = 2./(nAzi-1);

	nSun++;
	nEarth++; // field will include cosSun==0 and cosEarth==0!
	field.makeSize(nSun, nEarth, nAzi);
	
	{ // set factors for cosSun=0 by hand
		// UPDATE: set these factors to 1 instead of 0!!!
		for (int j=0; j<nEarth; j++)
		{
			for (int k=0; k<nAzi; k++)
				entry(0, j, k) = 1;
		};
	};
	{ // read in factors
		double dummy;
		for (int i=1; i<nSun; i++)
		{
			in>>dummy; // cosSun -- disregarded here
			for (int k=0; k<nAzi; k++)
				entry (i, 0, k) = 1; // set entry for cosEarth==0 by hand
			for (int j=1; j<nEarth; j++)
			{ // read in one set of values for cosEarth
				in>>dummy; // cosEarth -- disregarded
				for (int k=0; k<nAzi; k++)
					in>>dummy>>entry(i,j,k);
			}
		};
	};
	// done reading in -- check whether eof is reached
	if (!in.eof()) // maybe there is some blank left
	{
		// has some error occurred so far?
		if (!in)
		{
			stringstream dummy;
			dummy<<"Some error occurred while reading in "<<filename<<" -- must be malformatted!"<<ends;
			throw logic_error (dummy.str().c_str());
		};
		char c;
		in>>c; 
		if (! (!in)) 
		{ // there must be no char left to read in!
			stringstream dummy;
			dummy<<filename
				<<" is longer than it's supposed to be - comments are only allowed at the beginning of the file!"<<ends;
			throw logic_error (dummy.str().c_str());
		};// reading in a char didn't work, so eof is reached, indeed
	};
}; //XFile::ctor
	

	
