// XFileManager.cpp: implementation of the XFileManager class.
//
//////////////////////////////////////////////////////////////////////

#include "XFileManager.h"
#include "LagerrosLookedup.h"
#include <sstream>
#include <fstream>
#include <cmath>
using namespace std;


XFileManager::XFileManager (int gamma, unsigned int maxN)
:
files(),
maxN(maxN),
bluenessList(LagerrosLookedup::getBlueness()),
gammaList   (LagerrosLookedup::getGamma()),       // making lists static causes compiler errors - bummer!
gamma     (gamma),
directory (gammaList.dir(gamma))
{
	files.makeSize(maxN);
	currentN = 0;


}



void XFileManager::pop_first()
{
	if (currentN == 0)
		throw logic_error("Tried to delete more XFiles than there are!");
	delete files[0];
	for (int i=1; i<currentN; i++)
		files[i-1] = files[i];
	currentN--;
} // XFileManager::pop_first()



void XFileManager::push_back(unsigned int blueness)
{
	try
	{
		stringstream dummy;
		dummy<<directory<<'/'<<bluenessList.filename(blueness);
		{ // check it file exists at all
			ifstream in (dummy.str().c_str());
			if (!in)
			{
				stringstream dummy;
				dummy<<"Error in the internal file structure or non-existing blueness ("<<blueness<<") demanded!";
				throw logic_error (dummy.str().c_str());
			}
		};
		if (currentN == maxN)
			pop_first();
		files[currentN] = new XFile(dummy.str().c_str());
		currentN++;
	}
	catch (invalid_argument&) // file corrupted (thrown from XFile::ctor)
	{throw;}
	catch (logic_error&) // file doesn't exist
	{throw;}
	catch(...) // suppose: no memory available
	{
		if (currentN == 0) // no file left to release -- nothing I can do
			throw bad_alloc();
		pop_first(); // release memory
		push_back(blueness); // try again
	};
} // XFileManager::push_back



XFileManager::twoFilesPlusRemainder XFileManager::getFiles(double blueness)
{
	if (blueness < bluenessList.getFirst())
		throw invalid_argument("Too low a blueness chosen!");
	if (blueness > bluenessList.getLast())
		throw invalid_argument("Too high a blueness chosen!");
	double remainder = (blueness - bluenessList.getFirst()) / bluenessList.getStep();
	unsigned int index = floor(remainder);
	remainder -= index;

	const XFile& Xlow = getXFile (bluenessList.getFirst() + index     *bluenessList.getStep());
	const XFile& Xhigh= getXFile (bluenessList.getFirst() + (index+1) *bluenessList.getStep());

	return twoFilesPlusRemainder (remainder, Xlow, Xhigh);
} // XFileManager::getFiles()



const XFile& XFileManager::getXFile(int blueness)
{
	// check if corresponding file is already loaded
	for (int i=0; i<currentN; i++)
	{
		if (files[i]->getBlueness() == blueness) 
			return *files[i];
	};
	push_back(blueness); // load in file as last file
	return *files[currentN-1]; 
} // XFileManager::getXFile



void XFileManager::pop_all()
{
	for (int i=0; i<currentN; i++)
	{
		delete files[i];
		files[i]=0;
	};
	currentN=0;
}
