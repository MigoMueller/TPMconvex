// XFileManager.h: interface for the XFileManager class.
//
// Stores (and reads in) XFiles for later use in LagerrosLookedup
// 
//////////////////////////////////////////////////////////////////////

// Not supposed to be a base-class for anything, non-virtual dtor!!!


#if !defined(AFX_XFILEMANAGER_H__570D919A_9F4B_4BB7_AC70_FDF23371773D__INCLUDED_)
#define AFX_XFILEMANAGER_H__570D919A_9F4B_4BB7_AC70_FDF23371773D__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


//#include<vector>
#include "XFile.h"

class bluenessini;
class gammaini;


// lil dummy class returned by XFileManager.getFiles(double blueness) -- defined below
// contains all necessary information for interpolating in blueness
//class twoFilesPlusRemainder;

class XFileManager  
{
public:
	// lil dummy class returned by XFileManager.getFiles(double blueness) -- defined below
	// contains all necessary information for interpolating in blueness
	class twoFilesPlusRemainder;
protected:
	arrays::pointerArray<XFile*> files;
	const unsigned int maxN; // maximum number of XFiles to be stored (save on RAM!)
	unsigned int currentN;   // number of XFiles stored	
	const bluenessini& bluenessList;
	const gammaini&    gammaList;
	const int gamma;
	const std::string  directory;
public:
	// destroys all entries (save on RAM!)
	void pop_all();
	// returns two XFiles and a 0<'remainder'<1
	twoFilesPlusRemainder getFiles (double blueness);
	XFileManager (int gamma, unsigned int maxN=DefaultMaxN);
	~XFileManager(){};
protected:
	// gets XFile, either from array or calls push_back
	const XFile& getXFile(int blueness);
	// adds new XFile to end of list
	void push_back(unsigned int blueness);
	// deletes first entry in list, shifts element pointers by one position
	// throws logic_error, if there are no entries to the list
	void pop_first();
private: // prohibitive
	XFileManager();
	XFileManager (const XFileManager&);
	XFileManager& operator= (const XFileManager&);

	static const unsigned int DefaultMaxN;
public:
	class twoFilesPlusRemainder
	{
	public:
		const XFile& Xlow; // lower  blueness
		const XFile& Xhigh; // higher blueness
		//0<=interpol<=1
		const double interpol; // how far inbetween we are
		
		twoFilesPlusRemainder(double interpol, const XFile& Xlow, const XFile& Xhigh)
			: Xlow(Xlow), Xhigh(Xhigh), interpol(interpol)
		{};
	private: // prohibitive
		twoFilesPlusRemainder();
		//twoFilesPlusRemainder(const twoFilesPlusRemainder&); // copy-ctor IS necessary!!!
		twoFilesPlusRemainder& operator= (const twoFilesPlusRemainder&);
	}; // class XFileManager::twoFilesPlusRemainder

}; // class XFileManager




#endif // !defined(AFX_XFILEMANAGER_H__570D919A_9F4B_4BB7_AC70_FDF23371773D__INCLUDED_)
