
// read in files in the 'convex'-format as produced by OBJ2convex

// Author: michael.mueller@dlr.de
// Date:   2004 May 25

// Note: the std::vector dA lifes exactly as long as the object itself -
// so, if you want it to last longer, you'll need to copy it!

#if !defined(AFX_CONVEXFILE_H__265E4E6B_6B94_4FF8_9D9A_1132E737D236__INCLUDED_)
#define AFX_CONVEXFILE_H__265E4E6B_6B94_4FF8_9D9A_1132E737D236__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


#include<vector>
#include "vector3.h"
//#include<fstream>


class ConvexFile  
{
public:
	ConvexFile (const char filename[]);
	virtual ~ConvexFile(){};
	const std::vector<vector3> dA;
	const double intrinsicDiameter;
private:
	ConvexFile();
	ConvexFile(const ConvexFile&);
	ConvexFile& operator= (const ConvexFile&);
};

#endif // !defined(AFX_CONVEXFILE_H__265E4E6B_6B94_4FF8_9D9A_1132E737D236__INCLUDED_)

