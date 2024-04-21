// Utilities for handling of files and strings
// Author: Michael.Mueller@DLR.de
// Date:   2006/11/12

#if !defined(__MM_TEXTUTILS_H_INCLUDED)
#define __MM_TEXTUTILS_H_INCLUDED
#include<string>
#include<sstream>
#include<algorithm>

namespace MM_textutils
{
    // skips lines starting with any of the delimiters passed
    inline std::istream& skipCommentLines(std::istream& in, const std::string& delims);
	// Next lines starting with ONE delim will be read out and skipped
	// A line starting with TWO delims will be kept, but the two delims will be removed
	inline std::istream&  skipCommentLines12(std::istream& in, const std::string& delims);

    // cuts off trailing blanks (as specified in string 'delims')
    inline void trim_right(std::string& in, const std::string& delims = " \t\n\r");
    // cuts off leading blanks (as specified in string 'delims')
    inline void trim_left(std::string& in, const std::string& delims = " \t\n\r");
    // cuts off after the first comment delimiter as specified
    inline void comments_out(std::string& in, const std::string& delims = "#");
    // makes lower case
    inline void to_lower_case(std::string& in);
    // makes upper case
    inline void to_upper_case(std::string& in);

	// extracts filename-extension (including path) and extension from in-string
	// returns the delimiter found; 0 if none was found.
	inline char name_extension (const std::string& in, std::string& name, std::string& extension, 
		const std::string& delims = ".");
	// extracts path (without trailing delim) and proper file name from in-string
	// returns the delimiter found; 0 if none was found.
	inline char path_file       (const std::string& in, std::string& path, std::string& filename,
		const std::string& delims = "/\\");
}; // namespace MM_textutils


// inline implementations:
inline std::istream& MM_textutils::skipCommentLines(std::istream& in, const std::string& delims)
{
    static std::string dummystr;
    while(std::find(delims.begin(), delims.end(), in.peek()) != delims.end()) 
        getline(in, dummystr);
	return in;
}; //MM_textutils::skipComments

inline std::istream& MM_textutils::skipCommentLines12(std::istream& in, const std::string& delims)
{
	static std::string dummystr;
	static char c;
	while(std::find(delims.begin(), delims.end(), in.peek()) != delims.end())
	{
		in>>c; // read out delim
		if (std::find(delims.begin(), delims.end(), in.peek()) != delims.end())
		{ // second delim present: eat 2nd delim and return
			in>>c; 
			return in;
		}
		else
			getline(in, dummystr); // eat whole line + repeat loop
	}
	return in;
}; //MM_textutils::skipCommentLines12


inline void MM_textutils::trim_right(std::string& in, const std::string& delims)
{
    static std::string::size_type pos;
    pos = in.find_last_not_of(delims);
    if (pos == in.npos)  // if there are only blanks in the string
        in.erase();      // erase the whole lot
    else
        in.erase(pos+1); // delete the trailing blanks
}; // MM_textutils::trim_right

inline void MM_textutils::trim_left(std::string& in, const std::string& delims)
{
    static std::string::size_type pos;
    pos = in.find_first_not_of(delims);
    in.erase(0,pos);  
}; // MM_textutils::trim_left

inline void MM_textutils::comments_out(std::string& in, const std::string& delims)
{
    static std::string::size_type pos;
    pos = in.find_first_of(delims);
    if (pos != in.npos)  // if there are comments
        in.erase(pos); 
}; // MM_textutils::comments_out

inline void MM_textutils::to_lower_case(std::string& in)
    {std::transform(in.begin(), in.end(), in.begin(), tolower);};

inline void MM_textutils::to_upper_case(std::string& in)
    {std::transform(in.begin(), in.end(), in.begin(), toupper);};

inline char MM_textutils::name_extension (const std::string& in, std::string& name, std::string& extension, 
		const std::string& delims)
{
	static std::string::size_type pos;
	pos = in.find_last_of(delims);
	if (pos == in.npos) // delim not found
	{
		name = in;
		extension = "";
		return 0;
	}
	name = in.substr(0,pos);
	extension = in.substr(pos+1);
	return in[pos];
}; // MM_textutils::name_extension

inline char MM_textutils::path_file (const std::string& in, std::string& path, std::string& filename,
		const std::string& delims)
{
	char c = name_extension(in, path, filename, delims);
	if (c == 0) // no delim found --> filename comes back empty, but path should
		path.swap(filename);
	return c;
} //MM_textutils::path_file



#endif // #if !defined(__MM_TEXTUTILS_H_INCLUDED)
