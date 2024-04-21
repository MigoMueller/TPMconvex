// bfstream.h

// Binary file streams

// In particular: bifstream and bofstream,
// binary sibblings of ifstream and ofstream, 
// where the <<>>-operators are overwritten to work in binary format.

// Author: michael.mueller@dlr.de
// Date:   2004 June 1


#if !defined(bfstream_h_INCLUDED)
#define bfstream_h_INCLUDED

#include<fstream>

namespace std
{
	template <typename stream>
		class binary_stream : public stream
	{
	public:
		binary_stream(const char* const filename, ios::openmode mode = ios::binary)
			: stream (filename, mode | ios::binary)
		{};
		binary_stream() : stream() {};

		virtual void open(const char* const filename, ios::openmode mode = ios::binary)
			{stream::open(filename, ios::binary | mode);};
		virtual ~binary_stream(){};

		// Note: Only a template -- doesn't harm if stream doesn't support that operator
		template <class type>
			binary_stream<stream>& operator<< (const type& content)
		{
		  this->write ( (const char*) &content, sizeof(content));
		  return *this;
		};

		// Note: Only a template -- doesn't harm if stream doesn't support that operator
		template <class type>
			binary_stream<stream>& operator>> (type& content)
		{
		  this->read ( (char*) &content, sizeof(content));
		  return *this;
		};
	}; // class_template binary_stream<stream>

	typedef binary_stream<ifstream> bifstream;
	typedef binary_stream<ofstream> bofstream;
	typedef binary_stream<fstream>  bfstream;
}; // namespace std

#endif // !defined (bfstream_h_INCLUDED)
