#pragma once
#include<deque>

template<typename TYPE>
class deque_no_copy : 
	protected std::deque<TYPE*>
{
public:
	deque_no_copy<TYPE>(void) : std::deque<TYPE*> {};
	~deque_no_copy(void){};
	void push_back (const TYPE& val)
	{
		std::deque<TYPE*>.push_back(&val);
	};
	TYPE& back()
		{return *std::deque<TYPE*>.back()};
};
