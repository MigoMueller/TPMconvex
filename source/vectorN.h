// vectorT<T> - adds basic algebraic capabilities to vector<T>

// Author: michael.mueller@dlr.de
// Date:   Jan 20 2005

#if !defined (__vectorN_defined)
#define __vectorN_defined

#include<vector>
#include<stdexcept>

template<typename T>
class vectorT :
	public std::vector<T>
{
public:
	vectorT<T>(unsigned int N) 
		: std::vector<T>(N,0)
	{};

	vectorT<T>(const vectorT<T>& rhs)
		: std::vector<T>(rhs)
	{};

	vectorT<T>(const std::vector<T> rhs)
		: std::vector<T>(rhs)
	{};

	vectorT<T>& operator+= (const vectorT<T>& rhs) throw(std::invalid_argument)
	{
		if (this->size()!=rhs.size())
			throw std::invalid_argument("vectorT<T>: vectors added must be of equal size!");
		typename vectorT<T>::iterator my_it=this->begin(), ende=this->end();
		typename vectorT<T>::const_iterator rhs_it = rhs.begin();
		for ( ; my_it!=ende ; my_it++, rhs_it++)
			*my_it += *rhs_it;
		return *this;
	};
	
	vectorT<T>& operator-= (const vectorT<T>& rhs) throw(std::invalid_argument)
	{
		if (this->size()!=rhs.size())
			throw std::invalid_argument("vectorT<T>: vectors subtracted must be of equal size!");
		typename vectorT<T>::iterator my_it=this->begin(), ende=this->end();
		typename vectorT<T>::const_iterator rhs_it = rhs.begin();
		for ( ; my_it!=ende ; my_it++, rhs_it++)
			*my_it -= *rhs_it;
		return *this;
	};

	vectorT<T>& operator*= (T factor) throw()
	{
	  typename vectorT<T>::iterator it=this->begin(), ende=this->end();
	  for ( ; it!=ende; it++)
	    *it *= factor;
	  return *this;
	};

	vectorT<T>& operator/= (T factor) throw (std::invalid_argument)
	{
		if (factor==0)
			throw std::invalid_argument("vectorT<T>: attempted division by zero!");
		typename vectorT<T>::iterator it=this->begin(), ende=this->end();
		for ( ; it!=ende; it++)
			*it /= factor;
		return *this;
	};

	bool isNull() const throw()
	{
	  typename vectorT<T>::const_iterator it=this->begin();
	  while (it!=this->end())
	    {
	      if (*it != 0)
		return false;
	      it++;
	    }
	  return true;
	};

	virtual ~vectorT<T>(void){};
private:
	vectorT<T>();
}; // template <class T> vectorT

typedef vectorT<double> vectorN;

#endif //#if !defined...
