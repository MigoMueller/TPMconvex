

#if !defined(_CHECK_NEATM_INCLUDED)
#define _CHECK_NEATM_INCLUDED

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include<map>
#include<list>
#include<deque>
//#include"NEATMFile.h"


class checkNEATM  
{// define sub-class first
	class nameserver
	{// define sub-classes first
	public:
		class alphaIni
		{
		public:
			typedef std::map<int, std::string>::iterator       iterator;
			typedef std::map<int, std::string>::const_iterator const_iterator;
			typedef std::pair<iterator, bool>                  resultType; 
			typedef std::pair<int, std::string>                dataType;
		public:
			alphaIni(const std::string& fileName);
			//		const string& dirName(int alpha) const;
			const std::string& dirName (int alpha) const;
			// copies alphas into a list and returns it
			std::list<int> getAlphas() const;
			double getH () const {return H;};
			double getG () const {return G;};
			double getPV() const {return pV;};
			double getR () const {return rAU;};
			double getDelta() const {return deltaAU;};
		protected:
			double H, G, pV, rAU, deltaAU;
			std::map<int,std::string> data;
			// resultType data.insert (const dataType&):
		private:
			alphaIni();
			alphaIni(const alphaIni&);
			alphaIni operator= (const alphaIni&);		
		}; // class nameserver::alphaIni

		class cratersIni
		{
		public:
			typedef std::map<int, std::string>  innerMap;
			typedef innerMap::const_iterator    const_iteratorI;
			typedef innerMap::iterator          iteratorI;
			typedef std::pair<int, std::string> dataTypeI;
			typedef std::pair<iteratorI, bool>  resultTypeI;
			
			typedef std::map<int, innerMap>     outerMap;
			typedef outerMap::const_iterator    const_iteratorO;
			typedef outerMap::iterator          iteratorO;
			typedef std::pair<int, innerMap>    dataTypeO;
			typedef std::pair<iteratorO, bool>  resultTypeO;
		public:
			cratersIni(const std::string& fileName);
			// copies gammas and densities into a multimap and returns it
			std::multimap<int,int> getGammaDensity () const;
			const std::string& fileName (int openingAngleDeg, int craterDensity) const;
			//const outerMap& getData() const {return data;};
		protected:
			outerMap data;
		private:
			cratersIni();
			cratersIni(const cratersIni&);
			cratersIni operator=(const cratersIni&);
		}; // class nameserver::cratersIni
	// continue declaration of nameserver
	public:
		nameserver (const std::string& directoryRoot);

		// returns the (absolute) name of the file corresponding to the parameters
		const std::string operator () 
			(int phaseDeg, int openingAngleDeg, int craterDensityPercent) const;
	protected:
	private:
		nameserver();
		nameserver(const nameserver&);
		const nameserver operator=(const nameserver&);

		const std::string directoryRoot;
		alphaIni   alpha;
		cratersIni craters;
	public:
		const std::list<int>         alphaList;
		const std::multimap<int,int> GammaDensity;
	public:
		// returns rhs if it ends with '/', else rhs+'/'
		static inline const std::string verifyDir (const std::string& rhs);
	}; // class checkNEATM::nameserver

public:
	checkNEATM (const std::string& directoryRoot);
	~checkNEATM();
	
protected:
	//NEATMFile*   NEATM;
	nameserver   filenames;
	inline void CleanUp();
private: // prohibitive
	checkNEATM();
	checkNEATM (const checkNEATM&);
	checkNEATM operator= (const checkNEATM&);	
}; // class checkNEATM



#endif // !defined(_CHECK_NEATM_INCLUDED)
