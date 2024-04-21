#include "checkNEATM.h"
#include<iostream>
#include<string>
#include<sstream>
using namespace std;

int main(int nargs, const char*const argsv[])
{
	const string directoryRoot = "P:/NewThermalModel/ClassHierarchyFortran/spheres+ellipsoids/";

	cerr<<"This is checkNEATM."<<endl;

	{ // elli1.25
		cerr<<"Go for elli1.25:"<<endl;
		const string elliRoot=directoryRoot+"elli1.25/";
		{
			// latitude90 first:
			const string lat90=elliRoot+"lat90/";
			try { checkNEATM checker(lat90); }
			catch (exception& exc)
			{
				cerr<<"Caught an exception while fitting to latitude 90:"<<endl;
				cerr<<exc.what()<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			}
			catch (...)
			{
				cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude 90:"<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			};
			cerr<<"Finished latitude 90"<<endl;
		};
		for (int latitude=0; latitude<90; latitude+=10)
		{
			ostringstream dummy;
			dummy<<elliRoot<<"lat"<<latitude<<'/';
			const string latRoot=dummy.str();
			for (int azi=0; azi<360; azi+=20)
			{
				ostringstream dummy;
				dummy<<latRoot<<"azi"<<azi<<'/';
				try
				{
					checkNEATM anonymous(dummy.str());
				}
				catch (exception& exc)
				{
					cerr<<"Caught an exception while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<exc.what()<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				}
				catch(...)
				{
					cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				};
			};
			cerr<<"Finished latitude "<<latitude<<endl;
		}
		cerr<<endl;
	}; // elli1.25



	{ // elli1.5
		cerr<<"Go for elli1.5:"<<endl;
		const string elliRoot=directoryRoot+"elli1.5/";
		{
			// latitude90 first:
			const string lat90=elliRoot+"lat90/";
			try { checkNEATM checker(lat90); }
			catch (exception& exc)
			{
				cerr<<"Caught an exception while fitting to latitude 90:"<<endl;
				cerr<<exc.what()<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			}
			catch (...)
			{
				cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude 90:"<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			};
			cerr<<"Finished latitude 90"<<endl;
		};
		for (int latitude=0; latitude<90; latitude+=10)
		{
			ostringstream dummy;
			dummy<<elliRoot<<"lat"<<latitude<<'/';
			const string latRoot=dummy.str();
			for (int azi=0; azi<360; azi+=20)
			{
				ostringstream dummy;
				dummy<<latRoot<<"azi"<<azi<<'/';
				try
				{
					checkNEATM anonymous(dummy.str());
				}
				catch (exception& exc)
				{
					cerr<<"Caught an exception while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<exc.what()<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				}
				catch(...)
				{
					cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				};
			};
			cerr<<"Finished latitude "<<latitude<<endl;
		}
		cerr<<endl;
	}; // elli1.5


	{ // elli1.75
		cerr<<"Go for elli1.75:"<<endl;
		const string elliRoot=directoryRoot+"elli1.75/";
		{
			// latitude90 first:
			const string lat90=elliRoot+"lat90/";
			try { checkNEATM checker(lat90); }
			catch (exception& exc)
			{
				cerr<<"Caught an exception while fitting to latitude 90:"<<endl;
				cerr<<exc.what()<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			}
			catch (...)
			{
				cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude 90:"<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			};
			cerr<<"Finished latitude 90"<<endl;
		};
		for (int latitude=0; latitude<90; latitude+=10)
		{
			ostringstream dummy;
			dummy<<elliRoot<<"lat"<<latitude<<'/';
			const string latRoot=dummy.str();
			for (int azi=0; azi<360; azi+=20)
			{
				ostringstream dummy;
				dummy<<latRoot<<"azi"<<azi<<'/';
				try
				{
					checkNEATM anonymous(dummy.str());
				}
				catch (exception& exc)
				{
					cerr<<"Caught an exception while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<exc.what()<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				}
				catch(...)
				{
					cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				};
			};
			cerr<<"Finished latitude "<<latitude<<endl;
		}
		cerr<<endl;
	}; // elli1.75



	{ // elli2.0
		cerr<<"Go for elli2.0:"<<endl;
		const string elliRoot=directoryRoot+"elli2.0/";
		{
			// latitude90 first:
			const string lat90=elliRoot+"lat90/";
			try { checkNEATM checker(lat90); }
			catch (exception& exc)
			{
				cerr<<"Caught an exception while fitting to latitude 90:"<<endl;
				cerr<<exc.what()<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			}
			catch (...)
			{
				cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude 90:"<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			};
			cerr<<"Finished latitude 90"<<endl;
		};
		for (int latitude=0; latitude<90; latitude+=10)
		{
			ostringstream dummy;
			dummy<<elliRoot<<"lat"<<latitude<<'/';
			const string latRoot=dummy.str();
			for (int azi=0; azi<360; azi+=20)
			{
				ostringstream dummy;
				dummy<<latRoot<<"azi"<<azi<<'/';
				try
				{
					checkNEATM anonymous(dummy.str());
				}
				catch (exception& exc)
				{
					cerr<<"Caught an exception while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<exc.what()<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				}
				catch(...)
				{
					cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				};
			};
			cerr<<"Finished latitude "<<latitude<<endl;
		}
		cerr<<endl;
	}; // elli2.0


	{ // elli2.25
		cerr<<"Go for elli2.25:"<<endl;
		const string elliRoot=directoryRoot+"elli2.25/";
		{
			// latitude90 first:
			const string lat90=elliRoot+"lat90/";
			try { checkNEATM checker(lat90); }
			catch (exception& exc)
			{
				cerr<<"Caught an exception while fitting to latitude 90:"<<endl;
				cerr<<exc.what()<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			}
			catch (...)
			{
				cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude 90:"<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			};
			cerr<<"Finished latitude 90"<<endl;
		};
		for (int latitude=0; latitude<90; latitude+=10)
		{
			ostringstream dummy;
			dummy<<elliRoot<<"lat"<<latitude<<'/';
			const string latRoot=dummy.str();
			for (int azi=0; azi<360; azi+=20)
			{
				ostringstream dummy;
				dummy<<latRoot<<"azi"<<azi<<'/';
				try
				{
					checkNEATM anonymous(dummy.str());
				}
				catch (exception& exc)
				{
					cerr<<"Caught an exception while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<exc.what()<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				}
				catch(...)
				{
					cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				};
			};
			cerr<<"Finished latitude "<<latitude<<endl;
		}
		cerr<<endl;
	}; // elli2.25



	{ // elli2.5
		cerr<<"Go for elli2.5:"<<endl;
		const string elliRoot=directoryRoot+"elli2.5/";
		{
			// latitude90 first:
			const string lat90=elliRoot+"lat90/";
			try { checkNEATM checker(lat90); }
			catch (exception& exc)
			{
				cerr<<"Caught an exception while fitting to latitude 90:"<<endl;
				cerr<<exc.what()<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			}
			catch (...)
			{
				cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude 90:"<<endl;
				cerr<<"Proceed with next item on list:"<<endl;
			};
			cerr<<"Finished latitude 90"<<endl;
		};
		for (int latitude=0; latitude<90; latitude+=10)
		{
			ostringstream dummy;
			dummy<<elliRoot<<"lat"<<latitude<<'/';
			const string latRoot=dummy.str();
			for (int azi=0; azi<360; azi+=20)
			{
				ostringstream dummy;
				dummy<<latRoot<<"azi"<<azi<<'/';
				try
				{
					checkNEATM anonymous(dummy.str());
				}
				catch (exception& exc)
				{
					cerr<<"Caught an exception while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<exc.what()<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				}
				catch(...)
				{
					cerr<<"Caught a strange exception (out-of-memory?) while fitting to latitude "<<latitude<<", azimuth "<<azi<<':'<<endl;
					cerr<<"Proceed with next item on list:"<<endl;
				};
			};
			cerr<<"Finished latitude "<<latitude<<endl;
		}
		cerr<<endl;
	}; // elli2.5











/*	switch (nargs)
	{
	case 2:
		directoryRoot=argsv[1];
		break;
	default:
		cerr<<"Usage: "<<endl;
		cerr<<"checkNEATM directoryroot,"<<endl;
		cerr<<"in directoryroot two files named alpha.ini and craters.ini are expected..."<<endl;
		return 0;
	}; // switch (nargs)
	try
	{
		//const string directoryRoot = "N:/ClassHierarchyFortran/spheres+ellipsoids/newspheres/";
		checkNEATM checker (directoryRoot);


	}
	catch (exception& exc)
	{
		cerr<<exc.what()<<endl;
		return -1;
	}
	catch (...)
	{
		cerr<<"Caught some strange exception!"<<endl;
		return -2;
	};
	cerr<<"Successfully created NEATMfit.dat in "<<directoryRoot<<'.'<<endl;
*/
	return 0;

};