#ifndef root_H
#define root_H 1

// ROOT Headers
#include "TFile.h"
#include "TTree.h"

// C++ headers
#include <string>
#include <map>
#include <vector>
using namespace std;



class rTree
{
	public:
		rTree(){;}
		rTree(string, string, double verb);
		~rTree(){;}

		double verbosity;
	
		TTree *tree;
	
		// vType contains "d", "i", "s" respectively for:
		map<string, vector<double>* > vDoubles;
		map<string, vector<int>*>     vInts;
		map<string, vector<string>*>  vStrings;

		void addVariable(string vname, string vType);

		void insertVariable(string vname, string vType, double v);
		void insertVariable(string vname, string vType, int v);
		void insertVariable(string vname, string vType, string v);

		bool weGotSomething;
		void init();
		void fill();
		void filluser();
		void fillVoid();

};



#endif

