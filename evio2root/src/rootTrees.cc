// banks header
#include "rootTrees.h"

// C++ headers
#include <iostream>


rTree::rTree(string name, string description, double verb)
{
	tree = new TTree(name.c_str(), description.c_str());
	verbosity = verb;
}


void rTree::addVariable(string vname, string vType)
{
	
	if(vType.find("i") == 1)
	{
		vInts[vname] = new vector<int>;
		tree->Branch(vname.c_str(), &vInts[vname]);
	}
	else if(vType.find("d") == 1)
	{
		vDoubles[vname] = new vector<double>;
		tree->Branch(vname.c_str(), &vDoubles[vname]);
	}
	else if(vType.find("s") == 1)
	{
		vStrings[vname] = new vector<string>;
		tree->Branch(vname.c_str(), &vStrings[vname]);
	}
}

void rTree::insertVariable(string vname, string vType, double v)
{
	// variables that are not defined in the tree will not be written out
	if(vDoubles.find(vname) == vDoubles.end()) {
		if(verbosity > 2 && vType.find("d") == 1) {
			cout << " Warning: <double> variable >" << vname << "< is either not defined in the tree or not selected by user to be written out" << endl;
		}
		
	} else {
		if(vType.find("d") == 1) {
			vDoubles[vname]->push_back(v);
			weGotSomething = 1;
		}
	}
}
	   
void rTree::insertVariable(string vname, string vType, int v)
{
	if(verbosity > 2)
		if(vInts.find(vname) == vInts.end() && vType.find("i") == 1)
			cout << " Error: <int> variable >" << vname << "< not defined." << endl;

	if(vType.find("i") == 1)
		vInts[vname]->push_back(v);

	weGotSomething = 1;
}

void rTree::insertVariable(string vname, string vType, string v)
{
	if(verbosity > 2)
		if(vStrings.find(vname) == vStrings.end()  && vType.find("s") == 1)
			cout << " Error: <string> variable >" << vname << "< not defined." << endl;

	if(vType.find("s") == 1 )
		vStrings[vname]->push_back(v);

	weGotSomething = 1;
}


void rTree::fill()
{
// commenting the condition out cause ROOT
// does the right thing and add an empty entry
// to the tree if nothing there.

	if(weGotSomething)
		tree->Fill();
		
//	else
//		fillVoid();
}

void rTree::filluser()
{
// commenting the condition out cause ROOT
// does the right thing and add an empty entry
// to the tree if nothing there.

//  	if(weGotSomething){
 		tree->Fill();

}

void rTree::init()
{
	weGotSomething = 0;
	
	for(map<string, vector<double>* >::iterator it = vDoubles.begin(); it !=vDoubles.end(); it++)
		it->second->clear();
	
	for(map<string, vector<int>* >::iterator it = vInts.begin(); it !=vInts.end(); it++)
		it->second->clear();
	
	for(map<string, vector<string>* >::iterator it = vStrings.begin(); it !=vStrings.end(); it++)
		it->second->clear();
	
}

void rTree::fillVoid()
{
	for(map<string, vector<double>* >::iterator it = vDoubles.begin(); it !=vDoubles.end(); it++)
		it->second->push_back(-9999.0);
	
	for(map<string, vector<int>* >::iterator it = vInts.begin(); it !=vInts.end(); it++)
		it->second->push_back(-9999);
	
	for(map<string, vector<string>* >::iterator it = vStrings.begin(); it !=vStrings.end(); it++)
		it->second->push_back("na");

	tree->Fill();

}











