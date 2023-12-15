// C++ headers
#include <string>
#include <map>
#include <iostream>
using namespace std;

// banks header
#include "banks.h"
#include "string_utilities.h"


map<string, string> readConditions(evioDOMTree& EDT, double verbosity)
{
	map<string, string> simConditions;

	evioDOMNodeListP simCondBank = EDT.getNodeList(tagNumEquals(SIMULATION_CONDITIONS_BANK_TAG, 1));
	
	for(evioDOMNodeList::const_iterator iter=simCondBank->begin(); iter!=simCondBank->end(); iter++)
	{
		const vector<string> *vec = (*iter)->getVector<string>();
		
		for(unsigned i=0; i<vec->size(); i++)
		{
			vector<string> data = getCondition((*vec)[i]);
			simConditions[data[0]] = data[1];
			
		}
	}

	if(verbosity > 1 && simConditions.size())
	{
		cout << endl << "  Simulation Conditions Map: " << endl;
		for(map<string, string>::iterator it = simConditions.begin(); it != simConditions.end(); it++)
		{
			cout << "  > " << it->first << ": " << it->second << endl;
		}
	}
	
	return simConditions;
}

// condition key is string before ":" - value is after ":"
vector<string> getCondition(string input)
{
	vector<string> output;
	
	string thisInput = input;
	
	size_t keySeparator = thisInput.find(":");
	size_t optionStart  = 0;
		
	if(input.find("option") == 0)
	{
		optionStart = 7;
		thisInput = input.substr(optionStart, input.length());
		keySeparator = thisInput.find(":");
	}
	else if(input.find("parameter") == 0)
	{
		optionStart = 10;
		thisInput = input.substr(optionStart, input.length());
		keySeparator = thisInput.find(":");
	}
	
	output.push_back(input.substr(optionStart, keySeparator));
	output.push_back(input.substr(optionStart + keySeparator+1, input.length()));
	
	return output;
}


map<string, double> getHeaderBank(evioDOMTree& EDT, gBank bank, double verbosity)
{
	map<string, double> thisBank;
	unsigned long sizeOfBank = bank.name.size();
	
	evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(bank.idtag, 0));
	for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
	{
		const evioDOMNodeP node = *iter;
		if(node->isContainer())
		{
			evioDOMNodeList *variablesNodes = node->getChildList();
			for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
			{
				const evioDOMNodeP variable = *cIter;
				int vnum = variable->num;
				for(unsigned i=0; i<sizeOfBank; i++)
				{
					if(bank.gid[i] == vnum && variable->isLeaf())
					{
						// timestamp is special, and will be recorded on the key
						if(bank.type[i] == "Ns")
						{
							const vector<string> *vec = (*cIter)->getVector<string>();
							thisBank[(*vec)[0]] = -999;  // to identify it
						}
						
						if(bank.type[i] == "Ni")
						{
							const vector<int> *vec = (*cIter)->getVector<int>();
							thisBank[bank.name[i]] = (double) (*vec)[0];
						}
						
						if(bank.type[i] == "Nd")
						{
							const vector<double> *vec = (*cIter)->getVector<double>();
							thisBank[bank.name[i]] = (*vec)[0];
							
						}
					}
				}
			}
		}
	}
	
	if(verbosity > 0)
	{
		for(map<string, double>::iterator it = thisBank.begin(); it != thisBank.end(); it++)
		{
			cout << "  > " << it->first << ": " << it->second << endl;
		}
		cout << endl;
	}
	


	return thisBank;
}

vector<string> getUserHeaderBankNames(evioDOMTree& EDT, gBank bank, double verbosity)
{
  // Examine data to find out how many variables in userHeader bank and construct list of names

  vector<string> names;
  unsigned long sizeOfBank = 0;

  evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(bank.idtag, 0));
  
  for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
    {
      const evioDOMNodeP node = *iter;
      if(node->isContainer())
	{
	  evioDOMNodeList *variablesNodes = node->getChildList();
	  for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
	    {
	      const evioDOMNodeP variable = *cIter;
	      unsigned vnum = variable->num;
	      if (variable->isLeaf() && sizeOfBank < vnum)
		{
		  sizeOfBank = vnum;
		}
	    }
	}
    }

  for (unsigned i = 0; i < sizeOfBank; ++i)
    {
      string tmp = "userVar" ;
      if(i<9)        tmp +="00";
      else if (i<99) tmp +="0";
      
      tmp += to_string(i+1);
      names.push_back (tmp);
    }

  if (verbosity > 0)
    {
      cout << "getUserHeaderBankNames" << endl;
      for (vector<string>::iterator it = names.begin(); it != names.end(); it++)
	{
	  cout << "  > " << (*it) << endl;
	}
      cout << endl;
    }

  return names;
}

map<string, double> getUserHeaderBank(evioDOMTree& EDT, gBank bank, vector<string> names, double verbosity)
{
	map<string, double> thisBank;
	
	evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(bank.idtag, 0));

	for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
	{
	  const evioDOMNodeP node = *iter;
	  if(node->isContainer())
	    {
	      evioDOMNodeList *variablesNodes = node->getChildList();
	      for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
		{
		  const evioDOMNodeP variable = *cIter;
		  int vnum = variable->num;
		  if (variable->isLeaf())
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      thisBank[names[vnum-1]] = (*vec)[0];				  
		    }
		}
	    }
	}
	
	if(verbosity > 0)
	{
		for(map<string, double>::iterator it = thisBank.begin(); it != thisBank.end(); it++)
		{
			cout << "  > " << it->first << ": " << it->second << endl;
		}
		cout << endl;
	}
	


	return thisBank;
}

map<string, vector<hitOutput> > getRawIntDataBanks(evioDOMTree& EDT, vector<string> hitTypes, map<string, gBank> *banksMap, double verbosity)
{
	map<string, vector<hitOutput> > rawBanks;
	for(unsigned i=0; i<hitTypes.size(); i++)
		rawBanks[hitTypes[i]] = getRawIntDataBank(EDT, hitTypes[i], banksMap, verbosity);
	
	return rawBanks;
}



vector<hitOutput> getRawIntDataBank(evioDOMTree& EDT, string hitType, map<string, gBank> *banksMap, double verbosity)
{
	vector<hitOutput> thisDataBank;
	gBank thisHitBank = getBankFromMap(hitType, banksMap);	
	gBank rawBank = getBankFromMap("raws", banksMap);

	
	// getting raw data: bank tag, num is RAWINT_ID
	unsigned long sizeOfBank = rawBank.name.size();
	evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(thisHitBank.idtag + RAWINT_ID, 0));
	for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
	{
		const evioDOMNodeP node = *iter;
		if(node->isContainer())
		{
			evioDOMNodeList *variablesNodes = node->getChildList();
			for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
			{
				const evioDOMNodeP variable = *cIter;
				int vnum = variable->num;
				for(unsigned i=0; i<sizeOfBank; i++)
				{
					if(rawBank.gid[i] == vnum && variable->isLeaf())
					{
						
						if(rawBank.getVarType(rawBank.name[i]) == "i")
						{
							const vector<int> *vec = (*cIter)->getVector<int>();
							if(vec)
							{
								// resizing thisDataBank if necessary
								if(thisDataBank.size() == 0)
									thisDataBank.resize(vec->size());
								
								for(unsigned h=0; h<vec->size(); h++)
								{
									thisDataBank[h].setOneRaw(rawBank.name[i], (*vec)[h]);
								}
							}
						}
						if(rawBank.getVarType(rawBank.name[i]) == "d")
						{
							const vector<double> *vec = (*cIter)->getVector<double>();
							if(vec)
							{
								// resizing thisDataBank if necessary
								if(thisDataBank.size() == 0)
									thisDataBank.resize(vec->size());
								
								for(unsigned h=0; h<vec->size(); h++)
								{
									thisDataBank[h].setOneRaw(rawBank.name[i], (*vec)[h]);
								}
							}
						}
					}
				}
			}
		}
	}

	
	
	
	// this needs to be re-org a bit?
	if(verbosity > 0)
	{
		if(thisDataBank.size())
				cout << " >> Bank: " << thisHitBank.bankName <<  " raws: " << endl;
		for(unsigned h=0; h<thisDataBank.size(); h++)
		{
			map<string, double>  raws = thisDataBank[h].getRaws();
			if(raws.size())
			{
				for(map<string, double>::iterator it=raws.begin(); it != raws.end(); it++)
				{
					cout << "   > " << it->first << ": " << it->second << endl;
				}
				cout << endl;

			}
		}
		
	}
		
	
	return thisDataBank;
}


map<string, vector<hitOutput> > getDgtIntDataBanks(evioDOMTree& EDT, vector<string> hitTypes, map<string, gBank> *banksMap, double verbosity)
{
	map<string, vector<hitOutput> > dgtBanks;
	for(unsigned i=0; i<hitTypes.size(); i++)
		dgtBanks[hitTypes[i]] = getDgtIntDataBank(EDT, hitTypes[i], banksMap, verbosity);
	
	return dgtBanks;
}



vector<hitOutput> getDgtIntDataBank(evioDOMTree& EDT, string hitType, map<string, gBank> *banksMap, double verbosity)
{
	vector<hitOutput> thisDataBank;
	gBank thisHitBank = getBankFromMap(hitType, banksMap);
	gBank dgtBank = getDgtBankFromMap(hitType, banksMap);
	
	
	// getting raw data: bank tag, num is RAWINT_ID
	unsigned long sizeOfBank = dgtBank.name.size();
	evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(thisHitBank.idtag + DGTINT_ID, 0));
	for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
	{
		const evioDOMNodeP node = *iter;
		if(node->isContainer())
		{
			evioDOMNodeList *variablesNodes = node->getChildList();
			for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
			{
				const evioDOMNodeP variable = *cIter;
				int vnum = variable->num;
				for(unsigned i=0; i<sizeOfBank; i++)
				{
					if(dgtBank.gid[i] == vnum && variable->isLeaf())
					{
						if(dgtBank.getVarType(dgtBank.name[i]) == "i")
						{
							const vector<int> *vec = (*cIter)->getVector<int>();
							if(vec)
							{
								// resizing thisDataBank if necessary
								if(thisDataBank.size() == 0)
									thisDataBank.resize(vec->size());
								
								for(unsigned h=0; h<vec->size(); h++)
								{
									thisDataBank[h].setOneDgt(dgtBank.name[i], (*vec)[h]);
								}
							}
						}
						if(dgtBank.getVarType(dgtBank.name[i]) == "d")
						{
							const vector<double> *vec = (*cIter)->getVector<double>();
							if(vec)
							{
								// resizing thisDataBank if necessary
								if(thisDataBank.size() == 0)
									thisDataBank.resize(vec->size());
								
								for(unsigned h=0; h<vec->size(); h++)
								{
									thisDataBank[h].setOneDgt(dgtBank.name[i], (*vec)[h]);
								}
							}
						}
					}
				}
			}
		}
	}
	
	
	
	
	// this needs to be re-org a bit?
	if(verbosity > 0)
	{
		if(thisDataBank.size())
			cout << " >> Bank: " << thisHitBank.bankName <<  " dgt: " << endl;
		for(unsigned h=0; h<thisDataBank.size(); h++)
		{
			map<string, double>  dgts = thisDataBank[h].getDgtz();
			if(dgts.size())
			{
				for(map<string, double>::iterator it=dgts.begin(); it != dgts.end(); it++)
				{
					cout << "   > " << it->first << ": " << it->second << endl;
				}
				cout << endl;
				
			}
		}
		
	}
		
	return thisDataBank;
}



vector<generatedParticle> getGenerated(evioDOMTree& EDT, gBank bank, double verbosity)
{
	vector<generatedParticle> parts;
	vector<int> pid;
	vector<double> px;
	vector<double> py;
	vector<double> pz;
	vector<double> vx;
	vector<double> vy;
	vector<double> vz;
	vector<double> time;
	vector<int> multiplicity;

	
	
	unsigned long sizeOfBank = bank.name.size();

	evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(bank.idtag, 0));
	for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
	{
		const evioDOMNodeP node = *iter;
		if(node->isContainer())
		{
			evioDOMNodeList *variablesNodes = node->getChildList();
			for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
			{
				const evioDOMNodeP variable = *cIter;
				int vnum = variable->num;
				for(unsigned i=0; i<sizeOfBank; i++)
				{
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "pid")
					{
						const vector<int> *vec = (*cIter)->getVector<int>();
						for(unsigned int h=0; h<vec->size(); h++) pid.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "px")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) px.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "py")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) py.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "pz")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) pz.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "vx")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) vx.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "vy")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) vy.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "vz")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) vz.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "time")
					{
						const vector<double> *vec = (*cIter)->getVector<double>();
						for(unsigned int h=0; h<vec->size(); h++) time.push_back((*vec)[h]);
					}
					if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "multiplicity")
					{
						const vector<int> *vec = (*cIter)->getVector<int>();
						for(unsigned int h=0; h<vec->size(); h++) multiplicity.push_back((*vec)[h]);
					}
				}
			}
		}
	}
	
	
	
	for(unsigned p=0; p<pid.size(); p++)
	{
		generatedParticle part;
		
		part.PID = pid[p];
		part.momentum     = G4ThreeVector(px[p], py[p], pz[p]);
		part.vertex       = G4ThreeVector(vx[p], vy[p], vz[p]);

		part.time         = time[p];
		part.multiplicity = multiplicity[p];

		parts.push_back(part);
	}
	
	return parts;
}


vector<ancestorInfo> getAncestors(evioDOMTree& EDT, gBank bank, double verbosity)
{
  vector<ancestorInfo> trajs;
  vector<int> pid;
  vector<int> tid;
  vector<int> mtid;
  vector<double> trackE;
  vector<double> px;
  vector<double> py;
  vector<double> pz;
  vector<double> vx;
  vector<double> vy;
  vector<double> vz;
  
  unsigned long sizeOfBank = bank.name.size();
  
  evioDOMNodeListP thisBankNode = EDT.getNodeList(tagNumEquals(bank.idtag, 0));
  for(evioDOMNodeList::const_iterator iter=thisBankNode->begin(); iter!=thisBankNode->end(); iter++)
    {
      const evioDOMNodeP node = *iter;
      if(node->isContainer())
	{
	  evioDOMNodeList *variablesNodes = node->getChildList();
	  for(evioDOMNodeList::const_iterator cIter=variablesNodes->begin(); cIter!=variablesNodes->end(); cIter++)
	    {
	      const evioDOMNodeP variable = *cIter;
	      int vnum = variable->num;
	      for(unsigned i=0; i<sizeOfBank; i++)
		{
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "pid")
		    {
		      const vector<int> *vec = (*cIter)->getVector<int>();
		      for(unsigned int h=0; h<vec->size(); h++) pid.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "tid")
		    {
		      const vector<int> *vec = (*cIter)->getVector<int>();
		      for(unsigned int h=0; h<vec->size(); h++) tid.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "mtid")
		    {
		      const vector<int> *vec = (*cIter)->getVector<int>();
		      for(unsigned int h=0; h<vec->size(); h++) mtid.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "trackE")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) trackE.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "px")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) px.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "py")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) py.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "pz")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) pz.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "vx")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) vx.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "vy")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) vy.push_back((*vec)[h]);
		    }
		  if(bank.gid[i] == vnum && variable->isLeaf() && bank.name[i] == "vz")
		    {
		      const vector<double> *vec = (*cIter)->getVector<double>();
		      for(unsigned int h=0; h<vec->size(); h++) vz.push_back((*vec)[h]);
		    }
		}
	    }
	}
    }
		
  for(unsigned p=0; p<pid.size(); p++)
    {
      ancestorInfo part;
      
      part.pid = pid[p];
      part.tid = tid[p];
      part.mtid = mtid[p];
      part.trackE = trackE[p];
      part.p      = G4ThreeVector(px[p], py[p], pz[p]);
      part.vtx    = G4ThreeVector(vx[p], vy[p], vz[p]);
      
      trajs.push_back(part);
    }
  
  return trajs;
}







Mevent::Mevent(evioDOMTree& EDT, vector<string> hitTypes, map<string, gBank> *banksMap, double verbosity)
{
	headerBank = getHeaderBank(EDT, getBankFromMap("header", banksMap), verbosity);
	
	for(unsigned i=0; i<hitTypes.size(); i++)
	{
		rawBanks[hitTypes[i]] = getRawIntDataBank(EDT, hitTypes[i], banksMap, verbosity);
		dgtBanks[hitTypes[i]] = getDgtIntDataBank(EDT, hitTypes[i], banksMap, verbosity);
	}
	
}


