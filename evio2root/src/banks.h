#ifndef banks_H
#define banks_H 1


// EVIO headers
#include "evioUtil.hxx"
#include "evioFileChannel.hxx"
using namespace evio;

// gemc headers
#include "gbank.h"
#include "outputFactory.h"

map<string, string> readConditions(evioDOMTree& EDT, double verbosity);
vector<string> getCondition(string);


gBank getBank(gBank bank, double verbosity);

// normal banks are just containers of leafs
map<string, double> getHeaderBank(evioDOMTree& EDT, gBank bank, double verbosity);
vector<string> getUserHeaderBankNames(evioDOMTree& EDT, gBank bank, double verbosity);
map<string, double> getUserHeaderBank(evioDOMTree& EDT, gBank bank, vector<string> names, double verbosity);

// integrated raw and digitized banks
vector<hitOutput> getRawIntDataBank(evioDOMTree& EDT, string hitType, map<string, gBank> *banksMap, double verbosity);
vector<hitOutput> getDgtIntDataBank(evioDOMTree& EDT, string hitType, map<string, gBank> *banksMap, double verbosity);

// generated particle infos
vector<generatedParticle> getGenerated(evioDOMTree& EDT, gBank bank, double verbosity);

// ancestor info
vector<ancestorInfo> getAncestors (evioDOMTree& EDT, gBank bank, double verbosity);


map<string, vector<hitOutput> > getRawIntDataBanks(evioDOMTree& EDT, vector<string> hitTypes, map<string, gBank> *banksMap, double verbosity);
map<string, vector<hitOutput> > getDgtIntDataBanks(evioDOMTree& EDT, vector<string> hitTypes, map<string, gBank> *banksMap, double verbosity);


/// \class Mevent
/// <b>Mevent</b>\n\n
/// This class contains the event banks, organized by name, and the event number.\n
class Mevent
{
	public:
		 Mevent(){;}
		~Mevent(){;}

	
		// construct and fill the event based on the banks map
		Mevent(evioDOMTree& EDT, vector<string> hitTypes, map<string, gBank> *banksMap, double verbosity);

		map<string, double> headerBank;
	
		map<string, vector<hitOutput> > rawBanks;
		map<string, vector<hitOutput> > dgtBanks;
};



#endif
