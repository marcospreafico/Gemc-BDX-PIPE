
{
	
	
	
	TFile f("muon_bdx.root");
	
	
	vector<double> *hitn;
	vector<double> *layer;
	vector<double> *adcl;

	TTree *cormo = (TTree*) f.Get("cormo");
	
	
	
	cormo->SetBranchAddress("hitn",    &hitn);
	cormo->SetBranchAddress("layer",   &layer);
	cormo->SetBranchAddress("adcl",  &adcl);

	
	
	for(int e=0; e<40; e++)
	{
		hitn->clear();
		layer->clear();
		adcl->clear();

		cout << " Event: " << e+1 << endl;

		cormo->GetEvent(e);
	
	
		cout << " Hit Size: " << hitn->size() << endl;
	
	
		for(unsigned i=0; i<hitn->size(); i++)
		{
			cout << (*hitn)[i] << " " <<  (*adcl)[i] << " " << (*layer)[i] << endl;
		}
		
		
		cout << " end of event " << e+1 << endl;
	}
	
}