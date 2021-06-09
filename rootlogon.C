{
	cout<<"Loading Jixie's Style......done!"<<endl;
	if(!gSystem->AccessPathName("/home/jixie/MyStyle.h")) {
		gROOT->LoadMacro("/home/jixie/MyStyle.h");
		SetMyFitStyle();
	}
	else if(!gSystem->AccessPathName("/home/a-compton/companaEcal/script/MyStyle.h")) {
		gROOT->LoadMacro("/home/a-compton/companaEcal/script/MyStyle.h");
		SetMyFitStyle();
	}
	gSystem->Load("libXSTree");
}
