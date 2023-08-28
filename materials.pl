use strict;

our %configuration;

# One can change the "variation" here if one is desired different from the config.dat
# $configuration{"variation"} = "myvar";

sub define_bdx_materials
{
	# uploading the mat definition
	
	
	# Iron
	my %mat = init_mat();
	$mat{"name"}          = "BDX_Iron";
	$mat{"description"}   = "beam dump iron blocks";
	$mat{"density"}       = "7.874";
	$mat{"ncomponents"}   = "1";
	$mat{"components"}    = "G4_Fe 1.";
	print_mat(\%configuration, \%mat);
	

	# scintillator
	%mat = init_mat();
	$mat{"name"}          = "scintillator";
	$mat{"description"}   = "ft scintillator material C9H10 1.032 g/cm3";
	$mat{"density"}       = "1.032";
	$mat{"ncomponents"}   = "2";
	$mat{"components"}    = "C 9 H 10";
	print_mat(\%configuration, \%mat);


    
    # CsI(Tl)
    %mat = init_mat();
    $mat{"name"}          = "CsI_Tl";
    $mat{"description"}   = "Babar crystal CsI(Tl) 4.53 g/cm3";
    $mat{"density"}       = "4.53";
    $mat{"ncomponents"}   = "2";
    $mat{"components"}    = "G4_CESIUM_IODIDE 0.999 G4_Tl 0.001";
    print_mat(\%configuration, \%mat);
    
    
    
	# micromegas mylar
	%mat = init_mat();
	$mat{"name"}          = "bdx_mylar";
	$mat{"description"}   = "ft micromegas mylar 1.40g/cm3";
	$mat{"density"}       = "1.4";
	$mat{"ncomponents"}   = "3";
	$mat{"components"}    = "G4_H 0.041958 G4_C 0.625017 G4_O 0.333025";
	print_mat(\%configuration, \%mat);
    
    # BDX_Concrete
    %mat = init_mat();
    $mat{"name"}          = "BDX_Concrete";
    $mat{"description"}   = "Concrete";
    $mat{"density"}       = "2.7";
    $mat{"ncomponents"}   = "10";
    $mat{"components"}    = "G4_H 0.01 G4_C 0.001 G4_O 0.529107 G4_Na 0.016 G4_Mg 0.002 G4_Al 0.033872 G4_Si 0.337021 G4_K 0.013 G4_Ca 0.044 G4_Fe 0.014";
    print_mat(\%configuration, \%mat);
    

    # BDX_Pb_shots cponsidering rho=11.35/2=5.675 gr/cm3
    %mat = init_mat();
    $mat{"name"}          = "BDX_Pb_shots";
    $mat{"description"}   = "Small shots of lead";
    $mat{"density"}       = "5.675";
    $mat{"ncomponents"}   = "1";
    $mat{"components"}    = "G4_Pb 1.";
    print_mat(\%configuration, \%mat);
    # BDX_W_shots cponsidering  density -> 19.3/2=9.65 gr/cm3
    %mat = init_mat();
    $mat{"name"}          = "BDX_W_shots";
    $mat{"description"}   = "Small shots of tungstane";
    $mat{"density"}       = "9.65";
    $mat{"ncomponents"}   = "1";
    $mat{"components"}    = "G4_W 1.";
    print_mat(\%configuration, \%mat);
    
    # BDX_Dirt from A.s. Hoover et Al IEEE Symposium 2009
    %mat = init_mat();
    $mat{"name"}          = "BDX_Dirt";
    $mat{"description"}   = "Dirt";
    $mat{"density"}       = "1.9";
    $mat{"ncomponents"}   = "8";
    $mat{"components"}    = "G4_H 0.021 G4_C 0.016 G4_O 0.577 G4_Al 0.05 G4_Si 0.271 G4_K 0.013 G4_Ca 0.041 G4_Fe 0.011";
    print_mat(\%configuration, \%mat);
   
    # BDX_Dirt -CLAY- from L.Martin N.Mercier S.Incerti Geant4 simulations for sedimentary grains in infinite matrix conditions: The case of alpha dosimetry, Radiation Measurements 70 (2014) 39e47  55% in mass of SiO2, 35% of Al2O3 and 10% of Fe2O3.
    %mat = init_mat();
    $mat{"name"}          = "BDX_SiO2";
    $mat{"description"}   = "Silicon Dioxide";
    $mat{"density"}       = "2.65";
    $mat{"ncomponents"}   = "2";
    $mat{"components"}    = "O 2 Si 1";
    print_mat(\%configuration, \%mat);
    $mat{"name"}          = "BDX_Al2O3";
    $mat{"description"}   = "Al oxide";
    $mat{"density"}       = "3.95";
    $mat{"ncomponents"}   = "2";
    $mat{"components"}    = "Al 2 O 3";
    print_mat(\%configuration, \%mat);
    $mat{"name"}          = "BDX_Fe2O3";
    $mat{"description"}   = "Fe oxide";
    $mat{"density"}       = "5.24";
    $mat{"ncomponents"}   = "2";
    $mat{"components"}    = "Fe 2 O 3";
    print_mat(\%configuration, \%mat);

    
    
    %mat = init_mat();
    $mat{"name"}          = "BDX_Dirt_real";
    $mat{"description"}   = "Dirt";
    $mat{"density"}       = "1.9";
    $mat{"ncomponents"}   = "3";
    $mat{"components"}    = "BDX_SiO2 0.55 BDX_Al2O3 0.35 BDX_Fe2O3 0.1";
    print_mat(\%configuration, \%mat);
    
    
    # Epoxy
    %mat = init_mat();
    $mat{"name"}          = "Epoxy";
    $mat{"description"}   = "Epoxy";
    $mat{"density"}       = "1.16";
    $mat{"ncomponents"}   = "4";
    $mat{"components"}    = "H 32 N 2 O 4 C 15";
    print_mat(\%configuration, \%mat);
  
    
    
    # CarbonFiber
    %mat = init_mat();
    $mat{"name"}          = "CarbonFiber";
    $mat{"description"}   = "CarbonFiber";
    $mat{"density"}       = "1.750";
    $mat{"ncomponents"}   = "2";
    $mat{"components"}    = "G4_C 0.745 Epoxy 0.255";
    print_mat(\%configuration, \%mat);
   
    
    
    # Borotron 5% B
    %mat = init_mat();
    $mat{"name"}          = "Borotron";
    $mat{"description"}   = "Polyethylene + 5% B (weigth)";
    $mat{"density"}       = "0.94";
    $mat{"ncomponents"}   = "3";
    $mat{"components"}    = "G4_C 0.814 G4_H 0.136 G4_B 0.05";
    print_mat(\%configuration, \%mat);
    
   
	
}


