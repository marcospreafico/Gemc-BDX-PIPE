use strict;
use warnings;

our %configuration;

my $degrad = 0.01745329252;
my $cic    = 2.54;

my $a = $configuration{"a"}/2.;
my $A = $configuration{"A"}/2.;
my $wrap = $configuration{"wrap"}/2.;
my $l = $configuration{"l"}/2.;

my $L = $configuration{"L"};
my $n_upper = $configuration{"n_upper"};
my $n_lower = $configuration{"n_lower"};


my $flag_rotation = 0;  # if this flag is enabled, turn one row of crystal every two in the central module
if ($configuration{"variation"} eq "rotated") {$flag_rotation = 1;}

if ($configuration{"variation"} eq "ideal") {
    $a = 2.9/2; #2.9/2; #side of crystal in central module
    $A = 2.9/2; #side of crystal in outer layers
    $wrap = 0.1; #2; #wrap thickness
    $l = 20.0/2; #crystal lenght
    $L = 22.5; #0;  #side of central module
    $n_lower = 5; #number of crystal of the first outer layer
    $n_upper = 0; #number of crystal of the second outer layer
}

if($configuration{"variation"} eq "threeveto"){
     $a = 2.9/2; #2.9/2; #side of crystal in central module
     $A = 2.9/2; #side of crystal in outer layers
     $wrap = 0.2; #2; #wrap thickness
     $l = 20.0/2; #crystal lenght
     $L = 16.5; #0;  #side of central module
     $n_lower = 5; #number of crystal of the first outer layer
     $n_upper = 1; #number of crystal of the second outer layer
}

#Crystal parameter

#my $n_lower = 5; #number of crystal of the first outer layer
#my $n_upper = 1; #number of crystal of the second outer layer


#aggiungere flag con cui fare i cristalli tutti in verticale o alternando verticale orizzontale



###########################################################################################
# Build Crystal Volume and Assemble calorimeter
###########################################################################################

sub make_main
{
    my %detector = init_det();
    $detector{"name"}        = "main_volume";
    $detector{"mother"}      = "root";
    $detector{"description"} = "World";
    $detector{"color"}       = "666666";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    
    my $X = 0.;
    my $Y = 0.;
    my $Z = 0.;
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    my $par1 = 100.;
    my $par2 = 200.;
    my $par3 = 400.;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_Galactic";
    print_det(\%configuration, \%detector);
}

sub make_vessel
{
    my %detector = init_det();
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.; # depth in the pipe
    my $rotX=90.;
    my $rotY=0.;
    my $rotZ=0.;
    
    $detector{"mother"}      = "main_volume";   
    
    my $vs_top_tk=0.5/2;
    my $vs_side_tk=2.0;
    my $vs_lg=89.0/2;
    my $vs_od = 16*2.54; #vessel outer diameter 
    
    my $vs_or = $vs_od/2; #outer radius
    my $vs_ir = $vs_or - $vs_side_tk; #inner radius
    
    my $par4 = 0.;    my $par5 = 360.;
    
    $detector{"name"}        = "mutest_vessel_mother";
    $detector{"description"} = "Mother volume for the vessel";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$rotX*deg $rotY*deg $rotZ*deg";
    $detector{"dimensions"}  = "0*cm $vs_or*cm $vs_lg*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);

    $detector{"name"}        = "mutest_vessel";
    $detector{"mother"}      = "mutest_vessel_mother";
    $detector{"description"} = "Vessel";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$vs_ir*cm $vs_or*cm $vs_lg*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "mutest_vessel_air";
    $detector{"mother"}      = "mutest_vessel_mother";
    $detector{"description"} = "air in vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $vs_ir*cm $vs_lg*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_AIR";
    print_det(\%configuration, \%detector);
    
    $Z = $vs_lg-$vs_top_tk ;
    
    $detector{"name"}        = "mutest_vessel_top";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $vs_ir*cm $vs_top_tk*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);

    $Z= -$Z ;

    
    $detector{"name"}        = "mutest_vessel_bottom";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $vs_ir*cm $vs_top_tk*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    

}

sub make_veto
{
    my $or = $_[0]; 
    my $vetoid = $_[1]; 
    
    my %detector = init_det();
    my $X = 0. ;
    my $Y = 0. ;
    my $Z = 0.; # depth in the pipe
    my $rotX=0.;
    my $rotY=0.;
    
    my $rotZ=180.;
    
    $detector{"mother"}      = "mutest_vessel_air";   
    
    my $ov_top_tk=0.8;
    my $ov_side_tk=0.8*2;
    my $ov_lg= 53.;
    

    my $ov_or= $or; #(16*2.54-4);
    my $ov_ir = $ov_or-$ov_side_tk;  #OutRad-thick
    
    my $par1 = $ov_ir/2; #InRad
    my $par2 = $ov_or/2; #InRad
    my $par3  =$ov_lg/2.;#length
    my $par4 = 0.;
    my $par5 = 360.;
    ###
    
    my $NchOV = 8; # MAX 8 otherwise overlaps with OV_TOP/BOTTOM (8,9)
    my $DeltaAng = 360.0/$NchOV;
    my $DeltaToll = 0.;#Tolerance (in cm) between adiacent sectors
    my $DeltaAngShift = 0.; #fixed shift
    $DeltaAngShift = 90. -$DeltaAng/2; #to have #1 perfectly back
    ###
    for(my $ib=1; $ib<($NchOV+1); $ib++)
	#for(my $ib=1; $ib<(2); $ib++)
    {

        $par4=($ib-1)*$DeltaAng+$DeltaAngShift+$DeltaToll/2;
        $par5=$DeltaAng-$DeltaToll;# -0.1deg tolerance
	$detector{"name"}        = "mutest_OV_$ib"."_"."$vetoid";
	$detector{"description"} = "OV in the pipe"." "."$vetoid";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Tube";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	if($vetoid==1){
	    $detector{"rotation"}    = "0*deg 0*deg 22.5*deg";
	}
	$detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm $par4*deg $par5*deg";
	if($vetoid != 2){
	    $detector{"color"}       = "ff8000";
	    $detector{"material"}    = "ScintillatorB";
	    $detector{"sensitivity"} = "veto";
	    $detector{"hit_type"}    = "veto";
	    my $veto_channel = 100*$vetoid+$ib; 
	    $detector{"identifiers"} = "sector manual 0 veto manual 7 channel manual  $veto_channel";
	}
	elsif($vetoid == 2){
	    
	    $detector{"color"}       = "808080";
	    $detector{"material"} = "G4_W"; 
	}
	
	print_det(\%configuration, \%detector);
        #print  " OV segments = ",$par4,"-",$par5,"\n";
    }
    print  " OV INNER radius = ",$par1,"\n";
    print  " OV OUTER radius = ",$par2,"\n";
    print  " OV segmentations Nch = ",$NchOV,"\n";
    
    $Z = $ov_lg/2+$ov_top_tk/2 ;
    $par2 = $ov_or/2; #InRad
    $par3  =$ov_top_tk/2;#thickness
    $par4 = 0.;
    $par5 = 360.;


    $detector{"name"}        = "mutest_OV_top"."_"."$vetoid";
    $detector{"mother"}      = "mutest_vessel_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "cc6804";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    my $top_channel = 100*$vetoid+9; 
    $detector{"identifiers"} = "sector manual 0 veto manual 7 channel manual $top_channel";
    #print_det(\%configuration, \%detector);
    
    $Z = -($ov_lg/2-$ov_top_tk/2) ;
    $par2 = $ov_ir/2; #InRad
    
    $detector{"name"}        = "mutest_OV_bottom"."_"."$vetoid";
    $detector{"mother"}      = "mutest_OV_air";
    $detector{"description"} = "top vessel ";
    $detector{"color"}       = "cc6804";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "0*cm $par2*cm $par3*cm $par4*deg $par5*deg";
    $detector{"material"}    = "ScintillatorB";
    $detector{"sensitivity"} = "veto";
    $detector{"hit_type"}    = "veto";
    my $bottom_channel = 100*$vetoid+10; 
    $detector{"identifiers"} = "sector manual 0 veto manual 7 channel manual 10";
    #print_det(\%configuration, \%detector);
}

sub make_ecal
{
   
    #crystal parameters:
    #my $a = 2.9/2; #2.5/2;
    my $aa = $a + $wrap;
    #my $l = 20.0/2;
    #my $L = 19.0; #20.0; #side of central module  

    my $X = 0.; 
    my $Y = 0.; 
    my $Z = 0.; 

    my $module = $_[0];
    
    
    my $chcount = 0; 

    my $n_core = int($L/(2.*$aa)); #number of crystals inside the central core
    
    my $outer_pos = $n_core*$aa + $aa ;
    #my $n_lower = 5; #I set these by hand
    #my $n_upper = 0; 
    
    #central core 
    for(my $ix = 0; $ix< $n_core; $ix++){
	for(my $iy = 0; $iy< $n_core; $iy++){
	    my $X = -$n_core*$aa + $aa + 2*$ix*$aa;#-$L/2+ $a + $ix * 2 * $a;
	    my $Y = -$n_core*$aa + $aa + 2*$iy*$aa;#-$L/2+ $a + $iy * 2 * $a;
        #if detector rotation is enabled the distance between matrixes is determined by the central square side
        if($module == 0){
            $Z = -$l-$wrap;
            if ($flag_rotation == 1) {$Z = -$L/2.-$wrap}
        }elsif($module == 1) {
        $Z = +$l+$wrap;
            if ($flag_rotation == 1) {$Z = +$L/2.+$wrap}
        }
        
        #correct coordinates for rotated rows
        if ( $flag_rotation == 1 && $ix % 2 == 0){
            $Z = $Z-$n_core*$aa + $aa + 2*$iy*$aa;
            $Y = 0;
        }
        my %detector = init_det();
        my $thisx = 10*$ix;
        my $thisy = 10*$iy;
	    $detector{"name"}        = "crs_central_"."$thisx"."_"."$iy"."_"."$module";
	    $detector{"mother"}      = "mutest_vessel_air";
        $detector{"description"} = "Panda PbWO4 crystal central"."$thisx"."$iy"."$module";
	    $detector{"color"}       = "1c86ea"; #"000000";#"ffffff";
	    $detector{"style"}       = 0;
	    $detector{"visible"}     = 1;
	    $detector{"type"}        = "Box";
	    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        if($flag_rotation == 1 && $ix % 2 ==0 ){
            $detector{"rotation"} = "90*deg 0*deg 0*deg";
        }
	    $detector{"dimensions"}  = "$a*cm $a*cm $l*cm";
	    $detector{"material"}    = "G4_PbWO4";
	    $detector{"sensitivity"} = "crs";
	    $detector{"hit_type"}    = "crs";
        my $man=400+$chcount;
	    $chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	    print_det(\%configuration, \%detector);
	}
    }
    
    #reset Z value for outer crystals
    if($module == 0){
        $Z = -$l-$wrap;
        if ($flag_rotation == 1) {$Z = -$L/2.-$wrap}
    }elsif($module == 1) {
    $Z = +$l+$wrap;
        if ($flag_rotation == 1) {$Z = +$L/2.+$wrap}
    }
    
    #top side 
    #lower part
    #my $A = 2.9/2;
    my $AA = $A+$wrap;
    for(my $ix = 0; $ix< $n_lower; $ix++){
	#my $X = $L/2+ $AA ;
       $X = $outer_pos ;
    my $Y = -$n_lower*$AA + $AA+ 2*$ix*$AA;
	#-$L/2+ 3*$a + $ix * 2 * $a; 
	my %detector = init_det();
        my $thisy = $n_core*10+0;
        my $thisx = ($n_core-$n_lower)*5+10*$ix;
	$detector{"name"}        = "crs_top_"."$thisx"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal top"."$thisx"."$thisy"."$module";
	$detector{"color"}       = "1c86ea"; #"000000";#"ff0000";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
	$detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }
    #upper part
    for(my $ix = 0; $ix< $n_upper; $ix++){
	my $X = $L/2+ 3*$AA ; 
	my $Y = -$n_upper*$AA + $AA + $ix * 2 * $AA; 
	my %detector = init_det();
        my $thisy = $n_core*10+10;
        my $thisx = ($n_core-$n_upper)*5+10*$ix;
	$detector{"name"}        = "crs_top_"."$thisx"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal top"."$thisx"."$thisy"."$module";
	$detector{"color"}       = "1c86ea";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }
    
    #bottom side
    #lower part 
    for(my $ix = 0; $ix< $n_lower; $ix++){
	my $X = -$outer_pos;
	my $Y = -$n_lower*$AA + $AA+ 2*$ix*$AA; 
	my %detector = init_det();
        my $thisy = -10;
        my $thisx = ($n_core-$n_lower)*5+10*$ix;
	$detector{"name"}        = "crs_bottom_"."$thisx"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal bottom"."$thisx"."$thisy"."$module";
	$detector{"color"}       = "1c86ea"; #"000000";#"00ff00";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }
    #upper part
    for(my $ix = 0; $ix< $n_upper; $ix++){
	my $X = -$L/2- 3*$AA ; 
	my $Y = -$n_upper*$AA + $AA + $ix * 2 * $AA; 
	my %detector = init_det();
        my $thisy = -20;
        my $thisx = ($n_core-$n_upper)*5+10*$ix;
	$detector{"name"}        = "crs_bottom_"."$thisx"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal bottom"."$thisx"."$thisy"."$module";
	$detector{"color"}       = "1c86ea";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }

    #left side
    #lower part 
    for(my $ix = 0; $ix< $n_lower; $ix++){
	my $Y = -$outer_pos ;
	my $X = -$n_lower*$AA + $AA+ 2*$ix*$AA; 
	my %detector = init_det();
        my $thisx = -10;
        my $thisy = ($n_core-$n_lower)*5+10*$ix;
	$detector{"name"}        = "crs_left_"."$thisx"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal left"."$thisx"."$thisy"."$module";
	$detector{"color"}       = "1c86ea"; #"000000";#"0000ff";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
        print_det(\%configuration, \%detector);
    }
    #upper part
    for(my $ix = 0; $ix< $n_upper; $ix++){
	my $Y = -$L/2 - 3*$AA ;  
	my $X = -$n_upper*$AA + $AA + $ix * 2 * $AA; 
	my %detector = init_det();
        my $thisx = -20;
        my $thisy = ($n_core-$n_upper)*5+10*$ix;
	$detector{"name"}        = "crs_left_"."$thisx"."_"."$thisx"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal left"."$ix"."$thisy"."$module";
	$detector{"color"}       = "1c86ea";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }

    #right side 
    #lower part 
    for(my $ix = 0; $ix< $n_lower; $ix++){
	my $Y = $outer_pos;
	my $X = -$n_lower*$AA + $AA+ 2*$ix*$AA; 
	
	my %detector = init_det();
        my $thisx = $n_core*10+0;
        my $thisy = ($n_core-$n_lower)*5+10*$ix;
	$detector{"name"}        = "crs_right_"."$ix"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal right"."$ix"."$thisy"."$module";
	$detector{"color"}       = "1c86ea"; #"000000";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }
    #upper part
    for(my $ix = 0; $ix< $n_upper; $ix++){
	my $Y = $L/2+ 3*$AA ;
	my $X = -$n_upper*$AA + $AA + $ix * 2 * $AA; 
	my %detector = init_det();
        my $thisx = $n_core*10+10;
        my $thisy = ($n_core-$n_upper)*5+10*$ix;
	$detector{"name"}        = "crs_right_"."$thisx"."_"."$thisy"."_"."$module";
	$detector{"mother"}      = "mutest_vessel_air";
	$detector{"description"} = "Panda PbWO4 crystal right"."$ix"."$thisy"."$module";
	$detector{"color"}       = "1c86ea"; #"000000";
	$detector{"style"}       = 0;
	$detector{"visible"}     = 1;
	$detector{"type"}        = "Box";
	$detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
	$detector{"rotation"}    = "0*deg 0*deg 0*deg";
	$detector{"dimensions"}  = "$A*cm $A*cm $l*cm";
	$detector{"material"}    = "G4_PbWO4";
	$detector{"sensitivity"} = "crs";
	$detector{"hit_type"}    = "crs";
	my $man=400+$chcount;
	$chcount = $chcount + 1;
        $detector{"identifiers"} = "sector manual $man xch manual $thisx ych manual $thisy zch manual $module";
	print_det(\%configuration, \%detector);
    }
    if($module == 0){

	print "     %%%%%     \n";
	print "\n";
	
        my $Npipe = 15;
    print " number of PIPEs = $Npipe \n";
	my $n_module = 2.*$chcount;
	my $n_tot = $Npipe*$n_module;

	print "Total number of crystals per matrix = $chcount \n";
	print "Total number of crystals per module = $n_module \n";
	print "Total number of crystals in the full detector = $n_tot \n";
	
       # my  $V_core = 980;
       # my   $V_out  = 980;
	my $V_core = 8.31*8*$a*$a*$l;  #MASS crs in core
	my $V_out = 8.31*8*$A*$A*$l;   #MASS crs in external boundery
    print "Mass single PbWO crsytal = $V_core\n";
	
        my $V_tot = $Npipe*2*($n_core*$n_core*$V_core+4.*$n_lower*$V_out+4.*$n_upper*$V_out);
	
    print "Total mass in the full detector = $V_tot\n";
        
        
      ########### Griglia #############
        
        my $density_Pb = 11.3;
        my $density_W = 19.3;
    my $V_GRID_CRS = $density_W*4*($wrap/2)*2*$a*2*$l;
        my $V_GRID = $n_module*$V_GRID_CRS;
        my $V_GRID_tot = $Npipe*$n_module*$V_GRID_CRS;
        print "Total mass grid in 1 module  = $V_GRID\n";
        print "Total mass grid   = $V_GRID_tot\n";
    
    my $V_BDX = 4.53*800*4.7*5.4*32.5;
    my $frac =$V_tot/$V_BDX;
   
        my $frac2 =($V_tot+$V_GRID_tot)/$V_BDX;
    
    print "Total mass BDX detector = $V_BDX\n";
    
    print "BDX-PIPE corresponds to $frac BDX\n";
    print "BDX-PIPE with grid corresponds to $frac2 BDX\n";
        
	#my $n_2x2 = ($V_tot)/(2.0*2.0*20.0);
	#my $n_2x2_tot = (15.0*$V_tot)/(2.0*2.0*20.0);
	
	#my $n_bdxmini = $n_2x2 / 44.;
	#my $n_bdxmini_tot = $n_2x2_tot/44.;

	#print "Each module corresponds to $n_bdxmini BDX-MINI detectors \n";
	#print "The full detector corresponds to $n_bdxmini_tot detectors \n";

	#my $bdxpipe_frac_bdx = $n_2x2_tot/4334.0;

	#print "The full detector correpsonds to $bdxpipe_frac_bdx BDX (BDX = 125 BDX-MINI) \n";

	print "\n";
	print "     %%%%%     \n";
       
    }
    }





sub make_all
{
    make_main(); 
    make_vessel(); 
    my $r1 = (16*2.54-4); 
    make_veto($r1, 0); 
    make_veto($r1-1.6, 1); 
    make_veto($r1-3.2, 2);
    if($configuration{"variation"} eq "threeveto"){
        make_veto($r1-4.8, 3);
    }
    make_ecal(0);
    make_ecal(1);
}










1;

