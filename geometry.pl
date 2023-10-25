use strict;
use warnings;

#TO DO: fix alveolus dimension to value measured

our %configuration;

my $degrad = 0.01745329252;
my $cic    = 2.54;

my $crystal = "CsI";
$crystal = $configuration{"crs"};

# spacings
my $cr_mylar = $configuration{"cr_mylar"}/2.0; # Mylar wrapping thikness
my $cr_air = $configuration{"cr_air"} /2.0; # Air gap
my $cr_alv = $configuration{"cr_alv"} /2.0; # Alveolus thikness (0.03cm if Cfiber, 0.2 cm if Al)

# Effective distance between cyrstal blocks
my $crs_gap = $configuration{"crs_gap"};
my $veto_gap = $configuration{"veto_gap"}; # spacing between veto and layer inside

my $nmodules = int($configuration{"modules"}); # ECal+Veto                     #  <----- X Y ||
my $ncol = int($configuration{"columns"}); # Number of columns (X width)      #             ||
my $nrow = int($configuration{"rows"}); # Number of rows (Y width)         #             \/
my $ndep = int($configuration{"depth"}); # Number of ncol x nrow matrixes

my $max_w = $configuration{"paddle_max_width"}/2; #maximum width of paddle scintillator

my $pb_tk = $configuration{"pb_tk"}/2;
my $iv_tk = $configuration{"iv_tk"}/2;
my $ov_tk = $configuration{"ov_tk"}/2;

# Cal center in X, Y, Z
my $X0=0; my $Y0=0; my $Z0=0;

# crystal parameters -- all in cm -- DEFAULT = CSI
my $crs_x=4.7/2 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
my $crs_y=4.8/2 ; # Endcap: short side Y 4.7
my $crs_X=5.8/2 ; # Endcap: long side X (5+4.6)/2=4.8cm
my $crs_Y=6.0/2 ; # Endcap: long side Y 5.4
my $crs_l=31.6/2.; # Endcap: lenght side Y 32.5

my $crs_available = 0;

if($configuration{"crs"}eq "CsI"){
    $crs_available = 800;
    $crs_x=4.7/2 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
    $crs_y=4.8/2 ; # Endcap: short side Y 4.7
    $crs_X=5.8/2 ; # Endcap: long side X (5+4.6)/2=4.8cm
    $crs_Y=6.0/2 ; # Endcap: long side Y 5.4
    $crs_l=31.6/2.; # Endcap: lenght side Y 32.5

}
if($configuration{"crs"}eq "PbWO4"){
    $crs_available = 4000;
    $ncol = $ncol*2;
    $nrow = $nrow*2;
    $ndep = int($ndep*3.0/2.0);
    $crs_x=2.1/2 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
    $crs_y=2.1/2 ; # Endcap: short side Y 4.7
    $crs_X=2.9/2 ; # Endcap: long side X (5+4.6)/2=4.8cm
    $crs_Y=2.9/2 ; # Endcap: long side Y 5.4
    $crs_l=20.0/2.; # Endcap: lenght side Y 32.5
}



my $alv_dz = $crs_l*2.0+$crs_gap;
my $alv_dx = 2.0*($crs_X+$cr_mylar+$cr_alv+$cr_alv);
my $alv_dy = 2.0*($crs_Y+$cr_mylar+$cr_alv+$cr_alv);

if($configuration{"variation"} eq "pbwo_core"){
    $crs_available = 80000;
    # for this configuration the crystal dimension for reference is the size of BaBar crystals
    $crs_x=4.7/2 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
    $crs_y=4.8/2 ; # Endcap: short side Y 4.7
    $crs_X=5.8/2 ; # Endcap: long side X (5+4.6)/2=4.8cm
    $crs_Y=6.0/2 ; # Endcap: long side Y 5.4
    $crs_l=31.6/2.; # Endcap: lenght side Y 32.5
    
    $alv_dz = $crs_l*2.0+$crs_gap;
    $alv_dx = 2.0*($crs_X+$cr_mylar+$cr_alv+$cr_alv);
    $alv_dy = 2.0*($crs_Y+$cr_mylar+$cr_alv+$cr_alv);
    
    # for this configuration the module size is fixed:
    $ncol = 10;
    $nrow = 10;
}

my $core_side = $configuration{"core_side"};

# parameters for vertical crystal configuration (only considered if flag on)
my $nplane = int($configuration{"planes"});
my $plane_side = int($configuration{"plane_side"});
my $plane_depth = int(($plane_side*$alv_dx)/$alv_dz)+1;
my $matrix_side = ($plane_side*$alv_dx/2.0);
if($plane_depth * $alv_dz/2.0 > $matrix_side ){
    $matrix_side = $plane_depth * $alv_dz/2.0;
}

my $pipe_radius = 304.8/2; # cm pipe radius

my $crs_count = 0;

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
    my $par1 = 200;
    my $par2 = 200.;
    my $par3 = 400.;
    $detector{"dimensions"}  = "$par1*cm $par2*cm $par3*cm";
    $detector{"material"}    = "G4_AIR";
    $detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
}

sub make_pipes{
    my %detector = init_det();
    my $X = 0. ;   my $Y = 0. ;
    my $Z = 0.; # depth in the pipe
    my $rotX=90.;  my $rotY=0.;  my $rotZ=0.;
    
    my $vs_lg=100.0/2;
    my $vs_or = $pipe_radius; #outer radius
    my $vs_side_tk = 2.;
    my $vs_ir = $vs_or - $vs_side_tk; #inner radius
    
    my $par4 = 0.;    my $par5 = 360.;
    
    $detector{"mother"}      = "main_volume";
    
    
    $detector{"name"}        = "pipe_1";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "Pipe 1";
    $detector{"color"}       = "A05070";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Tube";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "$rotX*deg $rotY*deg 0*deg";
    $detector{"dimensions"}  = "$vs_ir*cm $vs_or*cm $vs_lg*cm $par4*deg $par5*deg";
    $detector{"material"}    = "G4_Fe";
    print_det(\%configuration, \%detector);
    
    
}

open(my $traslation_table, '>', 'crs_traslation_table.txt') or die "Can't open output file, sorry";

sub make_crystal{
    my $im = $_[0];
    my $ix = $_[1];
    my $iy = $_[2];
    my $iz = $_[3];
    
    my $rotx = 0;
    my $roty = 0;
    my $rotz = 0;
    
    my $rot=0;
    
    if($configuration{"variation"} eq "pbwo_core"){
        my $start = ($nrow-$core_side)/2;
        my $stop = ($nrow-$core_side)/2+$core_side;
        
        if($ix>$start-1 && $ix<$stop && $iy>$start-1 && $iy<$stop){
            # for this configuration I assume lead crystals to be the merge of 4
            $crs_x=2.1 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
            $crs_y=2.1 ; # Endcap: short side Y 4.7
            $crs_X=2.9 ; # Endcap: long side X (5+4.6)/2=4.8cm
            $crs_Y=2.9 ; # Endcap: long side Y 5.4
            $crs_l=20.0/2.; # Endcap: lenght side Y 32.5
            $crs_gap = (2*(31.6+$configuration{"crs_gap"})-3.0*2.0*$crs_l)/3.0;
        }else{
            $crs_x=4.7/2 ; # Endcap: short side X (4.3+3.9)/2=4.1cm
            $crs_y=4.8/2 ; # Endcap: short side Y 4.7
            $crs_X=5.8/2 ; # Endcap: long side X (5+4.6)/2=4.8cm
            $crs_Y=6.0/2 ; # Endcap: long side Y 5.4
            $crs_l=31.6/2.; # Endcap: lenght side Y 32.5
            $crs_gap = $configuration{"crs_gap"};
            
        }
        
        $alv_dz = $crs_l*2.0+$crs_gap;
    }
    
    $crs_count = $crs_count+1;
    
    # All parallelepiped are defined as
    #     _ _________
    #  y |_|_________|
    #     x     l
    
    # For crystals the capital letter refers to the larger value of the face side
    
    # Wrapped crystals
    my $wr_x = $crs_x+$cr_mylar;
    my $wr_X = $crs_X+$cr_mylar;
    my $wr_y = $crs_y+$cr_mylar;
    my $wr_Y = $crs_Y+$cr_mylar;
    my $wr_l = $crs_l+$cr_mylar;
    
    # Air around crystals
    my $air_X = $wr_X+$cr_air;
    my $air_Y = $wr_Y+$cr_air;
    my $air_l = $wr_l+$cr_air;
    
    # Crystal alveolus
    my $al_X = $air_X+$cr_alv;
    my $al_Y = $air_Y+$cr_alv;
    my $al_l = $air_l+$cr_alv;
    
    if($crs_count<$crs_available+1){
        my %detector = init_det();
        
        $detector{"name"}        = "cry_alveol_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"mother"}      = "module_$im"."_crs_mother";
        $detector{"description"} = "Carbon/Al container_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        # pos =  center of module - half size of block + alv. side
        my $X = -($ncol-1)*$al_X + 2*$al_X*$ix;
        my $Y = -($nrow-1)*$al_Y + 2*$al_Y*$iy;
        my $Z = $iz * $alv_dz- 0.5*$alv_dz*($ndep-1) ;
                
        
        
        if($configuration{"vertical_crystals"} eq 1){   #in this case ix is the number of the colunm/row, iy the position along the row and iz the plane
            if($iz % 2 == 0){
                $rotx = 90;
                $X = -($plane_side-1)*$al_X+ 2*$al_X*$ix;
                $Y = -($plane_depth-1)*0.5*$alv_dz + $alv_dz*$iy;
                $Z = -0.5*$alv_dy*($nplane-1) + $iz*$alv_dy;
            }elsif($iz % 2 != 0){
                $roty = 90;
                $X = -($plane_depth-1)*0.5*$alv_dz + $alv_dz*$iy;
                $Y = -($plane_side-1)*$al_X+ 2*$al_X*$ix;
                $Z = -0.5*$alv_dy*($nplane-1) + $iz*$alv_dy;
                print("$X $Y $Z $plane_depth \n");
                
            }    }
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "$rotx*deg $roty*deg $rotz*deg";
        $detector{"dimensions"}  = "$al_X*cm $al_Y*cm $al_l*cm";
        $detector{"material"}    = "G4_Al";
        print_det(\%configuration, \%detector);
        
        # Air layer
        %detector = init_det();
        $detector{"name"}        = "cry_air_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"mother"}      = "cry_alveol_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"description"} = "Air $ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"color"}       = "00fff1";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$air_X*cm $air_Y*cm $air_l*cm";
        $detector{"material"}    = "G4_AIR";
        $detector{"material"}   = "KryptoniteLight";
        print_det(\%configuration, \%detector);
        
        # Mylar wrapping
        %detector = init_det();
        $detector{"name"}        = "cry_mylar_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"mother"}      = "cry_air_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"description"} = "Mylar wrapping_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"color"}       = "00fff2";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$wr_x*cm $wr_X*cm $wr_y*cm $wr_Y*cm $wr_l*cm";
        $detector{"material"}    = "bdx_mylar";
        print_det(\%configuration, \%detector);
        
        # Crystals
        %detector = init_det();
        $detector{"name"}        = "crystal_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"mother"}      = "cry_mylar_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"description"} = "Crystal_$ix"."_"."$iy"."_"."$iz"."_"."$im";
        $detector{"color"}       = "00ffff";
        $detector{"style"}       = 1;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Trd";
        $detector{"pos"}         = "0*cm 0*cm 0*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$crs_x*cm $crs_X*cm $crs_y*cm $crs_Y*cm $crs_l*cm";
        $detector{"material"}    = "CsI_Tl";
        if($configuration{"crs"}eq "PbWO4"){$detector{"material"} = "G4_PbWO4";}
        if($configuration{"variation"} eq "pbwo_core"){
            my $start = ($nrow-$core_side)/2;
            my $stop = ($nrow-$core_side)/2+$core_side;
            if($ix>$start-1 && $ix<$stop && $iy>$start-1 && $iy<$stop){
                $detector{"color"}       = "ff00ff";
                $detector{"material"} = "G4_PbWO4";
            }else{
                $detector{"color"}       = "00ffff";
                $detector{"material"}    = "CsI_Tl";
            }
        }
        $detector{"sensitivity"} = "crs";
        $detector{"hit_type"}    = "crs";
        
        my $i_im = $im;
        my $i_ix = $ix+1;
        my $i_iy = $iy+1;
        my $i_iz = $iz+1;
        $detector{"identifiers"} = "sector manual $i_im xch manual $i_ix ych manual $i_iy zch manual $i_iz";
        #$detector{"identifiers"} = "sector manual $i_im xch manual $X ych manual $Y zch manual $Z";
        print_det(\%configuration, \%detector);
        
        
        print $traslation_table "$i_im $i_ix $i_iy $i_iz $X $Y $Z \n";
    }else{
        print("Exceeding number of crystal available");
        $crs_count = $crs_available;
    }
    
    
}

# BEGIN Porposal INNER LEAD
sub make_lead
{
    # Assuming fixed parameters for lead
    # Dimensions defined as follows
    #     _ _________
    #  y |_|_________|
    #     x     z
    
    my $im = $_[0];
    
    # Front face
    my $pb_f_x = $ncol*$alv_dx/2.0;
    my $pb_f_y = $nrow*$alv_dy/2.0;
    my $pb_f_z = $pb_tk;
    # Side face
    my $pb_s_x = $pb_tk;
    my $pb_s_y = $nrow*$alv_dy/2.0;
    my $pb_s_z = $ndep*$alv_dz/2.0+2*$pb_tk;
    # Top face
    my $pb_t_x = $ncol*$alv_dx/2.0+2.0*$pb_tk;
    my $pb_t_y = $pb_tk;
    my $pb_t_z = $ndep*$alv_dz/2.0+2.0*$pb_tk;
    
    if($configuration{"vertical_crystals"} eq 1){
       
        
        print("Inside lead rotated routine\n");
       
        $pb_f_x = $matrix_side;
        $pb_f_y = $matrix_side;
        
        $pb_s_y = $matrix_side;
        $pb_s_z = $nplane*$alv_dy/2 + 2.0 *  $pb_tk;
        print($pb_t_z);
        $pb_t_x = $matrix_side+2.0*$pb_tk;
        $pb_t_z = $nplane*$alv_dy/2 + 2.0 * $pb_tk;
        print(" $pb_t_z \n");
    }
    
    my %detector = init_det();
    
    $detector{"mother"}      = "module_$im"."_mother";
    $detector{"material"}    = "G4_Pb";
    
    # Top/bottom
    $detector{"name"}        = "lead_top_$im";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $Y = ($nrow)*$alv_dy/2.0+$pb_tk;  
    if($configuration{"vertical_crystals"} eq 1){
        $Y = $matrix_side+$pb_tk;
    }
    $detector{"pos"}         = "0*cm $Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_t_x*cm $pb_t_y*cm $pb_t_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_in_bottom_$im";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Y = -$Y;
    $detector{"pos"}         = "0*cm $Y*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_t_x*cm $pb_t_y*cm $pb_t_z*cm";
    print_det(\%configuration, \%detector);
    
    # Front/back
    $detector{"name"}        = "lead_in_upstream_$im";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $Z = ($ndep)*$alv_dz/2.0+$pb_tk;
    if($configuration{"vertical_crystals"} eq 1){
        $Z = ($nplane)*$alv_dy/2.0+$pb_tk;
    }
    $detector{"pos"}         = "0*cm 0*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_f_x*cm  $pb_f_y*cm $pb_f_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_in_downstream_$im";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $Z = -$Z;
    $detector{"pos"}         = "0*cm 0*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_f_x*cm  $pb_f_y*cm $pb_f_z*cm";
    print_det(\%configuration, \%detector);
   
    $detector{"name"}        = "lead_in_right_$im";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 1;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    my $X = ($ncol)*$alv_dx/2.0+$pb_tk;  
    if($configuration{"vertical_crystals"} eq 1){
        $X = $matrix_side+$pb_tk;
    }
    $detector{"pos"}         = "$X*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_s_x*cm  $pb_s_y*cm $pb_s_z*cm";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "lead_in_left_$im";
    $detector{"description"} = "lead shield";
    $detector{"color"}       = "A9D0F5";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $X = -$X;
    $detector{"pos"}         = "$X*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$pb_s_x*cm  $pb_s_y*cm $pb_s_z*cm";
    print_det(\%configuration, \%detector);
}

sub make_iveto{
    my %detector = init_det();
    my $im = $_[0];
    
    print("Making iveto for module $im \n");
    $detector{"mother"}      = "module_$im"."_mother";
    
    #dimension of the volume inside the veto
    my $pb_x = $ncol * $alv_dx/2+2*$pb_tk+$veto_gap; # pb dimension + 1cm
    my $pb_y = $nrow * $alv_dy/2+2*$pb_tk+$veto_gap;
    my $pb_z = $ndep * $alv_dz/2+2*$pb_tk+$veto_gap;
    if($configuration{"vertical_crystals"} eq 1){
        $pb_x =  $matrix_side+2*$pb_tk+$veto_gap;
        $pb_y = $matrix_side+2*$pb_tk+$veto_gap;
        $pb_z = $nplane*$alv_dy/2+2*$pb_tk+$veto_gap;
    }

    # Front face
    my $n_f = int($pb_x/$max_w)+1;
    my $iv_f_x = $pb_x/$n_f;
    my $iv_f_y = $pb_y;
    my $iv_f_z = $iv_tk;
    # Side face
    my $iv_s_x = $iv_tk;
    my $iv_s_y = $pb_y;
    my $n_s = int(($pb_z+2.0*$iv_tk)/$max_w)+1;
    my $iv_s_z = ($pb_z+2.0*$iv_tk)/$n_s;
    # Top face
    my $iv_t_x = $pb_x+2.0*$iv_tk;
    my $iv_t_y = $iv_tk;
    my $n_t = int(($pb_z+2.0*$iv_tk)/$max_w)+1;
    my $iv_t_z = ($pb_z+2.0*$iv_tk)/$n_t;
    
   
    
    print("Front: max w = $max_w - dimension $pb_x - number of paddles $n_f - side of paddles $iv_f_x \n");
    print("Side : max w = $max_w - dimension $pb_z - number of paddles $n_s - side of paddles $iv_s_z \n");
    print("Top  : max w = $max_w - dimension $pb_z - number of paddles $n_t - side of paddles $iv_t_z \n");
    
    for(my $if = 0; $if<$n_f; $if++){
        $detector{"name"}        = "iveto_front_$if"."module_$im";
        $detector{"description"} = "inner veto front";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = -$pb_x+$iv_f_x*(1.0+2.0*$if);#$cal_centx;
        my $Y = 0;
        my $Z = -($pb_z+$iv_tk);
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_f_x*cm $iv_f_y*cm $iv_f_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        my $ch_id = 1000*$im+000+00+$if;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 3";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "iveto_back_$if"."module_$im";
        $detector{"description"} = "inner veto back";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $Z = -$Z;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_f_x*cm $iv_f_y*cm $iv_f_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $ch_id = 1000*$im+000+50+$if;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 4";
        print_det(\%configuration, \%detector);
    }
    
    for(my $is = 0; $is<$n_s; $is++){
        $detector{"name"}        = "iveto_left_$is"."module_$im";
        $detector{"description"} = "inner veto left";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = -$pb_x-$iv_tk;
        my $Y = 0;
        my $Z = -($pb_z+2.0*$iv_tk)+$iv_s_z*(1.0+2.0*$is);
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_s_x*cm $iv_s_y*cm $iv_s_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        my $ch_id = 1000*$im+100+00+$is;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 6";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "iveto_right_$is"."module_$im";
        $detector{"description"} = "inner veto right";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = -$X;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_s_x*cm $iv_s_y*cm $iv_s_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $ch_id = 1000*$im+100+50+$is;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 5";
        print_det(\%configuration, \%detector);
    }
    
    for(my $it = 0; $it<$n_t; $it++){
        $detector{"name"}        = "iveto_top_$it"."module_$im";
        $detector{"description"} = "inner veto top";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = 0;
        my $Y = -($pb_y+$iv_tk);
        my $Z = -($pb_z+2.0*$iv_tk)+$iv_t_z*(1.0+2.0*$it);
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_t_x*cm $iv_t_y*cm $iv_t_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        my $ch_id = 1000*$im+200+00+$it;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 1";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "iveto_bottom_$it"."module_$im";
        $detector{"description"} = "inner veto bottom";
        $detector{"color"}       = "0000FF";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $Y = -$Y;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$iv_t_x*cm $iv_t_y*cm $iv_t_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $ch_id = 1000*$im+200+50+$it;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 4 channel manual 2";
        print_det(\%configuration, \%detector);
    }
    
    
    
}

sub make_oveto{
    my %detector = init_det();
    my $im = $_[0];
    print("Making oveto for module $im \n");
    $detector{"mother"}      = "module_$im"."_mother";;
    
    my $pb_x = $ncol * $alv_dx/2+2*$pb_tk+$veto_gap; # pb dimension + 1cm
    my $pb_y = $nrow * $alv_dy/2+2*$pb_tk+$veto_gap;
    my $pb_z = $ndep * $alv_dz/2+2*$pb_tk+$veto_gap;
    if($configuration{"vertical_crystals"} eq 1){
        $pb_x =  $matrix_side+2*$pb_tk+$veto_gap;
        $pb_y = $matrix_side+2*$pb_tk+$veto_gap;
        $pb_z = $nplane*$alv_dy/2+2*$pb_tk+$veto_gap;
    }
    
    #dimension of the volume inside the veto
    my $iv_x = $pb_x+2.0*$iv_tk+$veto_gap; # pb dimension + 1cm
    my $iv_y = $pb_y+2.0*$iv_tk+$veto_gap;
    my $iv_z = $pb_z+2.0*$iv_tk+$veto_gap;
    
    
    # Front face
    my $n_f = int($iv_x/$max_w)+1;
    my $ov_f_x = $iv_x/$n_f;
    my $ov_f_y = $iv_y;
    my $ov_f_z = $ov_tk;
    # Side face
    my $ov_s_x = $ov_tk;
    my $ov_s_y = $iv_y;
    my $n_s = int(($iv_z+2.0*$ov_tk)/$max_w)+1;
    my $ov_s_z = ($iv_z+2.0*$ov_tk)/$n_s;
    # Top face
    my $ov_t_x = $iv_x+2.0*$ov_tk;
    my $ov_t_y = $ov_tk;
    my $n_t = int(($iv_z+2.0*$ov_tk)/$max_w)+1;
    my $ov_t_z = ($iv_z+2.0*$ov_tk)/$n_t;
    
   
    
    print("Front: max w = $max_w - dimension $iv_x - number of paddles $n_f - side of paddles $ov_f_x \n");
    print("Side : max w = $max_w - dimension $iv_z - number of paddles $n_s - side of paddles $ov_s_z \n");
    print("Top  : max w = $max_w - dimension $iv_z - number of paddles $n_t - side of paddles $ov_t_z \n");
    
    for(my $if = 0; $if<$n_f; $if++){
        $detector{"name"}        = "oveto_front_$if"."module_$im";
        $detector{"description"} = "outer veto front";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = -$iv_x+$ov_f_x*(1.0+2.0*$if);#$cal_centx;
        my $Y = 0;
        my $Z = -($iv_z+$ov_tk);
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_f_x*cm $ov_f_y*cm $ov_f_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        my $ch_id = 1000*$im+500+00+$if;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 5 channel manual 3";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "oveto_back_$if"."module_$im";
        $detector{"description"} = "outer veto back";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $Z = -$Z;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_f_x*cm $ov_f_y*cm $ov_f_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $ch_id = 1000*$im+500+50+$if;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 5 channel manual 4";
        print_det(\%configuration, \%detector);
    }
    
    for(my $is = 0; $is<$n_s; $is++){
        $detector{"name"}        = "oveto_left_$is"."module_$im";
        $detector{"description"} = "outer veto left";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = -$iv_x-$ov_tk;
        my $Y = 0;
        my $Z = -($iv_z+2.0*$ov_tk)+$ov_s_z*(1.0+2.0*$is);
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_s_x*cm $ov_s_y*cm $ov_s_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        my $ch_id = 1000*$im+600+00+$is;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 5 channel manual 6";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "oveto_right_$is"."module_$im";
        $detector{"description"} = "outer veto right";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $X = -$X;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_s_x*cm $ov_s_y*cm $ov_s_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $ch_id = 1000*$im+600+50+$is;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 5 channel manual 5";
        print_det(\%configuration, \%detector);
    }
    
    for(my $it = 0; $it<$n_t; $it++){
        $detector{"name"}        = "oveto_top_$it"."module_$im";
        $detector{"description"} = "outer veto top";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        my $X = 0;
        my $Y = -($iv_y+$ov_tk);
        my $Z = -($iv_z+2.0*$ov_tk)+$ov_t_z*(1.0+2.0*$it);
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_t_x*cm $ov_t_y*cm $ov_t_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        my $ch_id = 1000*$im+700+00+$it;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 5 channel manual 1";
        print_det(\%configuration, \%detector);
        
        $detector{"name"}        = "oveto_bottom_$it"."module_$im";
        $detector{"description"} = "outer veto bottom";
        $detector{"color"}       = "088A4B";
        $detector{"style"}       = 0;
        $detector{"visible"}     = 1;
        $detector{"type"}        = "Box";
        $Y = -$Y;
        $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
        $detector{"rotation"}    = "0*deg 0*deg 0*deg";
        $detector{"dimensions"}  = "$ov_t_x*cm $ov_t_y*cm $ov_t_z*cm ";
        $detector{"material"}    = "ScintillatorB";
        $detector{"sensitivity"} = "veto";
        $detector{"hit_type"}    = "veto";
        $ch_id = 1000*$im+700+50+$it;
        $detector{"identifiers"} = "sector manual $ch_id veto manual 5 channel manual 2";
        print_det(\%configuration, \%detector);
    }
    
    
    
}

open(my $output, '>', 'geom.txt') or die "Can't open output file, sorry";

sub make_module
{    
    my $im = $_[0];
    
    print("making module $im \n");
        
    
    
    my $calorimeter_X = $ncol * $alv_dx/2;
    my $calorimeter_Y = $nrow * $alv_dy/2;
    my $calorimeter_Z = $ndep * $alv_dz/2;
    if($configuration{"vertical_crystals"} eq 1){
        $calorimeter_X = $matrix_side;
        $calorimeter_Y = $matrix_side;
        $calorimeter_Z = $nplane*$alv_dy/2;
    }
    
    # effective thickness of veto
    my $lead_thickness = 2*$pb_tk;
    my $iv_thickness = $veto_gap+2.0*$iv_tk;
    my $ov_thickness = $veto_gap+2.0*$ov_tk;
    
    
    my $module_X = $calorimeter_X + $lead_thickness + $iv_thickness + $ov_thickness;
    my $module_Y = $calorimeter_Y + $lead_thickness + $iv_thickness + $ov_thickness;
    my $module_l = $calorimeter_Z + $lead_thickness + $iv_thickness + $ov_thickness;
    
    my $X = $X0;
    my $Y = $Y0;
    my $Z = $Z0+$im*2*$module_l-$module_l*($nmodules-1);
    
    
    print $output "$X $Y $Z \n";
    print $output "$calorimeter_X $calorimeter_Y $calorimeter_Z \n";
    print $output "$lead_thickness \n";
    print $output "$iv_thickness \n";
    print $output "$ov_thickness \n";
    
    
    my %detector = init_det();
    $detector{"name"}        = "module_$im"."_mother";
    $detector{"mother"}      = "main_volume";
    $detector{"description"} = "module_$im"."_mother";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 1;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "$X*cm $Y*cm $Z*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$module_X*cm $module_Y*cm $module_l*cm";
    $detector{"material"}    = "G4_AIR";
    $detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    $detector{"name"}        = "module_$im"."_crs_mother";
    $detector{"mother"}      = "module_$im"."_mother";
    $detector{"description"} = "module_$im"."_crs_mother";
    $detector{"color"}       = "ffffff";
    $detector{"style"}       = 0;
    $detector{"visible"}     = 0;
    $detector{"type"}        = "Box";
    $detector{"pos"}         = "0*cm 0*cm 0*cm";
    $detector{"rotation"}    = "0*deg 0*deg 0*deg";
    $detector{"dimensions"}  = "$module_X*cm $module_Y*cm $module_l*cm";
    $detector{"material"}    = "G4_AIR";
    $detector{"material"}    = "KryptoniteLight";
    print_det(\%configuration, \%detector);
    
    
    if($configuration{"vertical_crystals"} eq 1){
        for(my $iz=0; $iz< $nplane; $iz++){
            for(my $iy=0; $iy<$plane_depth; $iy++){
            for(my $ix = 0; $ix< $plane_side; $ix++){
                    make_crystal($im, $ix, $iy, $iz, 0, 0, 0);
                }
            }
        }
    }elsif($configuration{"variation"} eq "pbwo_core"){
        my $start = ($nrow-$core_side)/2;
        my $stop = ($nrow-$core_side)/2+$core_side;
            for(my $ix = 0; $ix< $ncol; $ix++){
                for(my $iy=0; $iy<$nrow; $iy++){
                    if($ix>$start-1 && $ix<$stop && $iy>$start-1 && $iy<$stop){
                        $ndep = int($configuration{"dep_csi"}*3.0/2.0);;
                    }else{
                        $ndep = $configuration{"dep_csi"};
                    }
                    for(my $iz=0; $iz< $ndep; $iz++){
                    make_crystal($im, $ix, $iy, $iz, 0, 0, 0);
                }
            }
        }
    }else{
        
        for(my $iz=0; $iz< $ndep; $iz++){
            for(my $ix = 0; $ix< $ncol; $ix++){
                for(my $iy=0; $iy<$nrow; $iy++){
                    make_crystal($im, $ix, $iy, $iz, 0, 0, 0);
                }
            }
        }
    }
    make_lead($im);
    make_iveto($im);
    make_oveto($im);
}



sub make_cry_module
{
    print("making module \n");
    for(my $im=0; $im<$nmodules; $im++){
        make_module($im);
    }
    
    print("Crystal distanca: \n");
    print("X: $alv_dx \n");
    print("Y: $alv_dy \n");
    print("Z: $alv_dz \n");
    #my $im = 0;
    #for(my $iz=0; $iz< $ndep; $iz++){
    #    for(my $ix = 0; $ix< $ncol; $ix++){
    #        for(my $iy=0; $iy<$nrow; $iy++){
    #            make_crystal($im, $ix, $iy, $iz, 0, 0, 0);
    #        }
    #    }
    #}
}




sub make_bdx{
    make_cry_module;
    #make_iveto;
    #make_in_lead;
    #make_oveto;
}

sub make_all
{
    make_main;
    make_pipes;
    make_cry_module;
    #if($configuration{"variation"} eq "proposal"){
    #    make_bdx;
    #}
    print("Number of crystals used: $crs_count \n");
}


1;
