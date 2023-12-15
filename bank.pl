use strict;
use warnings;

our %configuration;

# Variable Type is two chars.
# The first char:
#  R for raw integrated variables
#  D for dgt integrated variables
#  S for raw step by step variables
#  M for digitized multi-hit variables
#  V for voltage(time) variables
#
# The second char:
# i for integers
# d for doubles

# The ft banks id are:
#
# Tracker (trk):
# Tracker (hodo):
# Tracker (cal):

sub define_bdx_bank
{
	# uploading the hit definition
	my $bankId   = 100;
	my $bankname = "cormo";
	
	insert_bank_variable(\%configuration, $bankname, "bankid", $bankId, "Di", "$bankname bank ID");
	insert_bank_variable(\%configuration, $bankname, "sector",       1, "Di", "sector number");
	insert_bank_variable(\%configuration, $bankname, "layer",        2, "Di", "layer number");
	insert_bank_variable(\%configuration, $bankname, "paddle",       3, "Di", "paddle number");
	insert_bank_variable(\%configuration, $bankname, "adcl",         4, "Di", "adcl");
	insert_bank_variable(\%configuration, $bankname, "adcr",         5, "Di", "adcr");
	insert_bank_variable(\%configuration, $bankname, "tdcl",         6, "Di", "tdcl");
	insert_bank_variable(\%configuration, $bankname, "tdcr",         7, "Di", "tdcr");
	insert_bank_variable(\%configuration, $bankname, "adcb",         8, "Di", "adcb");
	insert_bank_variable(\%configuration, $bankname, "adcf",         9, "Di", "adcf");
	insert_bank_variable(\%configuration, $bankname, "tdcb",        10, "Di", "tdcb");
	insert_bank_variable(\%configuration, $bankname, "tdcf",        11, "Di", "tdcf");
	insert_bank_variable(\%configuration, $bankname, "hitn",        99, "Di", "hit number");
	
	
	$bankId   = 200;
	$bankname = "veto";
	
	insert_bank_variable(\%configuration, $bankname, "bankid", $bankId, "Di", "$bankname bank ID");
	insert_bank_variable(\%configuration, $bankname, "sector",       1, "Di", "sector number");
	insert_bank_variable(\%configuration, $bankname, "veto",         2, "Di", "veto number");
	insert_bank_variable(\%configuration, $bankname, "channel",      3, "Di", "channel number");
    insert_bank_variable(\%configuration, $bankname, "module",       4, "Di", "module number");
	insert_bank_variable(\%configuration, $bankname, "adc1",         5, "Di", "adc1");
	insert_bank_variable(\%configuration, $bankname, "adc2",         6, "Di", "adc2");
    insert_bank_variable(\%configuration, $bankname, "adc3",         7, "Di", "adc3");
    insert_bank_variable(\%configuration, $bankname, "adc4",         8, "Di", "adc4");
	insert_bank_variable(\%configuration, $bankname, "tdc1",         9, "Di", "tdc1");
	insert_bank_variable(\%configuration, $bankname, "tdc2",        10, "Di", "tdc2");
    insert_bank_variable(\%configuration, $bankname, "tdc3",        11, "Di", "tdc3");
    insert_bank_variable(\%configuration, $bankname, "tdc4",        12, "Di", "tdc4");
    insert_bank_variable(\%configuration, $bankname, "adc5",        13, "Di", "adc5");
    insert_bank_variable(\%configuration, $bankname, "adc6",        14, "Di", "adc6");
    insert_bank_variable(\%configuration, $bankname, "adc7",        15, "Di", "adc7");
    insert_bank_variable(\%configuration, $bankname, "adc8",        16, "Di", "adc8");
    insert_bank_variable(\%configuration, $bankname, "tdc5",        17, "Di", "tdc5");
    insert_bank_variable(\%configuration, $bankname, "tdc6",        18, "Di", "tdc6");
    insert_bank_variable(\%configuration, $bankname, "tdc7",        19, "Di", "tdc7");
    insert_bank_variable(\%configuration, $bankname, "tdc8",        20, "Di", "tdc8");
    insert_bank_variable(\%configuration, $bankname, "hitn",        99, "Di", "hit number");

    $bankId   = 300;
    $bankname = "crs";
    
    insert_bank_variable(\%configuration, $bankname, "bankid", $bankId, "Di", "$bankname bank ID");
    insert_bank_variable(\%configuration, $bankname, "sector",       1, "Di", "sector number");
    insert_bank_variable(\%configuration, $bankname, "xch",          2, "Di", "xch number");
    insert_bank_variable(\%configuration, $bankname, "ych",          3, "Di", "ych number");
    insert_bank_variable(\%configuration, $bankname, "zch",          4, "Di", "zch number");
    insert_bank_variable(\%configuration, $bankname, "adcl",         5, "Di", "adcl");
    insert_bank_variable(\%configuration, $bankname, "adcr",         6, "Di", "adcr");
    insert_bank_variable(\%configuration, $bankname, "tdcl",         7, "Di", "tdcl");
    insert_bank_variable(\%configuration, $bankname, "tdcr",         8, "Di", "tdcr");
    insert_bank_variable(\%configuration, $bankname, "adcb",         9, "Di", "adcb");
    insert_bank_variable(\%configuration, $bankname, "adcf",         10, "Di", "adcf");
    insert_bank_variable(\%configuration, $bankname, "tdcb",        11, "Di", "tdcb");
    insert_bank_variable(\%configuration, $bankname, "tdcf",        12, "Di", "tdcf");
    insert_bank_variable(\%configuration, $bankname, "hitn",        99, "Di", "hit number"); 


}



sub define_banks
{
	define_bdx_bank();
}

1;










