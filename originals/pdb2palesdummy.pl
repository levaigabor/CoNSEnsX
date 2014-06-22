#!/usr/bin/perl

print STDERR "
###################################################
# pdb2palesdummy.pl < pdb_file > pales_dummy_file #
#                                                 #
# Creates a PALES dummy file based on the PDB     #
# structure supplied.                             #
# Trying to make it work on multiple-chain        #
# molecules also.                                 #
#                                                 #
# Version 0.1                         15.11.2004. #
###################################################

";

 %G_SHORTCODES =      ('ALA','A',
                       'ASP','D',
                       'ASN','N',
                       'ARG','R',
                       'CYS','C',
                       'GLY','G',
                       'GLU','E',
                       'GLN','Q',
                       'HIS','H',
                       'ILE','I',
                       'LEU','L',
                       'LYS','K',
                       'MET','M',
                       'PHE','F',
                       'PRO','P',
                       'SER','S',
                       'THR','T',
                       'TRP','W',
                       'TYR','Y',
                       'VAL','V',
                    );


# Reading in PDB file
while(<>){

    chop;
    if ($_ =~ /^ATOM +[0-9]+ +CA  /){
	$atomline=$_;
	$resname=substr($atomline,17,3);
	$resnum=substr($atomline,22,5);$resnum=~s/ //g;
	$chainID=substr($atomline,21,1);
	$chains{$chainID}=1;
	$seq{$chainID}.=" $resname";
    }

}

# Writing sequence
foreach $chainid (sort keys %chains){
    print "DATA SEQUENCE ";
    $seq{$chainid}=~s/^ //;
    $counter=0;$linec=0;
    foreach $aa (split(/ /,$seq{$chainid})){
	print $G_SHORTCODES{$aa};
	$counter++;
	if ($counter == 10){
	    $linec++;
	    if ($linec == 5){
		print "\nDATA SEQUENCE ";
		$linec=0;
	    }
	    else{
		print " ";
	    }
	    $counter=0;
	}
    }
    print "\n";
}
	     
#print "\n";


#Writing dummy dipoles

print"
VARS RESID_I RESNAME_I ATOMNAME_I RESID_J RESNAME_J ATOMNAME_J D DD W
FORMAT %5d  %6s  %6s  %5d  %6s  %6s  %9.3f  %9.3f  %.2f

";

foreach $chainid (sort keys %chains){
    $resnum=1;
    foreach $aa (split(/ /,$seq{$chainid})){
        $resid=sprintf("%4d$chainid",$resnum);
        if ($aa ne "PRO"){
            printf(" $resid  %6s  %6s  $resid  %6s  %6s  %9.3f  %9.3f  %.2f\n",$aa, "H",$aa,"N",0.000,1.000,1.00);
            printf(" $resid  %6s  %6s  $resid  %6s  %6s  %9.3f  %9.3f  %.2f\n",$aa, "N",$aa,"C",0.000,1.000,1.00);
            printf(" $resid  %6s  %6s  $resid  %6s  %6s  %9.3f  %9.3f  %.2f\n",$aa, "C",$aa,"CA",0.000,1.000,1.00);
        }
        $resnum++;
    }


}

