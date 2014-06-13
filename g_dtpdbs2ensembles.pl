#!/usr/bin/perl -I/home/szpari/programs/Perl/

print STDERR "
                    :-] G R O M A C S [-:
                        helpful tools

                 :-] g_dtpdbs2ensembles.pl  [-:

Option     Filename  Type          Description
------------------------------------------------------------
";


use Getopt::Std;

$opt_f="*_dt*pdb";
$opt_o="*_t#.pdb";
$opt_t="-1";
$opt_p=0;

getopts('f:o:t:p:hv') || die "ERROR in options, exiting\n";


if ($opt_h){$HO="yes"} else{$HO="no"}
if($opt_v){$V=1;$VO="yes"}else{$VO="no"}

printf STDERR ("  -f  %12s   Input         File mask for input pdbs\n",$opt_f);
printf STDERR ("  -o  %12s   Output        File mask for output pdbs\n",$opt_o);


print STDERR"
      Option   Type  Value  Description
------------------------------------------------------
";

printf STDERR ("          -h   bool %9s  Print help info and quit\n",$HO);
printf STDERR ("          -v   bool %9s  Be loud and noisy\n",$VO);
printf STDERR ("          -p   int  %9s  Population size for sub-ensembles (0=no sub-ensembles)\n",$opt_p
);
printf STDERR ("          -t   int  %9s  Only write ensemble for frame at this time\n",$opt_t);

if ($opt_h) {
die"


Writes concatenated PDB files from those created from individual trajectories at specific time intervals
(e.g. with getstructuresatdt.sh).

Version 1.5 with population calculation support (option -p)

";

}


@infiles=`ls $opt_f`;
#foreach (@infiles){print STDERR  "I -> $_\n";}


$filenum=0;
foreach $inpdb (@infiles){

 $filenum++;
 chop ($inpdb);
 print STDERR "Reading file $inpdb...\n";
 $replicaid=$inpdb;
 $replicaid=~s/^.*[^0-9]([0-9]+)\_dt.*$/$1/;
 if ($opt_p > 0){
    $popid=$replicaid/$opt_p;
    $popid=~s/\..*$//;
    print STDERR "replica $replicaid -> pop $popid\n";
 }
 if ($filenum == 1){
  $fileroot=$inpdb;
  $fileroot=~s/_dt[0-9].*$//;
 }


 open (IP, "$inpdb");
 while(<IP>){
     $time=-1;
     chop;
     $_=~s/\r//g;
     #print STDERR "$inpdb -> r: $_\n";
     if ($_=~/^TITLE/){
        if ($opened){
           print OF "ENDMDL\n";
           close (OF);
           $opened=0;
        }
        $time=$_;
        $time=~s/^.* t= +([0-9\.]+) *$/$1/;
        $time=~s/\..*$//;
        print STDERR "File $inpdb, time $time\n" if $V;
        if (($opt_t == $time) || ($opt_t < 0)){
            $outfile=$fileroot.$opt_o;
            $outfile=~s/\*//g;
            $outfile=~s/\#/$time/;

            $modelnum=$filenum;

            if ($opt_p > 0){
              $outfile=~s/\_t/\_pop$popid\_t/;
              $modelnum-=$popid*$opt_p;
            }

            if (!$outfiles{$outfile}){
               print STDERR "File 1! ($inpdb)\n";
               open(OF,">$outfile") || die "Cannot create $outfile\n";
               print OF "HEADER    ENSEMBLE CREATED FROM FILES $opt_f\n";
               print OF "REMARK    ENSEMBLE AT TIME $time\n"; #if ($opt_t > 0);
               print OF "MODEL     $modelnum\n";
               $opened=1;
            }
            else{
               open(OF,">>$outfile") || die "Cannot append to $outfile\n";
               print OF "MODEL     $modelnum\n";
               $opened=1;
            }
            $outfiles{$outfile}=1;

        }
     }
     elsif ($_ =~/^ATOM/){
       print OF "$_\n" if $opened;
     }

     elsif ($_ =~/^END/){
       print OF "ENDMDL\n" if $opened;
     }


 }#while IF
 close(IF);
 close(OF);
 $opened=0;

}#_foreach inpdb


# Closing files
#foreach $outfile (keys %outfiles){
#   system "echo \"ENDMDL\nEND\n\" >> $outfile";
#}

