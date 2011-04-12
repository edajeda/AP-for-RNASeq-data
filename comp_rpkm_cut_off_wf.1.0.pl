#!/usr/bin/perl -w

use strict;
use warnings;

#Script for comparing rpkm values, and removing the ones where all of the replicates does not surpass the cut-off.Outputs a file with ids in 1st column and pass/fail (1/0) per sample in subsequent columns and a file with a list id ids where all replicates passed in both conditions.
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
comp_rpkm_cut_off_wf.1.0.pl  -i [infile...n] -icf [infile...n] -s [sample ID...n] -o [outfile]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s) containing rpkm-values, comma sep

-icf/--infile Infile(s) containing cut-off value, comma sep

-s/--sampleid The sample ID(s)

-nrrep1/--nrrep1 number of replicates one (Must be supplied to run)

-nrrep2/--nrrep1 number of replicates two (Must be supplied to run)

-oall/--outfile Output file for all samples with ID and passed/failed cut-off(Defaults to comp_rpkm_cut_off_pfcf_wf.txt)

-opid/--outfile Output file for all samples ID that passed (Defaults to comp_rpkm_cut_off_passcf_wf.txt)

=head3 I/O

Input format

1. rpkmforgenes.txt

2. rpkmforgenes_cut-off.txt

Output format:

1. comp_rpkm_cut_off_pfcf_wf.txt

2. comp_rpkm_cut_off_passcf_wf.txt

=cut
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
        qq{comp_rpkm_cut_off_wf.1.0.pl  -i [infile...n] -s [sample ID...n] -o [outfile.txt]
               -i/--infile Infile(s), comma sep
	       -icf/--infile Infile(s) containing cut-off value, comma sep
	       -s/--sampleid The sample ID(s)
	       -nrrep1/--nrrep1 number of replicates one (Must be supplied to run)
	       -nrrep2/--nrrep1 number of replicates two (Must be supplied to run)
               -oall/--outfile Output file for all samples with ID and passed/failed cut-off(Defaults to comp_rpkm_cut_off_pfcf_wf.txt)
	       -opid/--outfile Output file for all samples ID that passed (Defaults to comp_rpkm_cut_off_wf.txt)
	   };
    
}

my ($oall, $opid, $nrrep1, $nrrep2, $help) = ("comp_rpkm_cut_off_pfcf_wf.txt", "comp_rpkm_cut_off_passcf_wf.txt", 0, 0); #Arguments for program

my (@infn,@sid, @icf, @allP);
my (%allRPKM,%allRC, %allCF, %allPF);

GetOptions('i|infile:s'  => \@infn, #Comma separeted list
	   'icf|infilecf:s'  => \@icf, #Comma separeted list of infiles containing cut-off
	   's|sampleid:s'  => \@sid, #Comma separeted list, one below outdirdata
	   'nrrep1|nrrepone:n'  => \$nrrep1, #nr of replicates 1
	   'nrrep2|nrreptwo:n'  => \$nrrep2, #nr of replicates 2
           'oall|outfileall:s'  => \$oall, #Outfile ID with pass fail per replicate
	   'opid|outfilepassid:s'  => \$opid, #Outfile for ID where all replicates passed (both conditions)
           'h|help' => \$help,
           );

die $USAGE if( $help );

if (@infn == 0) {
    my $verbosity = 2;
 print"\n";
    pod2usage({-message => "Must supply an infile(s) as comma separeted list.\n",
	       -verbose => $verbosity
	       });
}
if (@icf == 0) {
    my $verbosity = 2;
 print"\n";
    pod2usage({-message => "Must supply an infile(s) with cut-off as comma separeted list.\n",
	       -verbose => $verbosity
	       });
}
if ($nrrep1 eq 0 && $nrrep2 eq 0) {
    my $verbosity = 2;
 print"\n";
    pod2usage({-message => "Must supply at least one nr of replicates.\n",
	       -verbose => $verbosity
	       });
}

@infn = split(/,/,join(',',@infn)); #Enables comma separated infiles(s),rpkm
@sid = split(/,/,join(',',@sid)); #Enables comma separeted list of sample IDs
@icf = split(/,/,join(',',@icf)); #Enables comma separated infiles(s), cut-off

for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {

    Readrpkm($infn[$inputfiles], $inputfiles);    	
}

for (my $inputfiles=0;$inputfiles<scalar(@icf);$inputfiles++) {
    
    Readcutoff($icf[$inputfiles], $inputfiles);	
}

Compare();

WriteAllPF($oall);

WritePID($opid);

sub Compare {
    
my $passcounter=0;
my $failcounter=0;

    for my $id ( keys %{ $allRPKM{0} } ) { #IDs for first inputfile	
	$passcounter=0;
	$failcounter=0;
	for my $inputfile ( keys %allRPKM ) { #For all input files per id
	    
	    if ($allRPKM{$inputfile}{$id}[0] > 3*$allCF{$inputfile}) { #3x the cut-off
		
		$allPF{$id}[$inputfile] = 1; #1 means above cut-off
		$passcounter++; #Sum up all passed per ID and inputfile
	    }
	    else {
		$allPF{$id}[$inputfile] = 0; #0 means below cut-off
		$failcounter++;
	    }
	}
	if ($passcounter eq scalar(keys %allRPKM) ) { #If at the end and all replicates passed
	    
	    push(@allP, $id);
	    
	}
	
    }
}

sub WriteAllPF {

    my $filename = shift;
    
    open (ALLPF, ">$filename") or die "Can't write to $filename: $!\n";
    
    print ALLPF "ID\t";
    
    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {  
	print ALLPF $sid[$sampleid], "\t";
    }
    
    print ALLPF "\n";
    
    for my $id ( keys %allPF ) { #IDs
	
	print ALLPF $id,"\t";
	
	for (my $pf=0;$pf <scalar( @{$allPF{$id} } );$pf++) { #All entries below/above cut-off
	    
	    print ALLPF $allPF{$id}[$pf], "\t";
	    
	}
	print ALLPF "\n";
    }
    close(ALLPF);
    return;
}

sub WritePID {
    
    my $filename = shift;
    
    open (PID, ">$filename") or die "Can't write to $filename: $!\n";

    print PID "ID\n";

    for (my $pf=0;$pf <scalar( @allP );$pf++) { #All entries above cut-off
	
	print PID $allP[$pf], "\n";
	
    }
    close(PID);
    return;
}

sub Readrpkm {
    
#$_[0] = filename
#$_[1] = Nr of inputfile
    
    open(RRPKM, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<RRPKM>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
        }
	if (/#/) { #Avoid headers
	    next;
	}
	if (/(\S+)/) {
	    
	    my @temp = split("\t",$_);	    #Loads tab separated info per line
	    shift(@temp); #Remove chr position
	    my $id = shift(@temp); # id string
	    $allRPKM{$_[1]}{$id} = [@temp]; # Hash{inputfile} of Hash{ID} of array[rpkm, Read_Counts]
	}
    } 	
    close(RRPKM);
    print STDERR "Finished Reading infile $_[0]","\n";
    return;
}

sub Readcutoff {
    
#$_[0] = filename
#$_[1] = Nr of inputfile
    
    open(RCF, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<RCF>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
        }
	
	if (/^(Cut-off:)/) { #grep cut-off
	    
	    my @temp = split("\t",$_);	    #Loads tab separated info per line
	    shift(@temp); #Remove "Cut-off:"
	    my $cutoff = shift(@temp); # id string
	    $allCF{$_[1]} = $cutoff; # Hash{inputfile} string cut-off
	}
    } 	
    close(RCF);
    print STDERR "Finished Reading cut off for infile $_[0]","\n";
    return;
}


