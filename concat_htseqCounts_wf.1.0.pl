#!/usr/bin/perl -w

use strict;
use warnings;

#Script for concatinating HTSeqCount data
#Copyright 2011 Henrik Stranneheim

=head1 SYNOPSIS
    
concat_htseqCounts_wf.1.0.pl  -i [infile...n] -s [sample ID...n] -o [outfile]

=head2 COMMANDS AND OPTIONS

-i/--infile Infile(s), comma sep

-S/--sampleid The sample ID(s), comma sep

-o/--outfile Output file (Defaults to concat_htseqCount_wf.txt)

=head3 I/O

Input format ( HTSeqCount.txt )

Output format:

1. concat_samplid_HTSeqCount.txt

=cut
    
use Pod::Usage;
use Pod::Text;
use Getopt::Long;

use vars qw($USAGE);

BEGIN {
    $USAGE =
        qq{concat_htseqCounts_wf.1.0.pl  -i [infile...n] -s [sample ID...n] -o [outfile.txt]
               -i/--infile Infile(s), comma sep
               -s/--sampleid The sample ID,comma sep
               -o/--outfile Output file (Defaults to concat_htseqCount_wf.txt)
	   };
    
}

my ($o, $sid, $help) = ("concat_htseqCount_wf.txt"); #Arguments for program

my (@infn,@sid);
my (%allRC);

GetOptions('i|infile:s'  => \@infn, #Comma separeted list
           's|sampleid:s'  => \@sid, #Comma separeted list, one below outdirdata
           'o|outfile:s'  => \$o, #Outfile
           'h|help' => \$help,
           );

die $USAGE if( $help );

if (@infn == 0) {
    my $verbosity = 2;
 print"\n";
    pod2usage({-message => "Must supply an infile directory as comma separeted list.\n",
	       -verbose => $verbosity
	       });
}
if ( scalar(@sid) eq 0) {
    
    print STDERR "\n";
    print STDERR "Must supply a sample ID as a comma separated list", "\n\n";
    die $USAGE;
}

@infn = split(/,/,join(',',@infn)); #Enables comma separated indir(s)
@sid = split(/,/,join(',',@sid)); #Enables comma separeted list of sample IDs

for (my $inputfiles=0;$inputfiles<scalar(@infn);$inputfiles++) {
    
    ReadHTSEQC($infn[$inputfiles], $inputfiles);	
}

WriteHTSEQC($o);

sub WriteHTSEQC {
 
    my $filename = shift;
    
    open (WHTSEQC, ">$filename") or die "Can't write to $filename: $!\n";
    
    print WHTSEQC "GeneID", "\t"; #Prints header geneid

    for (my $sampleid=0;$sampleid<scalar(@sid);$sampleid++) {

	print WHTSEQC $sid[$sampleid], "\t";

    }
    print WHTSEQC "\n";
    
    for my $inputfile ( keys %allRC ) { #For all input files
	
	if ($inputfile eq 0) { #Selects only the first file to print IDs
	    
	    for my $id ( keys %{ $allRC{0} } ) { #Print IDs from first file
		
		print WHTSEQC "$id","\t";
		
		for (my $i=0;$i<scalar( keys %allRC );$i++) { #prints read counts for all files

		    print WHTSEQC $allRC{$i}{$id}[0], "\t";
		    
		}
		print WHTSEQC "\n";
	    }
	}
	
    }
    close (WHTSEQC);
    return;
}

sub ReadHTSEQC {
    
#$_[0] = filename
#$_[1] = Nr of inputfile
    
    open(RHTSEQC, "<$_[0]") or die "Can't open $_[0]:$!, \n";    
    
    while (<RHTSEQC>) {
        chomp $_;
	
	if (m/^\s+$/) {		# Avoid blank lines
	    next;
        }		
	if (/(\S+)/) {
	    
	    if (m/no_feature/) { #Removes statistics we do not want in DESeq
		close(RHTSEQC);
		print STDERR "Finished Reading infile $_[0]","\n";
		return;
	    }
	    my @temp = split("\t",$_);	    #Loads ensemble id and read counts
	    my $enid = shift(@temp); # chr string
	    $allRC{$_[1]}{$enid} = [@temp]; # Hash{inputfile} of Hash{ID} of array[Read_Counts]
	}
    } 	
    close(RHTSEQC);
    print STDERR "Finished Reading infile $_[0]","\n";
    return;
}
