# Identify SNPs that create restriction site differences

#!/usr/bin/perl

# arguments:  Restriction_site_sequence_list, genome_sequence, list_of_SNP_coords, flanking_sequence_length

open(RSITES, "$ARGV[0]");

open(GENOME_SEQUENCE, "$ARGV[1]");

open(SNPlist, "$ARGV[2]");

open(OUTPUT_REF, '>', "SNPs_restriction_sites_in_PH-1");

print OUTPUT_REF "Contig\tSNP_position\tSequence_around_SNP\tRestriction_site\tNucleotide_in_Gz3639(position_in_restriction_site)\n\nSequence\n\n";

open(OUTPUT_COMP, '>', "SNPs_restriction_sites_in_Gz3639");

print OUTPUT_COMP "Contig\tSNP_position\tSequence_around_SNP\tRestriction_site\tNucleotide_in_PH-1(position_in_restriction_site)\n\nSequence\n\n";

# Read in restriction site list

while(<RSITES>) {

    chomp($_);

    @List = split(/\t/, $_);

    $RSite = $List[0];

    $Sequence = $List[1];

    $Sequence =~ s /\W+//;

    $RS_hash{$Sequence} = $RSite

}


# Read in genome sequence

$EOL = $/;

undef $/;

$TheFile = <GENOME_SEQUENCE>;

@List = split(/>/, $TheFile);

shift(@List);

foreach $Value (@List ) {

    @SubList = split(/\n/, $Value, 2);

   ($Contig, $Sequence) = @SubList;

    if($Contig =~ /(\d\.\d+)/) {

        $Contig = $1;

    }

    $Sequence =~ s/\W+//g;
 
   $Genome_hash{$Contig} = $Sequence

}

$/ = $EOL;


#Read in SNP coordinates list


while(<SNPlist>) {

    $TheLine = $_;

    chomp($TheLine);

    @SNPlist = split(/\t/, $TheLine);

    $Contig = $SNPlist[0];

    if($Contig =~ /(\d\.\d+)/) {

        $Contig = $1;

    }

    $Position = $SNPlist[1];

    $Nucleotide = $SNPlist[3];

    $NewNucleotide = $SNPlist[4];

    $NewNucleotide =~ s/\W+//;

#    print "$Contig\t$Position\t$Query\n";

# Capture 11 bp of sequence surrounding the SNP (5 bp on each side)


    $RefQuery = substr($Genome_hash{$Contig}, $Position - 6, 11);


# Generate the mutant sequence


    $CompQuery1 = substr($RefQuery, 0, 5);

    $CompQuery2 = substr($RefQuery, 6, 5);

    $CompQuery = "$CompQuery1"."$NewNucleotide"."$CompQuery2";

#    print "$RefQuery\t$CompQuery\n";


# Scan a 6 bp window through both sequences and look for RS in each window

    & SITE_IN_REFERENCE;

    & SITE_IN_COMPARISON;


} 


#SUB-ROUTINES:   


sub SITE_IN_REFERENCE {

    for ($i = 0; $i <= 5; $i++) {

         $Six_base_region_ref = substr($RefQuery, $i, 6);

         if(exists($RS_hash{$Six_base_region_ref})) {

            $NucleotidePosition = 6 - $i;

            $RSite = $RS_hash{$Six_base_region_ref};

            $Sequence_whole = $Genome_hash{$Contig};

            $FlankLength = $ARGV[3];

            $SeqLength = ($FlankLength*2)+1;

            $Start = $Position - $FlankLength;

            $End = $Position + $FlankLength;

            $SeqSegment = substr($Sequence_whole, $Start-1, $SeqLength);

            $SegLen = length($SeqSegment);

            print "$SegLen\n";

            print OUTPUT_REF "$Contig\t$Position\t$RefQuery\t$RSite (in PH-1)\t$NewNucleotide($NucleotidePosition)\n\n";

            print OUTPUT_REF ">$Contig"."_$Start-$End\n$SeqSegment\n\n"

        }

        else {

#            print "$Contig\t$Position\t$Query\tunknown\t$Nucleotide($NucleotidePosition)\n";

        }

    }

}

sub SITE_IN_COMPARISON {

    for ($i = 0; $i <= 5; $i++) {

        $Six_base_region_comp = substr($CompQuery, $i, 6);

        if(exists($RS_hash{$Six_base_region_comp})) {

            $NucleotidePosition = 6 - $i;

            $RSite = $RS_hash{$Six_base_region_comp};

            $Sequence_whole = $Genome_hash{$Contig};

            $FlankLength = $ARGV[3];

            $SeqLength = ($FlankLength*2)+1;

            $Start = $Position - $FlankLength;

            $End = $Position + $FlankLength;

            $SeqSegment = substr($Sequence_whole, $Start-1, $SeqLength);

            print OUTPUT_COMP "$Contig\t$Position\t$CompQuery\t$RSite (in Gz3639)\t$Nucleotide($NucleotidePosition)\n\n";

            print OUTPUT_COMP ">$Contig"."_$Start-$End\n$SeqSegment\n\n"

        }

        else {

#            print "$Contig\t$Position\t$CompQuery\tunknown\t$Nucleotide($NucleotidePosition)\n";

        }

    }

}

