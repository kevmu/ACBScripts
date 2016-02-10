#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

use Getopt::Long;

#

my ($gmap_summary_infile, $gmap_gff_infile, $blastx_infile, $output_dir);
my @options = (
'i=s'    => \$gmap_summary_infile,
'g=s'    => \$gmap_gff_infile,
'b=s'    => \$blastx_infile,
'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
    defined $gmap_summary_infile
    and defined $gmap_gff_infile
    and defined $blastx_infile
    and $output_dir
);

sub usage {
    
die << "USAGE";
    
Usage: $0 -i gmap_summary_infile -g gmap_gff_infile -b blastx_infile -o output_dir
    
USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
    mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %gmap_summary = ();
open(INFILE, "<$gmap_summary_infile") or die "Couldn't open file $gmap_summary_infile for reading, $!";
my $i = 0;
while(<INFILE>){
    chomp $_;
    if(($i ne 0) and ($_ !~ m/^$/)){
        #warn $_ . "\n";
        my @split_row_entry = split(/\t/, $_);
        my ($alignment_name, $coverage, $percent_id, $matches, $mismatches, $indels, $unknowns, $query_id, $query_start, $query_end, $query_length, $query_strand, $target_id, $target_start, $target_end, $target_length, $target_strand, $amino_acid_start, $amino_acid_end, $amino_acid_length, $amino_acid_changes, $num_exons, $query_align_block, $target_align_block, $align_identity_block, $intron_length_block) = @split_row_entry;
            my ($query_name, $query_suffix) = split(/\s/, $query_id, 2);
        
        my $alignment_id = "";
        if($alignment_name =~ m/path(\d+)of(\d+)/){
            my $path_num = $1;
            my $total_paths = $2;
            $alignment_id = join(".", $query_name, "path$path_num");
            
            
        }else{
            
            die "no path $alignment_name";
        }
        
        if(($coverage >= 80) and ($percent_id >= 90) and ($query_length >= 500)){
            $gmap_summary{$target_id}{$alignment_id} = join("\t", $query_name, $query_id, $coverage, $percent_id, $matches, $mismatches, $indels, $unknowns, $query_start, $query_end, $query_length, $query_strand, $target_start, $target_end, $target_length, $target_strand, $amino_acid_start, $amino_acid_end, $amino_acid_length, $amino_acid_changes, $num_exons);
        
        }
    }
    $i++;
}
close(INFILE) or die "Couldn't close file $gmap_summary_infile";



my %gmap_gff = ();
open(INFILE, "<$gmap_gff_infile") or die "Couldn't open file $gmap_gff_infile for reading, $!";
my $current_query_path = "";
while(<INFILE>){
    chomp $_;
    if(($_ !~ m/^$/) and ($_ !~ m/^#/)){
        #warn $_ . "\n";
        my @split_row_entry = split(/\t/, $_);
        my ($target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = @split_row_entry;
    # JH246661.1, $GCA_000230575.1_canSat3_genome, $gene, $7826, $10258, $., $+, $., $ID=gi|351590686|gb|JP449145.1|.path1;Name=gi|351590686|gb|JP449145.1|

        warn join("\t", $target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) . "\n";
        my %attributes = ();
        my @split_attributes = split(/;/, $attribute);
        foreach my $attribute_entry (@split_attributes){
            my ($attribute_id, $attribute_value) = split(/=/, $attribute_entry);
            
            $attributes{$attribute_id} = $attribute_value;
        }
        if($feature eq "gene"){
            my $query_name = $attributes{"ID"};
            $current_query_path = $query_name;
            $gmap_gff{$target_id}{$query_name}{"gene"} = join("\t", $target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute);
        }
        if($feature eq "mRNA"){
            
            my $query_name = $attributes{"Parent"};
            $current_query_path = $query_name;
            $gmap_gff{$target_id}{$query_name}{"mRNA"} = join("\t", $target_id, $source, $feature, $start, $end, $score, $strand, $frame, join(";", $attributes{"ID"}, $attributes{"Name"}, $attributes{"Parent"}));
        }
        if($feature eq "exon"){
            push(@{$gmap_gff{$target_id}{$current_query_path}{"exon"}}, join("\t", $target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute));
        
        }
        if($feature eq "CDS"){
            
            push(@{$gmap_gff{$target_id}{$current_query_path}{"CDS"}}, join("\t", $target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute));
        
        }
    }
}
close(INFILE) or die "Couldn't close file $gmap_gff_infile";

my %blastx = ();
open(INFILE, "<$blastx_infile") or die "Couldn't open file $blastx_infile for reading, $!";
$i = 0;
while(<INFILE>){
    chomp $_;
    if(($i ne 0) and ($_ !~ m/^$/)){
        warn $_ . "\n";
        my @split_row_entry = split(/\t/, $_);
        my ($query_name, $target_name, $query_coverage, $protein_query_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @split_row_entry;
#        die join("\t", $query_name, $target_name) . "\n";
        push(@{$blastx{$query_name}}, join("\t", $query_name, $target_name, $query_coverage, $protein_query_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score));
    }
    $i++;
}
close(INFILE) or die "Couldn't close file $blastx_infile";

my $filename = fileparse($gmap_gff_infile);
my $outfile = join('/', $output_dir, $filename);
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print "##gff-version   3" . "\n";
foreach my $target_id (sort keys %gmap_summary){
    warn $target_id . "\n";
    foreach my $alignment_id (sort keys $gmap_summary{$target_id}){
        
        my ($query_name, $query_id, $coverage, $percent_id, $matches, $mismatches, $indels, $unknowns, $query_start, $query_end, $query_length, $query_strand, $target_start, $target_end, $target_length, $target_strand, $amino_acid_start, $amino_acid_end, $amino_acid_length, $amino_acid_changes, $num_exons) = split(/\t/, $gmap_summary{$target_id}{$alignment_id});
#        warn join("===",$query_name, $query_id);
        if(defined($gmap_gff{$target_id}{$alignment_id}{"gene"})){
            warn $gmap_gff{$target_id}{$alignment_id}{"gene"} . "\n";
            
            warn $gmap_summary{$target_id}{$alignment_id} . "\n";
        }else{
            die "gene";
        }
        if(defined($gmap_gff{$target_id}{$alignment_id}{"mRNA"})){
            warn $gmap_gff{$target_id}{$alignment_id}{"mRNA"} . "\n";
            if(defined(@{$blastx{$target_id}{$query_name}})){
                warn join("\n", @{$blastx{$target_id}{$query_name}}) . "\n";
            }
            
        }else{
            die "mRNA";
        }
        if(defined(@{$gmap_gff{$target_id}{$alignment_id}{"exon"}})){
            
            warn join("\n", @{$gmap_gff{$target_id}{$alignment_id}{"exon"}}) . "\n";
            
        }else{
            die "exon";
        }
        
        if(defined(@{$gmap_gff{$target_id}{$alignment_id}{"CDS"}})){
            
            warn join("\n", @{$gmap_gff{$target_id}{$alignment_id}{"CDS"}}) . "\n";
            
        }
        
    }
    print "###" . "\n";
}
close(OUTFILE) or die "Couldn't close file $outfile";

(%gmap_summary, %gmap_gff, %blastx) = ();
