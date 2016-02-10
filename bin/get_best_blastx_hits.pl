#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

use Getopt::Long;

# perl get_best_blast_hits.pl -l /home/cookeadmin/workspace/adriana/Phloem/phloem-microarray-annotations/blastx_output-2014-07-03/phloem-blastx-file-list -f /home/cookeadmin/workspace/adriana/ncbi_nr_flat_files_2014-07-02/ncbi_nr_db_notes.txt -n phloem-master-blastx-file -o /home/cookeadmin/workspace/test

my ($blastx_infile, $output_dir);
my @options = (
      'i=s'    => \$blastx_infile,
      'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
      defined $blastx_infile
      and $output_dir
);

sub usage {
    
    die << "USAGE";
    
Usage: $0 -i blastx_infile -o output_dir
    
USAGE
}

## Create output directory if it doesn't already exist.
unless(-d $output_dir){
	mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %blast_annotations = ();
open(INFILE, "<$blastx_infile") or die "Couldn't open file $blastx_infile for reading, $!";
my $i = 0;
while(<INFILE>){
chomp $_;
if(($i ne 0) and ($_ !~ m/^$/)){
# 			warn $_ . "\n";
    #query_name	target_name	query_coverage	query_hsp_coverage	percent_identity	percent_positives	query_length	target_length	align_length	num_mismatch	num_gaps	query_start	query_end	target_start	target_end	e_value	bit_score
  my @split_row_entry = split(/\t/, $_);
  my $query_id = $split_row_entry[0];

  my $blast_entry = join("\t", @split_row_entry);
# 			warn $blast_entry . "\n";

  push(@{$blast_annotations{$query_id}}, [split(/\t/, $blast_entry)]);
}
$i++;
}
close(INFILE) or die "Couldn't close file $blastx_infile";

    

my $filename = fileparse($blastx_infile, qr/\.tsv\.txt/);
my ($blastx_filename, $blastx_dir) = fileparse($blastx_infile, qr/\.tsv\.txt/);

my $outfile = join('/', $output_dir, $blastx_filename . ".parsed.txt");
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "protein_query_coverage", "percent_identity", "percent_positives", "query_length", "target_length", "align_length", "num_mismatch", "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
foreach my $query_id (sort keys %blast_annotations){
    #warn $query_id . "\n";
      
	    my @blast_annotation_sorted = sort {$b->[16] <=> $a->[16]} @{$blast_annotations{$query_id}};
        my $i = 0;
	    foreach my $blast_annotation_entry (@blast_annotation_sorted){
            if($i < 5){
                my ($query_name, $target_name, $query_coverage, $query_hsp_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @$blast_annotation_entry;
                my $protein_query_coverage = ((($target_end - $target_start)/$target_length) * 100);
                #warn join("\t", $query_name, $target_name, $query_coverage, $query_hsp_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
                #die $protein_query_coverage;
                if(($protein_query_coverage >= 80) and ($percent_identity >= 80)){
                    print OUTFILE join("\t", $query_name, $target_name, $query_coverage, sprintf("%.2f", $protein_query_coverage), $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
                }
                
            }
           $i++;
	    }

}
close(OUTFILE) or die "Couldn't close file $outfile";


#my $filename = fileparse($blastx_infile, qr/\.tsv\.txt/);
#my ($blastx_filename, $blastx_dir) = fileparse($blastx_infile, qr/\.tsv\.txt/);
#
#my $outfile = join('/', $output_dir, $blastx_filename . ".tsv.txt");
#open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
#print OUTFILE join("\t", "query_name", "target_name", "query_coverage", "query_hsp_coverage", "percent_identity", "percent_positives", "query_length", "target_length", "align_length", "num_mismatch", "num_gaps", "query_start", "query_end", "target_start", "target_end", "e_value", "bit_score") . "\n";
#foreach my $blast_entry (sort {$b->[16] <=> $a->[16]} @blast_entries){
#      warn join("\t", @$blast_entry) . "\n";
#
#      
#      my ($query_name, $target_name, $query_coverage, $query_hsp_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = @$blast_entry;
#
#    print OUTFILE join("\t", $query_name, $target_name, $query_coverage, $query_hsp_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) . "\n";
#}
#close(OUTFILE) or die "Couldn't close file $outfile";
