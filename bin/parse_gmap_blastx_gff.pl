#!/usr/bin/perl
use warnings;
use strict;
use File::Basename;

use Getopt::Long;

#

my ($gmap_summary_infile, $gmap_gff_infile, $blastx_infile, $scaffolds_infile, $output_dir);
my @options = (
	'i=s'    => \$gmap_summary_infile,
	'g=s'    => \$gmap_gff_infile,
	'b=s'    => \$blastx_infile,
	's=s'    => \$scaffolds_infile,
	'o=s'    => \$output_dir,
);
&GetOptions(@options);


usage() unless (
    defined $gmap_summary_infile
    and defined $gmap_gff_infile
    and defined $blastx_infile
    and defined $scaffolds_infile
    and $output_dir
);

sub usage {
    
die << "USAGE";
    
Usage: $0 -i gmap_summary_infile -g gmap_gff_infile -b blastx_infile -s scaffolds_infile -o output_dir
    
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
        
        if(($coverage >= 80) and ($percent_id >= 90)){
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
            $gmap_gff{$target_id}{$query_name}{"mRNA"} = join("\t", $target_id, $source, $feature, $start, $end, $score, $strand, $frame, join(";", join("=", "ID", $attributes{"ID"}), join("=", "Name", $attributes{"Name"}), join("=", "Parent", $attributes{"Parent"})));
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

my %scaffolds_list = ();
open(INFILE, "<$scaffolds_infile") or die "Couldn't open file $scaffolds_infile for reading, $!";
while(<INFILE>){
    chomp $_;
    warn $_ . "\n";
	my ($target_id, $target_header) = split(/\t/, $_);
	$scaffolds_list{$target_id} = join("\t", $target_id, $target_header);
}
close(INFILE) or die "Couldn't close file $scaffolds_infile";

my $filename = fileparse($gmap_gff_infile, qr/\.gff/);
my $outfile = join('/', $output_dir, $filename . ".jbrowse.gff");
open(OUTFILE, ">$outfile") or die "Couldn't open file $outfile for writting, $!";
print OUTFILE "##gff-version   3" . "\n";
foreach my $target_id (sort keys %gmap_summary){
    warn $target_id . "\n";
    foreach my $alignment_id (sort keys %{$gmap_summary{$target_id}}){
        
        my ($query_name, $query_id, $coverage, $percent_id, $matches, $mismatches, $indels, $unknowns, $query_start, $query_end, $query_length, $query_strand, $target_start, $target_end, $target_length, $target_strand, $amino_acid_start, $amino_acid_end, $amino_acid_length, $amino_acid_changes, $num_exons) = split(/\t/, $gmap_summary{$target_id}{$alignment_id});
#        warn join("===",$query_name, $query_id);
#$target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, );
        if(defined($gmap_gff{$target_id}{$alignment_id}{"gene"})){
            #warn $gmap_gff{$target_id}{$alignment_id}{"gene"} . "\n";
            my ($target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $gmap_gff{$target_id}{$alignment_id}{"gene"});
			my %attributes = ();
			my @split_attributes = split(/;/, $attribute);
			foreach my $attribute_entry (@split_attributes){
				my ($attribute_id, $attribute_value) = split(/=/, $attribute_entry);
				
				$attributes{$attribute_id} = $attribute_value;
			}
			my $note_attribute = ""; 
			if(defined(@{$blastx{$query_name}})){
				my ($query_seq, $target_seq, $query_coverage, $protein_query_coverage, $percent_identity, $percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, $query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = split(/\t/, @{$blastx{$query_name}}[0]);
			
				my $target_seq_header = "";
				if($target_seq =~ m/ gi/){
					my ($target_header, $target_header_concat) = split(' gi', $target_seq, 2);
					$target_seq_header = $target_header;
				}
				else{
					$target_seq_header = $target_seq;
				}
				#die $target_header1;
				#	die join("\n", @{$blastx{$query_name}}) . "\n";
				$note_attribute = join("=", "Note", $target_seq_header);

			}
			
			my $name_attribute = join("=", "Name", $query_id);
			my $organism_attribute = join("=", "organism", "Cannabis sativa <a href%3D%22http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi%3Fmode%3DInfo%26id%3D3483%26lvl%3D3%26lin%3Df%26keep%3D1%26srchmode%3D1%26unlock%22%20target%3D%22_blank%22>[ Taxonomy Browser ]</a>");
			

			
			# ID=gi|351591582|gb|JP450028.1|.path1
			# ;
			# Name=gi|351591582|gb|JP450028.1| TSA: Cannabis sativa PK00237.2_1.CasaPuKu mRNA sequence
			# ;
			# organism=Cannabis sativa <a href%3D%22http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi%3Fmode%3DInfo%26id%3D3483%26lvl%3D3%26lin%3Df%26keep%3D1%26srchmode%3D1%26unlock%22%20target%3D%22_blank%22>[ Taxonomy Browser ]</a>
			# ;
			# Note=gi|645241791|ref|XP_008227245.1| PREDICTED: chromodomain-helicase-DNA-binding protein 1 [Prunus mume]
			# ;
			my ($target_toc_id, $target_toc_header) = split(/\t/, $scaffolds_list{$target_id});
			$target_toc_header =~ s/(scaffold\d+)/<b>$1<\/b>/g;
			$target_toc_header =~ s/,/%2C/g;
			my $gmap_target_id = join("%20", "<b>Target_ID:</b>", join(" ", $target_id, $target_toc_header));

			my $gmap_target_genebank_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $target_id, "%3Freport%3Dgenbank%22%20target%3D%22_blank%22>[ <b>GenBank</b> ]</a>");
			my $gmap_target_fasta_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $target_id, "%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>");
			my $gmap_target_ncbi_revisions_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $target_id, "%3Freport%3Dgirevhist%22%20target%3D%22_blank%22>[ <b>NCBI Revision History</b> ]</a>");

			my $scaffold_id = "";
			if($target_toc_header =~ m/(scaffold\d+)/){
				$scaffold_id = $1;
			}
			my $gmap_target_cannabis_unsc_gbrowser = join("", "<a href%3D%22http://genome.ccbr.utoronto.ca/cgi-bin/hgTracks%3Fdb%3DcanSat3%26position%3D", $scaffold_id, "%22%20target%3D%22_blank%22>[ <b>Cannabis UCSC Genome Browser</b> ]</a>");
									
			my $gmap_target_url_strings = join(" ", $gmap_target_genebank_url, $gmap_target_fasta_url, $gmap_target_ncbi_revisions_url, $gmap_target_cannabis_unsc_gbrowser);

			
			my ($jbrowse_query_id, $query_header) = split(/\s/, $query_id, 2);
			
			$query_header =~ s/(PK\d+\.\d+\_\d+)/<b>$1<\/b>/g;
			my @split_query_header = split(/\|/, $query_name);
			
			my $query_accession_num = $split_query_header[3];
			my $gmap_query_id = join(" ", "<b>Query_ID:</b>", $query_accession_num, $query_header);
			
			
			my $gmap_query_genebank_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $query_accession_num, "%3Freport%3Dgenbank%22%20target%3D%22_blank%22>[ <b>GenBank</b> ]</a>");
			my $gmap_query_fasta_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $query_accession_num, "%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>");
			my $gmap_query_ncbi_revisions_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $query_accession_num, "%3Freport%3Dgirevhist%22%20target%3D%22_blank%22>[ <b>NCBI Revision History</b> ]</a>");
			my $gmap_query_url_strings = join(" ", $gmap_query_genebank_url, $gmap_query_fasta_url, $gmap_query_ncbi_revisions_url);


			#my $gmap_alignment_overview = "<b>GMAP Alignment Overview</b>";
			my $gmap_query_coverage = join(" ", "<b>query_coverage:</b>", $coverage);
			my $gmap_percent_identity = join(" ", "<b>percent_identity:</b>", $percent_id);

			my $gmap_num_matches = join(" ", "<b>num_matches:</b>", $matches);
			my $gmap_num_mismatches = join(" ", "<b>num_mismatches:</b>", $mismatches);

			my $gmap_num_indels = join(" ", "<b>num_indels:</b>", $indels);
			my $gmap_num_unknowns = join(" ", "<b>num_unknowns:</b>", $unknowns);

			my $gmap_query_start = join(" ", "<b>query_start:</b>", $query_start);
			my $gmap_query_end = join(" ", "<b>query_end:</b>", $query_end);
			my $gmap_query_length = join(" ", "<b>query_length:</b>", $query_length, "bp");
			my $gmap_query_strand = join(" ", "<b>query_strand:</b>", $query_strand);
			
			my $gmap_target_start = join(" ", "<b>target_start:</b>", $target_start);
			my $gmap_target_end = join(" ", "<b>target_end:</b>", $target_end);
			my $gmap_target_length = join(" ", "<b>target_length:</b>", $target_length, "bp");
			my $gmap_target_strand = join(" ", "<b>target_strand:</b>", $target_strand); 
			
			my $gmap_amino_acid_start = join(" ", "<b>amino_acid_start:</b>", $amino_acid_start);
			my $gmap_amino_acid_end = join(" ", "<b>amino_acid_end:</b>", $amino_acid_end);
			my $gmap_amino_acid_length = join(" ", "<b>amino_acid_length:</b>", $amino_acid_length, "aa");
			my $gmap_amino_acid_changes = join(" ", "<b>amino_acid_changes:</b>", $amino_acid_changes);  
			my $gmap_num_exons = join(" ", "<b>num_exons:</b>", $num_exons);
			
			#my $gmap_alignment_summary = join("<br>", $gmap_alignment_overview, join("%3B ", $gmap_query_coverage, $gmap_percent_identity), join("%3B ", $gmap_num_matches, $gmap_num_mismatches), join("%3B ", $gmap_num_indels, $gmap_num_unknowns), join("%3B ", $gmap_query_start, $gmap_query_end, $gmap_query_length, $gmap_query_strand), join("%3B ", $gmap_target_start, $gmap_target_end, $gmap_target_length, $gmap_target_strand), join("%3B ", $gmap_amino_acid_start, $gmap_amino_acid_end, $gmap_amino_acid_length), $gmap_amino_acid_changes, $gmap_num_exons . " <br><br>");
			my $gmap_alignment_summary = join("<br>", join("%3B ", $gmap_query_coverage, $gmap_percent_identity), join("%3B ", $gmap_num_matches, $gmap_num_mismatches), join("%3B ", $gmap_num_indels, $gmap_num_unknowns), join("%3B ", $gmap_query_start, $gmap_query_end, $gmap_query_length, $gmap_query_strand), join("%3B ", $gmap_target_start, $gmap_target_end, $gmap_target_length, $gmap_target_strand), join("%3B ", $gmap_amino_acid_start, $gmap_amino_acid_end, $gmap_amino_acid_length), $gmap_amino_acid_changes, $gmap_num_exons . " <br><br>");

			my $gmap_concat_attributes = join("<br><br> ", join("<br>", $gmap_target_id, $gmap_target_url_strings), join("<br>", $gmap_query_id, $gmap_query_url_strings), $gmap_alignment_summary);

			## gmap_alignment=<b>Target_ID:</b>%20<b>JH243385.1</b> Cannabis sativa unplaced genomic scaffold <b>scaffold25337</b>%2C whole genome shotgun sequence
			 # <br> 
			## <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JH243385.1%3Freport%3Dgenbank%22%20target%3D%22_blank%22>[ <b>GenBank</b> ]</a>
			## <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JH243385.1%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>
			## <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JH243385.1%3Freport%3Dgirevhist%22%20target%3D%22_blank%22>[ <b>NCBI Revision History</b> ]</a>
			## <a href%3D%22http://genome.ccbr.utoronto.ca/cgi-bin/hgTracks%3Fhgsid%3D68876%26position%3Dscaffold25337%22%20target%3D%22_blank%22>[ <b>Cannabis UCSC Genome Browser</b> ]</a>
			 # <br><br>
			# <b>Query_ID:<b> <b>JP450028.1</b> TSA: Cannabis sativa <b>PK00237.2_1</b>.CasaPuKu mRNA sequence
			# <br> 
			# <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JP450028.1%3Freport%3Dgenbank%22%20target%3D%22_blank%22>[ <b>GenBank</b> ]</a>
			 # <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JP450028.1%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>
			 # <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JP450028.1%3Freport%3Dgirevhist%22%20target%3D%22_blank%22>[ <b>NCBI Revision History</b> ]</a>
			 # <br><br>
			# <a href%3D%22http://ec2-54-201-126-170.us-west-2.compute.amazonaws.com/jbrowse/data/gmap_summary/JH243385.1_JP450028.1_gmap_summary.html%22%20target%3D%22_blank%22>[ <b>GMAP Alignment Overview</b> ]</a>
			# <br>
			# <b>query_coverage:</b> 100.0%3B <b>percent_identity:</b> 100.0
			# <br>
			# <b>num_matches:</b> 4901%3B <b>num_mismatches:</b> 1
			# <br> 
			# <b>num_indels:</b> 0%3B <b>num_unknowns:</b> 0<br><br>
						#die join(";\n", $attributes{"ID"}, $name_attribute, $organism_attribute, $note_attribute);
			my $gmap_align_attribute = join("=", "gmap_alignment", $gmap_concat_attributes);
			my $gene_attributes = "";
			
			if($note_attribute ne ""){
				$gene_attributes = join(";", join("=", "ID", $attributes{"ID"}), $name_attribute, $organism_attribute, $note_attribute, $gmap_align_attribute);
			}else{

				$gene_attributes = join(";", join("=", "ID", $attributes{"ID"}), $name_attribute, $organism_attribute, $gmap_align_attribute);
			}
			print OUTFILE join("\t", $target_id, "GMAP%5CBLASTX", $feature, $start, $end, $score, $strand, $frame, $gene_attributes) . "\n";
        }else{
            die "gene";
        }
        if(defined($gmap_gff{$target_id}{$alignment_id}{"mRNA"})){
            #print OUTFILE $gmap_gff{$target_id}{$alignment_id}{"mRNA"} . "\n";
			my ($target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $gmap_gff{$target_id}{$alignment_id}{"mRNA"});
			
			#print OUTFILE join("\t", $target_id, "GMAP%5CBLASTX", $feature, $start, $end, $score, $strand, $frame, $attribute) . "\n";
			
			my %attributes = ();
			my @split_attributes = split(/;/, $attribute);
			foreach my $attribute_entry (@split_attributes){
				my ($attribute_id, $attribute_value) = split(/=/, $attribute_entry);
				
				$attributes{$attribute_id} = $attribute_value;
			}
			
			my $blastx_annotation_attributes = "";
			if(defined(@{$blastx{$query_name}})){
				my %blastx_hits = ();
				for(my $i = 0; $i < scalar(@{$blastx{$query_name}}); $i++){
					my ($query_seq, $target_seq, $query_coverage, $protein_query_coverage, $percent_identity, 
						$percent_positives, $query_length, $target_length, $align_length, $num_mismatch, $num_gaps, 
						$query_start, $query_end, $target_start, $target_end, $e_value, $bit_score) = split(/\t/, @{$blastx{$query_name}}[$i]);
				
					my $query_seq_header = "";
					if($query_seq =~ m/ gi/){
						my ($query_header, $query_header_concat) = split(' gi', $query_seq, 2);
						$query_seq_header = $query_header;
					}
					else{
						$query_seq_header = $query_seq;
					}
					
					
					my ($query_seq_id, $query_sequence_header) = split(/\s/, $query_id, 2);
					
					my @split_query_header = split(/\|/, $query_seq_id);
					
					
					my $query_gi_num = "<b>$split_query_header[1]</b>";
					my $query_accession_num = "<b>$split_query_header[3]</b>";
                    my $query_accession_id = $split_query_header[3];
                    $query_sequence_header =~ s/(PK\d+\.\d+\_\d+)/<b>$1<\/b>/g;
                    
					$query_seq_header = join("| ", $split_query_header[0], $query_gi_num, $split_query_header[2], $query_accession_num, $query_sequence_header);
                    #die $query_seq_header;
					my $target_seq_header = "";
					if($target_seq =~ m/ gi/){
						my ($target_header, $target_header_concat) = split(' gi', $target_seq, 2);
						$target_seq_header = $target_header;
					}
					else{
						$target_seq_header = $target_seq;
					}
					
					my ($target_seq_id, $target_sequence_header) = split(/\s/, $target_seq_header, 2);
					
					my @split_target_header = split(/\|/, $target_seq_id);
					my $target_accession_num = $split_target_header[3];
					
					my $target_annotation = $target_sequence_header;
					
					my $blastx_counter = ($i + 1);
					my $blast_hit_name = join("", "Blastx_hit", $blastx_counter);
					
					#die $query_seq;
					my $blastx_query_genebank_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $query_accession_id, "%3Freport%3Dgenbank%22%20target%3D%22_blank%22>[ <b>GenBank</b> ]</a>");
					my $blastx_query_fasta_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $query_accession_id, "%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>");
					my $blastx_query_ncbi_revisions_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/", $query_accession_id, "%3Freport%3Dgirevhist%22%20target%3D%22_blank%22>[ <b>NCBI Revision History</b> ]</a>");
					my $blastx_query_url_strings = join(" ", $blastx_query_genebank_url, $blastx_query_fasta_url, $blastx_query_ncbi_revisions_url);
					
					my $blastx_target_genepept_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/protein/", $target_accession_num, "%3Freport%3Dgenpept%22%20target%3D%22_blank%22>[ <b>GenPept</b> ]</a>");
					my $blastx_target_fasta_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/protein/", $target_accession_num, "%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>");
					my $blastx_target_ident_protiens_url = join("", "<a href%3D%22http://www.ncbi.nlm.nih.gov/protein/", $target_accession_num, "%3Freport%3Dipg%22%20target%3D%22_blank%22>[ <b>Identical Proteins</b> ]</a>");
					my $blastx_target_url_strings = join(" ", $blastx_target_genepept_url, $blastx_target_fasta_url, $blastx_target_ident_protiens_url);
					
					my $blastx_annotation = join(" ", "Target_ID:", "<b>$target_accession_num</b>", $target_annotation);
                    my $blastx_query_id = join(" ", "Query_ID:", "<b>$query_accession_num</b>", $query_sequence_header);
					my $blastx_query_coverage = join(" ", "<b>query_transcript_coverage:</b>", $query_coverage);
					my $blastx_protein_query_coverage = join(" ", "<b>query_protein_coverage:</b>", $protein_query_coverage);
					
					my $blastx_percent_identity = join(" ", "<b>percent_identity:</b>", $percent_identity);
					my $blastx_percent_positives = join(" ", "<b>percent_positives:</b>", $percent_positives);
					
					my $blastx_num_mismatches = join(" ", "<b>num_mismatches:</b>", $num_mismatch);
					my $blastx_num_gaps = join(" ", "<b>num_gaps:</b>", $num_gaps);
				
					my $blastx_query_length = join(" ", "<b>query_length:</b>", $query_length, "bp");
					my $blastx_target_length = join(" ", "<b>target_length:</b>", $target_length, "aa");
					my $blastx_align_length = join(" ", "<b>alignment_length:</b>", $align_length, "aa");
					
					my $blastx_query_start = join(" ", "<b>query_start:</b>", $query_start);
					my $blastx_query_end = join(" ", "<b>query_end:</b>", $query_end);
					
					my $blastx_target_start = join(" ", "<b>target_start:</b>", $target_start);
					my $blastx_target_end = join(" ", "<b>target_end:</b>", $target_end);

					my $blastx_e_value = join(" ", "<b>e_value:</b>", $e_value);
					my $blastx_bit_score = join(" ", "<b>bit_score:</b>", $bit_score);

					
					my $blastx_alignment_summary = join("<br>", join("%3B ", $blastx_query_coverage, $blastx_protein_query_coverage), join("%3B ", $blastx_percent_identity, $blastx_percent_positives),  join("%3B ", $blastx_num_mismatches, $blastx_num_gaps), join("%3B ", $blastx_query_start, $blastx_query_end, $blastx_query_length), join("%3B ", $blastx_target_start, $blastx_target_end, $blastx_target_length, $blastx_align_length), join("%3B ", $blastx_e_value, $blastx_bit_score));

					my $blastx_concat_attributes = join("<br><br> ", join("<br>", $blastx_query_id, $blastx_query_url_strings), join("<br>", $blastx_annotation, $blastx_target_url_strings), $blastx_alignment_summary);
					
					$blastx_hits{$blast_hit_name} = join("=", $blast_hit_name, $blastx_concat_attributes);
					
# Blastx_hit1=gi|<b>351591582</b>|gb|<b>JP450028.1</b>| TSA: Cannabis sativa <b>PK00237.2_1</b>.CasaPuKu mRNA sequence
# <br>
# <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JP450028.1%3Freport%3Dgenbank%22%20target%3D%22_blank%22>[ <b>GenBank</b> ]</a>
# <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JP450028.1%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>
# <a href%3D%22http://www.ncbi.nlm.nih.gov/nuccore/JP450028.1%3Freport%3Dgirevhist%22%20target%3D%22_blank%22>[ <b>NCBI Revision History</b> ]</a>
# <br><br><b>XP_008227245.1</b>| PREDICTED: chromodomain-helicase-DNA-binding protein 1 [Prunus mume]
# <br>
# <a href%3D%22http://www.ncbi.nlm.nih.gov/protein/XP_008227245.1%3Freport%3Dgenpept%22%20target%3D%22_blank%22>[ <b>GenPept</b> ]</a>
# <a href%3D%22http://www.ncbi.nlm.nih.gov/protein/XP_008227245.1%3Freport%3Dfasta%22%20target%3D%22_blank%22>[ <b>FASTA</b> ]</a>
# <a href%3D%22http://www.ncbi.nlm.nih.gov/protein/XP_008227245.1%3Freport%3Dipg%22%20target%3D%22_blank%22>[ <b>Identical Proteins</b> ]</a>
# <br><br>
# <a href%3D%22http://ec2-54-201-126-170.us-west-2.compute.amazonaws.com/jbrowse/data/blastx_summaries/JP450028.1-XP_008227245.1_blastx_summary.html%22 target%3D%22_blank%22>[ <b>BLASTx Alignment Summary</b> ]</a>
# <br>
# <b>coverage:</b> 90%3B <b>percent_identity:</b> 80.32%3B <b>percent_positives:</b> 88.24
# <br>
# <b>query_length:</b> 4903 bp%3B <b>target_length:</b> 1760aa%3B <b>alignment_length:</b> 1479aa
# <br>
# <b>mismatches:</b> 272%3B  <b>gaps:</b> 9
# <br>
# <b>query_start:</b> 489%3B  <b>query end:</b> 4892
# <br>
# <b>target_start:</b> 1%3B <b>target_end:</b> 1471
# <br>
# <b>e_value:</b> < 1e-179%3B <b>bit_score:</b> 2324
				}
				
				my @blastx_attributes = ();
				foreach my $blast_hit_name (sort keys %blastx_hits){
					push(@blastx_attributes, $blastx_hits{$blast_hit_name});
				}
				$blastx_annotation_attributes = join(";", @blastx_attributes);
			}
            my @split_query_name = split(/\|/, $query_name);
            
            my $query_accession_name = $split_query_name[3];
            
            my $blastx_query_search = join("=", "alignment_search", "<a href%3D%22http://blast.ncbi.nlm.nih.gov/Blast.cgi%3FPROGRAM%3Dblastx%26PAGE_TYPE%3DBlastSearch%26BLAST_SPEC%3D%26LINK_LOC%3Dblasttab%26LAST_PAGE%3Dblastx%26QUERY%3D$query_accession_name%22 target%3D%22_blank%22>$query_accession_name (Run NCBI BLASTX)</a>");
			my $mnrna_attributes = "";
			if($blastx_annotation_attributes ne ""){
				$mnrna_attributes = join(";", $attribute, $blastx_query_search, $blastx_annotation_attributes);
			}else{

				$mnrna_attributes = join(";", $attribute, $blastx_query_search);
			}
			
			print OUTFILE join("\t", $target_id, "GMAP%5CBLASTX", $feature, $start, $end, $score, $strand, $frame, $mnrna_attributes) . "\n";

			
        }else{
            die "mRNA";
        }
        if(defined(@{$gmap_gff{$target_id}{$alignment_id}{"exon"}})){
            
			foreach my $exon_entry (@{$gmap_gff{$target_id}{$alignment_id}{"exon"}}){
				my ($target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $exon_entry);
				print OUTFILE join("\t", $target_id, "GMAP%5CBLASTX", $feature, $start, $end, $score, $strand, $frame, $attribute) . "\n";
            }
        }else{
            die "exon";
        }
        
        if(defined(@{$gmap_gff{$target_id}{$alignment_id}{"CDS"}})){
            
			foreach my $cds_entry (@{$gmap_gff{$target_id}{$alignment_id}{"CDS"}}){
				my ($target_id, $source, $feature, $start, $end, $score, $strand, $frame, $attribute) = split(/\t/, $cds_entry);
				print OUTFILE join("\t", $target_id, "GMAP%5CBLASTX", $feature, $start, $end, $score, $strand, $frame, $attribute) . "\n";
            }
        }
        
    }
    print OUTFILE "###" . "\n";
}
close(OUTFILE) or die "Couldn't close file $outfile";

(%gmap_summary, %gmap_gff, %blastx) = ();
