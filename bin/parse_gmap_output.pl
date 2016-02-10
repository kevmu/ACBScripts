#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

use Bio::SeqIO;
use File::Basename;
use Switch;

# perl process_gmap_output-new.pl -i /home/cookeadmin/workspace/cathy/AllSNPs_output_GMAP -o /home/cookeadmin/cathy_snps

my ($gmap_infile, $output_dir);
GetOptions(
      'i=s'    => \$gmap_infile,
      'o=s'    => \$output_dir,
);

usage() unless (
      defined $gmap_infile
      and defined $output_dir
);

sub usage {
    
die <<"USAGE";
    
Usage: $0 -i gmap_infile -o output_dir
    
Description - 
    
OPTIONS:
     
      -i gmap_infile -
      -o output_dir -
    
USAGE
}

# Create output directory if it doesn't already exist.
unless(-d $output_dir){
      mkdir($output_dir, 0777) or die "Can't make directory: $!";
}

my %parsed_gmap = ();
my $gmap_query_id = "";
open(INFILE, "<$gmap_infile") or die "Couldn't open file $gmap_infile for reading, $!";
while(<INFILE>){
      chomp $_;

      if($_ =~ m/^>/){
	    $gmap_query_id = $_;
          
	    $gmap_query_id =~ s/>//g;
		
          #die join("\t", ($gmap_query_id, $gmap_query_title));
		  push(@{$parsed_gmap{$gmap_query_id}}, $gmap_query_id);
      }else{
          warn $_ . "\n";
	    push(@{$parsed_gmap{$gmap_query_id}}, $_);
      }
      
}
close(INFILE) or die "Couldn't close file $gmap_infile";

my $outfilename = fileparse($gmap_infile);
my $gmap_outfile = join('/', $output_dir, $outfilename . ".parsed.txt");
open(OUTFILE, ">$gmap_outfile") or die "Couldn't open file $gmap_outfile for writting, $!";                          
print OUTFILE join("\t", "alignment_name","coverage","percent_id","matches","mismatches","indels", "unknowns","query_id","query_start","query_end","query_length","query_strand","target_id","target_start",
		  "target_end","target_length","target_strand","amino_acid_start","amino_acid_end","amino_acid_length",
		  "amino_acid_changes","num_exons","query_align_block","target_align_block", "align_identity_block", 
		  "intron_length_block","target_gap_positions","query_gap_positions","mismatch_positions") . "\n";
my %gmap_output_final = ();
foreach my $gmap_query_id (keys %parsed_gmap){
      my %processed_gmap = ();
      my $gmap_entry = join("\t", @{$parsed_gmap{$gmap_query_id}});
    
      my ($alignment_paths,$path_num,$coverage,$percent_identity,$matches,$mismatches,$indels,
	    $unknowns,$query_start,$query_end,$query_length,$query_strand,$target_id,
	    $target_start,$target_end,$target_length,$target_strand,$amino_acid_start,
	    $amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons,$intron_blocks);

      my ($current_query_start, $current_query_end, $current_query_length, $current_target_id, 
      $current_target_start, $current_target_end, $current_target_length);
      my $current_path = 1;
      my $counter = 1;
      foreach my $field (@{$parsed_gmap{$gmap_query_id}}){
        
	    if($field =~ m/^Paths \((\d+)\):/){
		  $alignment_paths = $1;
		  if($alignment_paths eq 0){ # Exit out of loop if $alignment_paths equals 0 indicating no alignment info to parse.
			last;
		  }
	    }
	    
	    if($field =~ m/Path (\d+): query ([\d,]+)\.\.([\d,]+) \(([\d,]+) bp\) => genome ([\w\_\-\.\d\|\:\s]+):([\d,]+)\.\.([\d,]+) \(-*([\d,]+) bp\)/){
	    
		  ($path_num, $query_start, $query_end, $query_length, $target_id, $target_start, $target_end, $target_length) = ($1, $2, $3, $4, $5, $6, $7, $8);
#die join("\t",$path_num, $query_start, $query_end, $query_length, $target_id, $target_start, $target_end, $target_length);
		  $target_start =~ s/,//g;
		  $target_end =~ s/,//g;
            ;
		  my $target_temp = -1;
		  if($target_start > $target_end){
			$target_temp = $target_start;
			$target_start = $target_end;
			$target_end = $target_temp;
		  }
	    }

	    
	    if(defined($path_num)){
		  warn "path: $path_num eq $current_path";
		  if($path_num eq $current_path){
			if($field =~ m/cDNA direction: (sense|antisense|indeterminate)/){
			      
			      my $strandedness = $1; 
			      switch ($strandedness) {
				    case "sense" {
					  $query_strand = "+"; 
				    }
				    case "antisense" { 
					  $query_strand = "-";
				    }
				    case "indeterminate" {
					  $query_strand = "?"; 
				    }
				    else { 
					  die "Strandedness: $strandedness" 
				    }
			      }
			}
			if($field =~ m/Genomic pos: [\w\_\-\.\d\|\:\s]+:[\d,]+,[\d,]+\.\.[\d,]+,[\d,]+ \(([+-]) strand\)/){
			      $target_strand = $1;
			}
			if($field =~ m/Number of exons: (\d+)/){
			      $num_exons = $1;
			}
			if($field =~ m/Coverage: (\d+\.\d+) \(query length: \d+ bp\)/){
			      $coverage = $1;
			}
			if($field =~ m/Percent identity: (\d+\.\d+) \((\d+) matches, (\d+) mismatches, (\d+) indels, (\d+) unknowns\)/){

			      ($percent_identity,$matches,$mismatches,$indels,$unknowns) = ($1,$2,$3,$4,$5);

			}
			if($field =~ m/Translation: (\d+)\.\.(\d+) \((\d+) aa\)/){

			      ($amino_acid_start,$amino_acid_end,$amino_acid_length) = ($1,$2,$3);

			      my $amino_acid_temp = -1;
			      if($amino_acid_start > $amino_acid_end){
				    $amino_acid_temp = $amino_acid_start;
				    $amino_acid_start = $amino_acid_end;
				    $amino_acid_end = $amino_acid_temp;
			      }

			}
			if($field =~ m/Amino acid changes: ((.+)|)/){
			      
			      my $has_aa_change = $2;
			      if(defined($has_aa_change)){
				    my $amino_acid_changes_list = $2;
				    my @amino_acid_change_entries = split(/, /, $amino_acid_changes_list);
				    my @amino_acid_changes_final = ();
				    foreach my $amino_acid_change_entry (@amino_acid_change_entries){
					  warn $amino_acid_change_entry;
					  my ($aa_change, $aa_position);
					  if($amino_acid_change_entry =~ m/(.+) \[(\d+)\]/){
						($aa_change, $aa_position) = ($1, $2);
					  }
					  push(@amino_acid_changes_final, join(";", $aa_change, $aa_position));
				    }

				    $amino_acid_changes = join(",", @amino_acid_changes_final);

			      }else{
				    $amino_acid_changes = "N/A";
			      }
			      ($current_query_start, $current_query_end, $current_query_length, $current_target_id, $current_target_start, $current_target_end, $current_target_length) = ($query_start, $query_end, $query_length, $target_id, $target_start, $target_end, $target_length);
			}
		  }
		  if($path_num > $current_path){

      # 		  print OUTFILE join("\t", join("-", $gmap_query_id, "cmp", $current_target_id, join("_", "path", $current_path, "of", $alignment_paths)),$coverage,$percent_identity,$matches,$mismatches,$indels,
      # 		  $unknowns,$gmap_query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
      # 		  $target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons) . "\n";
			$amino_acid_start = "N/A" unless(defined($amino_acid_start));
			$amino_acid_end = "N/A" unless(defined($amino_acid_end));
			$amino_acid_length = "N/A" unless(defined($amino_acid_length));

			$processed_gmap{$current_path} = join("\t", join("-", $gmap_query_id, "cmp", $current_target_id, "path" . $current_path . "of" . $alignment_paths),$coverage,$percent_identity,$matches,$mismatches,$indels,
			$unknowns,$gmap_query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
			$target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons);
			$current_path++;

			warn "path: $path_num > $current_path";
		  }
	    }

	    
      }
      next if ($alignment_paths eq 0); # Exit out of loop if $alignment_paths equals 0 indicating no alignment info to parse.
#       print OUTFILE join("\t", join("-", $gmap_query_id, "cmp", $current_target_id, join("_", "path", $current_path, "of", $alignment_paths)),$coverage,$percent_identity,$matches,$mismatches,$indels,
#       $unknowns,$gmap_query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
#       $target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons) . "\n";

      $amino_acid_start = "N/A" unless(defined($amino_acid_start));
      $amino_acid_end = "N/A" unless(defined($amino_acid_end));
      $amino_acid_length = "N/A" unless(defined($amino_acid_length));

      $processed_gmap{$current_path} = join("\t", join("-", $gmap_query_id, "cmp", $current_target_id, "path" . $current_path . "of" . $alignment_paths),$coverage,$percent_identity,$matches,$mismatches,$indels,
      $unknowns,$gmap_query_id,$current_query_start,$current_query_end,$current_query_length,$query_strand,$current_target_id,$current_target_start,$current_target_end,$current_target_length,
      $target_strand,$amino_acid_start,$amino_acid_end,$amino_acid_length,$amino_acid_changes,$num_exons);
			
      my $found_align_tag = "false";
      my $found_target_align = "false";
      my (@query_align_blocks,@target_align_blocks,@intron_length_blocks,@align_identity_blocks,@query_protein_fragments, @target_protein_fragments, @query_alignment,@target_alignment,@consensus_alignment) = ();
      
      my ($block_path_num,$query_align_block,$target_align_block,$intron_length_block,$align_identity_block,$query_align_sequence,$target_align_sequence,$consensus_align_sequence);

      my $current_block_path = 1;
      foreach my $field (@{$parsed_gmap{$gmap_query_id}}){
	    if($field =~ m/Alignments:/){
		  $found_align_tag = "true";
# 		  die $found_align_tag;
	    }

	    if($found_align_tag eq "true"){
		  
		  if($field =~ m/Alignment for path (\d+):/){
			$block_path_num = $1;
		  }

		  if(defined($block_path_num)){
			warn "Block path: $block_path_num eq $current_block_path";
			if($block_path_num eq $current_block_path){
			      if($field =~ m/^\s+[\+-][\w\_\-\.\d\|\:\s]+:(\d+-\d+)  \((\d+-\d+)\)   (\d+)%( (<-|\(-|==|-\)|->)   ...(\d+)...  (\d+\.\d+), (\d+\.\d+)|)/){
				    my ($target_align, $query_align, $align_identity, $has_intron) = ($1,$2,$3,$4);
		  # 			warn "$target_align, $query_align, $align_identity, $has_intron";

				    push(@target_align_blocks, $target_align);

				    push(@query_align_blocks, $query_align);
				    
				    push(@align_identity_blocks, $align_identity);
				    if(defined($has_intron) and $has_intron ne ""){
					  my ($intron_length, $intron_coverage, $intron_identity) = ($6,$7,$8);
					  $intron_blocks = join(";", $intron_length, $intron_coverage, $intron_identity);
				    }else{

					  $intron_blocks = "N/A";
				    }
				    push(@intron_length_blocks, $intron_blocks);
                  }
#                  else{
#
#                $field =~ s/^\s+//;
#                warn "$field and $counter\n";
#                if(($field =~ m/^\d+\s[\s\.\:]+$/)){
#                    warn $field;
#                    $counter = 1;
#                }
#                if(($field =~ m/^aa\.g\s+\d+\s([ABCDEFGHIKLMNPQRSTVWXYZ\s\*]+)$/) and ($counter eq 2)){
#                    warn $field;
#                    my $target_protein_partial = $1;
#                    $target_protein_partial =~ s/\s|\*//g;
#                    push(@target_protein_fragments, $target_protein_partial);
#                }
#                if(($field =~ m/^\s+$/) and ($counter eq 2)){
#                    warn $field;
#                }
#                #  +tscaffold7014:134446 TAGCTTGACAATGTAAAAAAGCCTCAGAATTTCATGCTTTTGTCTATTAC
#			    if(($field =~ m/^[\+-][\w\d\.\-\s\_]+:\d+\s([AGCTN\.\s]+)$/i) and ($counter eq 3)){
#                      
#                    warn $field;
#                    my $target_align_partial = $1;
#				    push(@target_alignment, $target_align_partial);
#			      
#                }
#                if(($field =~ m/^[\|\=\s\.\>\<\-]+$/) and ($counter eq 4)){
#                    
#                    warn $field;
#                }
#                # 46 TAGCTTGACAATGTAAAAAAGCCTCAGAATTTCATGCTTTTGTCTATTAC
#                if(($field =~ m/^\d+\s([AGCTN\s\d]+)$/i) and ($counter eq 5)){
#                      warn $field;
#                    my $query_align_partial = $1;
#				    push(@query_alignment, $query_align_partial);
#                }
#                if(($field =~ m/^aa\.c\s+\d+\s([ABCDEFGHIKLMNPQRSTVWXYZ\s\*]+)$/) and ($counter eq 6)){
#                    warn $field;
#                    my $query_protein_partial = $1;
#                    $query_protein_partial =~ s/\s|\*//g;
#                    push(@query_protein_fragments, $query_protein_partial);
#                    
#                }
#                      
#                      if(($field =~ m/^\s+$/) and ($counter > 6)){
#                          warn $field;
#                      }
#                    $counter++;
#                  }
			}
			if($block_path_num > $current_block_path){
       
			      $query_align_block = join(",", @query_align_blocks);
			      $target_align_block = join(",", @target_align_blocks);
			      $intron_length_block = join(",", @intron_length_blocks);
			      $align_identity_block = join(",", @align_identity_blocks);
			      
#			      $query_align_sequence = join("", @query_alignment);
#			      $target_align_sequence = join("", @target_alignment);
#                warn "Target Sequence: $target_align_sequence";
#               warn "query Sequence: $query_align_sequence";
#                warn $query_align_sequence . "\n";
#                warn $target_align_sequence . "\n";
#                $target_align_sequence =~ s/[AGCTN][AGCTN][AGCTN]\.\.\.[AGCTN][AGCTN][AGCTN]//gi;
#                $target_align_sequence =~ s/\s/-/g;
#                
#                $query_align_sequence =~ s/\s\s\s\d+\s\s\s//g;
#                $query_align_sequence =~ s/\s/-/g;
#
#			      my ($query_align_seq_length,$target_align_seq_length) = (length($query_align_sequence),length($target_align_sequence));
#			      warn join("\t",$processed_gmap{$current_block_path});
#                warn "Target Sequence: $target_align_sequence";
#                warn "query Sequence: $query_align_sequence";
#                die "Error: query_align_seq_length=$query_align_seq_length and target_align_seq_length=$target_align_seq_length are not equal" unless($target_align_seq_length eq $query_align_seq_length);
#
#			      warn join("\t", $target_align_seq_length, $query_align_seq_length);
#
#			      my @target_align_seq_chars = split('', $target_align_sequence);
#			      my @query_align_seq_chars = split('', $query_align_sequence);
#
#			      my (@target_gap_positions,@query_gap_positions,@mismatch_positions) = ();
#			      for(my $i = 0; $i < $query_align_seq_length; $i++){
#				    my $position = ($i + 1);
#				    if($target_align_seq_chars[$i] ne $query_align_seq_chars[$i]){
#					  warn join("\t", $target_align_seq_chars[$i], $query_align_seq_chars[$i]);
#					  if($target_align_seq_chars[$i] eq "-"){
#						push(@target_gap_positions, $position);
#					  }elsif($query_align_seq_chars[$i] eq "-"){
#						push(@query_gap_positions, $position);
#					  }else{
#						push(@mismatch_positions, $position);
#					  }
#				    }
#			      }
#
#			      my ($target_gap_position_list,$query_gap_position_list,$mismatch_position_list) = "";
#			      $target_gap_position_list = join(",", @target_gap_positions);
#			      $query_gap_position_list = join(",", @query_gap_positions);
#			      $mismatch_position_list = join(",", @mismatch_positions);
#
#			      $target_gap_position_list = "N/A" if($target_gap_position_list eq "");
#			      $query_gap_position_list = "N/A" if($query_gap_position_list eq "");
#			      $mismatch_position_list = "N/A" if($mismatch_position_list eq "");

#			      $processed_gmap{$current_block_path} = join("\t",$processed_gmap{$current_block_path},$query_align_block,$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,$query_gap_position_list,$mismatch_position_list);
                $processed_gmap{$current_block_path} = join("\t",$processed_gmap{$current_block_path},$query_align_block,$target_align_block,$align_identity_block,$intron_length_block);
			      (@query_align_blocks,@target_align_blocks,@intron_length_blocks,@align_identity_blocks,@query_alignment,@target_alignment,@consensus_alignment) = ();

			      warn $processed_gmap{$current_block_path};
#			      warn $target_align_sequence;
#			      warn $query_align_sequence;

			      warn "Block path: $block_path_num > $current_block_path";
			      $current_block_path++;
			}
		  }
	    }
      }

      $query_align_block = join(",", @query_align_blocks);
      $target_align_block = join(",", @target_align_blocks);
      $intron_length_block = join(",", @intron_length_blocks);
      $align_identity_block = join(",", @align_identity_blocks);

#      $query_align_sequence = join("", @query_alignment);
#      $target_align_sequence = join("", @target_alignment);
#
#    
#    warn $query_align_sequence . "\n";
#    warn $target_align_sequence . "\n";
#    
#    $target_align_sequence =~ s/[AGCTN][AGCTN][AGCTN]\.\.\.[AGCTN][AGCTN][AGCTN]//gi;
#    $target_align_sequence =~ s/\s/-/g;
#    
#    $query_align_sequence =~ s/\s\s\s\d+\s\s\s//g;
#    $query_align_sequence =~ s/\s/-/g;

    #die $query_align_sequence;

#      my ($query_align_seq_length,$target_align_seq_length) = (length($query_align_sequence),length($target_align_sequence));
#      warn join("\t",$processed_gmap{$current_block_path});
#    
#    warn $query_align_sequence . "\n";
#    warn $target_align_sequence . "\n";
#
#    die "Error: query_align_seq_length=$query_align_seq_length and target_align_seq_length=$target_align_seq_length are not equal" unless($target_align_seq_length eq $query_align_seq_length);
#
#      warn join("\t", $target_align_seq_length, $query_align_seq_length);
#
#      my @target_align_seq_chars = split('', $target_align_sequence);
#      my @query_align_seq_chars = split('', $query_align_sequence);
#
#      my (@target_gap_positions,@query_gap_positions,@mismatch_positions) = ();
#      for(my $i = 0; $i < $query_align_seq_length; $i++){
#	    my $position = ($i + 1);
#	    if($target_align_seq_chars[$i] ne $query_align_seq_chars[$i]){
#		  warn join("\t", $target_align_seq_chars[$i], $query_align_seq_chars[$i]);
#		  if($target_align_seq_chars[$i] eq "-"){
#			push(@target_gap_positions, $position);
#		  }elsif($query_align_seq_chars[$i] eq "-"){
#			push(@query_gap_positions, $position);
#		  }else{
#			push(@mismatch_positions, $position);
#		  }
#	    }
#      }
#
#      my ($target_gap_position_list,$query_gap_position_list,$mismatch_position_list) = "";
#      $target_gap_position_list = join(",", @target_gap_positions);
#      $query_gap_position_list = join(",", @query_gap_positions);
#      $mismatch_position_list = join(",", @mismatch_positions);
#
#      $target_gap_position_list = "N/A" if($target_gap_position_list eq "");
#      $query_gap_position_list = "N/A" if($query_gap_position_list eq "");
#      $mismatch_position_list = "N/A" if($mismatch_position_list eq "");

#      $processed_gmap{$current_block_path} = join("\t",$processed_gmap{$current_block_path},$query_align_block,$target_align_block,$align_identity_block,$intron_length_block,$target_gap_position_list,$query_gap_position_list,$mismatch_position_list);
          $processed_gmap{$current_block_path} = join("\t",$processed_gmap{$current_block_path},$query_align_block,$target_align_block,$align_identity_block,$intron_length_block);
      warn $processed_gmap{$current_block_path};
#      warn $target_align_sequence;
#      warn $query_align_sequence;

      foreach my $processed_path_num (keys %processed_gmap){
	    my @processed_gmap_entry = split(/\t/, $processed_gmap{$processed_path_num});
	    my $query_name = $processed_gmap_entry[7];
	    $gmap_output_final{$query_name}{$processed_path_num} = $processed_gmap{$processed_path_num};
      }


}

foreach my $query_name (sort {$a cmp $b} keys %gmap_output_final){
      foreach my $processed_path_num (sort %{$gmap_output_final{$query_name}}){
	    print OUTFILE $gmap_output_final{$query_name}{$processed_path_num} . "\n" if(defined($gmap_output_final{$query_name}{$processed_path_num}));
      }
}
close(OUTFILE) or die "Couldn't close file $gmap_outfile";
