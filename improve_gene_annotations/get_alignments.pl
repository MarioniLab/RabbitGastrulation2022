use strict;
use warnings;

use Text::CSV;
use Bio::EnsEMBL::Registry;
use List::Util;


my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
  -host => 'ensembldb.ensembl.org',
  -user => 'anonymous'
);

my $query_species = 'human';
my @species_names = @{ $registry->get_all_species() };

my $human_slice_adaptor = $registry->get_adaptor( 'human', 'core', 'Slice' );
my $rabbit_slice_adaptor = $registry->get_adaptor( 'oryctolagus_cuniculus', 'core', 'Slice' );

# Get the GenomeDB adaptor
my $genome_db_adaptor = $registry->get_adaptor( 'Multi', 'compara', 'GenomeDB' );

# Fetch GenomeDB objects for mouse and rabbit:
my $human_genome_db = $genome_db_adaptor->fetch_by_name_assembly('homo_sapiens');
my $rabbit_genome_db = $genome_db_adaptor->fetch_by_name_assembly('oryctolagus_cuniculus');

# Get the MethodLinkSpeciesSetAdaptor
my $method_link_species_set_adaptor = $registry->get_adaptor( 'Multi', 'compara', 'MethodLinkSpeciesSet');

# Fetch the MethodLinkSpeciesSet object corresponding to LASTZ_NET alignments between human and mouse genomic sequences
my $human_rabbit_lastz_net_mlss = $method_link_species_set_adaptor->fetch_by_method_link_type_GenomeDBs( "LASTZ_NET", [$human_genome_db, $rabbit_genome_db] );

# Get the GenomicAlignBlockAdaptor
my $genomic_align_block_adaptor = $registry->get_adaptor( 'Multi', 'compara', 'GenomicAlignBlock' );




# Input file
my $csv = Text::CSV->new({sep_char => "\t"});
my $infi = 'data-in/hs_mm_candidates.tsv';
open(my $fh, "<:encoding(UTF-8)", $infi);
$csv->getline($fh);


#Output file
my $out_tsv = Text::CSV->new({sep_char => "\t",eol=>$/});
my $outfi = 'data-out/add_unannotated_genes/alignments.tsv';
my @out_headers = ("human_gene_id","human_gene_name","rabbit_align_seq","rabbit_align_start","rabbit_align_end","rabbit_align_strand","rabbit_align_length");
open(my $out_fh, ">:encoding(UTF-8)", $outfi) or die "Failed to create $outfi: $!";
$out_tsv->print($out_fh,\@out_headers);


sub exportAlignment {
	my @out_data = @_;
	$out_tsv->print($out_fh,\@out_data);
}





sub getThreePrimeEnd {
	my ($slice_adaptor, $query_transcript) = @_;

	my $strand = $query_transcript->strand();
	my $chrom = $query_transcript->seq_region_name();
	my $start_pos;
	my $end_pos;

	my $three_prime_utr = $query_transcript->three_prime_utr_Feature();
	my $three_prime_exon = @{$query_transcript->get_all_Exons()}[-1];
	if(! defined $three_prime_utr) {$three_prime_utr = $three_prime_exon;}

	if($query_transcript->strand()==1) {
		$start_pos = $three_prime_exon->seq_region_start();
		$end_pos = $three_prime_utr->seq_region_end();
		
	} elsif($query_transcript->strand==-1) {
		$start_pos = $three_prime_utr->seq_region_start();
		$end_pos = $three_prime_exon->seq_region_end();
	} else {
		return;
	}

	
	my $query_slice = $slice_adaptor->fetch_by_region( 'toplevel', $chrom, $start_pos, $end_pos );
	
	return $query_slice;

}



sub filterSliceAlignments {
	my ($query_species,@genomic_aligns) = @_;

	my @output_slices;
	my $seq_seen;
	
	foreach my $genomic_align (@genomic_aligns) {
		my $cur_species = $genomic_align->genome_db->get_scientific_name;
		if($cur_species ne $query_species) {next;}
		
		my $slice = $genomic_align->get_Slice;
		if(!$seq_seen) {$seq_seen = $slice->seq_region_name();}
		
		my @slice_genes = @{$slice->get_all_Genes};
		my @slice_exons = @{$slice->get_all_Exons};

		my $alignment_length = $genomic_align->length;
		
		#check it's not overlaying existing gene
		my $slice_gene = "NA";
		if(@slice_genes) { 
			my $slice_gene = $slice_genes[0]->external_name();
		} 
		
		if($slice_gene ne "NA") {
			print "Skipping alignment - slice overlaps existing annotation: ", $slice_gene , ".\n";
			next;
		}

		if($slice->seq_region_name() ne $seq_seen) {
			print "Skipping alignment - max scoring alignments of the same slice are mapping to different chromosomes. Keeping the first alignment only. \n "; 
			next;
		}
		
		print "Accepted alignment: ", $cur_species, "\t", $slice->seq_region_name, ":", $slice->seq_region_start, "-", $slice->seq_region_end,"\t",
			$slice->strand(), "\t", $alignment_length, "\n";

		push(@output_slices,$slice);

   	
	}


	return @output_slices;

} 


sub getSliceAlignments {
	my ($query_species,$query_slice) = @_;

	# Fetch all the GenomicAlignBlocks corresponding to this Slice from the pairwise alignments (LASTZ_NET) between human and rabbit
	my @genomic_align_blocks = @{ $genomic_align_block_adaptor->fetch_all_by_MethodLinkSpeciesSet_Slice( $human_rabbit_lastz_net_mlss, $query_slice ) };
	
	if(!@genomic_align_blocks) {print "No alignment found.\n"; return;}

	my @align_block_scores = map { $_->score() } @genomic_align_blocks;
	my @sorted_align_block_scores = sort { $align_block_scores[$b] <=> $align_block_scores[$a]} 0..$#align_block_scores;

	my $max_score = $align_block_scores[$sorted_align_block_scores[0]];

	my @slice_alignments;
	foreach my $genomic_align_block( @genomic_align_blocks ) {

		# Only choose alignments with maximum scores
		my $alignment_score = $genomic_align_block->score();
		
		if($alignment_score < $max_score) {
			print "Skipping alignment block - alignment block score (",$alignment_score, ") < maximum found (", $max_score, "). \n";
			next;
		} 

		my $restricted_gab = $genomic_align_block->restrict_between_reference_positions($query_slice->seq_region_start(), $query_slice->seq_region_end());
	
		# fetch the GenomicAligns and move through
		my @genomic_aligns = @ { $restricted_gab->get_all_GenomicAligns };
		if(!@genomic_aligns) {print "No alignment found.\n"; next; }

		my @gb_alignments = filterSliceAlignments($query_species,@genomic_aligns);
		
		push(@slice_alignments,@gb_alignments);

	}

	return (\@slice_alignments,$max_score);
}



sub getThreePrimeAlignments {
	my ($query_species,$slice_adaptor,$query_transcript) = @_;

	my $query_slice = getThreePrimeEnd($slice_adaptor,$query_transcript);
	my ($alignments,$score) = getSliceAlignments($query_species,$query_slice);
	
	return ($alignments,$score);
}


#TODO: Export print statements to a log file. 
sub getAlignments {
	my ($query_species,$slice_adaptor,$gene_id) = @_;

	# Get the SliceAdaptor and fetch a slice
	my $gene_slice = $slice_adaptor->fetch_by_gene_stable_id($gene_id);
	#my $gene_slice = $slice_adaptor->fetch_by_region( 'toplevel', $seq_region, $seq_region_start, $seq_region_end );


	print "\nFinding ", $query_species, " 3' alignments for gene: ", $gene_id, "... \n";


	my @query_transcripts = @{$gene_slice->get_all_Transcripts_by_type('protein_coding')};
	if(! @query_transcripts) {
		print "No transcripts found.\n";
		return;
	}


	my @output_alignments;
	my $seq_seen;
	my $max_score_seen = 0; 

	foreach my $query_transcript (@query_transcripts) {
		print "\n------------------------------------------------------\n";
		print "Transcript - ", $query_transcript->stable_id(), "\t", $query_transcript->seq_region_name(), "\t", $query_transcript->seq_region_start(),
		"\t", $query_transcript->seq_region_end(), "\n\n";

		my ($alignments,$score) = getThreePrimeAlignments($query_species,$slice_adaptor,$query_transcript);

		my @alignments;
		if($alignments) {@alignments = @$alignments;}
		else {next;}
		
		# Store chromosome and score of first non-empty alignments
		if(!$seq_seen && @alignments) {
			$seq_seen = $alignments[0]->seq_region_name();
		}


		# If transcript aligns to a different chromosome as previously seen, keep the alignments with the highest score. 
		if($alignments[0]->seq_region_name() ne $seq_seen) {
			print "Transcript is aligning to a different chromosome (chrom ", $alignments[0]->seq_region_name(),", score: ", $score, ").  Keeping alignments with maximum score (chrom: ";
			if($score > $max_score_seen) {
				print $alignments[0]->seq_region_name(), ", score: ", $score, "\n";
				@output_alignments=@alignments;
				$max_score_seen = $score;
				$seq_seen = $alignments[0]->seq_region_name();
				next;
			}
			else {
				print $seq_seen, " , score: ", $max_score_seen, "\n";
				next;
			}	
		}

		if($score > $max_score_seen) {$max_score_seen = $score;}


		push(@output_alignments,@alignments);

	}
	
	print "\n------------------------------------------------------\n\n\n\n";

	return @output_alignments;

};


sub main {

	while(my $line = <$fh>) {
		if($csv->parse($line)) {
			my @columns = $csv->fields();
			print "$columns[6]\t$columns[7]\t$columns[8]\t$columns[9]\t$columns[10]\n\n";
		
			my $ref_gene_id = $columns[6];
			my $seq_region = $columns[8];
			my $seq_region_start = $columns[9];
			my $seq_region_end   = $columns[10];

			my @align_slices = getAlignments("Oryctolagus cuniculus",$human_slice_adaptor,$ref_gene_id);
			foreach my $align_slice (@align_slices) {
				exportAlignment($columns[6],$columns[7],$align_slice->seq_region_name,
						$align_slice->seq_region_start,$align_slice->seq_region_end,
						$align_slice->strand,$align_slice->length);
			} 
		

		}
	}

	close $fh;
	close $out_fh;

}

main();

#getAlignments("Oryctolagus cuniculus",$human_slice_adaptor,"ENSG00000170500");
# Example aligning rabbit region with human
#getAlignments("Homo sapiens",$rabbit_slice_adaptor,10,10003000,10005000);









