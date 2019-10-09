#! /usr/bin/perl

#######################################
# Parse GO files 
#######################################
#
# Script to parse GO onntology files 
#
#
# Created by: Sara Ballouz
# Date: 12th Feb 2013
#
#######################################
# Modified: 
#
#######################################
#
#
#
#
#######################################
# Use: 
#
# perl parse_GO.pl <gene_ontology> <output file>

#######################################
# Main
#######################################

#if($#ARGV =! 2) {
#	die "$0: Please specify gene ontology file, gene association file and an output file name\n"; 

#} 

$gene_ontol_file = $ARGV[0];

parse_ontology($gene_ontol_file);


sub parse_ontology($){
	my $FILE = $_[0];

	open( IN, "< $FILE");
	open (OUT, "> $FILE.rel");
	open( VOC, "> $FILE.voc"); 
	while(<IN>){
		chomp;
		my $line = $_; 
		if ($line =~ m/^id: (GO:.+)/){
			$GOID = $1;
			$HASH{$GOID} =1;
		} 
		elsif( $line =~ m/^name: (.+)/){
			$GOname = $1;
		}
	        elsif( $line =~ m/^namespace: (.+)/){
	        	$GOnspace = $1;
		}
		elsif( $line =~ m/^is_a: (GO:.+) !/){
			print OUT "$GOID\t$1\tis_a\n";
			push(@{$IS_A{$1}}, $GOID );
		}
		elsif( $line =~ m/^relationship: (.+) (GO:.+) !/){
                        my $rel_id = $1;
			my $rel_go = $2; 
			print OUT "$GOID\t$rel_go\t$rel_id\n";
			if( $rel_id =~ m/part_of/){
				push(@{$PART_OF{$rel_go}}, $GOID );
			}
			if( $rel_id =~ m/^regulates/){
                                push(@{$REG{$rel_go}}, $GOID );
                        }
			if( $rel_id =~ m/positive/){
                                push(@{$POS_REG{$rel_go}}, $GOID );
                        }
			if( $rel_id =~ m/negative/){
                                push(@{$NEG_REG{$rel_go}}, $GOID );
                        }
	
		}

		elsif( $line =~ m/^\[Term\]/){
			print VOC "$GOID\t$GOname\t$GOnspace\n";
			
		}
	}
	close OUT; 
	close IN; 
	close VOC; 
	
	open( MAP, "> $FILE.map");

	foreach $GOID (keys %HASH){
		my @is_a = @{$IS_A{$GOID}}    ? @{$IS_A{$GOID}}    : "-";
		my @part = @{$PART_OF{$GOID}} ? @{$PART_OF{$GOID}} : "-"; 
		my @reg  = @{$REG{$GOID}}     ? @{$REG{$GOID}}     : "-";
		my @pos  = @{$POS_REG{$GOID}} ? @{$POS_REG{$GOID}} : "-";
		my @neg  = @{$NEG_REG{$GOID}} ? @{$NEG_REG{$GOID}} : "-";
		print MAP "$GOID\t@is_a\t@part\t@reg\t@pos\t@neg\n"; 	
	}	
	close MAP;

}




