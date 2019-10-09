#! /usr/bin/perl


#
#
# Author: Sara Ballouz <s.ballouz@victorchang.edu.au>
#
# Created: 04/12/2007
# Usage


#################################################################################
# 				Modules 					#
#################################################################################
use strict;
my $file   = $ARGV[0]; # gene_ontology.1_2.obo.rel
my $out    = $ARGV[1]; # output? 

my %PHEN_ONT;
my $i = 0;
my @col;
my $range;
my %DONE; 

sub create_hash(){
	open(ONTOLOGY,"< $file ") or die;
	while(<ONTOLOGY>){

		if(m/^GO:(.+)\tGO:(.+)\t(.+)\n/){
			my $id   = $2;
			my $is_a = $1;
			my $rel  = $3;  
			if( $rel =~ m/is_a/ || $rel =~ m/part_of/){
				push (@{$PHEN_ONT{$id}}, $is_a);	
			}		
		}

	}
	close ONTOLOGY;
	my @ke = keys %PHEN_ONT;
	print "$#ke Done\n";	
	
}



sub get_node($$){
	my $value = $_[0];
	my $level = $_[1]; 
	$level++;
	
	# So it doesn't go in an infinite loop de loooooooop 
	# Not needed
	#if($DONE{$value}){
	#	next;
	#} else {
	#	$DONE{$value}="done";
	#}
	
	if (! exists $PHEN_ONT{$value}){
		#print $value;
		return;	
	} else {
		my @temp = @{$PHEN_ONT{$value}};
		my $n = $#temp;
		my $i = 0;
		print TREE "";
		foreach (@temp) {			
			$value = $_;			
			get_node($value,$level);
			print TREE "$value";	
			
			if($i != $n){
				print TREE " ";
			}
			$i++;		
		}
		print TREE " ";	
	
	}
	
	
}


sub plot(){
	open(FILE, "cut -f1,2 $file | tr '\t' '\n' | sort | uniq |");
	open(TREE,"> $out.temp"); #with "loop" condition removed 
	while(<FILE>){
		
		chomp;
		my %DONE="";
		if (m/(GO:.+)/){
			print TREE "$1\t";	
			create_trees($1,"");	
		} else {
			next;
		}

	}
	close TREE;
	close FILE;
}

sub create_trees($$){

	my $id = $_[0];
	my $descr = $_[1];
	my @d = split("\s", $descr);
	$id =~ s/GO://g;
	
	print TREE "";
	get_node($id,1);
	print TREE "\n";
	
}

sub clean(){

	open (IN, "< $out.temp");
	open (OUT, "> $out"); 
	while(<IN>){
		if(m/(GO:.+)\t(.+)\n/){
			my $id = $1;
			my $desc = $2;
				
			my @goids = split(/\s/, $desc);
			 
			@goids = uniq(@goids);   
			print OUT "$id\tGO:";
			print OUT join " GO:", @goids;
			print OUT "\n"; 
		} else {
			print OUT $_; 
		}	
	}



}


sub uniq {
        my %seen;
        return grep { !$seen{$_}++ } @_;
}

create_hash();
plot();
clean();


