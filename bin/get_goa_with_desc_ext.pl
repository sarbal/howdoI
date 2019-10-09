#! /usr/bin/perl
#
########################################
## Parse gene association file
########################################


$GOA_FILE = $ARGV[0];
$PHENTREE = $ARGV[1]; 


sub uniq {
        my %seen;
        return grep { !$seen{$_}++ } @_;
}


open (IN, "< $GOA_FILE" );

while(<IN>){
	if(m/^(.+)\t(.+)\t(GO:.+)\t(.+)\n/){
		my $name = $1;
		my $uid  = $2;
		my $goid = $3;
		my $iia  = $4;  
		push(@{$GOA{$goid}}, "$uid\t$iia");
		$ID2NAME{$uid} = $name;
		#my $key = "$uid$goid";
		#push(@{$IIA{$key}},$iia);
	}
	

}

close IN; 



open (IN, "< $PHENTREE" );
while(<IN>){
        if(m/(GO:.+)\t(.+)\n/){
		my $goid = $1;
                my $desc = $2;
		my @genes = @{$GOA{$goid}};
		my @desc  = split(/\s/, $desc);
	

		foreach $id (@desc){
			push(@genes, @{$GOA{$id}});
			
		}
		@genes = uniq(@genes);

		foreach $key (@genes){
			@id = split(/\t/,$key); 
			my $uid = $id[0];
			my $iia = $id[1]; 
		 	print "$ID2NAME{$uid}\t$uid\t$goid\t$iia\n";
		}        
	} elsif(m/(.+)\t\n/){
		my $goid = $1;
		#print $goid; 
		my @genes = @{$GOA{$goid}};
		@genes = uniq(@genes);
		foreach $key (@genes){
			@id = split(/\t/,$key);
                        my $uid = $id[0];
                        my $iia = $id[1];

                        print "$ID2NAME{$uid}\t$uid\t$goid\t$iia\n";
                }
		#print ; 
	}
	
}
close IN;




