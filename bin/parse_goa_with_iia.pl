#! /usr/bin/perl
#
########################################
## Parse gene association file
########################################
## 
## Created by: Sara Ballouz
## Date: 5th Mar 2013



$ID2NAME  = $ARGV[0];
$GOA_FILE = $ARGV[1];

@ID = `cut -f1 $ID2NAME`; 
chomp(@ID);

@NAME = `cut -f2 $ID2NAME`;
chomp(@NAME);

@ID2NAME{@ID} = @NAME;
@NAME2ID{@NAME} = @ID;


# foreach $key (keys %ID2NAME){
#	print "$key\n";
#}


open (IN, "< $GOA_FILE");
open (LOG, "> $GOA_FILE.log");
open (OUT, "> $GOA_FILE.parsed");


while (<IN>){
	
	if(m/^(.+)\t(GO:.+)\t(.+)\t(.+)\n/){
		my $name = $1;
		my $goid = $2;
		my $iia  = $3;
		my $alt  = $4; 
		
		if($NAME2ID{$name}){
			print OUT "$name\t$NAME2ID{$name}\t$goid\t$iia\n";
		}
		else { 
			my @altnames = split(/\|/, $alt);
			my $appr;	
			foreach $altname (@altnames){
				 
				if($NAME2ID{$altname}){
					$appr = $altname;
				}

			}
			if($appr){	#print "$name missing trying $list[0]\n"; 
				print OUT "$appr\t$NAME2ID{$appr}\t$goid\t$iia\n"; #
				print LOG "Used alternate $appr, $name was missing\n";
			} else { 
				 print LOG "$name missing, tried @altnames\n";
			}	
		}
	}
	
}


close IN;
close OUT;
close LOG; 



