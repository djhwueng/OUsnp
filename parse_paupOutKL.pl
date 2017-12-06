#!/usr/bin/perl -w
use diagnostics;
use strict;
open(IN,"model.scores");
open(OUT,">summary.txt");
print OUT "model\tloglikelihood\tK\tAC\tAG\tAT\tCG\tCT\tGT\tfA\tfC\tfG\tfT\tP-Inv\tgamma shape";
my @modelNames=("JC","JC+G","F81","F81+G","K80","K80+G","HKY","HKY+G","TrNef","TrNef+G","TrN","TrN+G","K3P","K3P+G","K3Puf","K3Puf+G","TIMef","TIMef+G","TIM","TIM+G","TVMef","TVMef+G","TVM","TVM+G","SYM","SYM+G","GTR","GTR+G","GTR+G+I");
my @modelK=(1,2,4,5,2,3,5,6,3,4,6,7,3,4,6,7,4,5,7,8,5,6,8,9,6,7,9,10,11); #add the number of free parameters here
my @transitions=(1,1,1,1,1,1);
my $lnL=0;
my $K=0;
my $model="none";
my @basefreq=(0.25,0.25,0.25,0.25);
my $inv=0;
my $shape=0;
my @labelsplit=();
my @paramsplit=();
my $labelline="";
my $modelcount=0;
while(<IN>) {
	$shape=0;
	my $inline=$_;
	chomp $inline;
	if ($inline=~m/Tree/) { #reset everything
		@transitions=(1,1,1,1,1,1);
		$lnL=0;
		$model=$modelNames[$modelcount];
		$K=$modelK[$modelcount];
		$modelcount++;
		@basefreq=(0.25,0.25,0.25,0.25);
		@labelsplit=split(/\t+/,$inline);
		$labelline=$inline;
		@paramsplit=();
	}
	elsif ($inline=~m/1/) {
		@paramsplit=split(/\t+/,$inline);
		$lnL=$paramsplit[1]; #remember perl counts from 0
		if ($labelline=~m/freq/) { #freq was a free parameter
			my $basefreqindex=0;
			for (my $index=0; $index<scalar @labelsplit; $index++) {
				if ($labelsplit[$index]=~m/freq/) {
					$basefreq[$basefreqindex]=$paramsplit[$index];
					$basefreqindex++;
				}
			}
		}
		if ($labelline=~m/ratio/) { #ti/tv was a free parameter
			my $titv="";
			for (my $index=0; $index<scalar @labelsplit; $index++) {
				if ($labelsplit[$index]=~m/ratio/) {
					$titv=$paramsplit[$index];
				}
			}
			$transitions[1]=$titv;
			$transitions[4]=$titv;
		}		
		if ($labelline=~m/R/) { #gtr matrix was a free parameter
			my $transitionindex=0;
			for (my $index=0; $index<scalar @labelsplit; $index++) {
				if ($labelsplit[$index]=~m/R/) {
					$transitions[$transitionindex]=$paramsplit[$index];
					$transitionindex++;
				}
			}
		}	
               if ($labelline=~m/inv/) { #p-inv was a free parameter
			for (my $index=0; $index<scalar @labelsplit; $index++) {
				if ($labelsplit[$index]=~m/inv/) {
					$inv=$paramsplit[$index];
				}
			}
		}
 		if ($labelline=~m/shape/) { #gamma shape was a free parameter
			for (my $index=0; $index<scalar @labelsplit; $index++) {
				if ($labelsplit[$index]=~m/shape/) {
					$shape=$paramsplit[$index];
                                        if ($shape=="infinity"){$shape=300}
				}
			}
		}
                
		print OUT "\n$model\t$lnL\t$K\t$transitions[0]\t$transitions[1]\t$transitions[2]\t$transitions[3]\t$transitions[4]\t$transitions[5]\t$basefreq[0]\t$basefreq[1]\t$basefreq[2]\t$basefreq[3]\t$inv\t$shape";
	}
}
close IN;
close OUT;
if ($modelcount!=scalar @modelNames) {
	system("mv summary.txt WRONG_NAMES_summary.txt");
}
