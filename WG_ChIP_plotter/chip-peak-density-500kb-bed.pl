#!/usr/bin/perl 
use warnings;        
use strict;
#for sliding window=500KB
#usage: perl program input output
my ($in,$out)=@ARGV;
open IN, "<$in" or die($!);
open OU, ">$out" or die ($!); # Output goes exactly where I tell it to. This is very common $! die comment tells the program to shut off if it cant do this!.

my %value;
my (%start,%end);   # declared hashes but where are the keys assigned and how is the data split up?
my %fc;			# never used fold change
while(<IN>){     # data split up here, while there is an input file
	next if !/^C/;   # go to the next line if the line doesnt start with an upper case c (prevent wierdness in the beginning)  ! (not) /^C/ match C
	chomp;                                    	 #  removes all newlines /n the result is a flat sentence with no returns
	my ($chr,$start,$end,$fc)=(split/\t/)[0,1,2,4];   # 4 columns chromosome, start, end, fc  split by a tab 1st, 2nd, 3rd, 5th columns, I believe this is for the data structure...
	$value{$chr}{$start}=($end-$start)+1;      # setting up the hash {chr} is the keys for {start}
	$end{$chr}{$start}=$end;		# setting up hash end {chr} is again key for start}
	$fc{$chr}{$start}=$fc;			# never used fc fold change hash, its set up by keys {$chr}{$start} = $fc
}
# so what are we doing here ....

for my $chr(sort {$a cmp $b}keys %value){      # sort, compare a$ vs b$ smallest to largest logically alphabetically, compare $b vs $a largest to smallest, compare $a <=> $b  only if a and b are numbers
	my $len=
		($chr eq "Chr1")?60:    # loop declaring length of sliding windows chromosome 1 is 30 MB / 500 kb = 60 (sliding windows)?
		($chr eq "Chr2")?38:
		($chr eq "Chr3")?45:
		($chr eq "Chr4")?36:
		($chr eq "Chr5")?52:
	my $sum;
	for my $i(0..$len){     #  for my $i from o to length 60  
		for my $n(sort {$a<=>$b}keys %{$value{$chr}}){  			# $n is another variable, sorting $n large to small by keys 
			if (($n < $i*500000) && ($end{$chr}{$n} > $i*500000) && ($end{$chr}{$n} < ($i+1)*500000)){
				$sum+=($end{$chr}{$n}-($i*500000)+1)/$value{$chr}{$n};   # $sum+=   means sum plus whatever value given the above condition 
			}elsif($n >= $i*500000 && $end{$chr}{$n} <= ($i+1)*500000){
				$sum++;   						 #$sum++ = $sum + 1 in human language
			}elsif(($n > $i*500000) && ($n < ($i+1)*500000) && ($end{$chr}{$n} > ($i+1)*500000)){
				$sum+=((($i+1)*500000)-$n+1)/$value{$chr}{$n};    
			}
		} 			# this is just a loop && evalutes the truth of the statement then does something conditionally, in this case the && requires that all these conditions be met
		printf OU "$chr\t%.4f\n",($sum);   # print output $chr tab (answer to $sum) to 4 floating decimal places.
		$sum=0;
	}
}

close IN;
close OU;
