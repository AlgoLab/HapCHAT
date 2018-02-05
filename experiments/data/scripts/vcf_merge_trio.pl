#!/usr/bin/env perl
$mother_filename=$ARGV[0];
$father_filename=$ARGV[1];
$child_filename=$ARGV[2];

open(FP,$child_filename);
my %hash;
my %hash1;
my %hash2;
my @array;
while($line=<FP>)
{
	if($line!~/^\#/)
	{
		$line=~s/\n//g;
		@var=split(" ",$line);
		@var1=split("\:",$var[-1]);
		if($var[4]!~/\,/ && $var[3]!~/\,/ && length($var[4])>=1 && length($var[3])>=1) 
		{
			$hash{"$var[0]_$var[1]_$var[3]_$var[4]"} = "$var[0] $var[1] $var[2] $var[3] $var[4] $var[5] $var[6] $var[7] GT $var1[0]";
		}	
	}
}

open(FP1,$father_filename);
while($line=<FP1>)
{
        if($line!~/^\#/)
        {
		$line=~s/\n//g;
		@var=split(" ",$line);
		@var1=split("\:",$var[-1]);
		$hash1{"$var[0]_$var[1]_$var[3]_$var[4]"} = $var1[0];	
	}
}

open(FP2,$mother_filename);
while($line=<FP2>)
{
        if($line!~/^\#/)
        {
		$line=~s/\n//g;
		@var=split(" ",$line);
		@var1=split("\:",$var[-1]);
		$hash2{"$var[0]_$var[1]_$var[3]_$var[4]"} = $var1[0];
        }
}
close(FP);
$count=0;
open(FP,$child_filename);
while($line=<FP>){
	if($line=~/^\#/)
	{
		if($count==314)
		{
			$line=~s/\n/	HG003	HG004\n/g;
		}
		print $line;
	}
	else
	{
		@var=split(" ",$line);
		if($var[4]!~/\,/ && $var[3]!~/\,/ && length($var[4])>=1 && length($var[3])>=1){
			if($hash1{"$var[0]_$var[1]_$var[3]_$var[4]"}!="" && $hash2{"$var[0]_$var[1]_$var[3]_$var[4]"}!="")
			{
				push(@array,$hash{"$var[0]_$var[1]_$var[3]_$var[4]"}." ".$hash1{"$var[0]_$var[1]_$var[3]_$var[4]"}." ".$hash2{"$var[0]_$var[1]_$var[3]_$var[4]"}."\n");
			}
		}
	}
	$count++;
}
%seen=();
@unique = grep { ! $seen{$_} ++ } @array;
for($i=0;$i<@array;$i++)
{
	print $unique[$i];
}
close(FP);
