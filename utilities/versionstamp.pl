#!bin/perl

open HEADREF, "../../.git/HEAD ";
while(<HEADREF>)
	{
	if(m/ref:\s*([\w\/]*)/){
	$ref=$1;	
	
	}

}
close HEADREF;



open HEAD, "../../.git/".$ref;
read HEAD, $hash,20 ;
close HEAD;

open VERISIONINC,">", "../../src/version.inc";
print VERISIONINC "\#define GITHASH \"$hash\"\n";
close VERISIONINC;
