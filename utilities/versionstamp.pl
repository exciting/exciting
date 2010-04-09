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
read HEAD, $hasht,40 ;

$hash1=substr($hasht, 0, 20);
$hash2=substr($hasht, 20, 39);

close HEAD;

open VERISIONINC,">", "../../src/version.inc";
print VERISIONINC "\#define GITHASH \"$hash1\"\n";
print VERISIONINC "\#define GITHASH2 \"$hash2\"\n";

print VERISIONINC "\#define VERSIONFROMDATE ", `date "+/%y,%m,0/"`