#!bin/perl
use date::format;
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
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst)=localtime(time);
$year=$year-100;
$mon=$mon+1;
print VERISIONINC "\#define VERSIONFROMDATE /$year,$mon,0/\n";
close VERISIONINC;
