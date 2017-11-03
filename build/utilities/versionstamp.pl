#!/usr/bin/env perl

# Script generates src/version.inc from git commit hash and
# current date. In version.inc GITHASH, GITHASH2 and VERSIONFROMDATE 
# are defined.

use strict;
use warnings;

# Read git branch name from git HEAD (usual content of HEAD is "ref: refs/heads/master")
my $ref = '';
open my $HEADREF, "../../.git/HEAD ";
while (<$HEADREF>) {
  if (m/ref:\s*([\w\/-]*)/) {
    $ref = $1;
  }
}
close $HEADREF;


my $hash1 = '';
my $hash2 = '';

# if file exists (i.e. ../../.git/refs/heads/master)
if ( -e  "../../.git/" . $ref ) {

  # read hash of current commit and print it

  open my $HEAD, "../../.git/" . $ref;
  read $HEAD, my $hasht, 40;

  $hash1 = substr( $hasht, 0,  20 );
  $hash2 = substr( $hasht, 20, 39 );

  close $HEAD;

} else {

  open my $REFS ,"../../.git/packed-refs";

  while(<$REFS>){
    if(m/$ref/) {
      $hash1 = substr( $_, 0,  20 );
      $hash2 = substr( $_, 20, 20 );
    }
  }

  close $REFS;
}

# Write version number to file

open my $VERISIONINC, ">", "../../src/version.inc";
print $VERISIONINC "\#define GITHASH \"$hash1\"\n";
print $VERISIONINC "\#define GITHASH2 \"$hash2\"\n";

print $VERISIONINC "\#define VERSIONFROMDATE ", `date "+/%y,%m,%d/"`;

close $VERISIONINC;
