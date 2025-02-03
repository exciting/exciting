#!/usr/bin/env perl

# Script generates src/version.inc from git commit hash, compiler version
# and the current date.
# In version.inc GITHASH, GITHASH2, COMPILERVERSION and VERSIONFROMDATE
# get defined.

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

# For a given string containing the name of the compiler:
# call the shell to get the compiler version.
#
# Recognised compiler strings: 'ifort' or 'gcc'
#
# Returns string $compiler_version
sub GetCompiler {
    # Get compiler command F90
    my ($F90_line) = @_;
    my @F90_line = split(' ', $F90_line);
    my $F90 = $F90_line[2];

    # Find compiler version
    # This isn't declaring a string,
    # it's a call to the shell, which returns to $stdout
    my $stdout = `$F90 --version`;
    my @arr = split("\n", $stdout);
    my $compiler_version = $arr[0];

    # Check if compiler is supported.
    my $compiler_is_valid
        = (index($compiler_version, 'ifort') != -1) 
        + (index($compiler_version, 'GNU') != -1);

    if ( $compiler_is_valid == 0){
      $compiler_version = "";
    }

    return $compiler_version
}

# Match 'F90' in make.inc and use that to establish the compiler type
open ( my $make_inc, '<', '../make.inc' ) or die $!;

my $compiler = "";
while (my $line = <$make_inc> ) {
   if($line =~ /^F90 =/) {
    $compiler = GetCompiler($line);
   }
    # Stop checking file if compiler has been matched
    last if ($compiler ne "");
}

close ( $make_inc );

# Compiler choice is not recognised
if ($compiler eq "") {
    print "Compiler type could not be established from build/make.inc \n",
        "Please ensure FC = 'ifort' or 'gfortran', or update ",
        "build/utilities/versionstamp.pl \n";
    exit
}


# Write version number to file
open my $VERISIONINC, ">", "../../src/version.inc";
print $VERISIONINC "\#define GITHASH \"$hash1\"\n";
print $VERISIONINC "\#define GITHASH2 \"$hash2\"\n";
print $VERISIONINC "\#define COMPILERVERSION \"$compiler\"\n";
print $VERISIONINC "\#define VERSIONFROMDATE ", `date "+/%y,%m,%d/"`;
close $VERISIONINC;
