#!/usr/bin/env perl

# Script copies platform dependent *.inc make files
# from build/platforms to build and adds boolean flags
# for the usage of MPI and ScaLapack.

use strict;
use warnings;

# If nothing was passed 
if (@ARGV == 0) {

  print "---------------------------------------------------------\n";

  # open platforms subdirectory
  opendir(my $PDIR, "build/platforms") or die("Cannot open directory");

  # write filenames of *.inc files into array
  my @makeincfiles= sort(readdir($PDIR));

  # Counter for architectures or platforms
  my $count=1;

  my @fileslist=[];

  # Loop over platform specific make include files
  foreach my $file (@makeincfiles){

    my $platform="";

    # If filename matches the form "make.inc.<arch>"
    # The extension <arch> is saved in $1
    if( $file =~ m/make\.inc\.(.+$)/) {

      $platform=$1;

      print $count." ".$platform."\n";

      $count++;

      # Add matching file to file array
      push(@fileslist,$file);

      
      if ($count%20==0) {
        print "type enter for more";
        my $wait=<>;
      }

    }
  }

  print "\nEnter the number of the platform that suites your system best:  ";

  my $sel=<>;

  if ($sel > $count-1 || $sel < 1 || $sel =~ m/^$/ || $sel !~ m/^\d+$/) {
    print "\ntry again\n\n";
    exit;
  } else {
    print "\nYou use the makefile from:\n\n build/platforms/" . $fileslist[$sel];
    print "\n\nIf the compilation fails, edit \"build/make.inc\" and execute \"make\" again.\n"
  }

  my $filename="build/platforms/" . $fileslist[$sel];

  # Copy selected platform include file to build directory
  my @args=("cp", $filename, "build/make.inc");
  my $return= system(@args);


  # Option selections 


  my $usempi = 0;
  my $usescalapack = 0;

  # Ask for mpi usage

  my $optdone = 0;
  while($optdone == 0){

    print "\nIf you have MPI installed you can build exciting with k-point parallelization support.\n\n";
    print "Build MPI binary ? (yes/No)  ";
    my $MPI=<>;

    if($MPI =~ m/yes/i){

      system("echo \"BUILDMPI = true\">>build/make.inc"); 
      $optdone=1;
      $usempi=1;

    } elsif($MPI =~ m/no/i) {

      system("echo \"BUILDMPI = false\">>build/make.inc");
      $optdone=1;
      $usempi=0;

    } else {
      print "Please choose yes or no";
      $optdone=0;
    }

  }

  # Ask for scalapack usage

  $optdone = 0;
  if($usempi == 1) {

    while($optdone == 0){

      print "\nIf you have ScaLapack installed you can use it in parts of exciting.\n\n";
      print "Build with ScaLapack? (yes/No)  ";
      my $SCL=<>;

      if($SCL =~ m/yes/i){

        system("echo \'LIBS_MPI = \$(LIB_SCLPK) \'>>build/make.inc"); 
        system("echo \'MPIF90_CPP_OPTS += \$(CPP_SCLPK) \'>>build/make.inc"); 
        $optdone=1;
        $usescalapack=1;

      } elsif($SCL =~ m/no/i) {

        $optdone=1;
        $usescalapack=0;
        system("echo \'LIBS_MPI = '' \'>>build/make.inc"); 

      } else {
        print "Please choose yes or no";
        $optdone=0;
      }

    }
  }

}
