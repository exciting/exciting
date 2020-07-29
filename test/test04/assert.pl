use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");
#do things to assert test

#confirm that the calculation has successfully finished
open INFO, "runso/INFO.OUT";
$status=failed;
$didsv=failed;
while(<INFO>)
	{
	if (m/ EXCITING .+stopped/){
		$status="passed";
	}
	if (m/ Moments/){
		$didsv="passed";
	}
        if (m/ SCF iteration number :\s*(\d+)/) {
                $iterations = $1;
        }
        if (m/ Total energy /) {
                ($first, $toten)=split /:/, $_, 2;
                $toten =~ s/^\s+//;
        }
}
close INFO;


 Test::writetestreport({
 		directory=>"test04/runso",
 		name=>"SO run",
 		description=>"The test run with the spin-orbit interaction finished without errors",
 		status=>$status
 		}, $writer);
 
 		
#confirm that the second variation was used 

 Test::writetestreport({
 		directory=>"test04/runso",
 		name=>"second variation",
 		description=>"The second variation was used",
 		status=>$didsv
 		}, $writer);
 

$maxscl=15;
if ($iterations<$maxscl){
$status="passed";
}

 Test::writetestreport({
 		directory=>"test04/runso",
 		name=>"SCF iterations",
 		description=>"The number of SCF is $iterations (should be below $maxscl)",
 		status=>$status
 		}, $writer);

 
#compare total energies
$tol=1e-6;
$refen=-578.073720387;
$status="failed";
$err=$toten-$refen;
if (abs($err)<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test04/runso",
 		name=>"totalenergy",
 		description=>"The test is passed if the total energy differs\
 		 from the reference by less than $tol. 
 		 Obtained difference is $err Ha ($toten Ha vs. $refen Ha).",
 		status=>$status
 		}, $writer);


#compare eigenvalues

open EIGVAL, "runso/EIGVAL.OUT";
while(<EIGVAL>)
	{
	if (m/ 4 /){
                $_ =~ s/^\s+//;
                ($first, $EVAL4, $third)=split / \s/, $_, 3;
	}
	if (m/ 5 /){
                $_ =~ s/^\s+//;
                ($first, $EVAL5, $third)=split / \s/, $_, 3;
	}
	if (m/ 8 /){
                $_ =~ s/^\s+//;
                ($first, $EVAL8, $third)=split / \s/, $_, 3;
	}
	if (m/ 9 /){
                $_ =~ s/^\s+//;
                ($first, $EVAL9, $third)=split / \s/, $_, 3;
                last;
	}
}
close EIGVAL;


#test SO splitting

$status="failed";
$soref=49.7874; #meV
$sosplit=($EVAL5-$EVAL4)*1000*27.211399; #meV
$tol=1e-3;
if (abs($sosplit-$soref)<=$tol) {
    $status="passed";
}

 Test::writetestreport({
 		"directory"=>"test04/runso",
 		"name"=>"SO splitting is right",
 		"description"=>"The test is passed if the SO splitting of VBM ($sosplit meV) differs 
 		 from the reference ($soref meV) by less than $tol meV.",
 		status=> $status
 		}, $writer);


#test degeneracy
$status="failed";
$endiff=($EVAL8-$EVAL5); #Ha
$tol=1e-8;
if ($endiff<=$tol) {
    $status="passed";
}

 Test::writetestreport({
 		"directory"=>"test04/runso",
 		"name"=>"Test degeneracy of VBM",
 		"description"=>"The test is passed if VBM is degenerate. 
                 The difference in energies of bands 5 and 8 at Gamma is $endiff Ha 
                 compared to the allowed difference of $tol Ha.",
 		status=> $status
 		}, $writer);


#Gamma-Gamma gap
$status="failed";
$gapref=2.46929100; #eV
$gap=($EVAL9-$EVAL8)*27.211399; #eV
$tol=1e-6;
if ($endiff<=$tol) {
    $status="passed";
}

 Test::writetestreport({
 		"directory"=>"test04/runso",
 		"name"=>"Gamma-Gamma gap",
 		"description"=>"The test is passed if the Gamma-Gamma gap ($gap eV) 
                agrees with the reference ($gapref eV) within $tol eV.",
 		status=> $status
 		}, $writer);

 Test::closereport($writer);

