use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");
#do things to assert test

#confirm that the calculation has successfully finished
open INFO, "runcollinear/INFO.OUT";
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

        if (m/ total moment /) {
                ($first, $totmom)=split /:/, $_, 2;
                $totmom =~ s/^\s+//;
        }

}
close INFO;

 Test::writetestreport({
 		directory=>"test12/runcollinear",
 		name=>"SO run",
 		description=>"The test run with collinear magnetism finished without errors",
 		status=>$status
 		}, $writer);
 
 		
#confirm that the second variation was used 

 Test::writetestreport({
 		directory=>"test12/runcollinear",
 		name=>"second variation",
 		description=>"The second variation was used",
 		status=>$didsv
 		}, $writer);
 

$maxscl=25;
if ($iterations<$maxscl){
$status="passed";
}

 Test::writetestreport({
 		directory=>"test12/runcollinear",
 		name=>"SCF iterations",
 		description=>"The number of SCF is $iterations (should be below $maxscl)",
 		status=>$status
 		}, $writer);

 
#compare total energies
$tol=1e-6;
$refen=-1270.57692665;
$status="failed";
$err=$toten-$refen;
if (abs($err)<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test12/runcollinear",
 		name=>"totalenergy",
 		description=>"The test is passed if the total energy differs\
 		 from the reference by less than $tol. 
 		 Obtained difference is $err Ha ($toten Ha vs. $refen Ha).",
 		status=>$status
 		}, $writer);

#compare magnetic moments
$tol=3e-5;
$refmom=-2.00000003;
$status="failed";
$err=$totmom-$refmom;
if (abs($err)<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test12/runcollinear",
 		name=>"magnetic moments",
 		description=>"The test is passed if the total magnetic moment differs\
 		 from the reference by less than $tol. 
 		 Obtained difference is $err ($totmom vs. $refmom).",
 		status=>$status
 		}, $writer);


#confirm that the calculation has successfully finished
open INFO, "runnoncollinear/INFO.OUT";
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

        if (m/ total moment /) {
                ($first, $totmom)=split /:/, $_, 2;
                $totmom =~ s/^\s+//;
                ($totmom, $second, $third)=split /   /, $totmom, 3;
                $totmom =~ s/^\s+//;
        }

}
close INFO;

 Test::writetestreport({
 		directory=>"test12/runnoncollinear",
 		name=>"SO run",
 		description=>"The test run with magnetism and the spin-orbit interaction finished without errors",
 		status=>$status
 		}, $writer);
 
 		
#confirm that the second variation was used 

 Test::writetestreport({
 		directory=>"test12/nonruncollinear",
 		name=>"second variation",
 		description=>"The second variation was used",
 		status=>$didsv
 		}, $writer);
 

$maxscl=25;
if ($iterations<$maxscl){
$status="passed";
}

 Test::writetestreport({
 		directory=>"test12/nonruncollinear",
 		name=>"SCF iterations",
 		description=>"The number of SCF is $iterations (should be below $maxscl)",
 		status=>$status
 		}, $writer);

 
#compare total energies
$tol=3e-5;
$refen=-1270.57748365;
$status="failed";
$err=$toten-$refen;
if (abs($err)<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test12/nonruncollinear",
 		name=>"totalenergy",
 		description=>"The test is passed if the total energy differs\
 		 from the reference by less than $tol. 
 		 Obtained difference is $err Ha ($toten Ha vs. $refen Ha).",
 		status=>$status
 		}, $writer);

#compare magnetic moments
$tol=1e-6;
$refmom=-2.00181244;
$status="failed";
$err=$totmom-$refmom;
if (abs($err)<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test12/nonruncollinear",
 		name=>"magnetic moment",
 		description=>"The test is passed if the total magnetic moment differs\
 		 from the reference by less than $tol. 
 		 Obtained difference is $err ($totmom vs. $refmom).",
 		status=>$status
 		}, $writer);

 Test::closereport($writer);

