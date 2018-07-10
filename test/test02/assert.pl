use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");
#do things to assert test

open INFO, "runarp/INFO.OUT";
$status=failed;
$didarpack=failed;
while(<INFO>)
	{
#	if (m/\| EXCITING .+stopped/){
	if (m/ EXCITING .+stopped/){
		$status="passed";
	}
		if (m/ARPACK iterations/){
		$didarpack="passed";
	}
}
close INFO;


 Test::writetestreport({
 		directory=>"test02/runarp",
 		name=>"arpack run",
 		description=>"The test run  using arpack finished without errors",
 		status=>$status
 		}, $writer);
 
  Test::writetestreport({
 		directory=>"test02/runarp",
 		name=>"arpack run",
 		description=>"The arpack solver was correctly invoked",
 		status=>$didarpack
 		}, $writer);
 		
 open INFO, "runlapack/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/ EXCITING .+stopped/){
		$status="passed";
	}
}
close INFO;


 Test::writetestreport({
 		directory=>"test02/runlapack",
 		name=>"lapack_run",
 		description=>"The test run  using lapack finished without errors",
 		status=>$status
 		}, $writer);
 



#compare total energies
$tol=1e-6;
open TOTEARP, "runarp/TOTENERGY.OUT";
open TOTELAP, "runlapack/TOTENERGY.OUT";
$status="failed";
$totearp= <TOTEARP>;
$totelap= <TOTELAP>;
$err=$totearp-$totelap;
if (abs($totearp-$totelap)<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test02/ ",
 		name=>"totalenergy_compare_LAPACK_ARPACK ",
 		description=>"The test is passed if the total energy differs\
 		 less than $tol between lapack and Arpack 
 		 difference=$err",
 		status=>$status
 		}, $writer);

#compare eigenvalues
$tol=1e-5;
 

%statuseigvalr=Test::assert_file_same_within( "./reference/EIGVAL.REF",
	"runlapack/EIGVAL.OUT",$tol);
 Test::writetestreport({
 		"directory"=>"test02/ ",
 		"name"=>"Eigenvalues_same_as_reference_file",
 		"description"=>"The test is passed if the eigenvalues  differ
 		 less than $tol between runlapack/EIGVAL.OUT and reference file reference/EIGVAL.REF
 		 difference=" .  %statuseigvalr->{maxerror},
 		"status"=> %statuseigvalr->{status}
 		}, $writer);

 Test::closereport($writer);
