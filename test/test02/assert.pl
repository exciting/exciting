use lib "../perl/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");
#do things to assert test

open INFO, "runarp/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/\| EXCITING version.+stopped/){
		$status="passed";
	}
}
close INFO;


 Test::writetestreport({
 		directory=>"test02/runarp",
 		name=>"arpackrun",
 		description=>"The test run  using arpack finished without errors",
 		status=>$status
 		}, $writer);
 
 open INFO, "runlapack/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/\| EXCITING version.+stopped/){
		$status="passed";
	}
}
close INFO;


 Test::writetestreport({
 		directory=>"test02/runlapack",
 		name=>"laackrun",
 		description=>"The test run  using lapack finished without errors",
 		status=>$status
 		}, $writer);
 



#compare total energies
$tol=1e-9;
open TOTEARP, "runarp/TOTENERGY.OUT";
open TOTELAP, "runlapack/TOTENERGY.OUT";
$status="failed";
$totearp= <TOTEARP>;
$totelap= <TOTELAP>;
$err=$totearp-$totelap;
if ($totearp-$totelap<=$tol){
$status="passed";
}
 Test::writetestreport({
 		directory=>"test02/ ",
 		name=>"totalenergy compare LAPACK ARPACK ",
 		description=>"The test is passed if the total energy differs\
 		 less than $tol between lapack and Arpack 
 		 difference=$err",
 		status=>$status
 		}, $writer);

#compare eigenvalues
$tol=1e-8;
%statuseigval=Test::assert_file_same_within( "runarp/EIGVAL.OUT",
	"runlapack/EIGVAL.OUT",$tol);
 Test::writetestreport({
 		"directory"=>"test02/ ",
 		"name"=>"eigenvalue comparison LAPACK ARPACK ",
 		"description"=>"The test is passed if the eigenvalues  differ
 		 less than $tol between lapack and Arpack 
 		 difference=" .  %statuseigval->{maxerror},
 		"status"=> %statuseigval->{status}
 		}, $writer);


 Test::closereport($writer);
