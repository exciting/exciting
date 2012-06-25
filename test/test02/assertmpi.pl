use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("reportmpi.xml");
#do things to assert test

open INFO, "runlapackmpi/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/\| EXCITING.+stopped/){
		$status="passed";
	}
}
close INFO;


 Test::writetestreport({
 		directory=>"test02/runlapackmpi",
 		name=>"mpi run",
 		description=>"The test run  using arpack finished without errors",
 		status=>$status
 		}, $writer);
 
 open INFO, "runlapack/INFO.OUT";
$status=failed;



#compare total energies
$tol=3e-8;
open TOTEARP, "runlapackmpi/TOTENERGY.OUT";
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
 		name=>"totalenergy_compare_LAPACK_LAPACKMPI ",
 		description=>"The test is passed if the total energy differs\
 		 less than $tol between lapack and MPI run  
 		 difference=$err",
 		status=>$status
 		}, $writer);

#compare eigenvalues
$tol=1.e-5;
%statuseigval=Test::assert_file_same_within( "runlapackmpi/EIGVAL.OUT",
	"runlapack/EIGVAL.OUT",$tol);
 Test::writetestreport({
 		"directory"=>"test02/ ",
 		"name"=>"eigenvalue_comparison_LAPACK ,MPI ",
 		"description"=>"The test is passed if the eigenvalues  differ
 		 less than $tol between lapack and MPI run 
 		 difference=" .  %statuseigval->{maxerror},
 		"status"=> %statuseigval->{status}
 		}, $writer);

 Test::closereport($writer);
