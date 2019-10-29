use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");

# Test 1
open INFO, "runlapack/INFO.OUT";
$status=failed;
while(<INFO>)
	{
        if (m/\| EXCITING.+stopped/){
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

# Test 2
open INFO, "runarp/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/\| EXCITING.+stopped/){
		$status="passed";
	}
}
close INFO;

Test::writetestreport({
 		directory=>"test02/runarp",
 		name=>"arpack run",
 		description=>"The test ARPACK finished without errors",
 		status=>$status
 		}, $writer);

#compare total energies
# $tol=1e-2;
# open TOTEARP, "runarp/TOTENERGY.OUT";
# open TOTELAP, "runlapack/TOTENERGY.OUT";
# $status="failed";
# $totearp= <TOTEARP>;
# $totelap= <TOTELAP>;
# $err=$totearp-$totelap;
# if (abs($totearp-$totelap)<=$tol){
# $status="passed";
# }
#  Test::writetestreport({
#  		directory=>"test02/ ",
#  		name=>"totalenergy_compare_LAPACK_ARPACK ",
#  		description=>"The test is passed if the total energy differs\
#  		 less than $tol between lapack and Arpack
#  		 difference=$err",
#  		status=>$status
#  		}, $writer);

$tol = 0.0001;

%statusse=Test::assert_file_same_within( "./reference/EIGVAL.REF", "runlapack/EIGVAL.OUT", $tol);

Test::writetestreport({
 		"directory"=>"test02/ ",
 		"name"=>"Eigenvalues are the same as in the reference_file",
 		"description"=>"The test is passed if the eigenvalues differ
 		 less than ".sprintf("%.6f",$tol)." between EIGVAL.OUT and reference file.
 		 maxerror=".sprintf("%.6f",$statusse{"maxerror"}),
 		"status"=>$statusse{"status"}}, $writer);


 Test::closereport($writer);
