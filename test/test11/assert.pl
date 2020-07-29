use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$tol = 1.0;

$writer= Test::initreport("report.xml");

%statusse=Test::assert_file_same_within("./reference/EIGVAL-pbe0.REF", "run-pbe0/EIGVAL.OUT", $tol);

Test::writetestreport({
    "directory"=>"test11/",
    "name"=>"Hybrid functional: PBE0",
    "description"=>"The test is passed if the eigenvalues differ
     less than ".sprintf("%.2f",$tol)." between run-pbe0/EIGVAL.OUT and reference file
     reference/EIGVAL-pbe0.REF difference = ".sprintf("%.4f",$statusse{"maxerror"}),
    "status"=>$statusse{"status"}}, $writer);

%statusse=Test::assert_file_same_within("./reference/EIGVAL-hse.REF", "run-hse/EIGVAL.OUT", $tol);

Test::writetestreport({
    "directory"=>"test11/",
    "name"=>"Hybrid functional: HSE",
    "description"=>"The test is passed if the eigenvalues differ
     less than ".sprintf("%.2f",$tol)." between run-hse/EIGVAL.OUT and reference file 
     reference/EIGVAL-hse.REF difference = ".sprintf("%.4f",$statusse{"maxerror"}),
    "status"=>$statusse{"status"}}, $writer);

 Test::closereport($writer);
