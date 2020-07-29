use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");

$tol = 1.0;
%statusse=Test::assert_file_same_within( "./reference/ldos.REF", "runldos/ldos.out", $tol);
Test::writetestreport({
 		"directory"=>"test13/ ",
 		"name"=>"Local Density of States (ldos)",
 		"description"=>"The test is passed if the computed data differs
 		 less than ".sprintf("%.4f",$tol)." from the reference one:
 		 difference=".sprintf("%.4f",$statusse{"maxerror"}),
 		"status"=>$statusse{"status"}}, $writer);

$tol = 0.1;
%statusse=Test::assert_file_same_within( "./reference/band_edges.REF", "runldos/band_edges.out", $tol);
Test::writetestreport({
 		"directory"=>"test13/ ",
 		"name"=>"Local Density of States (band_edges)",
 		"description"=>"The test is passed if the computed data differs
 		 less than ".sprintf("%.4f",$tol)." from the reference one:
 		 difference=".sprintf("%.4f",$statusse{"maxerror"}),
 		"status"=>$statusse{"status"}}, $writer);

Test::closereport($writer);
