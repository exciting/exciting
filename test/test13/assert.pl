use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$tol = 0.1;

$writer= Test::initreport("report.xml");

%statusse=Test::assert_file_same_within( "./reference/band_edges.REF", "runldos/band_edges.out", $tol);
Test::writetestreport({
 		"directory"=>"test13/ ",
 		"name"=>"Local Density of States (band_edges)",
 		"description"=>"The test is passed if the computed data differs
 		 less than ".sprintf("%.4f",$tol)." from the reference one:
 		 difference=".sprintf("%.4f",$statusse{"maxerror"}),
 		"status"=>$statusse{"status"}}, $writer);

Test::closereport($writer);
