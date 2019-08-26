use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$tol = 1.0;

$writer= Test::initreport("report.xml");

%statusse=Test::assert_file_same_within( "./reference/EVALQP.REF", "rungw/EVALQP.DAT", $tol);

Test::writetestreport({
 		"directory"=>"test08/ ",
 		"name"=>"GW quasiparticle banstructure",
 		"description"=>"The test is passed if the quasiparticle banstructure differs
 		 less than ".sprintf("%.2f",$tol)." between rungw/EVALQP.OUT and reference file 
 		 reference/EVALQP.REF difference=".sprintf("%.4f",$statusse{"maxerror"}),
 		"status"=>%statusse->{"status"}}, $writer);

Test::closereport($writer);
