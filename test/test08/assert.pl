use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");




%statusse=Test::assert_file_same_within( "./reference/EVALQP.REF",
	"rungw/EVALQP.DAT",1.0);
 Test::writetestreport({
 		"directory"=>"test08/ ",
 		"name"=>"GW quasiparticle banstructure",
 		"description"=>"The test is passed if the quasiparticle banstructure differs
 		 less than $tol between runlgw/EVALQP.OUT and reference file 
 		 reference/EVALQP.REF difference=".  %statusse->{maxerror},
 		"status"=> %statusse->{status}}, $writer);

 Test::closereport($writer);
