use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");

%statusse=Test::assert_file_same_within( "./reference/EIGVAL-pbe0-hf.REF",
	"run-pbe0-hf/EIGVAL.OUT",1.0);
 Test::writetestreport({
 		"directory"=>"test11/ ",
 		"name"=>"Hybrid functional: PBE0-HF",
 		"description"=>"The test is passed if the quasiparticle banstructure differs
 		 less than $tol between run-pbe0-hf/EIGVAL.OUT and reference file 
 		 reference/EIGVAL-pbe0-hf.REF difference=".  %statusse->{maxerror},
 		"status"=> %statusse->{status}}, $writer);

%statusse=Test::assert_file_same_within( "./reference/EIGVAL-pbe0-oep.REF",
  "run-pbe0-oep/EIGVAL.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test11/ ",
    "name"=>"Hybrid functional: PBE0-OEP",
    "description"=>"The test is passed if the quasiparticle banstructure differs
     less than $tol between run-pbe0-oep/EIGVAL.OUT and reference file 
     reference/EIGVAL-pbe0-oep.REF difference=".  %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

 Test::closereport($writer);
