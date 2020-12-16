use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");

$tol=1.0;

%statusse=Test::assert_file_same_within( "./reference/spintext.xml",
 "run/spintext.xml",$tol);
Test::writetestreport({                                                   
   "directory"=>"test14/ ",
   "name"=>"Properties: spintexture",
   "description"=>"The test is passed if the output file differs
    less than $tol from the reference one ". $statusse{maxerror},
   "status"=> $statusse{status}}, $writer); 

Test::closereport($writer);
                        
