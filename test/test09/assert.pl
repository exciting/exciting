use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;
use File::Compare;

$writer = Test::initreport("report.xml");

if ( compare( "err", "errref" ) == 0 ) {
	$statusse = "passed";

}
else {
	$statusse = "failed";
}
Test::writetestreport(
	{
		"directory" => "test09/ ",
		"name"      => "Schema validation",
		"description" =>
"The test is passed if input.xsd and all linked schemas are valid schema and validate ../test02/runlapack/input.xml",
		"status" => $statusse
	},
	$writer
);

Test::closereport($writer);
