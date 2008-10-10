use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");
#do things to assert test

open INFO, "run/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/\| EXCITING version.+stopped/){
		$status="passed";
	}
	if (m/Iteration number :   (\d+)/){
	$iterations=$1;
}
}
close INFO;


 Test::writetestreport({
 		directory=>"test03/run",
 		name=>"debug binary run",
 		description=>"The test run  using arpack finished without errors",
 		status=>$status
 		}, $writer);
 		
 		$iterationsref=13;
 		if ($iterations=13){$passed="true";}else{$passed="false";}
 		
 		Test::writetestreport({
 		directory=>"test03/run",
 		name=>"check iteration number of Al",
 		description=>"iterations is $iterations reference is $iterationsref",
 		status=>$status
 		}, $writer);
 		
 		
 		open INFO, "runmixer2/INFO.OUT";
$status=failed;
while(<INFO>)
	{
	if (m/\| EXCITING version.+stopped/){
		$status="passed";
	}
	if (m/Iteration number :   (\d+)/){
	$iterations=$1;
}
}
Test::writetestreport({
 		directory=>"test03/runmixer2",
 		name=>"debug binary run mixer2",
 		description=>"The test run  using arpack finished without errors",
 		status=>$status
 		}, $writer);
 		
 		$iterationsref=13;
 		if ($iterations=13){$passed="true";}else{$passed="false";}
 		
 		Test::writetestreport({
 		directory=>"test03/runmixer2",
 		name=>"check iteration number of Al for mixer2",
 		description=>"iterations is $iterations reference is $iterationsref",
 		status=>"unspecified"
 		}, $writer);
 		
 Test::closereport($writer);
