use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer = Test::initreport("report.xml");

#do things to assert test

open INFO, "run/INFO.OUT";
$status     = 'failed';
$iterations = 0;
$rightmixer = 'failed';
while (<INFO>) {
	if (m/\| EXCITING.+stopped/) {
		$status = "passed";
	}
	if (m/ SCF iteration number :\s*(\d+)/) {
		$iterations = $1;
	}
	if (m/Using adaptive step size linear potential mixing/) {
		$rightmixer = 'passed';
	}
}
close INFO;
Test::writetestreport(
	{
		directory   => "test03/run",
		name        => "is_it_linadapt",
		description => "linadapt mixer was invoked",
		status      => $rightmixer
	},
	$writer
);

Test::writetestreport(
	{
		directory   => "test03/run",
		name        => "debug binary run",
		description => "The test run standard mixing finished without errors",
		status      => $status
	},
	$writer
);

$iterationsref = 38;
if   ( $iterations < $iterationsref+5 &&  
	 $iterations > $iterationsref-5) { $status = "passed"; }
else                                   { $status = "failed"; }

Test::writetestreport(
	{
		directory   => "test03/run",
		name        => "check iteration number of Al",
		description => "iterations is $iterations reference is $iterationsref",
		status      => $status
	},
	$writer
);

open INFO, "runmixer2/INFO.OUT";
$status     = "failed";
$iterations = 0;
$rightmixer = 'failed';
while (<INFO>) {
	if (m/\| EXCITING .+stopped/) {
		$status = "passed";
	}
	if (m/Iteration number :\s*(\d+)/) {
		$iterations = $1;
	}
	if (m/Using multisecant Broyden potential mixing/) {
		$rightmixer = 'passed';
	}
}
Test::writetestreport(
	{
		directory   => "test03/runmixer2",
		name        => "is_it_msec",
		description => "used msec",
		status      => $rightmixer
	},
	$writer
);
Test::writetestreport(
	{
		directory => "test03/runmixer2",
		name      => "debug binary run mixer2",
		description =>
		  "The test run  using multicecant broyden finished without errors",
		status => $status
	},
	$writer
);

$iterationsref = 12;
if   ( $iterations <= $iterationsref ) { $status = 'passed'; }
else                                   { $status = "failed"; }

Test::writetestreport(
	{
		directory   => "test03/runmixer2",
		name        => "check iteration number of Al for mixer2",
		description => "iterations is $iterations reference is $iterationsref",
		status      => $status
	},
	$writer
);

open INFO, "runmixer3/INFO.OUT";
$status     = 'failed';
$iterations = 0;
$rightmixer = 'failed';
while (<INFO>) {
	if (m/\| EXCITING .+stopped/) {
		$status = "passed";
	}
	if (m/Iteration number :\s*(\d+)/) {
		$iterations = $1;
	}

	if (m/Using Pulay potential/) {
		$rightmixer = 'passed';
	}
}

Test::writetestreport(
	{
		directory   => "test03/runmixer3",
		name        => "is_it_Pulay",
		description => "Pulay Mixing is used",
		status      => $rightmixer
	},
	$writer
);

Test::writetestreport(
	{
		directory   => "test03/runmixer3",
		name        => "Pulay mixing works",
		description => "The test run  using pulay (3) finished without errors",
		status      => $status
	},
	$writer
);

$iterationsref = 12;
if   ( $iterations <= $iterationsref ) { $status = 'passed'; }
else                                   { $status = "failed"; }

Test::writetestreport(
	{
		directory   => "test03/runmixer2",
		name        => "check iteration number of Al for Pulay mixing (3) ",
		description => "iterations is $iterations reference is $iterationsref",
		status      => $status
	},
	$writer
);

Test::closereport($writer);
