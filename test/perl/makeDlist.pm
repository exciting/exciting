package Summary;
use Data::Dumper;
use XML::Writer;

sub makeDlist {
	my ( $dir, $file );
	 
	opendir TESTROOT, "./";
	my @thefiles = readdir(TESTROOT);
	closedir(TESTROOT);
	my $output = new IO::File(">directories.xml");
	my $writer = new XML::Writer(
		OUTPUT      => $output,
		DATA_MODE   => 'true',
		DATA_INDENT => 2
	);
	$writer->xmlDecl('UTF-8');
	$writer->startTag("testdirectories");
	$writer->emptyTag( "githash", hash => get_git_hash("../src/version.inc") );

	foreach $dir (@thefiles) {

		if ( $dir =~ m/test.*/ ) {
			$writer->startTag("dir");
			$writer->characters($dir);
			$writer->endTag("dir");
		}
	}
	$writer->endTag("testdirectories");
	$output = $writer->getOutput();
	$output->close();
}

sub get_git_hash {

	open VERS, $_[0];
	while (<VERS>) {
		if (m/GITHASH\s+\"([\w|\d]+)\"/) {
			return $1;
		}
		else {
			return "no git";
		}

	}
}

return (1);
