use lib "../test/perl/";
use lib "../test/perl/lib";
use lib "../build/utilities/lib";
use XML::Simple;
use XML::Writer;
use IO::File;
use List::Util qw[min max];
use Data::Dumper;
use Text::Wrap;

foreach $sourcedir (@ARGV) {
	if ( opendir( DIR, "$sourcedir" ) ) {
		@filelist = readdir(DIR);
		closedir DIR;

		foreach $file (@filelist) {
			if ( $file =~ m/.+\.[fF]90$/ ) {
				if ( not( $file =~ m/modtetra/ ) ) {

					push( @fileslist, "$sourcedir$file" );
				}
			}
		}
	}
	elsif ( open( FILE, $sourcedir ) ) {
		push( @fileslist, $sourcedir );
		close FILE;
	}
	else {
		die "canot open file $/ $! $sourcedir";
	}
}

#open xmlvaldb
# create object
$xml = new XML::Simple;

# read XML file
system('xsltproc attributelist.xsl excitinginput.xsd >data.xml');
$data = $xml->XMLin("data.xml");

# print output
# print Dumper($data);

FILE:
foreach $file (@fileslist) {

	local $/ = undef;
	print "open file ", $file, "\n";
	open FILE, $file or die "Couldn't open file: $!";
	binmode FILE;
	$string = <FILE>;

  VAR:
	for my $att ( @{ $data->{attribute} } ) {
		if ( not( $att->{xpath} =~ m/^input/ ) ) { next VAR; }
		if ( $string =~ s/(.*::(.*[\s,])?\s*$att->{vname}(\(|\s|,|\n))/!replaced by inputstructure$1/g ) {
			print $att->{vname}, "\n";
		}
	}

	open FILE, ">", $file or die "Couldn't open file: $!";
	print FILE $string, "\n";
	close FILE;
}
