use XML::eXistDB;
use XML::eXistDB::RPC;
use Digest::MD5 qw(md5 md5_hex md5_base64);

my ($db) = XML::eXistDB::RPC->new(
	user        => "chm",
	password    => "test",
	destination => "http://localhost:8080/exist/xmlrpc",
	chunk_size=>1000000
);
my ( $rc, $exists ) = $db->hasCollection("calculations");

use strict;
my @inputs;
my @calcdir;

sub recursefiles($) {
	my ($path) = @_;

	## append a trailing / if it's not there
	$path .= '/' if ( $path !~ /\/$/ );

	## print the directory being searched
	##print $path,"\n";

	## loop through the files contained in the directory
	for my $eachFile ( glob( $path . '*' ) ) {

		## if the file is a directory
		if ( -d $eachFile ) {
			## pass the directory to the routine ( recursion )
			recursefiles($eachFile);
		}
		else {

			## print the file ... tabbed for readability
			if ( $eachFile =~ m/input.xml/ ) {
				push( @inputs,  $eachFile );
				push( @calcdir, $path );
			}
		}
	}
}

## initial call ... $ARGV[0] is the first command line argument
recursefiles( $ARGV[0] );

foreach (@calcdir) {
	my ($dir) = $_;
	open( FILE, $dir . "input.xml" );

	#print $dir. "input.xml\n";
	local $/ = undef;
	my ($data) = <FILE>;
	close(FILE);
	my ($digest) = md5_hex($data);
	print($digest);
	print "\n";

	my ( $rc, my $ok ) = $db->createCollection( "calculations/" . $digest );
	opendir( DIR, "$dir" ) || die("Cannot open directory $!\n");
	my (@filelist) = readdir(DIR);
	closedir DIR;
	for my $eachFile (@filelist) {
		if ( $eachFile =~ m/.*.xml/ ) {

			#print $eachFile. "\n";

			open( FILE, $dir . $eachFile );
			local $/ = undef;
			my ($xmldata) = <FILE>;
			close(FILE);
			print "uploading "
			  . "calculations/"
			  . $digest . "/"
			  . $eachFile . "\n";
			my ( $rc, my $ok ) =
			  $db->uploadDocument( "calculations/" . $digest . "/" . $eachFile,
				$xmldata,   replace => "<true>" )
			  ;
		}
	}
}
