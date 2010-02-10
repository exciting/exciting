use XML::eXistDB;
use XML::eXistDB::RPC;
use XML::Simple;
use Data::Dumper;
use Digest::MD5 qw(md5 md5_hex md5_base64);
use Term::ReadPassword;
print "Enter DB user: ";
$user = <STDIN>;
$password = read_password('password: ');
$user  =~s/\n//;
my ($db) = XML::eXistDB::RPC->new(
	user        => $user,
	password    => $password,
	destination => "http://g44228:8080/exist/xmlrpc",
	chunk_size  => 1000000
);
my ( $rc, $exists ) = $db->hasCollection("calculations");
$rc == 0 or die "ERROR:cannot connect': $exists ($rc)\n";

use strict;
my @inputs;
my @calcdir;
my @filelist;
my $dir;

my $setfile = $ARGV[0];
my $xms     = XML::Simple->new();
my $doc     = $xms->XMLin($setfile);

foreach my $set ( @{ $doc->{set} } ) {
	my ($dir) = $set->{path};

	open( FILE, $dir . "input.xml" );

	print $dir. "input.xml\n";
	local $/ = undef;
	my ($data) = <FILE>;
	close(FILE);
	my ($digest) = md5_hex($data);
	$set->{hash} = $digest;

	my ( $rc, my $ok ) = $db->createCollection( "calculations/" . $digest );
	opendir( DIR, "$dir" ) || die("Cannot open directory $!\n");
	my (@filelist) = readdir(DIR);
	closedir DIR;
	for my $eachFile (@filelist) {
		if ( $eachFile =~ m/.*.xml/ ) {

			print $eachFile. "\n";
			if ( ( -s $dir . $eachFile ) < 5000000 ) {
				open( FILE, $dir . $eachFile );
				local $/ = undef;
				my ($xmldata) = <FILE>;
				close(FILE);
				print "uploading "
				  . "calculations/"
				  . $digest . "/"
				  . $eachFile . "\n";

					my ( $rc, my $ok ) = $db->uploadDocument(
						"calculations/" . $digest . "/" . $eachFile,
						$xmldata, replace => "<true>" );
						 $rc==0  or die "ERROR:cannot upload': $ok ($rc)\n";
			}
			else {
				print $dir
				  . $eachFile
				  . " Is too big please upload manually if you insist\n";
				my ($filesize) = -s $dir . $eachFile;
				print "Filesize: " . $filesize . "\n";
			}
		}
	}
}

my $data = $xms->XMLout( $doc, rootname => "experiment" );
my ($digest) = md5_hex($data);
my ( $rc, my $ok ) = $db->uploadDocument( "experiments/" . $digest . ".xml",
	$data, replace => "<true>" );
$rc == 0 or die "ERROR:cannot upload': $ok ($rc)\n";

