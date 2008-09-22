use lib "./perl/";
use lib "./perl/lib";
use Summary;
use XML::Simple;
use XML::Writer;
use IO::File;
use List::Util qw[min max];
use Data::Dumper;

my @allreports=Summary::collectreports();
Summary::make_summary (@allreports);
Summary::write_history("");