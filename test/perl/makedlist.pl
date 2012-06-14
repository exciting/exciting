use lib "./perl/";
use lib "./perl/lib";
use makeDlist;
use XML::Writer;
use IO::File;
use List::Util qw[min max];
use Data::Dumper;

Summary::makeDlist();
#$error=Summary::make_summary (@allreports);
#Summary::write_history("");
#print $error,"\n";
#exit ($error);