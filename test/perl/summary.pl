use lib "./perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use List::Util qw[min max];
use Data::Dumper;

my @allreports=collectreports();
make_summary (@allreports);



sub collectreports{
my $dir ,$file;
@reports=();
    opendir TESTROOT , "./";
   my @thefiles= readdir(TESTROOT);
    closedir(TESTROOT);
    
    foreach $dir (@thefiles)
    {
	if($dir=~m/test.*/)
	{
	    
	    opendir TESTDIR ,$dir;
	   my @thefiles2= readdir(TESTDIR);
	    closedir(TESTDIR);
	    foreach $file (@thefiles2)
	    {
		if($file=~m/.*\.xml/)
		{
		    push (@reports, readreport($dir . "/" . $file));
		}
	    }
	}
    }
    return @reports;
}



sub make_summary{
%merged=();
    $reports=@_[0];
    foreach $report (@reports)
    {
		$tests=$report->{test};
		while ( ($k,$v) = each(%$tests) ) {
	    	$merged{$k} = $v;
		}
	
    }
   	print Dumper %merged;
   	$xml = new XML::Simple(NoAttr=>1, RootName=>'report');
    $data{test}=\%merged;
open ALL,"all.xml";
    print ALL $xml->XMLout(\%data)
    close ALL;

}


sub readreport($file)
{
    $xml = new XML::Simple(NoAttr=>1, RootName=>'report');
    
# read XML file
#print @_[0];
    $data = $xml->XMLin(@_[0]);
    
    return $data;
}



