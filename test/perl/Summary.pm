package Summary;
use Data::Dumper;
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
    $reports=@_[0];
    foreach $report (@reports)
    {
	$tests=$report->{test};
	while ( ($k,$v) = each(%$tests) ) 
		{
		    $merged{$k} = $v;
		    if($v->{status} eq "failed")
		    {
			$failedtests{$k}=$v;
		    } 
		    elsif($v->{status} eq "passed")
		    {
			$passedtests{$k}=$v;
	    		}
		    else
		    {
	    			$unspecifiedtests{$k}=$v;
	    			print "$k\n";
			    }
		    
		}
	 return keys %$failedtests
    }

    $xml = new XML::Simple(NoAttr=>1, RootName=>'report',XMLDecl=>"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<?xml-stylesheet href=\"./report.xsl\" type=\"text/xsl\"?>");
    $data{test}=\%merged;
     age_files("all.xml");
    unless (open (ALL,">./report/all.xml")){
	die "Sorry, I couldn't create all.xml: $!";
    };
    print ALL $xml->XMLout(\%data);
    close ALL;
    
    age_files("failed.xml");
    unless (open (FAILED ,">./report/failed.xml")){
	die "Sorry, I couldn't create failed.xml: $!";
    };
    $data{test}=\%failedtests;
    print FAILED $xml->XMLout(\%data);
    close FAILED;
      
      age_files("passed.xml");
    unless (open (PASSED ,">./report/passed.xml")){
	die "Sorry, I couldn't create passed.xml: $!";
    };
    
    $data{test}=\%passedtests;
    print PASSED $xml->XMLout(\%data);
    close PASSED;
     
      age_files("unspecified.xml");
    unless (open (UNSPEC ,">./report/unspecified.xml")){
	die "Sorry, I couldn't create unsecified.xml: $!";
    };
    $data{test}=\%unspecifiedtests;
    print UNSPEC $xml->XMLout(\%data);
    close UNSPEC;
    
    
    
    $npassed=keys %passedtests;
    $nfailed=keys %failedtests;
    $nall=keys %merged;
    $nunspec=keys %unspecifiedtests;
    $ppassed=$npassed/$nall*100;
    $punspec=$nunsped/$nall*100;
    $pfailed=$nfailed/$nall*100;
    age_files("stats.xml");
    {
	my $output = new IO::File(">report/stats.xml");
  	my $writer = new XML::Writer(OUTPUT => $output,DATA_MODE => 'true', DATA_INDENT => 2);
  	$writer->xmlDecl( 'UTF-8' );
  	$writer->pi('xml-stylesheet', 'href="./stats.xsl" type="text/xsl"');
  	$writer->startTag("statistics");
	$writer->startTag("run");
	$writer->emptyTag("githash",hash=>get_git_hash("../src/version.inc"));
 	$writer->emptyTag("passed",count=>$npassed,percent=>$npassed/$nall*100.0);
  	$writer->emptyTag("failed",count=>$nfailed,percent=>$nfailed/$nall*100.0);
   	$writer->emptyTag("all",count=>$nall);
	$writer->emptyTag("unspecified",count=>$nunspec,percent=>$nunspec/$nall*100.0);
	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  	$mon=$mon+1;
	$writer->emptyTag("timestamp",time=>time(),timestring=>"$mday.$mon $hour:$min");
	$writer->endTag("run");
	$writer->endTag("statistics");
	
	$output=$writer-> getOutput();
  	$output->close();
    }
}


sub readreport
{
    $xml = new XML::Simple( RootName=>'report');
    
# read XML file
#print @_[0];
    $data = $xml->XMLin(@_[0]);
    
    return $data;
}


sub get_git_hash{
    
    open VERS ,@_[0];
    while(<VERS>){
	if (m/GITHASH\s+\"([\w|\d]+)\"/){
	    return $1;
	}
	else
	{return "no git"
	 }
	
    }
}

sub age_files{
    my $filetoage=@_[0];
    my $dir="./report";
    if (!opendir REPORTROOT , $dir)
    {die "error open $dir:$!";}
    my @thefiles= readdir(REPORTROOT);
    closedir (REPORTROOT);
    @thefiles=reverse(sort(@thefiles));
    foreach $file (@thefiles)
	
    {
	if ($file=~m/^((\d+)\.$filetoage)/)
	{
	    $n=$2+1;
	    if($n>=100){
		rename("$dir/$1","notmorethan100$filetoage")
		}else{
		    if($n<10){  rename("$dir/$1","$dir/0$n.$filetoage");}
		    else
		    {rename("$dir/$1","$dir/$n.$filetoage");}
	}
	} 	
    }
    foreach $file (@thefiles)
	
    {
	if($file=~m/^$filetoage/)
	{
	    rename("$dir/$file","$dir/01.$filetoage");
	}   
    }

}

sub write_history
{
    my $dir ="report";
    my $file;
    @stats=();
    if (!opendir REPORTROOT , $dir)
    {die "error open $dir:$! ";}
    my @thefiles= readdir(REPORTROOT);
    closedir (REPORTROOT);
    @thefiles=sort(@thefiles);
    foreach $file (@thefiles)
    {
		if($file=~m/((.*)stats.xml)/)
		{
		
	 	 $stat= readreport($dir . "/" . $file) ;
	 	 $stats{run}{$2}=$stat->{run};
		}
    }

$xml = new XML::Simple( RootName=>'statistics',XMLDecl=>"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<?xml-stylesheet  href=\"./stats.xsl\" type=\"text/xsl\"?>");
unless (open (ALL,">./report/index.xml")){
	die "Sorry, I couldn't create index.xml: $!";
    };
    print ALL $xml->XMLout(\%stats);
    close ALL;
return 1;
}

return 1;
