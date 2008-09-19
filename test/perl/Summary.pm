package Summary;
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
	    		}
	    	
		}
	
    }

   	$xml = new XML::Simple(NoAttr=>1, RootName=>'report',XMLDecl=>"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<?xml-stylesheet href=\"./report.xsl\" type=\"text/xsl\"?>");
    $data{test}=\%merged;
    unless (open (ALL,">./report/all.xml")){
  		die "Sorry, I couldn't create all.xml: $!";
    };
    print ALL $xml->XMLout(\%data);
    close ALL;
    
    unless (open (FAILED ,">./report/failed.xml")){
  		die "Sorry, I couldn't create failed.xml: $!";
    };
    $data{test}=\%failedtests;
    print FAILED $xml->XMLout(\%data);
    close FAILED;
    
    unless (open (PASSED ,">./report/passed.xml")){
  		die "Sorry, I couldn't create passed.xml: $!";
    };
    $data{test}=\%passedtests;
    print PASSED $xml->XMLout(\%data);
    close PASSED;
    
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
    {
    my $output = new IO::File(">report/stats.xml");
  	my $writer = new XML::Writer(OUTPUT => $output,DATA_MODE => 'true', DATA_INDENT => 2);
  	$writer->xmlDecl( 'UTF-8' );
  	$writer->pi('xml-stylesheet', 'href="./stats.xsl" type="text/xsl"');
  	$writer->startTag("statistics");
  		$writer->startTag("run");
 	$writer->emptyTag("passed",count=>$npassed,percent=>$npassed/$nall*100.0);
  	$writer->emptyTag("failed",count=>$nfailed,percent=>$nfailed/$nall*100.0);
   	$writer->emptyTag("all",count=>$nall);
    $writer->emptyTag("unspecified",count=>$nunspec,percent=>$nunspec/$nall*100.0);
    $writer->emptyTag("githash",hash=>get_git_hash("../src/version.inc"));
    ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  	$mon=$mon+1;
    $writer->emptyTag("timestamp",time=>time(),timestring=>"$mday.$mon $hour:$min");
   $writer->endTag("run");
    $writer->endTag("statistics");
 
    $output=$writer-> getOutput();
  	$output->close();
  	}
  	{
  	my $output = new IO::File(">report/index.html");
  	my $writer = new XML::Writer(OUTPUT => $output,DATA_MODE => 'true', DATA_INDENT => 2);
  	$writer->xmlDecl( 'UTF-8' );
  	$writer->startTag(html);$writer->startTag(body);
  	$writer->startTag(h1);
  	$writer->characters("Test Result");
  	$writer->endTag(h1);
  	$writer->startTag(p);
  	($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
  	$mon=$mon+1;
  	$writer->characters("Test suite run from"." $mday.$mon $hour:$min :  ");
  	$writer->endTag(p);$writer->startTag(p);
  	$url= "http://chart.apis.google.com/chart?cht=bhs".
  		"&chd=t:$ppassed|$punspec|$pfailed&".
  		"chs=500x80&".
	  	"chdl=passed|unspecified|failed&chdlp=t&".
  	 	"chl=passed:$npassed|unspecified:$nunspec|failed:$nfailed&".
  	 	"chco=006600,f0f000,cc0033";
 	$writer->emptyTag("img",src=>"$url");
    $writer->endTag(p);$writer->startTag(p);
  	$writer->startTag(a, href=>"passed.xml");
  	$writer->characters("passed");
  	$writer->endTag(a);
  	$writer->startTag(a, href=>"unspecified.xml");
  	$writer->characters("unspecified");
  	$writer->endTag(a);
  	$writer->startTag(a, href=>"failed.xml");
  	$writer->characters("failed");
  	$writer->endTag(a);
    $writer->endTag(p);
  	$writer->endTag(body);
  	$writer->endTag(html);
	$output=$writer-> getOutput();
  	$output->close();
}
}


sub readreport($file)
{
    $xml = new XML::Simple(NoAttr=>1, RootName=>'report');
    
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

return 1;