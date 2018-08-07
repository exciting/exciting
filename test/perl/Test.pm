package Test;
use lib "./lib/XML";
use XML::Simple;
use XML::Writer;
use IO::File;
use List::Util qw[min max];

sub assert_file_same_within {
	$file1=$_[0];
	$file2=$_[1];
	$tol=$_[2];
	@numbers1=[];
	@numbers2=[];
	$error=[];
	open FILE1,$file1 or die "cannot open :$file1 $!";
	open FILE2, $file2 or die "cannot open :$file2 $!";

	while(<FILE1>)
	{
		$linefile1=$_;
		$linefile2=<FILE2>;
			while ($linefile1=~s/([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)//){
			$n1=$1;
			push(@numbers1,$n1)	;
			$linefile2=~s/([-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)//;
			$n2=$1;
			push(@numbers2,$n2);
			push(@error,abs($n1-$n2));
			
 		}
		
	}
        $maxerror=max(@error);
	print "Maximal deviation = $maxerror\n";
	$status='failed';
        
        if (max(@error) < $tol){
            $status='passed';
 	}
 	close FILE1;
 	close FILE2;

%test=(status => $status,
	line => $linenr,
	column => $collumn, 
	maxerror=> $maxerror,
	averageerror=> $averageerr
	);
return %test;
}

sub initreport{ #call with filename
  my $output = new IO::File(">@_[0]");
  my $writer = new XML::Writer(OUTPUT => $output,DATA_MODE => 'true', DATA_INDENT => 2);
  $writer->xmlDecl( 'UTF-8' );
  $writer->pi('xml-stylesheet', 'href="./report.xsl" type="text/xsl"');
  $writer->startTag("report");
  return $writer
}
  
sub closereport{
  $writer=@_[0] ;
  $writer->endTag("report");
  $writer->end();
  $output=$writer-> getOutput();
  $output->close();
}
  
sub writetestreport(%$) { #hash with values and writer object
  $writer=@_[1];
  $elements=@_[0];
  $writer->startTag("test");
  #,"directory" => @_[0]->{"directory"}, 
  #                  "name"=>@_[0]->{"name"});
                   
  foreach $key (keys %$elements ) { # once for each key of @_

	#if($key ne "name" && $key ne "directory"){
		$writer->startTag($key);
		$writer->characters($elements->{$key});
		$writer->endTag($key);
	#	}
			
  } 
  $writer->endTag("test");  
  return 1;
}

return 1
