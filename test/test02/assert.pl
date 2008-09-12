
use XML::Simple
my @report;
#do things to assert test


push (@report,{
	"name"->"",
	"status"=>"",
	"description"=>""
	});	

#repeat
my $xsimple = XML::Simple->new();

print $xsimple->XMLout(\%report,
                       noattr => 1,
                       xmldecl => '<?xml version="1.0">');
