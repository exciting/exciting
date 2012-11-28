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
	$useinput   = 0;
	$moddefined = 0;
	$inmodue    = 0;
	local $/ = undef;
	print "open file ", $file, "\n";
	open FILE, $file or die "Couldn't open file: $!";
	binmode FILE;
	$string     = <FILE>;
	$string=~s/\n\s+(subroutine|function)/\n$1/gi;
	@routines   = split( /\n((?=(subroutine))|(?=(function)))/i, $string );
	$partnumber = 0;
part:
	foreach $string (@routines) {

	$inmodule    = 0;
		if($string eq "subroutine" || $string eq "function"||$string eq "\n") { 
			$string="";
			$partnumber++;
			next part;
		}
		print "NEWPART",$partnumber,"\n", $string;

    $useinput   = 0;
	$moddefined = 0;
		
		$string =~ s/\\\n/<linebreak>/g;
		$string =~ s/&\n/<linebreak>/g;
		@stringar = "";
		@stringar = split( /\n/, $string );
		close FILE;
		local $/ = "\n";
		$i = 0;

		for ( $i = 0 ; $i < scalar(@stringar) ; ++$i ) {
			if (   not( $stringar[$i] =~ m/^\s*!/ )
				&& not( $stringar[$i] =~ m/function/i ) )
			{
				$stringar[$i] =~ s/(^\s*logical)\s*(?!,)(?=\w)/$1::/gi;
				$stringar[$i] =~ s/(^\s*real(\(\d\))*)\s*(?!,)(?=\w)/$1::/gi;
				$stringar[$i] =~ s/(^\s*integer(\(\d\))*)\s*(?!,)(?=\w)/$1::/gi;
				$stringar[$i] =~
				  s/(character(\((\d+|\*)\))?)\s*(?!,)(?=\w)/$1::/gi;
			}
			if ( $stringar[$i] =~ m/use modinput/ ) {
				$moddefined = 1;
			}
		}
	  VAR:
		for my $att ( @{ $data->{attribute} } ) {
			if(not($att->{xpath}=~m/^input/))
			{ next VAR;}
		  LINE:
			for ( $i = 0 ; $i < scalar(@stringar) ; ++$i ) {
				if (m/^\s*module/)       { $inmodule = 1 }
				if (m/^\s*end\s+module/) { $inmodule = 0 }
				$nstrcon        = 0;
				@stringconstant = "";
					if ( $stringar[$i] =~ m/^\s*!/ ||$stringar[$i] =~ m/^\s*\#/ ) {
					next LINE;
				}
				while ( $stringar[$i] =~ s/(\".*?\")/<stringconstant>/ ) {
					$stringconstant[$nstrcon] = $1;
					$nstrcon++;
				}
			
				if ( $stringar[$i] =~ m/subroutine.*($att->{vname})(?!\w)/mi ) {
					print "alarm!", $file, " var ", $att->{vname}, " ",
					  $stringar[$i], "\n";
					$nstrcon = 0;
					while ( $stringar[$i] =~
						s/<stringconstant>/$stringconstant[$nstrcon]/ )
					{
						$nstrcon++;
					}
					next VAR;
				}
				if ( $stringar[$i] =~ m/::/ ) {
					if ( $stringar[$i] =~ m/(?!\)).*?(\s|,|:)$att->{vname}(\s|,|\(|$|=).*?(?!\))/ ) {
						$stringar[$i] =~ s/<stringconstant>/$stringconstant/;
						if ($inmodule) { $stringar[$i] =~ s/(.*)/!$1/; }
						next VAR;
					}

				}
				if (   $att->{vname} eq 'position'
					&& $stringar[$i] =~ m/^\s*open/ )
				{
					$nstrcon = 0;
					while ( $stringar[$i] =~
						s/<stringconstant>/$stringconstant[$nstrcon]/ )
					{
						$nstrcon++;
					}
					next LINE;
				}

				if ( $stringar[$i] =~ m/only:/ ) {
					$stringar[$i] =~ s/(?<![\w|%])$att->{vname}(?!\w)//i;

				}

				if ( $stringar[$i] =~
					s/(?<![\w|%])$att->{vname}(?!\w)/$att->{fortranname}/gi )
				{
					print "replace ", $att->{vname}, "by: ",
					  $att->{fortranname}, "\n";
					$useinput = 1;
				}
				if ( $stringar[$i] =~
m/input%structure%speciesarray\(speciesindex\)%species%atomarray\(atomindex\)%atom%(\w+?)\(/
				  )
				{
					$var = $1;
					$stringar[$i] =~
s/%$var\(\s*([\w\d:]+)\s*,\s*([\w\d:]+)\s*,\s*([\w\d:]+)\s*\)/%$var\($1\)/i;
					$sp = $3;
					$at = $2;
					$stringar[$i] =~ s/speciesindex/$sp/;
					$stringar[$i] =~ s/atomindex/$at/;

				}
				elsif ( $stringar[$i] =~
m/input%structure%speciesarray\(speciesindex\)%species%(\w+?)(?!%)\(/gi
				  )
				{
					$var = $1;
					if ( $stringar[$i] =~ s/%$var\(\s*([\w\d:]+)\s*\)/%$var/i )
					{
						$sp = $1;
						$stringar[$i] =~ s/speciesindex/$sp/;
					}
				}

				$stringar[$i] =~
s/(^\s*if(?!\w).*%fixspin(?!\w).*\)\s+then)/!$1 fixxxmlinput \n if(.false.) then/gi
				  || $stringar[$i] =~
s/(^\s*if(?!\w).*%fixspin(?!\w).*\).*(?!(then)))/!$1 fixxxmlinput \n  /gi;
$stringar[$i] =~s/input%xs%tddft%acont = input%xs%tddft%acont/acont = input%xs%tddft%acont/;
$stringar[$i] =~s/input%xs%tddft%fxctypenumber\s*=\s*input%xs%tddft%fxctypenumber/ fxctype = input%xs%tddft%fxctypenumber/;
$stringar[$i] =~s/input%properties%bandstructure%character\((\d+)\)/character($1)/;
 $stringar[$i] =~s/input%xs%BSE%bsetypenumber = input%xs%BSE%bsetypenumber/bsetype = input%xs%BSE%bsetypenumber/;


		#	$stringar[$i] =~
		#s/(^\s*if(?!\w).*%xctype.*\)(.*))/!$1 fixxxmlinput \n if(.false.)$2/gi;

	   #$stringar[$i] =~ s/abs\((.*?%xctype(?!\w))\)/$1/;
	   #if (   $stringar[$i] =~ s/(^\s*call\s+xcifc(?!\w).*)/!$1 !fixxxmlinput /
	   #		|| $stringar[$i] =~
	   #		s/(^\s*call\s+getxcdata.*)/!$1 !fixxxmlinput / )
	   #{
	   #		$stringar[$i] =~ s/<linebreak>/\n!/g;
	   #	}
				$nstrcon = 0;
				if ( length( $stringar[$i] ) > 100 ) {
					$stringar[$i] =~ s/(?<![\sd])([-\+=\*]+)(?!\s)/ $1 /g;
				}
				while ( $stringar[$i] =~
					s/<stringconstant>/$stringconstant[$nstrcon]/ )
				{
					$nstrcon++;
				}
			}

		}

		for ( $i = 0 ; $i < scalar(@stringar) ; ++$i ) {
			if ( $stringar[$i] =~ m/only:/ ) {
				$stringar[$i] =~ s/(,\s*)+/,/g;
				$stringar[$i] =~ s/,\s*$//g;
			}
			
		}
		$string = join( "\n", @stringar );
		$string =~ s/<linebreak>/&\n/g;

		@stringar = split( /\n/, $string );

		for ( $i = 0 ; $i < scalar(@stringar) ; ++$i ) {
			if ( not( $stringar[$i] =~ m/^\s*!/ ) ) {
				$Text::Wrap::columns   = 120;
				$Text::Wrap::separator = "&\n";
				$stringar[$i] =~ s/,(?!\s)/, /g;
				$stringar[$i] = wrap( "", "    &", $stringar[$i] );
			}
		
		if (    $useinput
				and $i==2
				and not $moddefined)
			{
				$stringar[$i] = "use modinput\n".$stringar[$i];
					$moddefined = 1;
				
				if (   $stringar[$i] =~ m/end\s+(subroutine|module|function)/
					&& $moddefined eq 1 )
				{
					$moddefined = 0;
				}
			}
		}
		$string = join( "\n", @stringar );
		$string =~ s/&&/&/g;
		$string =~ s/\*\s+\*/\*\*/g;
		while ( $string =~ s/&\n\s+&\n\s+/&\n   /g ) { }
		$string="\n".$string;
	
		$routines[$partnumber] = $string;
		
		$partnumber++;
	}
		
	$string = join( "", @routines );
	




	open FILE, ">", $file or die "Couldn't open file: $!";
	print FILE $string, "\n";
	close FILE;
}
