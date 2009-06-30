#!perl

@fortran= split(/\n/,`xsltproc ../../xml/schematofortran.xsl $ARGV[0] `);

foreach $line(@fortran){
#	$line=~s/tns\://g;
	if($line =~m/\(\/(.*?)\/\)/){
		$string=$1;
		$string=~s/\s+/ /g;
		
	@array=split(/\s/,$string);
	$string=join(",",@array);
	$line =~s/\(\/(.*?)\/\)/\(\/$string\/\)/;
	}
	print $line,"\n";
}