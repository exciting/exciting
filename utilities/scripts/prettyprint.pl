foreach $sourcedir (@ARGV) {
	opendir( DIR, "$sourcedir" ) || die("Cannot open directory $!\n");
	@filelist = readdir(DIR);
	closedir DIR;

	foreach $file (@filelist) {
		if ( $file =~ m/.+\.[fF]90/ ) {
			push( @fileslist, "$sourcedir$file" );
		}
	}
}

print @fileslist;

foreach $file (@fileslist) {
	$return=-1;
	$command=  "f90ppr < $file > ./tmp.f90";
	print $command ,"\n";
	$return= system  $command;
	print "return value ",$return,"\n";
	if ($return==0){
		$mvcommand="mv ./tmp.f90 $file";
		print $mvcommand ,"\n";
		system $mvcommand;
		
	}
}
