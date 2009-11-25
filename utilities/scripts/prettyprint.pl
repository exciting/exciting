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
	
	$command=  "f90ppr < $file > ./tmp.f90";
	print $command ,"\n";
	sleep(1);
	system  $command
}
