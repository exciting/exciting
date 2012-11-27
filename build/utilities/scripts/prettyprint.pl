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
	$return  = -1;
	$command = "emacs -batch $file  -l $ENV{'PWD'}/../../build/utilities/emacs-format-file  -f emacs-format-function  ";
	print $command , "\n";
	$return = system $command;
	print "return value ", $return, "\n";
	 
	}
	

