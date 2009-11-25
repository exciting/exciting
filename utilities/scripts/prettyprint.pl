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
	$command=  "f90ppr < $file > ./pptmp";
	print $command ,"\n";
	$return= system  $command;
	print "return value ",$return,"\n";
	if ($return==0){
		$diffcommand="diff $file  ./pptmp ";
		$diff=`$diffcommand`;
		print "diff is:", $diff,"\n";
		if(not ($diff eq "")){
		$mvcommand="mv ./pptmp $file";
		print $mvcommand ,"\n";
		system $mvcommand;
		}else{
		system "rm pptmp"	
		}
	}
}
