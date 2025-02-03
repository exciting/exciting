foreach $sourcedir (@ARGV) {
	#if ( $sourcedir != "" ) 
	{
		opendir( DIR, "$sourcedir" )
		  || die("Cannot open directory $sourcedir :$!\n");
		@filelist = readdir(DIR);
		closedir DIR;

		foreach $file (@filelist) {
			if ( $file =~ m/.+\.[fF]90/ ) {
				push( @fileslist, "$sourcedir$file" );
			}
		}
	}
}

#print @fileslist;
$command = "../../build/utilities/scripts/protex -s " . join( " ", @fileslist );
open TEX, ">", "doc.tex";
print TEX `$command`;
close TEX;
system "pdflatex", "doc.tex";
system "pdflatex", "doc.tex";
