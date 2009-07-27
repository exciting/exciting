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

#print @fileslist;
$command = "../../utilities/scripts/protex -s " . join( " ", @fileslist );
open TEX, ">", "doc.tex";
print TEX `$command`;
close TEX;
system "pdflatex", "doc.tex";
system "pdflatex", "doc.tex";
