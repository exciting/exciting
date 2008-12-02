#generate makefile

$extractinterfacecommand = "perl ../../utilities/extractinterface.pl";
$interfacesdir="../../interfaces";
open IFCMAKEFILE, ">","../ifcmakefile";
print IFCMAKEFILE "all:interfaces\n\n";
print "@ARGV\n";
foreach $sourcedir (@ARGV) {
	opendir( DIR, "$sourcedir" )|| die("Cannot open directory $!\n");
	@filelist = readdir(DIR);
	closedir DIR;

	foreach $file (@filelist) {
		if ( $file =~ m/.+\.[fF]90/ ) {
			push( @interfacelist, "$interfacesdir/ifc_$file " );
			print IFCMAKEFILE
"$interfacesdir/ifc_$file:$sourcedir$file\n\t$extractinterfacecommand $sourcedir$file >$interfacesdir/ifc_$file\n";

		}
	}
}

print IFCMAKEFILE "\ninterfaces:@interfacelist\n\n";
close IFCMAKEFILE;