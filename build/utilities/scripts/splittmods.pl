while(<>){
if(m/^\s*module\s+(.+)/){
open FILE,">$1.F90";
print FILE  "#include \"maxdefinitions.inc\"\n";
}
print FILE;
if(m/^\s*end\s+module/){
close FILE,
}
}
