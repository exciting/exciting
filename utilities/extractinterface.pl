#!

while(<>){
if(m/^\s*subroutine(.*)/gim||m/^\s*function(.*)/gim) 
{
$interface.=$_;
$name=$1;
$name=~s/\(.*\)//;
$name=~s/\s*//g;
}
if(m/intent/i) 
{
$interface.=$_;
}

}
print "module ifc_$name\n";
print "begin interface\n";
print $interface;
print "end interface\n";
print "end module\n";