
    use lib "../test/perl/"; use lib "../test/perl/lib"; use lib "../utilities/lib"; use XML::Writer; use IO::File; use Switch; $npt="[-+]?[0-9]*\.?[0-9]*";
    local $/ = undef; open BLOCKINPUT, $ARGV[0] or die "Couldn't open file:$ARGV[0] $!"; binmode BLOCKINPUT; $blockinput = <BLOCKINPUT>."\n\n";
$blockinput=~s/[!:](.*)\n/\n/g;
$/="\n";
$enumhashfixspin{"0"}="none";
$enumhashfixspin{"1"}="total FSM";
$enumhashfixspin{"2"}="localmt FSM";
$enumhashfixspin{"3"}="both";
$enumhashstype{"0"}="Gaussian";
$enumhashstype{"1"}="Methfessel-Paxton 1";
$enumhashstype{"2"}="Methfessel-Paxton 2";
$enumhashstype{"3"}="Fermi Dirac";
$enumhashsolver{"1"}="Lapack";
$enumhashsolver{"2"}="Arpack";
$enumhashsolver{"3"}="DIIS";
$enumhashmixer{"1"}="lin";
$enumhashmixer{"2"}="msec";
$enumhashmixer{"3"}="pulay";
$enumhashxctype{"2"}="LDAPerdew-Zunger";
$enumhashxctype{"3"}="LSDAPerdew-Wang";
$enumhashxctype{"4"}="LDA-X-alpha";
$enumhashxctype{"5"}="LSDA-Barth-Hedin";
$enumhashxctype{"20"}="GGAPerdew-Burke-Ernzerhof";
$enumhashxctype{"21"}="GGArevPBE";
$enumhashxctype{"22"}="GGAPBEsol";
$enumhashxctype{"26"}="GGA-Wu-Cohen";
$enumhashxctype{"30"}="GGAArmiento-Mattsson";
$enumhashxctype{""}="EXX";
$enumhashxctype{"0"}="none";

sub getblock {
  $blockname=@_[0];
  if($blockinput=~m/\n\s*$blockname\s*\n\s*(.*?)\n\s*\n/s){
  return($1."\n");
  }
  else{
    return "";
  }
}
sub getvalue{
  $myblock=getblock(@_[0]);
  $myblock=~s/.true./true/i;
  $myblock=~s/.false./false/i;
   $myblock=~s/\n$//m;
    $myblock=~s/\'(.*)\'/$1/g;
  return($myblock);
}

my %atthashorigin=();
   
if(getvalue("coord"))
  {
   $atthashorigin{"coord"}=getvalue("coord");
 
    }
   
my %atthashpoint=();
   
if(getvalue("coord"))
  {
   $atthashpoint{"coord"}=getvalue("coord");
 
    }
   
if(getvalue("label"))
  {
   $atthashpoint{"label"}=getvalue("label");
 
    }
   
my %atthashpath=();
   
if(getvalue("steps"))
  {
   $atthashpath{"steps"}=getvalue("steps");
 
    }
   
if(getvalue("outfileprefix"))
  {
   $atthashpath{"outfileprefix"}=getvalue("outfileprefix");
 
    }
   
my %atthashparallelogram=();
   
if(getvalue("grid"))
  {
   $atthashparallelogram{"grid"}=getvalue("grid");
 
    }
   
if(getvalue("outfileprefix"))
  {
   $atthashparallelogram{"outfileprefix"}=getvalue("outfileprefix");
 
    }
   
my %atthashbox=();
   
if(getvalue("grid"))
  {
   $atthashbox{"grid"}=getvalue("grid");
 
    }
   
if(getvalue("outfileprefix"))
  {
   $atthashbox{"outfileprefix"}=getvalue("outfileprefix");
 
    }
   
my %atthashinput=();
   
if(getvalue("xsltpath"))
  {
   $atthashinput{"xsltpath"}=getvalue("xsltpath");
 
    }
   
if(getvalue("scratchpath"))
  {
   $atthashinput{"scratchpath"}=getvalue("scratchpath");
 
    }
   
if(getvalue("id"))
  {
   $atthashinput{"id"}=getvalue("id");
 
    }
   
if(getvalue("depends"))
  {
   $atthashinput{"depends"}=getvalue("depends");
 
    }
   
my %atthashstructure=();
   
if(getvalue("sppath"))
  {
   $atthashstructure{"speciespath"}=getvalue("sppath");
 
    }
   
if(getvalue("molecule"))
  {
   $atthashstructure{"molecule"}=getvalue("molecule");
 
    }
   
if(getvalue("epslat"))
  {
   $atthashstructure{"epslat"}=getvalue("epslat");
 
    }
   
if(getvalue("autormt"))
  {
   $atthashstructure{"autormt"}=getvalue("autormt");
 
    }
   
if(getvalue("vacuum"))
  {
   $atthashstructure{"vacuum"}=getvalue("vacuum");
 
    }
   
if(getvalue("primcell"))
  {
   $atthashstructure{"primcell"}=getvalue("primcell");
 
    }
   
my %atthashsymmetries=();
   
if(getvalue("hrmg"))
  {
   $atthashsymmetries{"HermannMauguinSymbol"}=getvalue("hrmg");
 
    }
   
if(getvalue("HallSymbol"))
  {
   $atthashsymmetries{"HallSymbol"}=getvalue("HallSymbol");
 
    }
   
if(getvalue("SchoenfliesSymbol"))
  {
   $atthashsymmetries{"SchoenfliesSymbol"}=getvalue("SchoenfliesSymbol");
 
    }
   
if(getvalue("spaceGroupNumber"))
  {
   $atthashsymmetries{"spaceGroupNumber"}=getvalue("spaceGroupNumber");
 
    }
   
my %atthashlattice=();
   
if(getvalue("a"))
  {
   $atthashlattice{"a"}=getvalue("a");
 
    }
   
if(getvalue("b"))
  {
   $atthashlattice{"b"}=getvalue("b");
 
    }
   
if(getvalue("c"))
  {
   $atthashlattice{"c"}=getvalue("c");
 
    }
   
if(getvalue("ab"))
  {
   $atthashlattice{"ab"}=getvalue("ab");
 
    }
   
if(getvalue("ac"))
  {
   $atthashlattice{"ac"}=getvalue("ac");
 
    }
   
if(getvalue("bc"))
  {
   $atthashlattice{"bc"}=getvalue("bc");
 
    }
   
if(getvalue("ncell"))
  {
   $atthashlattice{"ncell"}=getvalue("ncell");
 
    }
   
my %atthashwspecies=();
   
if(getvalue("speciesfile"))
  {
   $atthashwspecies{"speciesfile"}=getvalue("speciesfile");
 
    }
   
my %atthashwpos=();
   
if(getvalue("coord"))
  {
   $atthashwpos{"coord"}=getvalue("coord");
 
    }
   
my %atthashcrystal=();
   
if(getvalue("scale"))
  {
   $atthashcrystal{"scale"}=getvalue("scale");
 
    }
   
if(getvalue("(sc1|sc2|sc3)"))
  {
   $atthashcrystal{"stretch"}=getvalue("(sc1|sc2|sc3)");
 
    }
   
my %atthashspecies=();
   
if(getvalue("spfname"))
  {
   $atthashspecies{"speciesfile"}=getvalue("spfname");
 
    }
   
if(getvalue("spsymb"))
  {
   $atthashspecies{"chemicalSymbol"}=getvalue("spsymb");
 
    }
   
if(getvalue("atomicNumber"))
  {
   $atthashspecies{"atomicNumber"}=getvalue("atomicNumber");
 
    }
   
my %atthashatom=();
   
if(getvalue("atposl"))
  {
   $atthashatom{"coord"}=getvalue("atposl");
 
    }
   
if(getvalue("bfcmt"))
  {
   $atthashatom{"bfcmt"}=getvalue("bfcmt");
 
    }
   
my %atthashLDAplusu=();
   
if(getvalue("notaname"))
  {
   $atthashLDAplusu{"L"}=getvalue("notaname");
 
    }
   
if(getvalue("notaname"))
  {
   $atthashLDAplusu{"U"}=getvalue("notaname");
 
    }
   
if(getvalue("notaname"))
  {
   $atthashLDAplusu{"J"}=getvalue("notaname");
 
    }
   
my %atthashgroundstate=();
   
if(getvalue("ngridk"))
  {
   $atthashgroundstate{"ngkgrid"}=getvalue("ngridk");
 
    }
   
if(getvalue("rgkmax"))
  {
   $atthashgroundstate{"rgkmax"}=getvalue("rgkmax");
 
    }
   
if(getvalue("epspot"))
  {
   $atthashgroundstate{"epspot"}=getvalue("epspot");
 
    }
   
if(getvalue("rmtapm"))
  {
   $atthashgroundstate{"rmtapm"}=getvalue("rmtapm");
 
    }
   
if(getvalue("swidth"))
  {
   $atthashgroundstate{"swidth"}=getvalue("swidth");
 
    }
   
if(getvalue("stype"))
  {
  $atthashgroundstate{"stype"}=$enumhashstype{getvalue("stype")};

    }
   
if(getvalue("isgkmax"))
  {
   $atthashgroundstate{"isgkmax"}=getvalue("isgkmax");
 
    }
   
if(getvalue("gmaxvr"))
  {
   $atthashgroundstate{"gmaxvr"}=getvalue("gmaxvr");
 
    }
   
if(getvalue("nempty"))
  {
   $atthashgroundstate{"nempty"}=getvalue("nempty");
 
    }
   
if(getvalue("nosym"))
  {
   $atthashgroundstate{"nosym"}=getvalue("nosym");
 
    }
   
if(getvalue("autokpt"))
  {
   $atthashgroundstate{"autokpt"}=getvalue("autokpt");
 
    }
   
if(getvalue("radkpt"))
  {
   $atthashgroundstate{"radkpt"}=getvalue("radkpt");
 
    }
   
if(getvalue("reducek"))
  {
   $atthashgroundstate{"reducek"}=getvalue("reducek");
 
    }
   
if(getvalue("tfibs"))
  {
   $atthashgroundstate{"tfibs"}=getvalue("tfibs");
 
    }
   
if(getvalue("tforce"))
  {
   $atthashgroundstate{"tforce"}=getvalue("tforce");
 
    }
   
if(getvalue("lmaxapw"))
  {
   $atthashgroundstate{"lmaxapw"}=getvalue("lmaxapw");
 
    }
   
if(getvalue("maxscl"))
  {
   $atthashgroundstate{"maxscl"}=getvalue("maxscl");
 
    }
   
if(getvalue("chgexs"))
  {
   $atthashgroundstate{"chgexs"}=getvalue("chgexs");
 
    }
   
if(getvalue("deband"))
  {
   $atthashgroundstate{"deband"}=getvalue("deband");
 
    }
   
if(getvalue("epschg"))
  {
   $atthashgroundstate{"epschg"}=getvalue("epschg");
 
    }
   
if(getvalue("epsocc"))
  {
   $atthashgroundstate{"epsocc"}=getvalue("epsocc");
 
    }
   
if(getvalue("solver"))
  {
  $atthashgroundstate{"solver"}=$enumhashsolver{getvalue("solver")};

    }
   
if(getvalue("mixtype"))
  {
  $atthashgroundstate{"mixer"}=$enumhashmixer{getvalue("mixtype")};

    }
   
if(getvalue("fromscratch"))
  {
   $atthashgroundstate{"fromscratch"}=getvalue("fromscratch");
 
    }
   
if(getvalue("lradstp"))
  {
   $atthashgroundstate{"lradstep"}=getvalue("lradstp");
 
    }
   
if(getvalue("nprad"))
  {
   $atthashgroundstate{"nprad"}=getvalue("nprad");
 
    }
   
if(getvalue("xctype"))
  {
  $atthashgroundstate{"xctype"}=$enumhashxctype{getvalue("xctype")};

    }
   
if(getvalue("evalmin"))
  {
   $atthashgroundstate{"evalmin"}=getvalue("evalmin");
 
    }
   
if(getvalue("lmaxvr"))
  {
   $atthashgroundstate{"lmaxvr"}=getvalue("lmaxvr");
 
    }
   
if(getvalue("fracinr"))
  {
   $atthashgroundstate{"fracinr"}=getvalue("fracinr");
 
    }
   
if(getvalue("lmaxinr"))
  {
   $atthashgroundstate{"lmaxinr"}=getvalue("lmaxinr");
 
    }
   
if(getvalue("lmaxmat"))
  {
   $atthashgroundstate{"lmaxmat"}=getvalue("lmaxmat");
 
    }
   
if(getvalue("kdotpgrid"))
  {
   $atthashgroundstate{"kdotpgrid"}=getvalue("kdotpgrid");
 
    }
   
if(getvalue("vkloff"))
  {
   $atthashgroundstate{"vkloff"}=getvalue("vkloff");
 
    }
   
if(getvalue("npsden"))
  {
   $atthashgroundstate{"npsden"}=getvalue("npsden");
 
    }
   
if(getvalue("packedmatrixstorage"))
  {
   $atthashgroundstate{"packedmatrixstorage"}=getvalue("packedmatrixstorage");
 
    }
   
my %atthashspin=();
   
if(getvalue("bfieldc"))
  {
   $atthashspin{"bfieldc"}=getvalue("bfieldc");
 
    }
   
if(getvalue("momfix"))
  {
   $atthashspin{"momfix"}=getvalue("momfix");
 
    }
   
if(getvalue("spinorb"))
  {
   $atthashspin{"spinorb"}=getvalue("spinorb");
 
    }
   
if(getvalue("spinsprl"))
  {
   $atthashspin{"spinsprl"}=getvalue("spinsprl");
 
    }
   
if(getvalue("vqlss"))
  {
   $atthashspin{"vqlss"}=getvalue("vqlss");
 
    }
   
if(getvalue("taufsm"))
  {
   $atthashspin{"taufsm"}=getvalue("taufsm");
 
    }
   
if(getvalue("reducebf"))
  {
   $atthashspin{"reducebf"}=getvalue("reducebf");
 
    }
   
if(getvalue("fixspin"))
  {
  $atthashspin{"fixspin"}=$enumhashfixspin{getvalue("fixspin")};

    }
   
my %atthashstructureoptimization=();
   
if(getvalue("epsforce"))
  {
   $atthashstructureoptimization{"epsforce"}=getvalue("epsforce");
 
    }
   
if(getvalue("tau0atm"))
  {
   $atthashstructureoptimization{"tau0atm"}=getvalue("tau0atm");
 
    }
   
if(getvalue("resume"))
  {
   $atthashstructureoptimization{"resume"}=getvalue("resume");
 
    }
   
my %atthashbandstructure=();
   
if(getvalue("scissor"))
  {
   $atthashbandstructure{"scissor"}=getvalue("scissor");
 
    }
   
if(getvalue("character"))
  {
   $atthashbandstructure{"character"}=getvalue("character");
 
    }
   
my %atthashdos=();
   
if(getvalue("lmirep"))
  {
   $atthashdos{"lmirep"}=getvalue("lmirep");
 
    }
   
if(getvalue("nwdos"))
  {
   $atthashdos{"nwdos"}=getvalue("nwdos");
 
    }
   
if(getvalue("ngrdos"))
  {
   $atthashdos{"ngrdos"}=getvalue("ngrdos");
 
    }
   
if(getvalue("scissor"))
  {
   $atthashdos{"scissor"}=getvalue("scissor");
 
    }
   
if(getvalue("nsmdos"))
  {
   $atthashdos{"nsmdos"}=getvalue("nsmdos");
 
    }
   
if(getvalue("wintdos"))
  {
   $atthashdos{"wintdos"}=getvalue("wintdos");
 
    }
   
my %atthashmasstensor=();
   
if(getvalue("deltaem"))
  {
   $atthashmasstensor{"deltaem"}=getvalue("deltaem");
 
    }
   
if(getvalue("ndspem"))
  {
   $atthashmasstensor{"ndspem"}=getvalue("ndspem");
 
    }
   
if(getvalue("vklem"))
  {
   $atthashmasstensor{"vklem"}=getvalue("vklem");
 
    }
   
my %atthashfermisurfaceplot=();
   
if(getvalue("nstfsp"))
  {
   $atthashfermisurfaceplot{"nstfsp"}=getvalue("nstfsp");
 
    }
   
if(getvalue("separate"))
  {
   $atthashfermisurfaceplot{"separate"}=getvalue("separate");
 
    }
   
my %atthashlinresponsetensor=();
   
if(getvalue("scissor"))
  {
   $atthashlinresponsetensor{"scissor"}=getvalue("scissor");
 
    }
   
my %atthashexcitedstates=();
   
if(getvalue("dummya"))
  {
   $atthashexcitedstates{"dummya"}=getvalue("dummya");
 
    }
   
my $output = new IO::File(">input.xml");
my $writer = new XML::Writer(OUTPUT => $output,DATA_MODE => 1);
$writer->xmlDecl("UTF-8");
$writer->pi('xml-stylesheet', 'href="inputtohtml.xsl" type="text/xsl"');

$writer->startTag("input", 
    "xsi:noNamespaceSchemaLocation"=>"excitinginput.xsd",
	"xmlns:xsi"=>"http://www.w3.org/2001/XMLSchema-instance",
	"xsltpath"=>"../../../xml/","scratchpath"=>"/tmp/chm/1"
	);
$writer->startTag("title");
$writer->endTag("title");
$writer->startTag("structure",%atthashstructure);
$writer->startTag("crystal",%atthashcrystal);
$string=getblock("avec")."\n";
while($string=~s/(.+)\n//){
  $writer->startTag("basevect");
  $writer->characters($1);
  $writer->endTag("basevect");
  }
  $writer->endTag("crystal");
  $string=getblock("atoms")."\n"; 
  $string=~s/(.*)\n//;
  $nspecies=$1;
  for(my $is=1;$is<=$nspecies;$is++)
    {
    $string=~s/(.*)\n//;
    $spfname=$1;
    $spfname=~s/['\s]//g;
    $spfname=~s/\.in/\.xml/g;
    $writer->startTag("species","speciesfile"=>$spfname);
    $string=~s/(.*)\n//;
    $natoms=$1;
    for(my $ia=1;$ia<=$natoms;$ia++)
      {
      $string=~s/\s*($npt\s+$npt\s+$npt)\s+//e;
      $coord=$1;
      $string=~s/($npt\s+$npt\s+$npt)\s*\n//e;
      $bfcmt=$1;
      $writer->startTag("atom","coord"=>$coord,"bfcmt"=>$bfcmt);
      $writer->endTag("atom");
    }
    $writer->endTag("species");
  }
$writer->endTag("structure");
$string=getblock("tasks")."\n";

while($string=~s/\s*(\d+)\s*\n//){
switch ($1) {
	case [0,1,2,3]{
	  if($1==0||$1==2)
	    {
	    $atthashgroundstate{"fromscratch"}="true";
	    }
	  if($1==1||$1==3)
	    {
	    $atthashgroundstate{"fromscratch"}="false";
	    }
      if(getvalue("tarpack")==".true."){
       $atthashgroundstate{"solver"}="Arpack";
      }

  
	  $writer->startTag("groundstate",%atthashgroundstate);
	  if(getvalue("spinpol")||getvalue("spinsprl"))
	    {
	    $writer->startTag("spin",%atthashspin);
	    $writer->endTag("spin");
	    }
	  $writer->endTag("groundstate");
	  if($1==2||$1==3)
	    {
	    $writer->startTag("structureoptimization",%atthashstructureoptimization);
	    $writer->endTag("structureoptimization");
	    }
	  }
    case [200, 201]
    {
     $writer->startTag("phonons",%atthashphonons);
     $writer->endTag("phonons");
    }
	}
}
   
$string=getblock("tasks")."\n";
$writer->startTag("properties");
while($string=~s/\s*(\d+)\s*\n//)
   {
   switch ($1) 
      {
        case[10]{
        $writer->startTag("dos",%atthashdos);
        $writer->endTag("dos");
      } 
      case[100]
      {
        $writer->startTag("fermisurfaceplot",%atthashfermisurfaceplot);
        $writer->endTag("fermisurfaceplot");
      } 
      case[20,21]
      {
       if($1==21){$atthashbandstructure{"character"}="true";}
       $writer->startTag("bandstructure",%atthashbandstructure);
       addplot1delement();
       $writer->endTag("bandstructure");
      }
      case[25]
      {
       $writer->startTag("masstensor",%atthashmasstensor);
       $writer->endTag("masstensor");
      }
      case[31, 32, 33]
      {
       $writer->startTag("chargedesityplot",%atthashhargedesityplot);
        if($1==31){addplot1delement;}
        if($1==32){addplot2delement;}
        if($1==33){addplot3delement;}
       $writer->endTag("chargedesityplot");
      }
      case[41, 42, 43]
      {
       $writer->startTag("exccplot",%atthashexccplot);
        if($1==41){addplot1delement;}
        if($1==42){addplot2delement;}
        if($1==43){addplot3delement;}
       $writer->endTag("exccplot");
      }
      case[51, 52, 53]
      {
       $writer->startTag("elfplot",%atthashelfplot);
        if($1==51){addplot1delement;}
        if($1==52){addplot2delement;}
        if($1==53){addplot3delement;}
       $writer->endTag("elfplot");
      }
        case[61, 62, 63]
      {
       $writer->startTag("wfplot",%atthashwfplot);
        if($1==61){addplot1delement;}
        if($1==62){addplot2delement;}
        if($1==63){addplot3delement;}
       $writer->endTag("wfplot");
      }
       case[ 72, 73]
      {
       $writer->startTag("mvecfield",%atthashmvecfield);
        if($1==72){addplot2delement;}
        if($1==73){addplot3delement;}
       $writer->endTag("mvecfield");
      }
      case[ 82, 83]
      {
       $writer->startTag("xcmvecfield",%atthashxcmvecfield);
        if($1==82){addplot2delement;}
        if($1==83){addplot3delement;}
       $writer->endTag("xcmvecfield");
      }
     case[ 115]
     {
      $writer->startTag("EFG",%atthashEFG);
      $writer->endTag("EFG");
     } 
     
      case[100,101]
      {
      if($1==101){$atthashfermisurfaceplot{"separate"}="true";}
       $writer->startTag("fermisurfaceplot",%atthashfermisurfaceplot);
       $writer->endTag("fermisurfaceplot");
      } 
    }
  }
$writer->endTag("properties");
$writer->endTag("input");
$writer->end();
$output->close();
sub addplot1delement
  {
  $writer->startTag("plot1d");
  $string=getblock("plot1d");
  $string=~s/\s*(\d+)\s+(\d+)\s*\n//;
  $writer->startTag("path","steps"=>$2);
  while($string=~s/(.+?)\s*\n//)
    {
      $writer->startTag("point","coord"=>$1);
      $writer->endTag("point");
    }
    $writer->endTag("path"); 
    $writer->endTag("plot1d"); 
  }
sub addplot2delement
  {
  $writer->startTag("plot2d");
  $string=getblock("plot2d");

  $string=~s/\n\s*(\d+\s+\d+\s+\d+)[\n\s]*$/\n/m;
  $writer->startTag("parallelogram","grid"=>$1);
  $string=~s/(.+)\n//;
  $writer->startTag("origin","coord"=>$1);
  $writer->endTag("origin");
  while($string=~s/(.+?)\s*\n//)
    {
      $writer->startTag("point","coord"=>$1);
      $writer->endTag("point");
    }
  $writer->endTag("parallelogram");
  $writer->endTag("plot2d");
  }
sub addplot3delement
  {
  $writer->startTag("plot3d");
  $string=getblock("plot3d");
  $string=~s/\n\s*(\d+\s+\d+\s+\d+)[\n\s]*$/\n/m;
  $writer->startTag("box","grid"=>$1);
  $string=~s/(.+)\n//;
  $writer->startTag("origin","coord"=>$1);
  $writer->endTag("origin");
  while($string=~s/(.+?)\s*\n//)
    {
      $writer->startTag("point","coord"=>$1);
      $writer->endTag("point");
    }
  $writer->endTag("box");
  $writer->endTag("plot3d");
  }
