
    use lib "../test/perl/"; use lib "../test/perl/lib"; use lib "../build/utilities/lib"; use XML::Writer; use IO::File; use Switch; $npt="[-+]?[0-9]*\.?[0-9]*";
    local $/ = undef; open BLOCKINPUT, $ARGV[0] or die "Couldn't open file:$ARGV[0] $!"; binmode BLOCKINPUT; $blockinput = <BLOCKINPUT>."\n\n";
$blockinput=~s/[!:](.*)\n/\n/g;
$/="\n";
$enumhashfixspin{"0"}="none";
$enumhashfixspin{"1"}="total FSM";
$enumhashfixspin{"2"}="localmt FSM";
$enumhashfixspin{"3"}="both";
$enumhashtype{"1"}="Lapack";
$enumhashtype{"2"}="Arpack";
$enumhashtype{"3"}="DIIS";
$enumhashdo{""}="fromscratch";
$enumhashdo{""}="fromfile";
$enumhashdo{""}="skip";
$enumhashstype{"0"}="Gaussian";
$enumhashstype{"1"}="Methfessel-Paxton 1";
$enumhashstype{"2"}="Methfessel-Paxton 2";
$enumhashstype{"3"}="Fermi Dirac";
$enumhashstype{"4"}="Square-wave impulse";
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
$enumhashxctype{"-2"}="EXX";
$enumhashxctype{"0"}="none";
$enumhashfxctype{"0"}="RPA";
$enumhashfxctype{"1"}="LRCstatic_NLF";
$enumhashfxctype{"2"}="LRCstatic";
$enumhashfxctype{"3"}="LRCdyn_NLF";
$enumhashfxctype{"4"}="LRCdyn";
$enumhashfxctype{"5"}="ALDA";
$enumhashfxctype{"7"}="MB1_NLF";
$enumhashfxctype{"8"}="MB1";
$enumhashscreentype{""}="full";
$enumhashscreentype{""}="diag";
$enumhashscreentype{""}="noinvdiag";
$enumhashscreentype{""}="longrange";
$enumhashsciavtype{""}="spherical";
$enumhashsciavtype{""}="screendiag";
$enumhashsciavtype{""}="invscreendiag";
$enumhashbsetype{""}="ip";
$enumhashbsetype{""}="rpa";
$enumhashbsetype{""}="singlet";
$enumhashbsetype{""}="triplet";
$enumhashtask{"301"}="xsgeneigvec";
$enumhashtask{"310"}="tetcalccw";
$enumhashtask{"320"}="writepmatxs";
$enumhashtask{"330"}="writeemat";
$enumhashtask{"340"}="df";
$enumhashtask{"345"}="df2";
$enumhashtask{"350"}="idf";
$enumhashtask{"401"}="scrgeneigvec";
$enumhashtask{"410"}="scrtetcalccw";
$enumhashtask{"420"}="scrwritepmat";
$enumhashtask{"430"}="screen";
$enumhashtask{"440"}="scrcoulint";
$enumhashtask{"441"}="exccoulint";
$enumhashtask{"445"}="BSE";
$enumhashtask{"450"}="kernxc_bse";
$enumhashtask{"23"}="writebandgapgrid";
$enumhashtask{"120"}="writepmat";
$enumhashtask{"121"}="dielectric";
$enumhashtask{"321"}="writepmatasc";
$enumhashtask{"322"}="pmatxs2orig";
$enumhashtask{"331"}="writeematasc";
$enumhashtask{"335"}="writepwmat";
$enumhashtask{"339"}="emattest";
$enumhashtask{"341"}="x0toasc";
$enumhashtask{"342"}="x0tobin";
$enumhashtask{"396"}="epsconv";
$enumhashtask{"398"}="fxc_alda_check";
$enumhashtask{"451"}="kernxc_bse3";
$enumhashtask{"499"}="testxs";
$enumhashtask{"700"}="xsestimate";
$enumhashtask{"701"}="xstiming";
$enumhashtask{"999"}="testmain";
$enumhashtask{"900"}="portstate(1)";
$enumhashtask{"901"}="portstate(2)";
$enumhashtask{"910"}="portstate(-1)";
$enumhashtask{"911"}="portstate(-2)";
$enumhashxstype{""}="TDDFT";
$enumhashxstype{""}="BSE";

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
   
if(getvalue("tshift"))
  {
   $atthashstructure{"tshift"}=getvalue("tshift");
 
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
   
if(getvalue("noreplace"))
  {
   $atthashlattice{"a"}=getvalue("noreplace");
 
    }
   
if(getvalue("noreplace"))
  {
   $atthashlattice{"b"}=getvalue("noreplace");
 
    }
   
if(getvalue("noreplace"))
  {
   $atthashlattice{"c"}=getvalue("noreplace");
 
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
   
if(getvalue("noreplace>"))
  {
   $atthashspecies{"rmt"}=getvalue("noreplace>");
 
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
   
if(getvalue("mommtfix"))
  {
   $atthashatom{"mommtfix"}=getvalue("mommtfix");
 
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
   
if(getvalue("nonreplace"))
  {
  $atthashgroundstate{"do"}=$enumhashdo{getvalue("nonreplace")};

    }
   
if(getvalue("ngridk"))
  {
   $atthashgroundstate{"ngridk"}=getvalue("ngridk");
 
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
   
if(getvalue("frozencore"))
  {
   $atthashgroundstate{"frozencore"}=getvalue("frozencore");
 
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
   
if(getvalue("mixtype"))
  {
  $atthashgroundstate{"mixer"}=$enumhashmixer{getvalue("mixtype")};

    }
   
if(getvalue("beta0"))
  {
   $atthashgroundstate{"beta0"}=getvalue("beta0");
 
    }
   
if(getvalue("betainc"))
  {
   $atthashgroundstate{"betainc"}=getvalue("betainc");
 
    }
   
if(getvalue("betadec"))
  {
   $atthashgroundstate{"betadec"}=getvalue("betadec");
 
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
   
if(getvalue("cfdamp"))
  {
   $atthashgroundstate{"cfdamp"}=getvalue("cfdamp");
 
    }
   
if(getvalue("nosource"))
  {
   $atthashgroundstate{"nosource"}=getvalue("nosource");
 
    }
   
if(getvalue("tevecsv"))
  {
   $atthashgroundstate{"tevecsv"}=getvalue("tevecsv");
 
    }
   
if(getvalue("nwrite"))
  {
   $atthashgroundstate{"nwrite"}=getvalue("nwrite");
 
    }
   
if(getvalue("ptnucl"))
  {
   $atthashgroundstate{"ptnucl"}=getvalue("ptnucl");
 
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
   
my %atthashHartreeFock=();
   
if(getvalue("epsengy"))
  {
   $atthashHartreeFock{"epsengy"}=getvalue("epsengy");
 
    }
   
my %atthashsolver=();
   
if(getvalue("solvertype"))
  {
  $atthashsolver{"type"}=$enumhashtype{getvalue("solvertype")};

    }
   
if(getvalue("packedmatrixstorage"))
  {
   $atthashsolver{"packedmatrixstorage"}=getvalue("packedmatrixstorage");
 
    }
   
if(getvalue("epsarpack"))
  {
   $atthashsolver{"epsarpack"}=getvalue("epsarpack");
 
    }
   
if(getvalue("evaltol"))
  {
   $atthashsolver{"evaltol"}=getvalue("evaltol");
 
    }
   
my %atthashOEP=();
   
if(getvalue("maxitoep"))
  {
   $atthashOEP{"maxitoep"}=getvalue("maxitoep");
 
    }
   
if(getvalue("tauoep"))
  {
   $atthashOEP{"tauoep"}=getvalue("tauoep");
 
    }
   
my %atthashRDMFT=();
   
if(getvalue("rdmxctype"))
  {
   $atthashRDMFT{"rdmxctype"}=getvalue("rdmxctype");
 
    }
   
if(getvalue("rdmmaxscl"))
  {
   $atthashRDMFT{"rdmmaxscl"}=getvalue("rdmmaxscl");
 
    }
   
if(getvalue("maxitn"))
  {
   $atthashRDMFT{"maxitn"}=getvalue("maxitn");
 
    }
   
if(getvalue("maxitc"))
  {
   $atthashRDMFT{"maxitc"}=getvalue("maxitc");
 
    }
   
if(getvalue("taurdmn"))
  {
   $atthashRDMFT{"taurdmn"}=getvalue("taurdmn");
 
    }
   
if(getvalue("taurdmc"))
  {
   $atthashRDMFT{"taurdmc"}=getvalue("taurdmc");
 
    }
   
if(getvalue("rdmalpha"))
  {
   $atthashRDMFT{"rdmalpha"}=getvalue("rdmalpha");
 
    }
   
if(getvalue("rdmtemp"))
  {
   $atthashRDMFT{"rdmtemp"}=getvalue("rdmtemp");
 
    }
   
my %atthashrelax=();
   
if(getvalue("epsforce"))
  {
   $atthashrelax{"epsforce"}=getvalue("epsforce");
 
    }
   
if(getvalue("tau0atm"))
  {
   $atthashrelax{"tau0atm"}=getvalue("tau0atm");
 
    }
   
if(getvalue("resume"))
  {
   $atthashrelax{"resume"}=getvalue("resume");
 
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
   
if(getvalue("sqados"))
  {
   $atthashdos{"sqados"}=getvalue("sqados");
 
    }
   
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
   
if(getvalue("wdos"))
  {
   $atthashdos{"winddos"}=getvalue("wdos");
 
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
   
my %atthashelnes=();
   
if(getvalue("vecql"))
  {
   $atthashelnes{"vecql"}=getvalue("vecql");
 
    }
   
my %atthasheliashberg=();
   
if(getvalue("mustar"))
  {
   $atthasheliashberg{"mustar"}=getvalue("mustar");
 
    }
   
my %atthashphonons=();
   
if(getvalue("reduceq"))
  {
   $atthashphonons{"reduceq"}=getvalue("reduceq");
 
    }
   
if(getvalue("deltaph"))
  {
   $atthashphonons{"deltaph"}=getvalue("deltaph");
 
    }
   
my %atthashxs=();
   
if(getvalue("emattype"))
  {
   $atthashxs{"emattype"}=getvalue("emattype");
 
    }
   
if(getvalue("dfoffdiag"))
  {
   $atthashxs{"dfoffdiag"}=getvalue("dfoffdiag");
 
    }
   
if(getvalue("lmaxapwwf"))
  {
   $atthashxs{"lmaxapwwf"}=getvalue("lmaxapwwf");
 
    }
   
if(getvalue("lmaxemat"))
  {
   $atthashxs{"lmaxemat"}=getvalue("lmaxemat");
 
    }
   
if(getvalue("emaxdf"))
  {
   $atthashxs{"emaxdf"}=getvalue("emaxdf");
 
    }
   
if(getvalue("broad"))
  {
   $atthashxs{"broad"}=getvalue("broad");
 
    }
   
if(getvalue("tevout"))
  {
   $atthashxs{"tevout"}=getvalue("tevout");
 
    }
   
if(getvalue("xstype"))
  {
  $atthashxs{"xstype"}=$enumhashxstype{getvalue("xstype")};

    }
   
if(getvalue("symmorph"))
  {
   $atthashxs{"symmorph"}=getvalue("symmorph");
 
    }
   
if(getvalue("fastpmat"))
  {
   $atthashxs{"fastpmat"}=getvalue("fastpmat");
 
    }
   
if(getvalue("fastemat"))
  {
   $atthashxs{"fastemat"}=getvalue("fastemat");
 
    }
   
if(getvalue("gather"))
  {
   $atthashxs{"gather"}=getvalue("gather");
 
    }
   
if(getvalue("tappinfo"))
  {
   $atthashxs{"tappinfo"}=getvalue("tappinfo");
 
    }
   
if(getvalue("dbglev"))
  {
   $atthashxs{"dbglev"}=getvalue("dbglev");
 
    }
   
if(getvalue("usegdft"))
  {
   $atthashxs{"usegdft"}=getvalue("usegdft");
 
    }
   
if(getvalue("gqmax"))
  {
   $atthashxs{"gqmax"}=getvalue("gqmax");
 
    }
   
if(getvalue("nosymxs"))
  {
   $atthashxs{"nosym"}=getvalue("nosymxs");
 
    }
   
if(getvalue("ngridkxs"))
  {
   $atthashxs{"ngridk"}=getvalue("ngridkxs");
 
    }
   
if(getvalue("vkloffxs"))
  {
   $atthashxs{"vkloff"}=getvalue("vkloffxs");
 
    }
   
if(getvalue("reducekxs"))
  {
   $atthashxs{"reducek"}=getvalue("reducekxs");
 
    }
   
if(getvalue("ngridqxs"))
  {
   $atthashxs{"ngridq"}=getvalue("ngridqxs");
 
    }
   
if(getvalue("reduceqxs"))
  {
   $atthashxs{"reduceq"}=getvalue("reduceqxs");
 
    }
   
if(getvalue("rgkmaxxs"))
  {
   $atthashxs{"rgkmax"}=getvalue("rgkmaxxs");
 
    }
   
if(getvalue("swidthxs"))
  {
   $atthashxs{"swidth"}=getvalue("swidthxs");
 
    }
   
if(getvalue("lmaxapwxs"))
  {
   $atthashxs{"lmaxapw"}=getvalue("lmaxapwxs");
 
    }
   
if(getvalue("lmaxmatxs"))
  {
   $atthashxs{"lmaxmat"}=getvalue("lmaxmatxs");
 
    }
   
if(getvalue("nemptyxs"))
  {
   $atthashxs{"nempty"}=getvalue("nemptyxs");
 
    }
   
if(getvalue("scissor"))
  {
   $atthashxs{"scissor"}=getvalue("scissor");
 
    }
   
my %atthashtddft=();
   
if(getvalue("intraband"))
  {
   $atthashtddft{"intraband"}=getvalue("intraband");
 
    }
   
if(getvalue("torddf"))
  {
   $atthashtddft{"torddf"}=getvalue("torddf");
 
    }
   
if(getvalue("tordfxc"))
  {
   $atthashtddft{"tordfxc"}=getvalue("tordfxc");
 
    }
   
if(getvalue("aresdf"))
  {
   $atthashtddft{"aresdf"}=getvalue("aresdf");
 
    }
   
if(getvalue("aresfxc"))
  {
   $atthashtddft{"aresfxc"}=getvalue("aresfxc");
 
    }
   
if(getvalue("fxcbsesplit"))
  {
   $atthashtddft{"fxcbsesplit"}=getvalue("fxcbsesplit");
 
    }
   
if(getvalue("acont"))
  {
   $atthashtddft{"acont"}=getvalue("acont");
 
    }
   
if(getvalue("nwacont"))
  {
   $atthashtddft{"nwacont"}=getvalue("nwacont");
 
    }
   
if(getvalue("lindhard"))
  {
   $atthashtddft{"lindhard"}=getvalue("lindhard");
 
    }
   
if(getvalue("epsdfde"))
  {
   $atthashtddft{"epsdfde"}=getvalue("epsdfde");
 
    }
   
if(getvalue("kerndiag"))
  {
   $atthashtddft{"kerndiag"}=getvalue("kerndiag");
 
    }
   
if(getvalue("lmaxalda"))
  {
   $atthashtddft{"lmaxalda"}=getvalue("lmaxalda");
 
    }
   
if(getvalue("alphalrc"))
  {
   $atthashtddft{"alphalrc"}=getvalue("alphalrc");
 
    }
   
if(getvalue("alphalrcdyn"))
  {
   $atthashtddft{"alphalrcdyn"}=getvalue("alphalrcdyn");
 
    }
   
if(getvalue("betalrcdyn"))
  {
   $atthashtddft{"betalrcdyn"}=getvalue("betalrcdyn");
 
    }
   
if(getvalue("mdfqtype"))
  {
   $atthashtddft{"mdfqtype"}=getvalue("mdfqtype");
 
    }
   
if(getvalue("fxctype"))
  {
  $atthashtddft{"fxctype"}=$enumhashfxctype{getvalue("fxctype")};

    }
   
my %atthashscreening=();
   
if(getvalue("nosymscr"))
  {
   $atthashscreening{"nosym"}=getvalue("nosymscr");
 
    }
   
if(getvalue("ngridkscr"))
  {
   $atthashscreening{"ngridk"}=getvalue("ngridkscr");
 
    }
   
if(getvalue("reducekscr"))
  {
   $atthashscreening{"reducek"}=getvalue("reducekscr");
 
    }
   
if(getvalue("vkloffscr"))
  {
   $atthashscreening{"vkloff"}=getvalue("vkloffscr");
 
    }
   
if(getvalue("rgkmaxscr"))
  {
   $atthashscreening{"rgkmax"}=getvalue("rgkmaxscr");
 
    }
   
if(getvalue("nemptyscr"))
  {
   $atthashscreening{"nempty"}=getvalue("nemptyscr");
 
    }
   
if(getvalue("screentype"))
  {
  $atthashscreening{"screentype"}=$enumhashscreentype{getvalue("screentype")};

    }
   
my %atthashBSE=();
   
if(getvalue("nosymbse"))
  {
   $atthashBSE{"nosym"}=getvalue("nosymbse");
 
    }
   
if(getvalue("reducekbse"))
  {
   $atthashBSE{"reducek"}=getvalue("reducekbse");
 
    }
   
if(getvalue("vkloffbse"))
  {
   $atthashBSE{"vkloff"}=getvalue("vkloffbse");
 
    }
   
if(getvalue("rgkmaxbse"))
  {
   $atthashBSE{"rgkmax"}=getvalue("rgkmaxbse");
 
    }
   
if(getvalue("scrherm"))
  {
   $atthashBSE{"scrherm"}=getvalue("scrherm");
 
    }
   
if(getvalue("fbzq"))
  {
   $atthashBSE{"fbzq"}=getvalue("fbzq");
 
    }
   
if(getvalue("sciavtype"))
  {
  $atthashBSE{"sciavtype"}=$enumhashsciavtype{getvalue("sciavtype")};

    }
   
if(getvalue("sciavbd"))
  {
   $atthashBSE{"sciavbd"}=getvalue("sciavbd");
 
    }
   
if(getvalue("sciavqhd"))
  {
   $atthashBSE{"sciavqhd"}=getvalue("sciavqhd");
 
    }
   
if(getvalue("sciavqwg"))
  {
   $atthashBSE{"sciavqwg"}=getvalue("sciavqwg");
 
    }
   
if(getvalue("sciavqbd"))
  {
   $atthashBSE{"sciavqbd"}=getvalue("sciavqbd");
 
    }
   
if(getvalue("bsedirsing"))
  {
   $atthashBSE{"bsedirsing"}=getvalue("bsedirsing");
 
    }
   
if(getvalue("lmaxdielt"))
  {
   $atthashBSE{"lmaxdielt"}=getvalue("lmaxdielt");
 
    }
   
if(getvalue("nleblaik"))
  {
   $atthashBSE{"nleblaik"}=getvalue("nleblaik");
 
    }
   
if(getvalue("nexcitmax"))
  {
   $atthashBSE{"nexcitmax"}=getvalue("nexcitmax");
 
    }
   
if(getvalue("nbfbse,nafbse"))
  {
   $atthashBSE{"nstlbse"}=getvalue("nbfbse,nafbse");
 
    }
   
if(getvalue("nbfce,nafce"))
  {
   $atthashBSE{"nstlce"}=getvalue("nbfce,nafce");
 
    }
   
if(getvalue("bsetype"))
  {
  $atthashBSE{"bsetype"}=$enumhashbsetype{getvalue("bsetype")};

    }
   
my %atthashtetra=();
   
if(getvalue("tetraocc"))
  {
   $atthashtetra{"tetraocc"}=getvalue("tetraocc");
 
    }
   
if(getvalue("tetradf"))
  {
   $atthashtetra{"tetradf"}=getvalue("tetradf");
 
    }
   
if(getvalue("tetrakordexc"))
  {
   $atthashtetra{"kordexc"}=getvalue("tetrakordexc");
 
    }
   
if(getvalue("tetracw1k"))
  {
   $atthashtetra{"cw1k"}=getvalue("tetracw1k");
 
    }
   
if(getvalue("tetraqweights"))
  {
   $atthashtetra{"qweights"}=getvalue("tetraqweights");
 
    }
   
my %atthashdosWindow=();
   
if(getvalue("points"))
  {
   $atthashdosWindow{"points"}=getvalue("points");
 
    }
   
if(getvalue("intv"))
  {
   $atthashdosWindow{"intv"}=getvalue("intv");
 
    }
   
if(getvalue("noreplace"))
  {
   $atthashdosWindow{"nsmdos"}=getvalue("noreplace");
 
    }
   
my %atthashdoonly=();
   
if(getvalue("noreplace"))
  {
  $atthashdoonly{"task"}=$enumhashtask{getvalue("noreplace")};

    }
   
my $output = new IO::File(">input.xml");
my $writer = new XML::Writer(DATA_MODE => 1);
$writer->xmlDecl("UTF-8");
$writer->pi('xml-stylesheet', 'href="xmlinput2html.xsl" type="text/xsl"');

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
	    $atthashgroundstate{"do"}="fromscratch";
	    }
	  if($1==1||$1==3)
	    {
	    $atthashgroundstate{"do"}="fromfile";
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
	    $writer->startTag("relax",%atthashrelax);
	    $writer->endTag("relax");
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
