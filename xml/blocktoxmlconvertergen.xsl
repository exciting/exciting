<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema">
  <xsl:output method="text"></xsl:output>
  <xsl:template match="/">
    use lib "../test/perl/"; use lib "../test/perl/lib"; use lib "../build/utilities/lib"; use XML::Writer; use IO::File; use Switch; $npt="[-+]?[0-9]*\.?[0-9]*";
    local $/ = undef; open BLOCKINPUT, $ARGV[0] or die "Couldn't open file:$ARGV[0] $!"; binmode BLOCKINPUT; $blockinput = &lt;BLOCKINPUT&gt;."\n\n";
$blockinput=~s/[!:](.*)\n/\n/g;
$/="\n";
<xsl:for-each select="//xs:restriction[xs:enumeration]">
<xsl:for-each select="xs:enumeration">
<xsl:text>$enumhash</xsl:text><xsl:value-of select="../../../@name"/>
<xsl:text>{"</xsl:text><xsl:value-of select="xs:annotation/xs:appinfo/oldnr"/>
<xsl:text>"}="</xsl:text><xsl:value-of select="@value"/><xsl:text>";
</xsl:text>
</xsl:for-each>
</xsl:for-each>
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
<xsl:for-each select="//xs:element[*/xs:attribute]">
my %atthash<xsl:value-of select="@name"/>=();
   <xsl:for-each select="*/xs:attribute">
if(getvalue("<xsl:value-of select="(@name|@ref|*/*/oldname)[last()]"/>"))
  {
  <xsl:choose>
  <xsl:when test="*/xs:restriction/xs:enumeration">
   <xsl:text>$atthash</xsl:text>
   <xsl:value-of select="../../@name"/>
   <xsl:text>{"</xsl:text>
   <xsl:value-of select="@name|@ref"/>
   <xsl:text>"}=$enumhash</xsl:text>
   <xsl:value-of select="@name|@ref"/>
   <xsl:text>{getvalue("</xsl:text>
   <xsl:value-of select="(@name|@ref|*/*/oldname)[last()]"/>
   <xsl:text>")};
</xsl:text>
  </xsl:when>
  <xsl:otherwise>
 <xsl:text> $atthash</xsl:text>
 <xsl:value-of select="../../@name"/>
 <xsl:text>{"</xsl:text>
 <xsl:value-of select="@name|@ref"/>
 <xsl:text>"}=getvalue("</xsl:text>
 <xsl:value-of select="(@name|@ref|*/*/oldname)[last()]"/>
 <xsl:text>");
 </xsl:text>
 </xsl:otherwise>
  </xsl:choose>
    }
   </xsl:for-each>
</xsl:for-each>
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
  for(my $is=1;$is&lt;=$nspecies;$is++)
    {
    $string=~s/(.*)\n//;
    $spfname=$1;
    $spfname=~s/['\s]//g;
    $spfname=~s/\.in/\.xml/g;
    $writer->startTag("species","speciesfile"=>$spfname);
    $string=~s/(.*)\n//;
    $natoms=$1;
    for(my $ia=1;$ia&lt;=$natoms;$ia++)
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
</xsl:template>
</xsl:stylesheet>