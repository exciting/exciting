<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
 xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:str="http://exslt.org/strings"
 xmlns:ex="http://xml.exciting-code.org/inputschemaextentions.xsd">
 <xsl:output method="text"/>

 <xsl:param name="importancelevels">
  <xsl:text>essential</xsl:text>
  <xs:annotation>
   <xs:documentation> In order to select the importance levels that should be included list them in
    the parameter "importancelevels". example: xsltproc --stringparam importancelevels "essential
    expert" schematolatex.xsl excitinginput.xsd >doc.tex </xs:documentation>
  </xs:annotation>
 </xsl:param>

 <xsl:template match="/">


  <xsl:text>

\documentclass{article}

\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{underscore}

\bibliographystyle {plain}
\errorstopmode
\usepackage{hyperref}
\usepackage{color}
\usepackage{longtable}
\hypersetup{colorlinks=true}
\begin{document}
\setlength{\LTleft}{0pt}
\newcommand{\vect}[1]{\mathbf{ #1}} 
\newcommand{\op}[1]{\mathbf {#1}}
\newcommand{\bra}[1]{\ensuremath{\left\langle #1\right|}}
\newcommand{\ket}[1]{\ensuremath{\left|#1\right\rangle}}
\newcommand{\braket}[2]{\ensuremath{\left\langle #1\vphantom{#2}\right.\left|\vphantom{#1}#2\right\rangle}}
\newcommand{\scalapack}{SCALAPACK }
\newcommand{\blas}[0]{BLAS }
\newcommand{\lapack}{LAPACK }
\newcommand{\arpack}{ARPACK }
\newcommand{\subsubsubsection}[1]{\paragraph{#1} \paragraph*{} }
\newcommand{\attref}[2]{{\tt \hyperref[#2att#1]{\color{green} #1}}}
\newcommand{\elementref}[1]{{\tt  \hyperref[#1]{\color{blue}  #1}}}
\newcommand{\exciting}{ {\usefont{T1}{lmtt}{b}{n} exciting} }

\begin{titlepage}
\begin{center}
\includegraphics{../../xml/exciting.jpg}
\ \\
\ \\
\ \\
\ \\
\ \\
\Huge \textbf{</xsl:text><xsl:apply-templates select="/xs:schema/xs:annotation/xs:appinfo/title"/><xsl:text>} \\
\ \\
\huge \textbf{\exciting \texttt{carbon}} \\
\Large
\vfill
\ \\
December 2015
\end{center}
\end{titlepage}

\newpage
\definecolor{green}{rgb}{0,0.5,0}
\section*{About this Document}
    </xsl:text>
  <xsl:apply-templates select="/xs:schema/xs:annotation/xs:documentation"/> \part{Input Elements}
   <xsl:call-template name="elementToLatex">
   <xsl:with-param name="myelement"
    select="//xs:element[@name=/xs:schema/xs:annotation[last()]/xs:appinfo/root]"/>
   <xsl:with-param name="level" select="0"/>
  </xsl:call-template>
  <xsl:text>\part{Reused Elements}
  The following elements can occur more than once in the input file. Therefore they are listed separately.
  </xsl:text>
  <xsl:for-each  
   select="/*/xs:element[@name!=/xs:schema/xs:annotation/xs:appinfo/root and contains($importancelevels,@ex:importance)]">
   <xsl:variable name="name" select="@name"/>
   <xsl:if test="count(//xs:element[@ref=$name])>1">
    <xsl:call-template name="elementToLatex">
     <xsl:with-param name="myelement" select="."/>
     <xsl:with-param name="level" select="0"/>
    </xsl:call-template>
   </xsl:if>
  </xsl:for-each>
  <xsl:text>\section{Data Types}
 
 The Input definition uses derived data types. These are described here.
  </xsl:text>
  <xsl:for-each select="/*/xs:simpleType">
   <xsl:call-template name="typetoDoc">
    <xsl:with-param name="typenode" select="."/>
   </xsl:call-template>
  </xsl:for-each>
  <xsl:text>
    \end{document}
    </xsl:text>
 </xsl:template>
 <xsl:template match="displaymath">
  <xsl:text> 
\begin{equation}
</xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
 
  <xsl:text>
\end{equation}
</xsl:text>
 </xsl:template>
 <xsl:template match="inlinemath">
  <xsl:text>$ </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
  <xsl:text> $</xsl:text>
 </xsl:template>
 <xsl:template match="pre">
  <xsl:text>{\tt </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="str:replace(., '_','\_')"/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
 </xsl:template>

 <xsl:template match="pre-bf">
  <xsl:text>{\usefont{T1}{lmtt}{b}{n} </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
 </xsl:template>

 <xsl:template match="blue">
  <xsl:text>{\usefont{T1}{lmtt}{b}{n} \color{blue} </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
 </xsl:template>

 <xsl:template match="green">
  <xsl:text>{\usefont{T1}{lmtt}{b}{n} \color{green} </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
 </xsl:template>

 <xsl:template match="it">
  <xsl:text>{\it </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
 </xsl:template>
 <xsl:template match="bf">
  <xsl:text>{\bf </xsl:text>
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
  <xsl:text>}</xsl:text>
 </xsl:template>
 <xsl:template match="text()">
  <xsl:call-template name="normalizespace">
   <xsl:with-param name="a" select="."/>
  </xsl:call-template>
 </xsl:template>
 <xsl:template match="p">
  <xsl:text>
  
  </xsl:text>
  <xsl:apply-templates select="./*|text()"/>
 </xsl:template>
 <xsl:template match="xs:documentation">
  <xsl:apply-templates
   select="text()|inlinemath|displaymath|pre|pre-bf|bf|blue|green|pre_ns|it|it_ns|p|exciting|a|list|li|attref|filename|filename_ns|elementref|elementref_ns"
  />
 </xsl:template>
 <xsl:template match="elementref">
  <xsl:text>\elementref{</xsl:text>
  <xsl:value-of select="."/>
  <xsl:text>}</xsl:text>
 </xsl:template>
 

 <xsl:template match="filename">
  <xsl:text>{\usefont{T1}{lmtt}{b}{n}  </xsl:text>
  <xsl:value-of select="str:replace(.,'_','\_')"/>
  <xsl:text>}</xsl:text>
 </xsl:template>
<xsl:template name="attref">
 <xsl:param name="att"/>
 <xsl:param name="parent" select="''"/>
  <xsl:text>\attref{</xsl:text>
  <xsl:value-of select="$att"/>

  <xsl:text>}</xsl:text>
 
  <xsl:text>{</xsl:text>
  <xsl:value-of select="$parent"/>
  <xsl:text>:}</xsl:text>
  
 </xsl:template>


 <xsl:template match="attref">
 <xsl:call-template name="attref">
  <xsl:with-param name="att" select="."/>
  <xsl:with-param name="parent" select="@parent"/>
 </xsl:call-template>
 </xsl:template>
 <xsl:template match="list">
  <xsl:text>
  \begin{itemize}
</xsl:text>
  <xsl:apply-templates select="./*|text()"/>
  <xsl:text>
  \end{itemize}
</xsl:text>
 </xsl:template>
 <xsl:template match="li">
  <xsl:text>\item </xsl:text>
  <xsl:apply-templates select="./*|text()"/>
  <xsl:text>
</xsl:text>
 </xsl:template>
 
 <xsl:template match="a">
  <xsl:text>{\usefont{T1}{lmtt}{b}{n} \color{black} </xsl:text>
  <xsl:value-of select="."/>

  <xsl:text> }(\url{</xsl:text>
  <xsl:value-of select="@href"/>
  <xsl:text>})</xsl:text>

 </xsl:template>
 <xsl:template match="exciting">
  <xsl:text>\exciting{}</xsl:text>
 </xsl:template>
 <xsl:template name="elementToLatex">
  <xsl:param name="myelement"/>
  <xsl:param name="level"/>

  <xsl:text>\</xsl:text>
  <xsl:value-of select="str:padding(0,'sub')"/>
  <xsl:text>section{ Element: {\textcolor{blue}{</xsl:text>
  <xsl:text/>
  <xsl:value-of select="$myelement/@name "/>
  <xsl:text>}}}
     \label{</xsl:text>
  <xsl:value-of select="$myelement/@name"/>
  <xsl:text>}
  </xsl:text>

  <xsl:apply-templates select="$myelement/xs:annotation/xs:documentation"/>
  <xsl:call-template name="TypeToDoc">
   <xsl:with-param name="contentnode" select="$myelement | //xs:element[@name=$myelement/@ref and contains($importancelevels,@ex:importance) ]"/>
  </xsl:call-template>
  <xsl:if test="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]"> This element
   <xsl:text>allows for specification of the following attributes: \begin{quotation}</xsl:text>
   <xsl:for-each
    select="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]">
    <xsl:sort select="@use='required'" order="descending"/>
    <xsl:sort select="@name|@ref"/>
    
    <xsl:call-template name="attref">
     <xsl:with-param name="att" select="@name|@ref"/>
      <xsl:with-param name="parent" select="../../@name"/>
     
    </xsl:call-template>
    <xsl:if test="@use='required'">
     <xsl:text> \nolinebreak {\color{red}(required)}</xsl:text>
    </xsl:if>
    <xsl:if test="position()!=last()">
     <xsl:text>, </xsl:text>
    </xsl:if>
   </xsl:for-each>
   <xsl:text>\end{quotation}</xsl:text>
  </xsl:if>
  <xsl:for-each select="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]">
   <xsl:sort select="@name|@ref"/>
   <xsl:call-template name="attributetolatex">
    <xsl:with-param name="myattribute" select="."/>
    <xsl:with-param name="level" select="$level"/>
   </xsl:call-template>
  </xsl:for-each>
  <xsl:for-each
   select="$myelement/*/*/xs:element[contains($importancelevels,@ex:importance) and @name]">
   <xsl:call-template name="elementToLatex">
    <xsl:with-param name="myelement" select="."/>
    <xsl:with-param name="level" select="$level+1"/>
   </xsl:call-template>
  </xsl:for-each>
  <xsl:for-each
   select="$myelement/*/*/xs:element[contains($importancelevels,@ex:importance) and @ref]">
   <xsl:variable name="ref" select="@ref"/>
   <xsl:if test="count(//xs:element[@ref=$ref])=1">
    <xsl:call-template name="elementToLatex">
     <xsl:with-param name="myelement" select="//xs:element[@name=$ref]"/>
     <xsl:with-param name="level" select="$level+1"/>
    </xsl:call-template>
   </xsl:if>
  </xsl:for-each>
 </xsl:template>
 <xsl:template name="attributetolatex">
  <xsl:param name="myattribute"/>
  <xsl:param name="level"/>
  <xsl:text>\subsection{Attribute: {\textcolor{green}{{</xsl:text>}
  <xsl:value-of select="$myattribute/@name |$myattribute/@ref"/>
  <xsl:text>}}}  \label{:att</xsl:text>
  <xsl:value-of select="$myattribute/@name |$myattribute/@ref"/>
  <xsl:text>}
  \label{</xsl:text>
  <xsl:value-of select="$myattribute/../../@name"/>
  <xsl:text>:att</xsl:text>
  <xsl:value-of select="$myattribute/@name |$myattribute/@ref"/>
  <xsl:text>}
    </xsl:text>
  <xsl:apply-templates select="$myattribute/xs:annotation/xs:documentation"/>
  <xsl:call-template name="TypeToDoc">
   <xsl:with-param name="contentnode"
    select="$myattribute | //xs:attribute[@name=$myattribute/@ref]"/>
  </xsl:call-template>
 </xsl:template>
 <xsl:template name="TypeToDoc">
  <xsl:param name="contentnode"/>
  <xsl:text>

\begin{longtable}{ll}


 </xsl:text>
  <xsl:choose>
   <xsl:when test="$contentnode/@type">
    <xsl:text> \bf{Type:} &amp; </xsl:text>
    <xsl:value-of select="str:replace(($contentnode/@type),'xs:','')"/>
    <xsl:if test="not(contains($contentnode/@type,'xs:'))">
     <xsl:text> (\ref{</xsl:text>
     <xsl:value-of select="$contentnode/@type"/>
     <xsl:text>}) 
</xsl:text>
    </xsl:if>
    <xsl:text>\\</xsl:text>
   </xsl:when>
   <xsl:when test="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
    <xsl:text> \bf{Type:} &amp; \bf{choose from:}  \\
</xsl:text>
    <xsl:for-each
     select="
 $contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
     <xsl:text> &amp; {\tt </xsl:text>
     <xsl:value-of select="str:replace(@value, '_','\_')"/>
     <xsl:text/>
     <xsl:text>}  \\
</xsl:text>
    </xsl:for-each>
    <xsl:text/>
   </xsl:when>
   <xsl:when test="$contentnode/xs:complexType/*[xs:element] ">
    <xsl:text>\bf{Contains: }  </xsl:text>
    <xsl:text>  </xsl:text>
    <xsl:for-each
     select="$contentnode/xs:complexType/*/xs:element[contains($importancelevels,@ex:importance)]">
     <xsl:text> &amp; \elementref{</xsl:text>
     <xsl:value-of select="./@name|@ref"/><xsl:text>}</xsl:text>
     <xsl:if test="@minOccurs=0">
      <xsl:choose>
       <xsl:when test="@maxOccurs='unbounded'">
        <xsl:text> (zero or more)</xsl:text>

       </xsl:when>
       <xsl:otherwise>
        <xsl:text> (optional)</xsl:text>
       </xsl:otherwise>
      </xsl:choose>
     </xsl:if>

     <xsl:if test="@minOccurs&gt;0">
      <xsl:text> (</xsl:text>
      <xsl:value-of select="@minOccurs"/>
      <xsl:text> times</xsl:text>
      <xsl:if test="@maxOccurs='unbounded'">
       <xsl:text> or more</xsl:text>
      </xsl:if>
      <xsl:text>) </xsl:text>
     </xsl:if>

   
      
     <xsl:text>  \\
</xsl:text>
    </xsl:for-each>
   </xsl:when>
   <xsl:otherwise>
    <xsl:text> \bf{Type:} &amp; no content \\
</xsl:text>
   </xsl:otherwise>
  </xsl:choose>
  <xsl:if test="$contentnode/@default">
   <xsl:text>
 \bf{Default:} &amp; ''{\tt </xsl:text>
   <xsl:value-of select="str:replace($contentnode/@default,'_','\_')"/>
   <xsl:text>}''\\
 </xsl:text>
  </xsl:if>
  <xsl:if test="$contentnode/@use or local-name($contentnode)='attribute'">
   <xsl:text>
 \bf{Use:} &amp; </xsl:text>
   <xsl:value-of select="$contentnode/@use"/>
   <xsl:if test="not($contentnode/@use) and local-name($contentnode)='attribute'">
    <xsl:text>optional</xsl:text>
   </xsl:if>
   <xsl:text> \\
 </xsl:text>
  </xsl:if>
  <xsl:choose>
   <xsl:when test="$contentnode/@ex:unit!=''">

    <xsl:text>
 \bf{Unit:}&amp;</xsl:text>
    <xsl:value-of select="$contentnode/@ex:unit"/>
    <xsl:text>  \\
  </xsl:text>
   </xsl:when>
  </xsl:choose>

  <xsl:text>\bf{XPath:}&amp; {\tt</xsl:text>
  <xsl:call-template name="genxpath">
   <xsl:with-param name="node" select="$contentnode"/>
   <xsl:with-param name="xpath" select="''"/>
  </xsl:call-template>
  <xsl:text>} \\
  </xsl:text>
  <xsl:if test="$contentnode/@name and name(..)='xs:schema'">
   <xsl:variable name="name" select="$contentnode/@name"/>
   <xsl:text>\bf{Parent:}  </xsl:text>
   <xsl:for-each
    select="//xs:element[@ref=$name  and contains($importancelevels,@ex:importance)]">
    <xsl:text> &amp; {\tt </xsl:text>
    <xsl:call-template name="genxpath">
     <xsl:with-param name="node" select="."/>
     <xsl:with-param name="xpath" select="''"/>
    </xsl:call-template>
    <xsl:text> } \\
</xsl:text>
   </xsl:for-each>
   
  </xsl:if>
  
  <xsl:text>

\end{longtable}
  </xsl:text>
 </xsl:template>
 <xsl:template name="genxpath">
  <xsl:param name="xpath"/>
  <xsl:param name="node"/>
  <xsl:variable name="current_name">
   <xsl:if test="$node/@name">
    <xsl:text>/\hyperref[</xsl:text>
    <xsl:if test="name($node)='xs:attribute'">
     <xsl:value-of select="../../@name"/>
     <xsl:text>:att</xsl:text>
      
    </xsl:if>
    <xsl:value-of select="$node/@name|$node/@ref"/>
 <xsl:text>]{</xsl:text>
   <xsl:if test="name($node)='xs:attribute'">
    <xsl:text>@</xsl:text>
   </xsl:if>
   <xsl:value-of select="$node/@name|$node/@ref"/>
   <xsl:text>}</xsl:text>
   </xsl:if>
   <xsl:value-of select="$xpath"/>
  </xsl:variable>
  <xsl:for-each select="$node[last()]">
   <xsl:variable name="name" select="$node/@name"/>
    <xsl:variable name="rootel" select="/xs:schema/xs:annotation/xs:appinfo/root"/>
   <xsl:choose>
    <xsl:when test="name(..)='xs:schema' and count(//xs:element[@ref=$name])=1">
     <xsl:call-template name="genxpath">
      <xsl:with-param name="node" select="//xs:element[*/*/xs:element[@ref=$name ]]"/>
      <xsl:with-param name="xpath">
       <xsl:value-of select="$current_name"/>
      </xsl:with-param>
     </xsl:call-template>
    </xsl:when>
    <xsl:when test="parent::node() ">
     <xsl:for-each select="parent::node()">
      <xsl:call-template name="genxpath">
       <xsl:with-param name="node" select="."/>
       <xsl:with-param name="xpath">
        <xsl:value-of select="$current_name"/>
       </xsl:with-param>
      </xsl:call-template>
     </xsl:for-each>
    </xsl:when>
    <xsl:when test="contains($xpath,$rootel)">
     <xsl:value-of select="$xpath"/>
    </xsl:when>
     
    <xsl:otherwise>
     <xsl:text>.</xsl:text>
     <xsl:value-of select="$xpath"/>
    </xsl:otherwise>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>
 <xsl:template name="typetoDoc">
  <xsl:param name="typenode"/>
  <xsl:text>
  \subsection{Type   </xsl:text>
  <xsl:value-of select="$typenode/@name"/>
  <xsl:text>
  } 
  \label{</xsl:text>
  <xsl:value-of select="$typenode/@name"/>
  <xsl:text>}
  </xsl:text>
  <xsl:apply-templates select="xs:annotation/xs:documentation"/>
 </xsl:template>
 <xsl:template name="normalizespace">
  <xsl:param name="a"/>
  <xsl:if test="substring($a,1,1)=' ' or substring($a,1,1)=$newline">
   <xsl:text> </xsl:text>
  </xsl:if>
  <xsl:value-of select="normalize-space($a)"/>
  <xsl:if test="substring($a,string-length($a),1)=' ' or substring($a,string-length($a),1)=$newline">
   <xsl:text> </xsl:text>
  </xsl:if>

 </xsl:template>
 <xsl:variable name="newline">
  <xsl:text>
</xsl:text>
 </xsl:variable>
</xsl:stylesheet>
