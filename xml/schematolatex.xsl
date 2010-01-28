<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:str="http://exslt.org/strings"
  xmlns:ex="inputschemaextentions.xsd">
  <xsl:output method="text"></xsl:output>
  <xsl:variable name="statuslevels">
    <xsl:text>essential</xsl:text>
  </xsl:variable>
  <xsl:template match="/">
    <xsl:text>

    \documentclass{article}



\usepackage{amsmath}

\bibliographystyle {plain}
\errorstopmode
\usepackage{hyperref}
\hypersetup{colorlinks=false}
\begin{document}
\newcommand{\exciting}{EXC!T!`NG }
\title{</xsl:text>
    <xsl:value-of select="/xs:schema/xs:annotation/xs:appinfo/title" />
    <xsl:text>} 
\author{exciting devteam}


\maketitle \newpage \tableofcontents

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
\newpage
    </xsl:text>
    <xsl:call-template name="elementToLatex">
      <xsl:with-param name="myelement" select="//xs:element[@name=/xs:schema/xs:annotation/xs:appinfo/root]" />
      <xsl:with-param name="level" select="0" />
    </xsl:call-template>
       <xsl:for-each select="/*/xs:element[@name!=/xs:schema/xs:annotation/xs:appinfo/root and contains($statuslevels,@ex:status)]">
       <xsl:call-template name="elementToLatex">
      <xsl:with-param name="myelement" select="." />
      <xsl:with-param name="level" select="0" />
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
    <xsl:value-of select="normalize-space(.)"></xsl:value-of>
    <xsl:text>
\end{equation}
</xsl:text>
  </xsl:template>
  <xsl:template match="inlinemath">
    <xsl:text> $ </xsl:text>
    <xsl:value-of select="normalize-space(.)"></xsl:value-of>
    <xsl:text> $</xsl:text>
  </xsl:template>
  <xsl:template match="pre">
    <xsl:text> {\tt </xsl:text>
    <xsl:value-of select="normalize-space(.)"></xsl:value-of>
    <xsl:text> }</xsl:text>
  </xsl:template>
  <xsl:template match="it">
    <xsl:text> {\it </xsl:text>
    <xsl:value-of select="normalize-space(.)"></xsl:value-of>
    <xsl:text> }</xsl:text>
  </xsl:template>
  <xsl:template match="bf">
    <xsl:text> {\bf</xsl:text>
    <xsl:value-of select="normalize-space(.)"></xsl:value-of>
    <xsl:text>} </xsl:text>
  </xsl:template>
  <xsl:template match="text()">
    <xsl:value-of select="normalize-space(.)" />
  </xsl:template>
  <xsl:template match="xs:documentation">
    <xsl:apply-templates select="text()|inlinemath|displaymath|pre|it" />
  </xsl:template>
  <xsl:template name="elementToLatex">
    <xsl:param name="myelement"></xsl:param>
    <xsl:param name="level"></xsl:param>
    
    <xsl:text>\</xsl:text>
    <xsl:value-of select="str:padding(0,'sub')" />
    <xsl:text>section{ Element:</xsl:text>
    <xsl:text> </xsl:text>
    <xsl:value-of select="$myelement/@name " />
    <xsl:text>}
     \label{</xsl:text> <xsl:value-of select="$myelement/@name"/>
     <xsl:text>}
  </xsl:text>
 
    <xsl:apply-templates select="$myelement/xs:annotation/xs:documentation"></xsl:apply-templates>
   <xsl:call-template name="TypeToDoc">
      <xsl:with-param name="contentnode" select="$myelement | //xs:element[@name=$myelement/@ref]" />
    </xsl:call-template>
    <xsl:for-each select="$myelement/*/xs:attribute[contains($statuslevels,@ex:status)]">
      <xsl:call-template name="attributetolatex">
        <xsl:with-param name="myattribute" select="." />
        <xsl:with-param name="level" select="$level" />
      </xsl:call-template>
    </xsl:for-each>
    <xsl:for-each select="$myelement/*/*/xs:element[contains($statuslevels,@ex:status)and @name]">
      <xsl:call-template name="elementToLatex">
        <xsl:with-param name="myelement" select="." />
        <xsl:with-param name="level" select="$level+1" />
      </xsl:call-template>
    </xsl:for-each>
  </xsl:template>
  <xsl:template name="attributetolatex">
    <xsl:param name="myattribute"></xsl:param>
    <xsl:param name="level"></xsl:param>
    <xsl:text>\subsection{@</xsl:text>
    <xsl:value-of select="$myattribute/@name |$myattribute/@ref" />
    <xsl:text>}  
    </xsl:text>
    <xsl:apply-templates select="$myattribute/xs:annotation/xs:documentation"></xsl:apply-templates>
    <xsl:call-template name="TypeToDoc">
      <xsl:with-param name="contentnode" select="$myattribute | //xs:attribute[@name=$myattribute/@ref]" />
    </xsl:call-template>
  </xsl:template>
  <xsl:template name="TypeToDoc">
    <xsl:param name="contentnode" />
    <xsl:text>

  \begin{center}
\begin{tabular*}{\textwidth}{ll}

 \bf{Type:} &amp; </xsl:text>
 <xsl:choose>
 <xsl:when test="$contentnode/@type">
 <xsl:value-of select="$contentnode/@type"/>
 <xsl:text>\\
</xsl:text>
 </xsl:when>
 <xsl:when test="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
<xsl:text> \bf{enumeration:}  \\
</xsl:text>
 <xsl:for-each select="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
<xsl:text> &amp; </xsl:text>
<xsl:value-of select="@value"/><xsl:text> </xsl:text>
<xsl:text>  \\
</xsl:text>
 </xsl:for-each>
 <xsl:text> </xsl:text>
 </xsl:when>
 <xsl:when test="$contentnode/xs:complexType/*[xs:element] ">
 <xsl:text>\bf{contains: }</xsl:text> 
 <xsl:text>  \\
</xsl:text>
 <xsl:for-each select="$contentnode/xs:complexType/*/xs:element[contains($statuslevels,@ex:status)]">
 <xsl:text> &amp; </xsl:text>
 <xsl:value-of select="./@name|@ref"/>
<xsl:if test="@minOccurs=0">
<xsl:text>(optional)</xsl:text>
</xsl:if>
<xsl:if test="@minOccurs&gt;0">
<xsl:text>(min </xsl:text>
<xsl:value-of select="@minOccurs"></xsl:value-of>
<xsl:text> times) </xsl:text>
</xsl:if>
<xsl:text>
See: \ref{</xsl:text>
 <xsl:value-of select="./@name|@ref"/>
<xsl:text>}</xsl:text>
 <xsl:text>  \\
</xsl:text>
 </xsl:for-each>
 </xsl:when>
 <xsl:otherwise>
 <xsl:text> no content \\
</xsl:text>
 </xsl:otherwise>
 </xsl:choose>

 <xsl:choose>
 <xsl:when test="$contentnode/@ex:unit!=''">
 <xsl:text>
 \bf{Unit:}&amp;</xsl:text>
 <xsl:value-of select="$contentnode/@ex:unit"></xsl:value-of>
  <xsl:text>  \\
  </xsl:text>
 </xsl:when>
 </xsl:choose>

 <xsl:text>\bf{XPath:}&amp;</xsl:text>
   <xsl:call-template name="genxpath" >
  <xsl:with-param name="node" select="$contentnode"/>
  <xsl:with-param name="xpath" select="''"/>
  </xsl:call-template>
  <xsl:text>\\
  

\end{tabular*}
\end{center}
  </xsl:text>
  </xsl:template>
   <xsl:template name="genxpath">
    <xsl:param name="xpath"/>
    <xsl:param name="node"/>
    <xsl:variable name="current_name">
     <xsl:if test="$node/@name">
       <xsl:text>/</xsl:text>
            <xsl:if test="name($node)='xs:attribute'">
              <xsl:text>@</xsl:text>
            </xsl:if>
            <xsl:value-of select="$node/@name|$node/@ref" />
          </xsl:if>
          <xsl:value-of select="$xpath" />
    </xsl:variable>
    <xsl:for-each  select="$node[last()]">
    <xsl:choose>
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
         <xsl:when test="contains($xpath,'input')">
         <xsl:value-of select="$xpath"/>
         </xsl:when>
         <xsl:otherwise>
         <xsl:text>/</xsl:text>
         <xsl:value-of select="$xpath"/>
         </xsl:otherwise>
         </xsl:choose>
         </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>