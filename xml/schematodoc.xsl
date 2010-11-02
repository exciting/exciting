<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
	xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema"
    xmlns:str="http://exslt.org/strings"
     xmlns:ex="inputschemaextentions.xsd"
  >
	<xsl:output method="text" />
	<xsl:variable name="statuslevels">
  <xsl:text>essential</xsl:text>
  </xsl:variable>
  <xsl:variable name="newline">
		<xsl:text>
		</xsl:text>
	</xsl:variable>
<xsl:template match="displaymath">
<xsl:text> 
\begin{equation}
</xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text>
\end{equation}
</xsl:text>
  </xsl:template>
   <xsl:template match="inlinemath">
 <xsl:text> $ </xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text> $</xsl:text>
  </xsl:template>
    <xsl:template match="pre">
 <xsl:text> {\tt </xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text> }</xsl:text>
  </xsl:template>
    <xsl:template match="it">
 <xsl:text> {\it </xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text> }</xsl:text>
  </xsl:template>
   <xsl:template match="bf">
 <xsl:text> {\bf</xsl:text> <xsl:value-of select="normalize-space(.)"></xsl:value-of><xsl:text>} </xsl:text>
  </xsl:template>
   <xsl:template match="text()">
 <xsl:value-of select="normalize-space(.)"/>
  </xsl:template>
  
   <xsl:template match="xs:documentation">
<xsl:apply-templates select="text()|inlinemath|displaymath|pre|it"/>
  </xsl:template>
	<xsl:template name="atributdescriptions" match="none">
  <xsl:if test="./@ex:status and contains($statuslevels,./@ex:status)">
	<!-- <xsl:apply-templates select=" contains($statuslevels,xs:status) and ./xs:annotation/xs:documentation" />
	<xsl:apply-templates select="contains($statuslevels,xs:status) and ./xs:complexType/xs:annotation/xs:documentation" />
 -->
 <xsl:apply-templates select="./xs:annotation/xs:documentation" />
  <xsl:apply-templates select="./xs:complexType/xs:annotation/xs:documentation" />
		<xsl:value-of select="$newline" />
		<xsl:if test="./xs:complexType/xs:attribute|./xs:attribute">
			\paragraph{Attributes:}
			<xsl:text>\begin{description}</xsl:text>
			<xsl:for-each select="./xs:complexType/xs:attribute|./xs:attribute">
				<xsl:text>\item[</xsl:text>
				<xsl:value-of select="./@name" />
				<xsl:value-of select="./@ref" />
				<xsl:text>]</xsl:text>
				<xsl:if test="./@ref">
				<xsl:apply-templates select="//*[@name=./@ref]/*/xs:documentation"/>
					<xsl:text> See \ref{</xsl:text>
					<xsl:value-of select="./@ref" />
					<xsl:text>}.  </xsl:text>
					
				</xsl:if>
				<xsl:if test="./@name">
					<xsl:text>\label{</xsl:text>
					<xsl:value-of select="./@name" />
					<xsl:text>}</xsl:text>
				</xsl:if>

				<xsl:apply-templates select="./*/xs:documentation" />
				
				<xsl:value-of select="$newline" />
						<xsl:if test="./xs:simpleType/xs:restriction/xs:enumeration/@value|./@type">	
					<xsl:text>
	
	\begin{description}</xsl:text>

	<xsl:text>\item[type] </xsl:text>

<xsl:if test="./xs:simpleType/xs:restriction/xs:enumeration/@value">
<xsl:text> select:
\begin{itemize}
</xsl:text>
<xsl:for-each select="./xs:simpleType/xs:restriction/xs:enumeration/@value  ">
<xsl:text>\item </xsl:text><xsl:value-of select="str:replace(.,'_','\_')"/>
</xsl:for-each>
<xsl:text> \end{itemize}</xsl:text>
</xsl:if>
					
					<xsl:value-of select="./@type" />

					<xsl:if test="./@use">
						<xsl:text> \item[use] </xsl:text>
						<xsl:value-of select="./@use" />
					</xsl:if>
					<xsl:if test="./@default">
						<xsl:text>\item[default-value] </xsl:text>
						<xsl:value-of select="./@default" />
					</xsl:if>
					<xsl:text> \end{description} 
	</xsl:text>
			</xsl:if>
			</xsl:for-each>
			<xsl:for-each select="./*/xs:attributeGroup|./xs:attributeGroup">
				<xsl:if test="./@ref">
					<xsl:text> \item[attribute group </xsl:text>
					<xsl:value-of select="./@ref" />
				<xsl:apply-templates select="//*[@name=./@ref]/*/xs:documentation"/>
					<xsl:text>] See \ref{</xsl:text><xsl:value-of select="./@ref" />
					<xsl:text>} </xsl:text>
					
				</xsl:if>
			</xsl:for-each>
			<xsl:text>\end{description}
		</xsl:text>
		</xsl:if>
    </xsl:if>
	</xsl:template>


	<xsl:template name="section">
		<xsl:param name="depth" />
		<xsl:choose>
			<xsl:when test="$depth=0">
				<xsl:text>\section{</xsl:text>
			</xsl:when>
			<xsl:when test="$depth=1">
				<xsl:text>\subsection{</xsl:text>
			</xsl:when>
			<xsl:when test="$depth=2">
				<xsl:text>\subsubsection{</xsl:text>
			</xsl:when>
			<xsl:when test="$depth=3">
				<xsl:text>\paragraph{</xsl:text>
			</xsl:when>
			<xsl:when test="$depth=4">
				<xsl:text>\subparagraph{</xsl:text>
			</xsl:when>
		</xsl:choose>
		<xsl:value-of select="./@name" />
		<xsl:value-of select="./@ref" />
		<xsl:choose>
			<xsl:when test="name(.)='xs:attributeGroup'">
				<xsl:text> attribute group}</xsl:text>
			</xsl:when>
			<xsl:otherwise>
				<xsl:text> element}</xsl:text>
			</xsl:otherwise>
		</xsl:choose>
		<xsl:if test="./@ref">
			<xsl:text>See \ref{</xsl:text>
			<xsl:value-of select="./@ref" />
			<xsl:text>}</xsl:text>
		<xsl:apply-templates select="//*[@name=./@ref]/*/xs:documentation"/>
			
		</xsl:if>
		<xsl:if test="./@name">
			<xsl:text>\label{</xsl:text>
			<xsl:value-of select="./@name" />
			<xsl:text>}</xsl:text>
		</xsl:if>
		<xsl:call-template name="atributdescriptions" />
		<xsl:for-each select="./*/*/xs:element|./*/*/xs:group|./*/xs:group">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="$depth+1" />
			</xsl:call-template>
		</xsl:for-each>
	</xsl:template>
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
\title{</xsl:text><xsl:value-of select="/xs:schema/xs:annotation/xs:appinfo/title"/><xsl:text>} 
\author{\exciting developers team\\
(C. Ambrosch-Draxl, Zohreh Basirat, Thomas Dengg,\\
Christian Meisenbichler, Dmitrii Nabok, Weine Olovsson,\\
Pasquale Pavone, Stephan Sagmeister, J\"urgen Spitaler)}

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
\newpage
		</xsl:text>
		<xsl:text>\section{Root element </xsl:text>
		<xsl:value-of select="/xs:schema/xs:annotation/xs:appinfo/root"/>
		<xsl:text>}
		\label{input}
		</xsl:text>

		<xsl:for-each select="/*/xs:element[@name=/xs:schema/xs:annotation/xs:appinfo/root]">
			<xsl:call-template name="atributdescriptions" />
		</xsl:for-each>
		<xsl:for-each select="/*/xs:element[@name=/xs:schema/xs:annotation/xs:appinfo/root]/*/*/xs:element">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="0" />
			</xsl:call-template>
		</xsl:for-each>
		\section{reused Elements}
		These elements make sense in more than only one context. In this documentation there are references placed if one of these applies.
		<xsl:for-each select="/*/xs:element[@name!=/xs:schema/xs:annotation/xs:appinfo/root]|/*/xs:group">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="1" />
			</xsl:call-template>
		</xsl:for-each>
		<xsl:text>
	\section{reused attributes}</xsl:text>
These attributes make sense in more than only one context.  In this documentation there are references placed if one of these applies.
	<xsl:for-each select="/*">
			<xsl:call-template name="atributdescriptions" />
		</xsl:for-each>

		<xsl:for-each select="/*/xs:attributeGroup">
			<xsl:call-template name="section">
				<xsl:with-param name="depth" select="1" />
			</xsl:call-template>
		</xsl:for-each>
		<xsl:text>
	\end{document}</xsl:text>
	</xsl:template>
</xsl:stylesheet>
