<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:xs="http://www.w3.org/2001/XMLSchema"
 xmlns:str="http://exslt.org/strings" xmlns:ex="inputschemaextentions.xsd">
 <xsl:output method="text" />
 <xsl:variable name="importancelevels">
 <!-- In order to select the importance levels that should be included list them in the variable "importancelevels". 
 example:
 <xsl:text>essential</xsl:text>
 or
 <xsl:text>essential,expert</xsl:text>
 
 -->
  <xsl:text>essential</xsl:text>
 </xsl:variable>
 <xsl:template match="/">
 <xsl:call-template name="elementToLatex">
   <xsl:with-param name="myelement" select="//xs:element[@name=/xs:schema/xs:annotation/xs:appinfo/root]" />
   <xsl:with-param name="level" select="0" />
  </xsl:call-template>
  <xsl:for-each select="/*/xs:element[@name!=/xs:schema/xs:annotation/xs:appinfo/root and contains($importancelevels,@ex:importance)]">
   <xsl:call-template name="elementToLatex">
    <xsl:with-param name="myelement" select="." />
    <xsl:with-param name="level" select="0" />
   </xsl:call-template>
  </xsl:for-each>
  
   </xsl:template>
 <xsl:template match="displaymath">
    <xsl:text> 
[[math label]] 
</xsl:text>
    <xsl:value-of select="." />
    <xsl:text>
[[/math]]
</xsl:text>
  </xsl:template>
  <xsl:template match="inlinemath">
    <xsl:text> [[$ </xsl:text>
    <xsl:value-of select="normalize-space(.)" />
    <xsl:text> $]]</xsl:text>
  </xsl:template>
  <xsl:template match="pre">
    <xsl:text> {{ </xsl:text>
    <xsl:value-of select="normalize-space(.)" />
    <xsl:text> }}</xsl:text>
  </xsl:template>
  <xsl:template match="it">
    <xsl:text> // </xsl:text>
    <xsl:value-of select="normalize-space(.)" />
    <xsl:text> //</xsl:text>
  </xsl:template>
  <xsl:template match="bf">
    <xsl:text> **</xsl:text>
    <xsl:value-of select="normalize-space(.)" />
    <xsl:text>** </xsl:text>
  </xsl:template>
  <xsl:template match="text()">
    <xsl:value-of select="normalize-space(.)" />
  </xsl:template>
  <xsl:template match="xs:documentation">
    <xsl:apply-templates select="text()|inlinemath|displaymath|pre|it" />
  </xsl:template>
 <xsl:template name="elementToLatex">
  <xsl:param name="myelement" />
  <xsl:param name="level" />
<xsl:text>
[[# </xsl:text>
  <xsl:value-of select="$myelement/@name" />
  <xsl:text>]]
</xsl:text>
  <xsl:text>+ Element: </xsl:text>
  <xsl:value-of select="$myelement/@name " />
  <xsl:text>
  
  </xsl:text>
  

  <xsl:apply-templates select="$myelement/xs:annotation/xs:documentation" />
  <xsl:call-template name="TypeToDoc">
   <xsl:with-param name="contentnode" select="$myelement | //xs:element[@name=$myelement/@ref]" />
  </xsl:call-template>
  <xsl:for-each select="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]">
   <xsl:sort select="@name|@ref" />
   <xsl:call-template name="attributetolatex">
    <xsl:with-param name="myattribute" select="." />
    <xsl:with-param name="level" select="$level" />
   </xsl:call-template>
  </xsl:for-each>
  <xsl:for-each select="$myelement/*/*/xs:element[contains($importancelevels,@ex:importance)and @name]">
   <xsl:call-template name="elementToLatex">
    <xsl:with-param name="myelement" select="." />
    <xsl:with-param name="level" select="$level+1" />
   </xsl:call-template>
  </xsl:for-each>
 </xsl:template>
 <xsl:template name="attributetolatex">
  <xsl:param name="myattribute" />
  <xsl:param name="level" />
  <xsl:text>
++ Attribute: </xsl:text>
  <xsl:value-of select="$myattribute/@name |$myattribute/@ref" />
  <xsl:text>  
    </xsl:text>
  <xsl:apply-templates select="$myattribute/xs:annotation/xs:documentation" />
  <xsl:call-template name="TypeToDoc">
   <xsl:with-param name="contentnode" select="$myattribute | //xs:attribute[@name=$myattribute/@ref]" />
  </xsl:call-template>
 </xsl:template>
 <xsl:template name="TypeToDoc">
  <xsl:param name="contentnode" />
  <xsl:text>

[[table ]]
[[row]]
[[cell style=" vertical-align:top;" ]] **Type:** [[/cell]] [[cell]]</xsl:text>
 <xsl:choose>
 <xsl:when test="$contentnode/@type">
 <xsl:value-of select="$contentnode/@type"/>
 <xsl:text>
</xsl:text>
 </xsl:when>
 <xsl:when test="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
<xsl:text> **choose from:**  
</xsl:text>
 <xsl:for-each select="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
<xsl:text> </xsl:text>
<xsl:value-of select="@value"/><xsl:text/> 
<xsl:text>  
</xsl:text>
 </xsl:for-each>
 <xsl:text/> 
 </xsl:when>
 <xsl:when test="$contentnode/xs:complexType/*[xs:element] ">
 <xsl:text> **contains:** </xsl:text> 
 <xsl:text>  
</xsl:text>
 <xsl:for-each select="$contentnode/xs:complexType/*/xs:element[contains($importancelevels,@ex:importance)]">
 <xsl:text>  [#</xsl:text>
<xsl:value-of select="./@name|@ref"/><xsl:text>   </xsl:text>
<xsl:value-of select="./@name|@ref"/><xsl:text>]</xsl:text>
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

 <xsl:text>  
</xsl:text>
 </xsl:for-each>
 </xsl:when>
 <xsl:otherwise>
 <xsl:text> no content 
</xsl:text>
 </xsl:otherwise>
 </xsl:choose>
 <xsl:text> [[/cell]][[/row]]</xsl:text>
 <xsl:choose>
 

 <xsl:when test="$contentnode/@ex:unit!=''">
 <xsl:text>
[[row]] [[cell]] **Unit:** [[/cell]][[cell]]</xsl:text>
 <xsl:value-of select="$contentnode/@ex:unit"/>
  <xsl:text>  [[/cell]] [[/row]]
  </xsl:text>
 </xsl:when>
 </xsl:choose>

 <xsl:text>[[row]] [[cell]] **XPath:** [[/cell]][[cell]] {{</xsl:text>
   <xsl:call-template name="genxpath" >
  <xsl:with-param name="node" select="$contentnode"/>
  <xsl:with-param name="xpath" select="''"/>
  </xsl:call-template>
  <xsl:text> }}[[/cell]] [[/row]]
  

[[/table]]
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
     <xsl:text>.</xsl:text>
      <xsl:value-of select="$xpath"/>
    </xsl:otherwise>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>
</xsl:stylesheet>