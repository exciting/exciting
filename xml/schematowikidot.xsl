<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
  xmlns:xs="http://www.w3.org/2001/XMLSchema" xmlns:str="http://exslt.org/strings"
  xmlns:regex="http://www.exslt.org/regexp"
  xmlns:ex="http://xml.exciting-code.org/inputschemaextentions.xsd">
  <xsl:output method="text"/>
  <xsl:param name="index" select="'false'"/>
  <xsl:param name="importancelevels">
    <xsl:text>essential</xsl:text>
    <xs:annotation>
      <xs:documentation> In order to select the importance levels that should be included list them
        in the parameter "importancelevels". example: xsltproc --stringparam importancelevels
        "essential expert" schematowikidot.xsl excitinginput.xsd >iref.txt --stringparam tabs "true"
        :activates the use of tabs in the element descriptions. </xs:documentation>
    </xs:annotation>
  </xsl:param>
  <xsl:param name="tabs" select="false"/>
  <xsl:template match="/">
    <xsl:apply-templates select="/xs:schema/xs:annotation/xs:documentation"/>
    <xsl:if test="$index='true'">
      <xsl:text>
   [[collapsible show="+ Show alphabetical index" hide="- Hide alphabetical index"]]
   The @ sign indicates an attribute.

</xsl:text>
      <xsl:for-each
        select="//xs:attribute[@name and contains($importancelevels,@ex:importance)]
        |//xs:element[@name and contains($importancelevels,@ex:importance)]">
        <xsl:sort select="@name"/>
        <xsl:variable name="plevel">
          <xsl:call-template name="isincluded">
            <xsl:with-param name="node" select="."/>
          </xsl:call-template>
        </xsl:variable>
        <xsl:if test="$plevel='include'">
          <xsl:text>|| [#</xsl:text>
          <xsl:if test="name(.)='xs:attribute'">
            <xsl:text>att</xsl:text>
          </xsl:if>
          <xsl:value-of select="@name"/>
          <xsl:text> </xsl:text>
          <xsl:value-of select="@name"/>
          <xsl:text>]</xsl:text>
          <xsl:choose>
            <xsl:when test="name(.)='xs:attribute'">
              <xsl:text>||@</xsl:text>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>|| </xsl:text>
            </xsl:otherwise>
          </xsl:choose>
          <xsl:text>||</xsl:text>
          <xsl:call-template name="genxpath">
            <xsl:with-param name="node" select="."/>
          </xsl:call-template>
          <xsl:text>||
</xsl:text>
        </xsl:if>
      </xsl:for-each>
      <xsl:text>
[[/collapsible]]
    </xsl:text>
    </xsl:if>
    <xsl:choose>
      <xsl:when test="/xs:schema/xs:annotation[last()]/xs:appinfo/root='/'">
        <xsl:for-each select="/*/xs:element">
        <xsl:call-template name="elementToLatex">
          <xsl:with-param name="myelement"
            select="."/>
          <xsl:with-param name="level" select="0"/>
        </xsl:call-template>
        </xsl:for-each>
      </xsl:when>
      <xsl:otherwise>
    
    <xsl:call-template name="elementToLatex">
      <xsl:with-param name="myelement"
        select="//xs:element[@name=/xs:schema/xs:annotation[last()]/xs:appinfo/root]"/>
      <xsl:with-param name="level" select="0"/>
    </xsl:call-template>
      </xsl:otherwise>
    </xsl:choose>
    <xsl:text>
+ Reused Elements
    
The following elements can occur more than once in the input file. There for they are listed separately.
  </xsl:text>
    <xsl:for-each
      select="/*/xs:element[@name!=/xs:schema/xs:annotation[last()]/xs:appinfo/root
          and contains($importancelevels,@ex:importance)  ]">
      <xsl:variable name="name" select="@name"/>
      <xsl:if test="count(//xs:element[@ref=$name])>1">
        <xsl:call-template name="elementToLatex">
          <xsl:with-param name="myelement" select="."/>
          <xsl:with-param name="level" select="0"/>
        </xsl:call-template>
      </xsl:if>
    </xsl:for-each>
    <xsl:text>
+ Data Types
 
 The Input definition uses derived data types. These are described here.
  </xsl:text>
    <xsl:for-each select="/*/xs:simpleType">
      <xsl:call-template name="typetoDoc">
        <xsl:with-param name="typenode" select="."/>
      </xsl:call-template>
    </xsl:for-each>
  </xsl:template>
  <xsl:template match="displaymath">
    <xsl:text> 
[[math label]] 
</xsl:text>
    <xsl:value-of select="."/>
    <xsl:text>
[[/math]]
</xsl:text>
  </xsl:template>
  <xsl:template match="inlinemath">
    <xsl:text> [[$ </xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text> $]] </xsl:text>
  </xsl:template>
  <xsl:template match="inlinemath_ns">
    <xsl:text> [[$ </xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text> $]]</xsl:text>
    <xsl:if test="@space">
      <xsl:text> </xsl:text>
    </xsl:if>
  </xsl:template>
  <xsl:template match="pre">
    <xsl:text> {{</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>}} </xsl:text>
  </xsl:template>
  <xsl:template match="ns_pre">
    <xsl:text>{{</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>}} </xsl:text>
  </xsl:template>
  <xsl:template match="pre_ns">
    <xsl:text> {{</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>}}</xsl:text>
  </xsl:template>
  <xsl:template match="ns_pre_ns">
    <xsl:text>{{</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>}}</xsl:text>
  </xsl:template>
  <xsl:template match="it">
    <xsl:text> //</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>// </xsl:text>
  </xsl:template>
  <xsl:template match="it_ns">
    <xsl:text> //</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>//</xsl:text>
  </xsl:template>
  <xsl:template match="bf">
    <xsl:text> **</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>**</xsl:text>
    <xsl:if test="not(@nospace='true')">
      <xsl:text> </xsl:text>
    </xsl:if>
  </xsl:template>
  <xsl:template match="pre-bf">
    <xsl:text> {{**</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>**}} </xsl:text>
  </xsl:template>
  <xsl:template match="pre-bf_ns">
    <xsl:text> {{**</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>**}}</xsl:text>
  </xsl:template>
  <xsl:template match="filename">
    <xsl:text> {{**//</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>//**}} </xsl:text>
  </xsl:template>
  <xsl:template match="filename_ns">
    <xsl:text> {{**//</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>//**}}</xsl:text>
  </xsl:template>
  <xsl:template match="green">
    <xsl:text> ##green|</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>## </xsl:text>
  </xsl:template>
  <xsl:template match="blue">
    <xsl:text> ##blue|</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>## </xsl:text>
  </xsl:template>
  <xsl:template match="text()">
    <xsl:value-of select="normalize-space(.)"/>
  </xsl:template>
  <xsl:template match="xs:documentation">
    <xsl:apply-templates
      select="text()|inlinemath|inlinemath_ns|displaymath|pre|pre_ns|ns_pre|ns_pre_ns|pre-bf|pre-bf_ns|it|it_ns|p|exciting|a|list|li|attref|bf"
    />
  </xsl:template>
  <xsl:template match="attref">
    <xsl:text> [#att</xsl:text>
    <xsl:value-of select="."/>
    <xsl:text> </xsl:text>
    <xsl:value-of select="."/>
    <xsl:text>]</xsl:text>
  </xsl:template>
  <xsl:template match="list">
    <xsl:text>
</xsl:text>
    <xsl:apply-templates select="./*|text()"/>
  </xsl:template>
  <xsl:template match="li">
    <xsl:text>* </xsl:text>
    <xsl:apply-templates select="./*|text()"/>
    <xsl:text>
</xsl:text>
  </xsl:template>
  <xsl:template match="p">
    <xsl:text>
  
  </xsl:text>
    <xsl:apply-templates select="./*|text()"/>
  </xsl:template>
  <xsl:template match="a">
    <xsl:text> [</xsl:text>
    <xsl:value-of select="@href"/>
    <xsl:text>  </xsl:text>
    <xsl:value-of select="."/>
    <xsl:text>]</xsl:text>
    <xsl:if test="@space">
      <xsl:text> </xsl:text>
    </xsl:if>
  </xsl:template>
  <xsl:template match="exciting">
    <xsl:text> {{**exciting**}} </xsl:text>
    <xsl:if test="@space">
      <xsl:text> </xsl:text>
    </xsl:if>
  </xsl:template>
  <xsl:template name="elementToLatex">
    <xsl:param name="myelement"/>
    <xsl:param name="level"/>
    <xsl:text>
[[# </xsl:text>
    <xsl:value-of select="$myelement/@name"/>
    <xsl:text>]]
</xsl:text>
    <xsl:text>+ ##blue|Element:## </xsl:text>
    <xsl:value-of select="$myelement/@name "/>
    <xsl:text>
  
  </xsl:text>
    <xsl:if test="$tabs">
      <xsl:text>
[[tabview]]
[[tab Description]]
</xsl:text>
    </xsl:if>
    <xsl:apply-templates select="$myelement/xs:annotation/xs:documentation"/>
    <xsl:if test="$tabs"> </xsl:if>
    <xsl:call-template name="TypeToDoc">
      <xsl:with-param name="contentnode" select="$myelement | //xs:element[@name=$myelement/@ref]"/>
    </xsl:call-template>
    <xsl:if test="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]">

      <xsl:choose>
        <xsl:when test="$tabs">
          <xsl:text>
      **List of attributes:** </xsl:text>
          <xsl:for-each
            select="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]">

            <xsl:text>[[# att</xsl:text>
            <xsl:value-of select="@name|@ref"/>
            <xsl:text>]] ##green|</xsl:text>
            <xsl:value-of select="@name|@ref"/>
            <xsl:text>## </xsl:text>
            <xsl:if test="not(position()=last())">
              <xsl:text>, </xsl:text>
            </xsl:if>
          </xsl:for-each>
          <xsl:text>
      [[/tab]]
   [[tab Attributes]]
   </xsl:text>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>
This element allows for specification of the following attributes:  </xsl:text>
        </xsl:otherwise>
      </xsl:choose>

    </xsl:if>
    <xsl:for-each select="$myelement/*/xs:attribute[contains($importancelevels,@ex:importance)]">
      <xsl:sort select="@name|@ref"/>
      <xsl:call-template name="attributetolatex">
        <xsl:with-param name="myattribute" select="."/>
        <xsl:with-param name="level" select="$level"/>
      </xsl:call-template>
    </xsl:for-each>
    <xsl:if test="$tabs">
      <xsl:text> 
   [[/tab]]
    [[/tabview]]
    </xsl:text>
    </xsl:if>
    <xsl:for-each select="$myelement/*/*/xs:element[contains($importancelevels,@ex:importance)]">
      <xsl:if test="@name">
        <xsl:call-template name="elementToLatex">
          <xsl:with-param name="myelement" select="."/>
          <xsl:with-param name="level" select="$level+1"/>
        </xsl:call-template>
      </xsl:if>
      <xsl:if test="@ref">
        <xsl:variable name="ref" select="@ref"/>     
        <xsl:if test="count(//xs:element[@ref=$ref])=1">
          <xsl:call-template name="elementToLatex">
            <xsl:with-param name="myelement" select="//xs:element[@name=$ref]"/>
            <xsl:with-param name="level" select="$level+1"/>
          </xsl:call-template>
        </xsl:if>
      </xsl:if>
    </xsl:for-each>
  </xsl:template>
  <xsl:template name="attributetolatex">
    <xsl:param name="myattribute"/>
    <xsl:param name="level"/>
    <xsl:text>
  [[# att</xsl:text>
    <xsl:value-of select="$myattribute/@name |$myattribute/@ref"/>
    <xsl:text>]]
  
++ ##green|Attribute:## </xsl:text>
    <xsl:value-of select="$myattribute/@name |$myattribute/@ref"/>
    <xsl:text>  
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

[[table ]]
[[row]]
</xsl:text>
    <xsl:choose>
      <xsl:when test="$contentnode/@type">
        <xsl:text>[[cell style=" vertical-align:top;" ]] **Type:** [[/cell]] [[cell]]</xsl:text>
        <xsl:if test="not(contains($contentnode/@type,'xs:'))">
          <xsl:text>[#</xsl:text>
          <xsl:value-of select="$contentnode/@type"/>
          <xsl:text> </xsl:text>
        </xsl:if>
        <xsl:value-of select="str:replace(($contentnode/@type),'xs:','')"/>
        <xsl:if test="not(contains($contentnode/@type,'xs:'))">
          <xsl:text>]
 </xsl:text>
        </xsl:if>
        <xsl:text>
</xsl:text>
      </xsl:when>
      <xsl:when test="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
        <xsl:text> [[cell style=" vertical-align:top;" ]] **Type:** [[/cell]] [[cell]] **choose from:**  
</xsl:text>
        <xsl:for-each
          select="$contentnode/xs:simpleType/xs:restriction[@base='xs:string']/xs:enumeration">
          <xsl:text> </xsl:text>
          <xsl:value-of select="@value"/>
          <xsl:text/>
          <xsl:text>  
</xsl:text>
        </xsl:for-each>
        <xsl:text/>
      </xsl:when>
      <xsl:when test="$contentnode/xs:complexType/*[xs:element] ">
        <xsl:text> [[cell style=" vertical-align:top;" ]] **contains:** [[/cell]] [[cell]]</xsl:text>
        <xsl:text>  
</xsl:text>
        <xsl:for-each
          select="$contentnode/xs:complexType/*/xs:element[contains($importancelevels,@ex:importance)]">
          <xsl:text>  [#</xsl:text>
          <xsl:value-of select="./@name|@ref"/>
          <xsl:text>   </xsl:text>
          <xsl:value-of select="./@name|@ref"/>
          <xsl:text>]</xsl:text>
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
        <xsl:text>[[cell style=" vertical-align:top;" ]] **Type:** [[/cell]] [[cell]] no  content  
</xsl:text>
      </xsl:otherwise>
    </xsl:choose>
    <xsl:text> [[/cell]][[/row]]</xsl:text>
    <xsl:if test="$contentnode/@default">
      <xsl:text>
[[row]] [[cell]] **Default:** [[/cell]][[cell]] "{{</xsl:text>
      <xsl:value-of select="$contentnode/@default"/>
      <xsl:text>}}" [[/cell]][[/row]]
 </xsl:text>
    </xsl:if>
    <xsl:if test="$contentnode/@use or local-name($contentnode)='attribute'">
      <xsl:text>
[[row]] [[cell]] **Use:** [[/cell]][[cell]]  </xsl:text>
      <xsl:value-of select="$contentnode/@use"/>
      <xsl:if test="not($contentnode/@use) and local-name($contentnode)='attribute'">
        <xsl:text>optional</xsl:text>
      </xsl:if>
      <xsl:text> [[/cell]][[/row]]
 </xsl:text>
    </xsl:if>
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
    <xsl:call-template name="genxpath">
      <xsl:with-param name="node" select="$contentnode"/>
      
      <xsl:with-param name="xpath" select="''"/>
    </xsl:call-template>
    <xsl:text> }}[[/cell]] [[/row]]
  </xsl:text>
    <xsl:if test="$contentnode/@name and name(..)='xs:schema'">
      <xsl:variable name="name" select="$contentnode/@name"/>
      <xsl:text>[[row]] [[cell style=" vertical-align:top;" ]] **Parent:** [[/cell]][[cell]] {{</xsl:text>
      <xsl:for-each
        select="//xs:element[@ref=$name  and contains($importancelevels,@ex:importance)]">
        <xsl:call-template name="genxpath">
          <xsl:with-param name="node" select="."/>
          <xsl:with-param name="xpath" select="''"/>
        </xsl:call-template>
        <xsl:text>
</xsl:text>
      </xsl:for-each>
      <xsl:text> }}[[/cell]] [[/row]]</xsl:text>
    </xsl:if>
    <xsl:text>
[[/table]]
  </xsl:text>
  </xsl:template>
  <xsl:template name="genxpath">
    <xsl:param name="xpath"/>
    <xsl:param name="node"/>
    <xsl:variable name="name" select="$node/@name"/>
    <xsl:variable name="current_name">
      <xsl:if test="$node/@name">
        <xsl:text>/</xsl:text>
        <xsl:choose>
          <xsl:when test="name($node)='xs:attribute'">
            <xsl:text>[#att</xsl:text>
            <xsl:value-of select="$node/@name|$node/@ref"/>
            <xsl:text> @</xsl:text>
            <xsl:value-of select="$node/@name|$node/@ref"/>
            <xsl:text>]</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>[#</xsl:text>
            <xsl:value-of select="$node/@name|$node/@ref"/>
            <xsl:text> </xsl:text>
            <xsl:value-of select="$node/@name|$node/@ref"/>
            <xsl:text>]</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:if>
      <xsl:value-of select="$xpath"/>
    </xsl:variable>
    <xsl:for-each select="$node[last()]">
      <xsl:choose>
        <xsl:when test="name(..)='xs:schema' and count(//xs:element[@ref=$name])=1">
          <xsl:call-template name="genxpath">
            <xsl:with-param name="node" select="//xs:element[*/*/xs:element[@ref=$name ]]"/>
            <xsl:with-param name="xpath">
              <xsl:value-of select="$current_name"/>
            </xsl:with-param>
          </xsl:call-template>
        </xsl:when>
        <xsl:when test="parent::node()">
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
      <xsl:if test="not(last())">
        <xsl:text>, </xsl:text>
      </xsl:if>
    </xsl:for-each>
  </xsl:template>
  <xsl:template name="typetoDoc">
    <xsl:param name="typenode"/>
    <xsl:text>
  [[# </xsl:text>
    <xsl:value-of select="$typenode/@name"/>
    <xsl:text>]]
++ Type   </xsl:text>
    <xsl:value-of select="$typenode/@name"/>
    <xsl:text>
  
  </xsl:text>
    <xsl:apply-templates select="xs:annotation/xs:documentation"/>
  </xsl:template>
  <xsl:template name="isincluded">
    <xsl:param name="node"/>
    <xsl:for-each select="$node">
      <xsl:choose>
        <xsl:when test="contains($importancelevels,@ex:importance) or not(@ex:importance)">
          <xsl:choose>
            <xsl:when test="parent::node()">
              <xsl:call-template name="isincluded">
                <xsl:with-param name="node" select="parent::node()"/>
              </xsl:call-template>
            </xsl:when>
            <xsl:otherwise>
              <xsl:text>include</xsl:text>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:when>
        <xsl:otherwise>
          <xsl:text>not</xsl:text>
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
  </xsl:template>
</xsl:stylesheet>
