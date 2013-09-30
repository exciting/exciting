<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" 
xmlns:xs="http://www.w3.org/2001/XMLSchema"
 xmlns:str="http://exslt.org/strings"
 xmlns:ex="http://xml.exciting-code.org/inputschemaextentions.xsd">
  <xsl:output method="text"/>
  <xsl:variable name="newline">
<xsl:text>
</xsl:text>
  </xsl:variable> 
  <xsl:param name="tool" select="exciting"/>
  <xsl:variable name="root" select="/xs:schema/xs:annotation[last()]/xs:appinfo/root"></xsl:variable>
  <xsl:template name="xstypetofortrantype">
    <xsl:param name="type"/>
    <xsl:param name="varname"/>
    <xsl:param name="maxoccurs" select="1"/>
<xsl:param name="enumeration" select="''"/>
    <xsl:variable name="allocatable">
      <xsl:if test="$maxoccurs='unbounded' or $maxoccurs&gt;1">
        <xsl:text>,pointer</xsl:text>
      </xsl:if>
    </xsl:variable>

    <xsl:choose>
      <xsl:when test="$type='xs:string'">
        <xsl:text> character(512)::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:text>
</xsl:text>
<xsl:if test="$enumeration=true()">
 <xsl:text> integer::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:text>number
</xsl:text>
</xsl:if>
      </xsl:when>
      <xsl:when test="$type='xs:boolean'">
        <xsl:text> logical::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:text>
</xsl:text>

      </xsl:when>
      <xsl:when test="$type='integerpair'">
        <xsl:text> integer</xsl:text>
        <xsl:value-of select="$allocatable"/>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:choose>
          <xsl:when test="$allocatable!=''">
            <xsl:text>(:,:)
</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>(2)
</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>

      <xsl:when test="$type='booleantriple'">
        <xsl:text> logical</xsl:text>
        <xsl:value-of select="$allocatable"/>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:choose>
          <xsl:when test="$allocatable!=''">
            <xsl:text>(:,:)
</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>(3)
</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>

      <xsl:when test="$type='integertriple'">
        <xsl:text> integer</xsl:text>
        <xsl:value-of select="$allocatable"/>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:choose>
          <xsl:when test="$allocatable!=''">
            <xsl:text>(:,:)
</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>(3)
</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:when test="$type='integerquadrupel'">
        <xsl:text> integer</xsl:text>
        <xsl:value-of select="$allocatable"/>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:choose>
          <xsl:when test="$allocatable!=''">
            <xsl:text>(:,:)
</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>(4)
</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:when test="$type='vect2d'">
        <xsl:text> real(8)</xsl:text>
        <xsl:value-of select="$allocatable"/>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:choose>
          <xsl:when test="$allocatable!=''">
            <xsl:text>(:,:)
</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>(2)
</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:when test="$type='vect3d'">
        <xsl:text> real(8)</xsl:text>
        <xsl:value-of select="$allocatable"/>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:choose>
          <xsl:when test="$allocatable!=''">
            <xsl:text>(:,:)
</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:text>(3)
</xsl:text>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
      <xsl:when test="$type='fortrandouble'">
        <xsl:text> real(8)</xsl:text>

        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:text>
</xsl:text>

      </xsl:when>
      <xsl:when test="$type='xs:integer'">
        <xsl:text> integer</xsl:text>
        <xsl:text>::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:text>
</xsl:text>

      </xsl:when>
      <xsl:when test="$type='xs:anyURI' or $type='xs:ID'  or $type='xs:IDREFS' ">
        <xsl:text> character(1024)::</xsl:text>
        <xsl:value-of select="$varname"/>
        <xsl:text>
</xsl:text>

      </xsl:when>
      <xsl:when test="$type=''"></xsl:when>
      <xsl:otherwise>
        <xsl:text>! missing something
</xsl:text>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>
  <xsl:template name="contenttotype">
    <xsl:choose>
      <xsl:when test="./*/*/xs:extension or ./*/*/xs:element  or ./*/xs:attribute or xs:attribute or ./*/*/xs:attribute">
        <xsl:text>
type </xsl:text>
        <xsl:value-of select="@name"/>
        <xsl:text>_type
! Do not edit! This file is automatically generated from the XML Schema
</xsl:text>
        <xsl:for-each select=" ./*/xs:attribute|xs:attribute|./*/*/xs:attribute">
          <xsl:variable name="type">
            <xsl:choose>
              <xsl:when test="@ref">
                <xsl:variable name="ref" select="./@ref"/>
                <xsl:value-of select="//xs:attribute[@name=$ref]/@type|//xs:attribute[@name=$ref]/*/xs:restriction/@base"/>
              </xsl:when>
              <xsl:when test="./*/*/xs:extension/@base">
                <xsl:variable name="name" select="./@name"/>
                <xsl:value-of select="./*/*/xs:extension/@base"/>
              </xsl:when>
              <xsl:when test="@type">
                <xsl:value-of select="./@type"/>
              </xsl:when>
              <xsl:when test="./xs:simpleType/xs:restriction/@base='xs:string'">
                <xsl:value-of select="'xs:string'"/>
              </xsl:when>
              <xsl:otherwise/>


            </xsl:choose>
          </xsl:variable>

          <xsl:call-template name="xstypetofortrantype"  >
            <xsl:with-param name="type" select="$type"/>
            <xsl:with-param name="varname" select="@name|@ref"/>
            <xsl:with-param name="enumeration" select="./xs:simpleType/xs:restriction/xs:enumeration"/>
         
          </xsl:call-template>
        </xsl:for-each>
        <xsl:for-each select="./*/*/xs:element">
      <xsl:choose> <!--    
            <xsl:when test="(not(./@type) and not(./*/*/xs:element)) and (not(./*/xs:attribute) and not(./@ref))">
              <xsl:text>  logical::</xsl:text>
              <xsl:value-of select="@name|@ref"/>
              <xsl:text>
</xsl:text>
            </xsl:when> --> 
            <xsl:when test="@maxOccurs&gt;1 or @maxOccurs='unbounded'">
            <xsl:variable name="ref" select="@ref"/>
              <xsl:choose>
                <xsl:when test="(not(./*/*/xs:element) and (not(./*/xs:attribute) or ( //xs:element[@name=$ref]/@type))
                and not(//xs:element[@name=$ref]/xs:complexType)) ">

                 

                  <xsl:call-template name="xstypetofortrantype">
                    <xsl:with-param name="type" select="//xs:element[@name=$ref]/@type|@type"/>
                    <xsl:with-param name="varname" select="@name|@ref"/>
                    <xsl:with-param name="maxoccurs" select="@maxOccurs"/>
                  </xsl:call-template>
                </xsl:when>
                <xsl:otherwise>
                  <!-- handle pointerarrays -->
                  <xsl:text>  type(</xsl:text>
                  <xsl:value-of select="@name|@ref"/>
                  <xsl:text>_type_array),pointer::</xsl:text>
                  <xsl:value-of select="@name|@ref"/>
                  <xsl:text>array(:)
</xsl:text>
                </xsl:otherwise>
              </xsl:choose>
            </xsl:when>
            <xsl:otherwise>
              <xsl:variable name="ref">
                <xsl:value-of select="@ref"/>
              </xsl:variable>
              <xsl:choose>
                <xsl:when test="(not(./*/*/xs:element) 
                and not(./*/xs:attribute)) 
                and ( @type or //xs:element[@name=$ref]/@type)">



                  <xsl:call-template name="xstypetofortrantype">
                    <xsl:with-param name="type" select="//xs:element[@name=$ref]/@type|@type"/>
                    <xsl:with-param name="varname" select="@name|@ref"/>
                  </xsl:call-template>
                </xsl:when>
                <xsl:otherwise>
                  <xsl:text>  type(</xsl:text>
                  <xsl:value-of select="@name|@ref"/>
                  <xsl:text>_type),pointer::</xsl:text>
                  <xsl:value-of select="@name|@ref"/>
                  <xsl:text>
</xsl:text>
                </xsl:otherwise>
              </xsl:choose>
            </xsl:otherwise>
          </xsl:choose>
        </xsl:for-each>
        <xsl:text>end type
</xsl:text>
      </xsl:when>
    </xsl:choose>
    
    <xsl:if test="(@maxOccurs='unbounded' and not(@type))">
      <xsl:text>
type </xsl:text>
      <xsl:value-of select="@name|@ref"/>
      <xsl:text>_type_array
type(</xsl:text>
      <xsl:value-of select="@name|@ref"/>
      <xsl:text>_type),pointer::</xsl:text>
      <xsl:value-of select="@name|@ref"/>
      <xsl:text>
 end type
    </xsl:text>
    </xsl:if>
     <xsl:if test="(not(xs:complexType/*) and not(@type)) and (not(@ref) )">
      <xsl:text>
type </xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>_type
logical::exists
 end type
    </xsl:text>
    </xsl:if>
  </xsl:template>
  <xsl:template match="/">
    <xsl:text>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!! DO NOT EDIT                          !!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! This file is automatically generated from the XML Schema
! all changes will be overwritten when anything changes in the Schema



module mod</xsl:text><xsl:value-of select="$root"/> <xsl:text>
use inputdom
implicit none
</xsl:text>


    <xsl:for-each select="//xs:element[not(@ref)]">
      <xsl:call-template name="contenttotype"/>
      <xsl:variable name="name" select="@name"/>
      <xsl:if test="not(@type)">
      <xsl:for-each select="//xs:element[(@ref=$name and not(@type))and (@maxOccurs>1 or @maxOccurs='unbounded')]">
      <xsl:if test="position()=1">
        <xsl:call-template name="contenttotype"/>
      </xsl:if>
      </xsl:for-each>
       </xsl:if>
    </xsl:for-each>
   

  
    
    <xsl:text>
   type(</xsl:text><xsl:value-of select="$root"></xsl:value-of> 
   <xsl:text>_type)::</xsl:text>
   <xsl:value-of select="$root"></xsl:value-of> <xsl:text>
contains
</xsl:text>
    <!-- generate functions -->
    <xsl:for-each select="//xs:element[(@name and not(@type)) ]">
      <xsl:variable name="struct">
        <xsl:text>getstruct</xsl:text>
        <xsl:value-of select="@name"/>
      </xsl:variable>
      <xsl:text>
function getstruct</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>(thisnode)
! Do not edit! This file is automatically generated from the XML Schema
implicit none
type(Node),pointer::thisnode
type(</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>_type),pointer::getstruct</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:value-of select="$newline"/>
      <xsl:if test="./*/xs:attribute">
        <xsl:text>type(Node),pointer::np

</xsl:text>
      </xsl:if>
      <xsl:text>
integer::len=1,i=0
allocate(getstruct</xsl:text>
      <xsl:value-of select="@name"/>
      <xsl:text>)  
#ifdef INPUTDEBUG      
      write(*,*)"we are at </xsl:text><xsl:value-of select="@name"/><xsl:text>"
#endif
      </xsl:text>
<xsl:if test="not (@type) and not (xs:complexType)">
<xsl:text>getstruct</xsl:text>
      <xsl:value-of select="@name"/><xsl:text>%exists=.false.
      </xsl:text>
<xsl:text>if (associated(thisnode))  getstruct</xsl:text>
      <xsl:value-of select="@name"/><xsl:text>%exists=.true.
      </xsl:text></xsl:if>
      <xsl:for-each select="./*/xs:attribute">
        <xsl:text>
nullify(np)  
np=>getAttributeNode(thisnode,"</xsl:text>
        <xsl:value-of select="@name|@ref"/>
        <xsl:text>")
</xsl:text>
<xsl:if test="@default">
  <xsl:value-of select="$struct"/>
        <xsl:text>%</xsl:text>
        <xsl:value-of select="@name|@ref"/>  <xsl:text>=</xsl:text>
<xsl:call-template name="defaultvaltofortran">
<xsl:with-param name="default" select="@default"/>
<xsl:with-param name="type" >
<xsl:variable name="ref" select="@ref"/>
<xsl:value-of select="@type"/>
<xsl:value-of select="//xs:attribute[@name=$ref]/@type"/>
<xsl:value-of select="//xs:attribute[@name=$ref]/*/xs:restriction/@base"/>
<xsl:value-of select="./*/xs:restriction/@base"/>
</xsl:with-param>
</xsl:call-template>
 <xsl:text>
</xsl:text>
</xsl:if>
<xsl:text>if(associated(np)) then
       call extractDataAttribute(thisnode,"</xsl:text>
        <xsl:value-of select="@name|@ref"/>
        <xsl:text>",</xsl:text>
        <xsl:value-of select="$struct"/>
        <xsl:text>%</xsl:text>
        <xsl:value-of select="@name|@ref"/>
        <xsl:text>)
       call removeAttribute(thisnode,"</xsl:text>
        <xsl:value-of select="@name|@ref"/>
        <xsl:text>")  </xsl:text>
        <xsl:if test="@use='required'">
        <xsl:text>
        else
        write(*,*)"Parser ERROR: The element '</xsl:text>
        <xsl:value-of select="../../@name"/>
        <xsl:text>' requires the attribute '</xsl:text>
        <xsl:value-of select="@name|@ref"/> 
        <xsl:text>' to be defined."</xsl:text>
        write(*,*)"stopped"
        stop
        </xsl:if>    
<xsl:text>
endif
</xsl:text>
<xsl:if test="./*/xs:restriction/xs:enumeration"> 
<xsl:text>getstruct</xsl:text><xsl:value-of select="../../@name"/>
<xsl:text>%</xsl:text><xsl:value-of select="@name|@ref"/><xsl:text disable-output-escaping="yes">number= &amp;
  stringtonumber</xsl:text>
<xsl:value-of select="../../@name"/><xsl:value-of select="@name|@ref"/>
<xsl:text>(getstruct</xsl:text><xsl:value-of select="../../@name"/>
<xsl:text>%</xsl:text>
<xsl:value-of select="@name|@ref"/>)
</xsl:if>
      </xsl:for-each>

      <xsl:for-each select="*/*/xs:element[not(@type) ]">
        <xsl:variable name="ref" select="@ref"/>
        <xsl:choose>
          <xsl:when test="//xs:element[@name=$ref]/@type">

          </xsl:when>
          <xsl:otherwise>
            <xsl:text>
            len= countChildEmentsWithName(thisnode,"</xsl:text>
            <xsl:value-of select="@name|@ref"/>
            <xsl:text>")
</xsl:text>
   <xsl:if test="@minOccurs>=1">
        <xsl:text>
        if(len.eq.0) then
        write(*,*)"Parser ERROR: The </xsl:text><xsl:value-of select="../../../@name"/>
        <xsl:text> element must contain at least </xsl:text>
        <xsl:value-of select="@minOccurs"/><xsl:text> </xsl:text>
      
        <xsl:if test="@maxOccurs>1">
        <xsl:text> and maximum </xsl:text>
        <xsl:value-of select="@maxOccurs"/>
         <xsl:text> </xsl:text>
        </xsl:if>
        <xsl:value-of select="@name|@ref"/>
       
        <xsl:text> element</xsl:text>
                <xsl:if test="@maxOccurs>1"><xsl:text>s</xsl:text></xsl:if>
        
        <xsl:text>"
        endif
        </xsl:text>
        </xsl:if>

            <xsl:if test="@maxOccurs='unbounded' or @maxOccurs&gt;1">
              <xsl:text>     
allocate(</xsl:text>
              <xsl:value-of select="$struct"/>
              <xsl:text>%</xsl:text>
              <xsl:value-of select="@name|@ref"/>
              <xsl:text>array(len))
</xsl:text>
            </xsl:if>
            <xsl:if test="not(@maxOccurs='unbounded' or @maxOccurs>1)" >
            <xsl:value-of select="$struct"/>
            <xsl:text>%</xsl:text>
            <xsl:value-of select="@name|@ref"/>
           <xsl:text>=>null()
</xsl:text>
</xsl:if>
            <xsl:text>Do i=0,len-1
</xsl:text>
            <xsl:value-of select="$struct"/>
            <xsl:text>%</xsl:text>
            <xsl:value-of select="@name|@ref"/>
            <xsl:if test="@maxOccurs='unbounded' or @maxOccurs>1">
              <xsl:text>array(i+1)%</xsl:text><xsl:value-of select="@name|@ref"/>
              </xsl:if>
            <xsl:text>=>getstruct</xsl:text>
            <xsl:value-of select="@name|@ref"/>
            <xsl:text>(&amp;
</xsl:text>
            <xsl:text>removeChild(thisnode,item(getElementsByTagname(thisnode,&amp;
"</xsl:text> <xsl:value-of select="@name|@ref"/>
            <xsl:text>"),0)) ) 
</xsl:text>
       <xsl:text>enddo
</xsl:text>   
          </xsl:otherwise>
        </xsl:choose>
      </xsl:for-each>
      <xsl:for-each select="*/*/xs:element[@type or @ref]">
      <xsl:variable name="ref" select="@ref"/>
      <xsl:if test="//xs:element[@name=$ref]/@type or @type">
     
      <xsl:text>
      len= countChildEmentsWithName (thisnode,"</xsl:text>
              <xsl:value-of select="@name|@ref"/>
              <xsl:text>")</xsl:text>
  <xsl:if test="@maxOccurs='unbounded' or @maxOccurs>1"> 
   <xsl:text>           
allocate(</xsl:text> <xsl:value-of select="$struct"/>%<xsl:value-of select="@name|@ref"/>  
<xsl:text>(</xsl:text> 
<xsl:choose>
<xsl:when test="//xs:element[@name=$ref]/@type='integerpair' or @type='integerpair'" >2</xsl:when>
<xsl:otherwise>3</xsl:otherwise>
</xsl:choose><xsl:text>,len))</xsl:text>
</xsl:if>
<xsl:if test="@minOccurs>=1">
<xsl:text>
if (len .lt. </xsl:text><xsl:value-of select="@minOccurs" />
<xsl:text>) then
  write(*,*) "Parser ERROR: "
  Write (*,*)"The Element: </xsl:text>
  <xsl:value-of select="@name|@ref"/> 
  <xsl:text> must occur at least </xsl:text><xsl:value-of select="@minOccurs" /><xsl:text> times in the"
   Write (*,*) "</xsl:text><xsl:value-of select="../../../@name" />
   <xsl:text> element"
  stop
endif</xsl:text>
</xsl:if>
<xsl:text>
Do i=1,len
</xsl:text>
      
       <xsl:value-of select="$newline"/>
      <xsl:value-of select="$struct"/>
            <xsl:text>%</xsl:text>
            <xsl:value-of select="@name|@ref"/> 
            <xsl:if test="@maxOccurs='unbounded' or @maxOccurs&gt;1 ">
            <xsl:text>(:,i)</xsl:text>
            </xsl:if>
      <xsl:text>=getvalueof</xsl:text>
      <xsl:value-of select="@name|@ref"/>  <xsl:text>(&amp;
      removechild(thisnode,item(getElementsByTagname(thisnode,&amp;
      "</xsl:text>
              <xsl:value-of select="@name|@ref"/>
              <xsl:text>"),0)))</xsl:text>
              <xsl:text>
end do
</xsl:text>
            
      </xsl:if>
      
      </xsl:for-each>
      
      <xsl:text>
      i=0
      len=0
      </xsl:text>
      <xsl:if test="@ex:importance!='ignore'">
      call  handleunknownnodes(thisnode)</xsl:if>
      <xsl:text>
end function
</xsl:text>
    </xsl:for-each>
    <xsl:for-each select="//xs:element[@type]">
   
    <xsl:text> 
function getvalueof</xsl:text> <xsl:value-of select="@name"/>
<xsl:text>(thisnode)
implicit none
type(Node),pointer::thisnode
</xsl:text>
<xsl:call-template name="xstypetofortrantype" >
<xsl:with-param name="type" select="@type"/>
    <xsl:with-param name="varname">
     <xsl:text>getvalueof</xsl:text> <xsl:value-of select="@name"/>
    </xsl:with-param>
</xsl:call-template>
<xsl:text>
#ifdef INPUTDEBUG
  write(*,*)"we are at </xsl:text><xsl:value-of select="@name"/><xsl:text>"
#endif  
   call extractDataContent(thisnode,  getvalueof</xsl:text> <xsl:value-of select="@name"/>
<xsl:text>)
end function</xsl:text>
</xsl:for-each>
 <xsl:for-each select="//xs:restriction[xs:enumeration]">
 <xsl:text>
 integer function  stringtonumber</xsl:text>
<xsl:value-of select="../../../../@name"/><xsl:value-of select="../../@name"/>
<xsl:text>(string) 
! Do not edit this Fortran file! It is generated from the Schema.
 character(80),intent(in)::string
 select case(trim(adjustl(string)))
</xsl:text>
<xsl:for-each select="xs:enumeration">
<xsl:text>case('</xsl:text>
<xsl:value-of select="@value"/>
<xsl:text>')
 stringtonumber</xsl:text>
<xsl:value-of select="../../../../../@name"/><xsl:value-of select="../../../@name"/>
<xsl:text>=</xsl:text>
<xsl:choose>
<xsl:when test="xs:annotation/xs:appinfo/oldnr">
<xsl:value-of select="normalize-space( xs:annotation/xs:appinfo/oldnr)"/>
</xsl:when>
<xsl:otherwise>
<xsl:value-of select="-1"/>
</xsl:otherwise>
</xsl:choose>
<xsl:text>
</xsl:text>
</xsl:for-each>

<xsl:text>case('')
 stringtonumber</xsl:text>
<xsl:value-of select="../../../../@name"/><xsl:value-of select="../../@name"/>
<xsl:text>=0
case default
write(*,*) "Parser ERROR: '", string,"' is not valid selection for </xsl:text> <xsl:value-of select="../../../@name|../../@name"/> <xsl:text> "
stop 
end select
end function

 </xsl:text>
 
 </xsl:for-each>

function countChildEmentsWithName(nodep,name)
  implicit none
  integer::countChildEmentsWithName
  type(Node),pointer ::nodep
  character(len=*),intent(in)::name
  type(NodeList),pointer::children
  type(Node),pointer::child
  
  integer::i
  children=>getChildNodes(nodep)
  countChildEmentsWithName=0
  do i=0,getlength(children)-1
    child=>item(children,i)
    if(name.eq.getNodeName(child)) countChildEmentsWithName=countChildEmentsWithName+1
  end do

end function  
    <xsl:if test="$root='input'">
      <xsl:text>
! these are some transient helper functions to simplify the port (should not be used)
function isspinorb()
logical::isspinorb
isspinorb=.false.
if(associated(input%groundstate%spin))then
if (input%groundstate%spin%spinorb) then
isspinorb=.true.
endif
endif
end function
function isspinspiral() 
logical::isspinspiral
isspinspiral=.false.
if(associated(input%groundstate%spin))then
if (input%groundstate%spin%spinsprl) then
isspinspiral=.true.
endif
endif
end function

function getfixspinnumber()
implicit none
integer::getfixspinnumber
getfixspinnumber=0
if(associated(input%groundstate%spin))then
getfixspinnumber=input%groundstate%spin%fixspinnumber
endif
end function

function istetraocc()
  implicit none
  logical ::istetraocc
  istetraocc =.false.
  if(associated(input%xs)) then
    if(associated(input%xs%tetra)) then
      istetraocc=input%xs%tetra%tetraocc
    endif
  endif
end function
</xsl:text>
</xsl:if>
  <xsl:text>
end module

</xsl:text>
  </xsl:template>
  
  <xsl:template name="defaultvaltofortran">
    <xsl:param name="default"/>
    <xsl:param name="type"/>
    <xsl:choose>
      <xsl:when test="$type='vect3d' or $type='integerpair' or $type='vect2d' or $type='integertriple' or $type='integerquadrupel'">
        <xsl:text>(/</xsl:text>
        <xsl:for-each select="str:split($default,' ')">
          <xsl:value-of select="."/>
          <xsl:if test="not(position()=count(str:split($default,' ')) )"><xsl:text>,</xsl:text></xsl:if> 
        </xsl:for-each>
        <xsl:text>/)</xsl:text>
      </xsl:when>
      <xsl:when test="$type='xs:boolean' ">
        <xsl:text>.</xsl:text><xsl:value-of select="$default"/><xsl:text>.</xsl:text>
      </xsl:when>
      <xsl:when test="$type='booleantriple' ">
        <xsl:text>(/</xsl:text>
        <xsl:for-each select="str:split($default,' ')">
          <xsl:text>.</xsl:text>
          <xsl:value-of select="."/>
          <xsl:if test="not(position()=count(str:split($default,' ')) )"><xsl:text>.,</xsl:text></xsl:if> 
        </xsl:for-each>
        <xsl:text>./)</xsl:text>
      </xsl:when>
      <xsl:when test="$type='xs:string'or $type='xs:anyURI'or $type='xs:ID' or $type='xs:IDREFS'">
        <xsl:text> "</xsl:text><xsl:value-of select="$default"/><xsl:text>"</xsl:text>
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="$default"/>
      </xsl:otherwise>
    </xsl:choose>
  </xsl:template>

</xsl:stylesheet>
  
