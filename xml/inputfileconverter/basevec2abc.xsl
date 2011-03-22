<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:math="http://exslt.org/math" xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
    xmlns:str="http://exslt.org/strings" version="1.0">
    <xsl:template name="norm">
        <xsl:param name="vectorstring"/>
        <xsl:value-of
            select="math:sqrt(             math:power(str:tokenize($vectorstring)[1],2)             +math:power(str:tokenize($vectorstring)[2],2)             +math:power(str:tokenize($vectorstring)[3],2)             )"
        />
    </xsl:template>
    <xsl:variable name="bohr2angstr" select="0.529177"/>
    <xsl:variable name="scale">
        <xsl:choose>
            <xsl:when test="/input/structure/crystal/@scale">
                <xsl:value-of select="/input/structure/crystal/@scale"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="1"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>
    <xsl:variable name="stretch">
        <xsl:choose>
            <xsl:when test="/input/structure/crystal/@stretch">
                <xsl:value-of select="/input/structure/crystal/@stretch"/>
            </xsl:when>
            <xsl:otherwise>
                <xsl:value-of select="'1 1 1'"/>
            </xsl:otherwise>
        </xsl:choose>
    </xsl:variable>
    <xsl:variable name="a">
        <xsl:call-template name="norm">
            <xsl:with-param name="vectorstring">
                <xsl:call-template name="scaledbasevec">
                    <xsl:with-param name="component" select="1"/>
                </xsl:call-template>
            </xsl:with-param>
        </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="b">
        <xsl:call-template name="norm">
            <xsl:with-param name="vectorstring">
                <xsl:call-template name="scaledbasevec">
                    <xsl:with-param name="component" select="2"/>
                </xsl:call-template>
            </xsl:with-param>
        </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="c">
        <xsl:call-template name="norm">
            <xsl:with-param name="vectorstring">
                <xsl:call-template name="scaledbasevec">
                    <xsl:with-param name="component" select="3"/>
                </xsl:call-template>
            </xsl:with-param>
        </xsl:call-template>
    </xsl:variable>
    <xsl:template name="scaledbasevec">
        <xsl:param name="component"/>
        <xsl:value-of
            select="str:tokenize(/input/structure/crystal/basevect[$component])[1]* $scale*str:tokenize($stretch)[$component]*$bohr2angstr"/>
        <xsl:text> </xsl:text>
        <xsl:value-of
            select="str:tokenize(/input/structure/crystal/basevect[$component])[2]* $scale*str:tokenize($stretch)[$component]*$bohr2angstr"/>
        <xsl:text> </xsl:text>
        <xsl:value-of
            select="str:tokenize(/input/structure/crystal/basevect[$component])[3] * $scale*str:tokenize($stretch)[$component]*$bohr2angstr"/>
        <xsl:text> </xsl:text>
    </xsl:template>
    <xsl:template name="angle">
        <xsl:param name="vecs1"/>
        <xsl:param name="vecs2"/>
        <xsl:variable name="vec1" select="str:tokenize($vecs1)"/>
        <xsl:variable name="vec2" select="str:tokenize($vecs2)"/>
        <xsl:variable name="norm1">
            <xsl:call-template name="norm">
                <xsl:with-param name="vectorstring" select="$vecs1"/>
            </xsl:call-template>
        </xsl:variable>
        <xsl:variable name="norm2">
            <xsl:call-template name="norm">
                <xsl:with-param name="vectorstring" select="$vecs2"/>
            </xsl:call-template>
        </xsl:variable>
        <xsl:value-of
            select="57.2957795130823*math:acos(($vec1[1]*$vec2[1]+$vec1[2]*$vec2[2]+$vec1[3]*$vec2[3])div $norm1 div $norm2)"
        />
    </xsl:template>
    <xsl:variable name="alpha">
        <xsl:call-template name="angle">
            <xsl:with-param name="vecs1" select="/input/structure/crystal/basevect[2]"/>
            <xsl:with-param name="vecs2" select="/input/structure/crystal/basevect[3]"/>
        </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="beta">
        <xsl:call-template name="angle">
            <xsl:with-param name="vecs1" select="/input/structure/crystal/basevect[3]"/>
            <xsl:with-param name="vecs2" select="/input/structure/crystal/basevect[1]"/>
        </xsl:call-template>
    </xsl:variable>
    <xsl:variable name="gamma">
        <xsl:call-template name="angle">
            <xsl:with-param name="vecs1" select="/input/structure/crystal/basevect[1]"/>
            <xsl:with-param name="vecs2" select="/input/structure/crystal/basevect[2]"/>
        </xsl:call-template>
    </xsl:variable>
</xsl:stylesheet>
