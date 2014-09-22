<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform" xmlns:str="http://exslt.org/strings" xmlns:math="http://exslt.org/math">
  <xsl:output method="text" />
  <!-- usage: xsltproc xmlinput2xsf.xsl input.xml -->
  <xsl:template name="norm">
    <xsl:param name="vectorstring" />
    <xsl:value-of
      select="math:sqrt(
math:power(str:tokenize($vectorstring)[1]*$scale,2)
+math:power(str:tokenize($vectorstring)[2]*$scale,2)
+math:power(str:tokenize($vectorstring)[3]*$scale,2)
)*$bohr2angstr" />
  </xsl:template>
  <xsl:variable name="bohr2angstr" select="0.529177" />
  <xsl:variable name="scale">
    <xsl:choose>
      <xsl:when test="/input/structure/crystal/@scale">
        <xsl:value-of select="/input/structure/crystal/@scale" />
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="1" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>
   <xsl:variable name="stretch">
    <xsl:choose>
      <xsl:when test="/input/structure/crystal/@stretch">
        <xsl:value-of select="/input/structure/crystal/@stretch" />
      </xsl:when>
      <xsl:otherwise>
        <xsl:value-of select="'1 1 1'" />
      </xsl:otherwise>
    </xsl:choose>
  </xsl:variable>
  <xsl:variable name="a">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="/input/structure/crystal/basevect[1]" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="b">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="/input/structure/crystal/basevect[2]" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="c">
    <xsl:call-template name="norm">
      <xsl:with-param name="vectorstring" select="/input/structure/crystal/basevect[3]" />
    </xsl:call-template>
  </xsl:variable>
  <xsl:variable name="cartesian" as="xs:boolean" select="/input/structure/@cartesian" /> 
  <xsl:template match="/">
    <xsl:text>CRYSTAL
PRIMVEC
</xsl:text>
    <!-- Convert vectors "basevect" into Angstrom, and multiply with 
         scale and corresponing stretch to get the actual basevectors,
         which are written to the xsf file. -->
    <xsl:for-each select="/input/structure/crystal/basevect">
    <xsl:variable name="basevn" select="position()"/>
      <xsl:for-each select="str:tokenize(.)">
        <xsl:value-of select="$scale*$bohr2angstr *./.* str:tokenize($stretch)[$basevn]" />
        <xsl:text>   </xsl:text>
      </xsl:for-each>
      <xsl:text>
</xsl:text>
    </xsl:for-each>
    <!-- Once again: Convert vectors "basevect" into Angstrom, and 
         multiply with scale and corresponing stretch to get the 
         actual basevectors, which are now saved in the variables 
         (a1,a2,a3), (b1,b2,b3), (c1,c2,c3) -->
    <xsl:variable name="a1" select=
      "str:tokenize(/input/structure/crystal/basevect[1])[1]*
      $bohr2angstr*$scale*str:tokenize($stretch)[1]" />
    <xsl:variable name="a2" select=
      "str:tokenize(/input/structure/crystal/basevect[1])[2]*
      $bohr2angstr*$scale*str:tokenize($stretch)[1]" />
    <xsl:variable name="a3" select=
      "str:tokenize(/input/structure/crystal/basevect[1])[3]*
      $bohr2angstr*$scale*str:tokenize($stretch)[1]" />
    <xsl:variable name="b1" select=
      "str:tokenize(/input/structure/crystal/basevect[2])[1]*
      $bohr2angstr*$scale*str:tokenize($stretch)[2]" />
    <xsl:variable name="b2" select=
      "str:tokenize(/input/structure/crystal/basevect[2])[2]*
      $bohr2angstr*$scale*str:tokenize($stretch)[2]" />
    <xsl:variable name="b3" select=
      "str:tokenize(/input/structure/crystal/basevect[2])[3]*
      $bohr2angstr*$scale*str:tokenize($stretch)[2]" />
    <xsl:variable name="c1" select=
      "str:tokenize(/input/structure/crystal/basevect[3])[1]*
      $bohr2angstr*$scale*str:tokenize($stretch)[3]" />
    <xsl:variable name="c2" select=
      "str:tokenize(/input/structure/crystal/basevect[3])[2]*
      $bohr2angstr*$scale*str:tokenize($stretch)[3]" />
    <xsl:variable name="c3" select=
      "str:tokenize(/input/structure/crystal/basevect[3])[3]*
      $bohr2angstr*$scale*str:tokenize($stretch)[3]" />
    <!--  Write out coordinates of atoms in the primitive cell;
          coordinates are in Cartesian coordinates and in Angstrom -->
    <!--  Note redundant line breaks not supported by Vesta viewer! -->
    <xsl:text>PRIMCOORD
</xsl:text>
    <xsl:value-of select="count(input/structure/species/atom)" />
    <xsl:text> 1</xsl:text>
    <xsl:for-each select="input/structure/species/atom">
      <xsl:text>
</xsl:text>
    <!--  Replace atomic symbol with Z_nuc as Jmol viewer does not understand atomic symbols -->
    <xsl:variable name="element" select="substring-before(../@speciesfile, '.')"/>
    <xsl:choose>
    <xsl:when test='$element="X"'>0</xsl:when>
    <xsl:when test='$element="H"'>1</xsl:when>
    <xsl:when test='$element="He"'>2</xsl:when>
    <xsl:when test='$element="Li"'>3</xsl:when>
    <xsl:when test='$element="Be"'>4</xsl:when>
    <xsl:when test='$element="B"'>5</xsl:when>
    <xsl:when test='$element="C"'>6</xsl:when>
    <xsl:when test='$element="N"'>7</xsl:when>
    <xsl:when test='$element="O"'>8</xsl:when>
    <xsl:when test='$element="F"'>9</xsl:when>
    <xsl:when test='$element="Ne"'>10</xsl:when>
    <xsl:when test='$element="Na"'>11</xsl:when>
    <xsl:when test='$element="Mg"'>12</xsl:when>
    <xsl:when test='$element="Al"'>13</xsl:when>
    <xsl:when test='$element="Si"'>14</xsl:when>
    <xsl:when test='$element="P"'>15</xsl:when>
    <xsl:when test='$element="S"'>16</xsl:when>
    <xsl:when test='$element="Cl"'>17</xsl:when>
    <xsl:when test='$element="Ar"'>18</xsl:when>
    <xsl:when test='$element="K"'>19</xsl:when>
    <xsl:when test='$element="Ca"'>20</xsl:when>
    <xsl:when test='$element="Sc"'>21</xsl:when>
    <xsl:when test='$element="Ti"'>22</xsl:when>
    <xsl:when test='$element="V"'>23</xsl:when>
    <xsl:when test='$element="Cr"'>24</xsl:when>
    <xsl:when test='$element="Mn"'>25</xsl:when>
    <xsl:when test='$element="Fe"'>26</xsl:when>
    <xsl:when test='$element="Co"'>27</xsl:when>
    <xsl:when test='$element="Ni"'>28</xsl:when>
    <xsl:when test='$element="Cu"'>29</xsl:when>
    <xsl:when test='$element="Zn"'>30</xsl:when>
    <xsl:when test='$element="Ga"'>31</xsl:when>
    <xsl:when test='$element="Ge"'>32</xsl:when>
    <xsl:when test='$element="As"'>33</xsl:when>
    <xsl:when test='$element="Se"'>34</xsl:when>
    <xsl:when test='$element="Br"'>35</xsl:when>
    <xsl:when test='$element="Kr"'>36</xsl:when>
    <xsl:when test='$element="Rb"'>37</xsl:when>
    <xsl:when test='$element="Sr"'>38</xsl:when>
    <xsl:when test='$element="Y"'>39</xsl:when>
    <xsl:when test='$element="Zr"'>40</xsl:when>
    <xsl:when test='$element="Nb"'>41</xsl:when>
    <xsl:when test='$element="Mo"'>42</xsl:when>
    <xsl:when test='$element="Tc"'>43</xsl:when>
    <xsl:when test='$element="Ru"'>44</xsl:when>
    <xsl:when test='$element="Rh"'>45</xsl:when>
    <xsl:when test='$element="Pd"'>46</xsl:when>
    <xsl:when test='$element="Ag"'>47</xsl:when>
    <xsl:when test='$element="Cd"'>48</xsl:when>
    <xsl:when test='$element="In"'>49</xsl:when>
    <xsl:when test='$element="Sn"'>50</xsl:when>
    <xsl:when test='$element="Sb"'>51</xsl:when>
    <xsl:when test='$element="Te"'>52</xsl:when>
    <xsl:when test='$element="I"'>53</xsl:when>
    <xsl:when test='$element="Xe"'>54</xsl:when>
    <xsl:when test='$element="Cs"'>55</xsl:when>
    <xsl:when test='$element="Ba"'>56</xsl:when>
    <xsl:when test='$element="La"'>57</xsl:when>
    <xsl:when test='$element="Ce"'>58</xsl:when>
    <xsl:when test='$element="Pr"'>59</xsl:when>
    <xsl:when test='$element="Nd"'>60</xsl:when>
    <xsl:when test='$element="Pm"'>61</xsl:when>
    <xsl:when test='$element="Sm"'>62</xsl:when>
    <xsl:when test='$element="Eu"'>63</xsl:when>
    <xsl:when test='$element="Gd"'>64</xsl:when>
    <xsl:when test='$element="Tb"'>65</xsl:when>
    <xsl:when test='$element="Dy"'>66</xsl:when>
    <xsl:when test='$element="Ho"'>67</xsl:when>
    <xsl:when test='$element="Er"'>68</xsl:when>
    <xsl:when test='$element="Tm"'>69</xsl:when>
    <xsl:when test='$element="Yb"'>70</xsl:when>
    <xsl:when test='$element="Lu"'>71</xsl:when>
    <xsl:when test='$element="Hf"'>72</xsl:when>
    <xsl:when test='$element="Ta"'>73</xsl:when>
    <xsl:when test='$element="W"'>74</xsl:when>
    <xsl:when test='$element="Re"'>75</xsl:when>
    <xsl:when test='$element="Os"'>76</xsl:when>
    <xsl:when test='$element="Ir"'>77</xsl:when>
    <xsl:when test='$element="Pt"'>78</xsl:when>
    <xsl:when test='$element="Au"'>79</xsl:when>
    <xsl:when test='$element="Hg"'>80</xsl:when>
    <xsl:when test='$element="Tl"'>81</xsl:when>
    <xsl:when test='$element="Pb"'>82</xsl:when>
    <xsl:when test='$element="Bi"'>83</xsl:when>
    <xsl:when test='$element="Po"'>84</xsl:when>
    <xsl:when test='$element="At"'>85</xsl:when>
    <xsl:when test='$element="Rn"'>86</xsl:when>
    <xsl:when test='$element="Fr"'>87</xsl:when>
    <xsl:when test='$element="Ra"'>88</xsl:when>
    <xsl:when test='$element="Ac"'>89</xsl:when>
    <xsl:when test='$element="Th"'>90</xsl:when>
    <xsl:when test='$element="Pa"'>91</xsl:when>
    <xsl:when test='$element="U"'>92</xsl:when>
    <xsl:when test='$element="Np"'>93</xsl:when>
    <xsl:when test='$element="Pu"'>94</xsl:when>
    <xsl:when test='$element="Am"'>95</xsl:when>
    <xsl:when test='$element="Cm"'>96</xsl:when>
    <xsl:when test='$element="Bk"'>97</xsl:when>
    <xsl:when test='$element="Cf"'>98</xsl:when>
    <xsl:when test='$element="Es"'>99</xsl:when>
    <xsl:when test='$element="Fm"'>100</xsl:when>
    <xsl:when test='$element="Md"'>101</xsl:when>
    <xsl:when test='$element="No"'>102</xsl:when>
    <xsl:when test='$element="Lr"'>103</xsl:when>
    <xsl:otherwise><xsl:message terminate="yes">Z_nuc for chemical element <xsl:value-of select="$element" /> is unknown!</xsl:message></xsl:otherwise>
      </xsl:choose>
      <xsl:text>    </xsl:text>
      <xsl:choose>
        <xsl:when test="$cartesian">
          <xsl:value-of select="str:tokenize(@coord)[1]*$bohr2angstr" />
          <xsl:text>  </xsl:text>
          <xsl:value-of select="str:tokenize(@coord)[2]*$bohr2angstr" />
          <xsl:text>  </xsl:text>
          <xsl:value-of select="str:tokenize(@coord)[3]*$bohr2angstr" />
        </xsl:when>
        <xsl:otherwise>
          <xsl:value-of select="str:tokenize(@coord)[1]*$a1 + str:tokenize(@coord)[2]*$b1  + str:tokenize(@coord)[3]*$c1" />
          <xsl:text>  </xsl:text>
          <xsl:value-of select="str:tokenize(@coord)[1]*$a2 + str:tokenize(@coord)[2]*$b2  + str:tokenize(@coord)[3]*$c2" />
          <xsl:text>  </xsl:text>
          <xsl:value-of select="str:tokenize(@coord)[1]*$a3 + str:tokenize(@coord)[2]*$b3  + str:tokenize(@coord)[3]*$c3" />
        </xsl:otherwise>
      </xsl:choose>
    </xsl:for-each>
    <xsl:text>
    
   </xsl:text>
  </xsl:template>
</xsl:stylesheet>
