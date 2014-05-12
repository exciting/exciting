<?xml version="1.0" encoding="UTF-8"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform" version="1.0">

    <xsl:template  name="reports2html" match="/">
        <html>
            <head>
                <script language="javascript"> 
function toggle(showHideDiv, switchTextDiv) {
	var ele = document.getElementById(showHideDiv);
	var text = document.getElementById(switchTextDiv);
	if(ele.style.display == "block") {
    		ele.style.display = "none";
		text.innerHTML = "show";
  	}
	else {
		ele.style.display = "block";
		text.innerHTML = "hide";
	}
} 
            </script>


                <style type="text/css">
                    .testlist {
                        display : none;
                    }
                    h2 {
                        border-top-style : solid;
                    }
                    .test {
                        width : 600px;
                    }</style>
            </head>
            <body>


                <xsl:for-each select="/reports/report">
                    <xsl:variable name="all" select="number(count(test))"/>
                    <xsl:variable name="testrun" select="count(/reports/report) - position()+1"/>
                    <p>
                        <div class="testsrun">
                            <h2> Testrun: <xsl:value-of select="$testrun"/>
                            </h2>
                            <div>
                                <span> date:<xsl:value-of select="@date"/>
                                </span>
                                <span> githash:<a>
                                    <xsl:attribute name="href">
                                        <xsl:text>https://github.com/exciting/exciting/commit/</xsl:text>
                                            <xsl:value-of select="@githash"/>
                                    </xsl:attribute>
                                    <xsl:value-of select="@githash"/></a>
                                </span>
                            </div>
                            <img>
                                <xsl:attribute name="src">
                                    <xsl:text>http://chart.apis.google.com/chart?cht=bhs&amp;chd=t:</xsl:text>
                                    <xsl:value-of
                                        select="number(count(test[status='passed'] )div $all * 100)"/>
                                    <xsl:text>|</xsl:text>
                                    <xsl:value-of
                                        select="number(count(test[status='unspecified'])) div $all *100"/>
                                    <xsl:text>|</xsl:text>
                                    <xsl:value-of
                                        select="number(count(test[status='failed'])) div $all *100"/>
                                    <xsl:text>&amp;chs=</xsl:text>
                                    <xsl:value-of select="$all * 10 "/>
                                    <xsl:text>x80&amp;chdl=passed|unspecified|failed&amp;chdlp=t&amp;chl=passed:</xsl:text>
                                    <xsl:value-of select="count(test[status='passed'])"/>
                                    <xsl:text>|unspecified:</xsl:text>
                                    <xsl:value-of select="count(test[status='uncpecified'])"/>
                                    <xsl:text>|failed:</xsl:text>
                                    <xsl:value-of select="count(test[status='failed'])"/>
                                    <xsl:text>&amp;chco=006600,f0f000,cc0033</xsl:text>
                                </xsl:attribute>
                            </img>
                            <div>
                                <a>
                                    <xsl:attribute name="href">
                                        <xsl:text>javascript:toggle('</xsl:text>
                                        <xsl:value-of select="position()"/>
                                        <xsl:text>passed','')</xsl:text>
                                    </xsl:attribute>

                                    <xsl:text> passed: </xsl:text>
                                    <xsl:value-of select="count(test[status='passed'])"/>
                                </a>
                                <xsl:text> </xsl:text>
                                <a>
                                    <xsl:attribute name="href">
                                        <xsl:text>javascript:toggle('</xsl:text>
                                        <xsl:value-of select="position()"/>
                                        <xsl:text>unspecified','')</xsl:text>
                                    </xsl:attribute>
                                    <xsl:text>unspecified: </xsl:text>
                                    <xsl:value-of select="count(test[status='unspecified'])"/>
                                </a>
                                <xsl:text> </xsl:text>
                                <a>
                                    <xsl:attribute name="href">
                                        <xsl:text>javascript:toggle('</xsl:text>
                                        <xsl:value-of select="position()"/>
                                        <xsl:text>failed','')</xsl:text>
                                    </xsl:attribute>
                                    <xsl:text>failed: </xsl:text>
                                    <xsl:value-of select="count(test[status='failed'])"/>
                                </a>
                            </div>

                        </div>
                        <div class="testlist">
                            <xsl:attribute name="id">
                                <xsl:value-of select="position()"/>
                                <xsl:text>failed</xsl:text>
                            </xsl:attribute>
                            <xsl:for-each select="test[status='failed']">
                                <xsl:call-template name="showtest">
                                    <xsl:with-param name="test" select="."/>
                                    <xsl:with-param name="color">PeachPuff</xsl:with-param>
                                </xsl:call-template>
                            </xsl:for-each>

                        </div>
                        <div class="testlist">
                            <xsl:attribute name="id">
                                <xsl:value-of select="position()"/>
                                <xsl:text>unspecified</xsl:text>
                            </xsl:attribute>
                            <xsl:for-each select="test[status='unspecified']">
                                <xsl:call-template name="showtest">
                                    <xsl:with-param name="test" select="."/>
                                    <xsl:with-param name="color">LemonChiffon</xsl:with-param>
                                </xsl:call-template>
                            </xsl:for-each>

                        </div>
                        <div class="testlist">
                            <xsl:attribute name="id">
                                <xsl:value-of select="position()"/>
                                <xsl:text>passed</xsl:text>
                            </xsl:attribute>
                            <xsl:for-each select="test[status='passed']">
                                <xsl:call-template name="showtest">
                                    <xsl:with-param name="test" select="."/>
                                    <xsl:with-param name="color">LightGreen</xsl:with-param>
                                </xsl:call-template>
                            </xsl:for-each>

                        </div>
                    </p>
                </xsl:for-each>
            </body>
        </html>
    </xsl:template>
    <xsl:template name="showtest">
        <xsl:param name="test"/>
        <xsl:param name="color"/>
        <div class="test">
            <div>

                <span>
                    <h3>
                        <xsl:attribute name="style">
                            <xsl:text>background-color:</xsl:text>
                            <xsl:value-of select="$color"/>
                        </xsl:attribute>
                        <xsl:value-of select="position()"/>: <xsl:value-of select="$test/name"
                        /></h3>
                </span>
                <xsl:text> </xsl:text>
                <span>
                   dir: <xsl:value-of select="$test/directory"/>
                </span>
            </div>
            <div>
                <xsl:value-of select="$test/description"/>
            </div>
        </div>
    </xsl:template>
</xsl:stylesheet>
