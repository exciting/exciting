PREFIX=ref:
STRPARAM=  --stringparam importancelevels "essential expert" --stringparam prefix "$(PREFIX)"
wiki/*: *.xsd ../schematowikidot.xsl common species
	xsltproc $(STRPARAM) ../schematowikidot.xsl   input.xsd >wiki/$(PREFIX)input;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   structure.xsd >wiki/$(PREFIX)structure;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   groundstate.xsd >wiki/$(PREFIX)groundstate;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   phonons.xsd >wiki/$(PREFIX)phonons;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   properties.xsd >wiki/$(PREFIX)properties;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   relax.xsd > wiki/$(PREFIX)relax;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   xs.xsd >wiki/$(PREFIX)xs;\
	xsltproc $(STRPARAM) ../schematowikidot.xsl   gw.xsd >wiki/$(PREFIX)gw;\
	
species:expandedschema
	xsltproc --stringparam prefix "/" ../schematowikidot.xsl  ../species.xsd >wiki/$(PREFIX)species

common:expandedschema
	xsltproc $(STRPARAM) --stringparam common "true"  ../schematowikidot.xsl  ../excitinginput.xsd >wiki/$(PREFIX)common
	
expandedschema:: 
	cd ../../ && $(MAKE)  expandedschema
