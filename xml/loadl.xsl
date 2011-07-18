<?xml version="1.0" encoding="UTF-8" ?>
<xsl:stylesheet version="1.0"
        xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text" />

<xsl:template match="/">

# @ job_type = parallel

# @ class = lessthour
# @ node           = 1
# @ tasks_per_node = 1
# @ arguments= 
# @ executable = PATH/TO/EXECUTABLE
<xsl:for-each select = "//@path">

# @ initialdir =  /PATH/TO/WORKINGDIR/<xsl:value-of select="../path"/>

# @ job_name  = example_<xsl:value-of select="../@name"/>
# @ output = $(job_name).out
# @ error = $(job_name).err
# @ resources = ConsumableCpus(1)
# @ queue

</xsl:for-each>
</xsl:template>
</xsl:stylesheet>
