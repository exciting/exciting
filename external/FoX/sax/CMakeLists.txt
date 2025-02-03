set(fox_sax_srcs
  FoX_sax.f90
)

set(fox_sax_srcs_pp
  m_sax_operate.F90
  m_sax_parser.F90
  m_sax_reader.F90
  m_sax_tokenizer.F90
  m_sax_types.F90
  m_sax_xml_source.F90
)
  
MK_LIBS(sax common utils fsys)