
use XML::Simple;
use XML::Writer;
use IO::File;

  my $output = new IO::File(">output.xml");

  my $writer = new XML::Writer(OUTPUT => $output,DATA_MODE => 'true', DATA_INDENT => 2);
  $writer->xmlDecl( 'UTF-8' );
  $writer->startTag("report");


  $writer->startTag("test", 
                    "directory" => "test01", name=>"my name to test");
  $writer->startTag("description");     
  $writer->characters("Hello, world! ><Üeß");
  $writer->endTag("description");
  $writer->startTag("status");
  $writer->characters("passed");
  $writer->endTag("status");  
  $writer->endTag("test");  


  $writer->endTag("report");
  $writer->end();
  $output->close();
my @report;

#do things to assert test


