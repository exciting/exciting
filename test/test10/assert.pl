use lib "../perl/";
use lib "../perl/lib/";
use XML::Simple;
use XML::Writer;
use IO::File;
use Test;

$writer= Test::initreport("report.xml");

%statusse=Test::assert_file_same_within( "./reference/RHO3D.xml",
	"run_properties/RHO3D.xml",1.0);
 Test::writetestreport({
 		"directory"=>"test10/ ",
 		"name"=>"Properties: chargedensityplot",
 		"description"=>"The test is passed if the output file differs
 		 less than $tol from the reference one". %statusse->{maxerror},
 		"status"=> %statusse->{status}}, $writer);

 %statusse=Test::assert_file_same_within( "./reference/VCL3D.xml",
  "run_properties/VCL3D.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: exccplot, VCL3D",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

%statusse=Test::assert_file_same_within( "./reference/VXC3D.xml",
  "run_properties/VXC3D.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: exccplot, VXC3D",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/WF3D.xml",
  "run_properties/WF3D.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: wfplot",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/LSJ.REF",
  "run_properties/LSJ.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: LSJ",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/ELF3D.xml",
  "run_properties/ELF3D.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: elfplot",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/EF3D.xml",
  "run_properties/EF3D.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: electricfield",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/EFG.REF",
  "run_properties/EFG.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: EFG",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/MOSSBAUER.REF",
  "run_properties/MOSSBAUER.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: mossbauer",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/EXPIQR.REF",
  "run_properties/EXPIQR.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: expiqr",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/EPSILON_33.REF",
  "run_properties/EPSILON_33.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: dielmat",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/KERR.REF",
  "run_properties/KERR.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: moke",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

%statusse=Test::assert_file_same_within( "./reference/CHI_111.REF",
  "run_properties/CHI_111.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: shg",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer); 

%statusse=Test::assert_file_same_within( "./reference/EFFMASS.REF",
  "run_properties/EFFMASS.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: masstensor",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

%statusse=Test::assert_file_same_within( "./reference/ELNES.REF",
  "run_properties/ELNES.OUT",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: elnes",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

 %statusse=Test::assert_file_same_within( "./reference/bandstructure.xml",
  "run_properties/bandstructure.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: bandstructure",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

 %statusse=Test::assert_file_same_within( "./reference/dos.xml",
  "run_properties/dos.xml",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: dos",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

 %statusse=Test::assert_file_same_within( "./reference/FERMISURF.bxsf",
  "run_properties/FERMISURF.bxsf",1.0);
 Test::writetestreport({
    "directory"=>"test10/ ",
    "name"=>"Properties: fermisurfaceplot",
    "description"=>"The test is passed if the output file differs
     less than $tol from the reference one". %statusse->{maxerror},
    "status"=> %statusse->{status}}, $writer);

 Test::closereport($writer);
