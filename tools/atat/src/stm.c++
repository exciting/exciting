#include <fstream.h>
#include "parse.h"
#include "getvalue.h"
#include "version.h"

Real interpol(Real x,Real x0,Real x1,Real y0,Real y1){
  return (y0+(y1-y0)*(x-x0)/(x1-x0));
}

int main(int argc, char *argv[]) {
  Real density=0.;
  char *strfilename="str.out";
  char *regionfilename="region.in";
  char *hrangefilename="";
  char *afterlabel="";
  int pdensityrange=0;
  int beginz=0;
  int roll=0;
  AskStruct options[]={
    {"","STM simulator " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-s","Input file defining the structure geometry (default: str.out)",STRINGVAL,&strfilename},
    {"-r","Input file defining the region scanned (default: region.in)",STRINGVAL,&regionfilename},
    {"-h","Input file defining the height range (default: calculated from input)",STRINGVAL,&hrangefilename},
    {"-a","Start read after this label",STRINGVAL,&afterlabel},
    {"-dr","Print density range",BOOLVAL,&pdensityrange},
    {"-id","Iso-density level",REALVAL,&density},
    {"-bz","",INTVAL,&beginz},
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }
  if (file_exists("roll.in")) {
    ifstream file("roll.in");
    file >> roll;
  }
  for (int i=0; i<strlen(afterlabel); i++) {
    if (afterlabel[i]=='_') {afterlabel[i]=' ';}
  }
  istream &infile=cin;
  AutoString line;
  do {
    get_string(&line,infile,"\n");
    skip_delim(infile,"\n");
  } while (!infile.eof() && strcmp(line,afterlabel)!=0);
  int nx,ny,nz;
  infile >> nx >> ny >> nz;
  Real mind=MAXFLOAT;
  Real maxd=-MAXFLOAT;
  Array<Array<Array<Real> > > raw;
  raw.resize(nz);
  for (int iz=0; iz<nz; iz++) {
    int riz=(nz+iz+roll)%nz;
    raw(riz).resize(ny);
    for (int iy=0; iy<ny; iy++) {
      raw(riz)(iy).resize(nx);
      for (int ix=0; ix<nx; ix++) {
        infile >> raw(riz)(iy)(ix);
	mind=MIN(mind,raw(riz)(iy)(ix));
	maxd=MAX(maxd,raw(riz)(iy)(ix));
      }
    }
  }
  if (pdensityrange) {
    cout << mind << " " << maxd << endl;
    exit(0);
  }
  Array2d<Real> h(nx,ny);
  for (int ix=0; ix<nx; ix++) {
    for (int iy=0; iy<ny; iy++) {
      int iz=nz-1-beginz;
      while (iz>=0 && raw(iz)(iy)(ix)<density) {
        iz--;
      }
      if (iz==(nz-1)) {h(ix,iy)=iz;}
      else if (iz<0) {h(ix,iy)=0;}
      else {
        h(ix,iy)=interpol(density,raw(iz)(iy)(ix),raw(iz+1)(iy)(ix),(Real)(iz),(Real)(iz+1));
      }
    }
  }

  ifstream strfile(strfilename);
  if (!strfile) {ERRORQUIT("Unable to open structure file");}
  rMatrix3d axes;
  rMatrix3d cell;
  read_cell(&axes,strfile);
  read_cell(&cell,strfile);
  cell=axes*cell;
  rMatrix3d icell=!cell;
  
  ifstream regionfile(regionfilename);
  if (!regionfile) {ERRORQUIT("Unable to open region file");}
  rVector3d v0,v1,v2;
  regionfile >> v0; v0=axes*v0;
  regionfile >> v1; v1=axes*v1;
  regionfile >> v2; v2=axes*v2;
  Real step=sqrt(norm(cell.get_column(0)^cell.get_column(1))/(Real)(nx*ny))/2.;
  int n1=(int)(norm(v1)/step);
  int n2=(int)(norm(v2)/step);
  rVector3d dv1,dv2;
  dv1=v1/n1;
  dv2=v2/n2;
  Real zscale=norm(cell.get_column(2))/((Real)nz);
  Array2d<Real> stm(n1,n2);
  Real minh=MAXFLOAT;
  Real maxh=-MAXFLOAT;
  for (int i1=0; i1<n1; i1++) {
    for (int i2=0; i2<n2; i2++) {
      rVector3d v=mod1(icell*(v0+i1*dv1+i2*dv2));
      int ix=(int)floor(v(0)*nx);
      int iy=(int)floor(v(1)*ny);
      Real fx=v(0)*nx-(Real)ix;
      Real fy=v(1)*ny-(Real)iy;
      Real a=interpol(fx,0.0,1.0,h(ix%nx,iy%ny),h((ix+1)%nx,iy%ny));
      Real b=interpol(fx,0.0,1.0,h(ix%nx,(iy+1)%ny),h((ix+1)%nx,(iy+1)%ny));
      stm(i1,i2)=interpol(fy,0.0,1.0,a,b)*zscale;
      minh=MIN(minh,stm(i1,i2));
      maxh=MAX(maxh,stm(i1,i2));
    }
  }
  if (strlen(hrangefilename)>0) {
    ifstream hfile(hrangefilename);
    if (!hfile) {ERRORQUIT("Unable to open height file");}
    hfile >> minh >> maxh;
  }
  {
    ofstream hfile("hrange.out");
    hfile << minh << " " << maxh << endl;
  }
  ostream &outfile=cout;
  int linecutlim=6;
  Real maxgray=255.;
  outfile << "P2" << endl;
  outfile << "# STM image" << endl;
  outfile << n1 << " " << n2 << endl;
  outfile << maxgray << endl;
  int linecut=0;
  for (int i2=0; i2<n2; i2++) {
    for (int i1=0; i1<n1; i1++) {
      Real g=((stm(i1,i2)-minh)*maxgray/(maxh-minh));
      if (g<0.) {g=0.;}
      if (g>=maxgray) {g=maxgray;}
      outfile << (int)g;
      linecut++;
      if (linecut>=linecutlim) {
        outfile << endl;
        linecut=0;
      }
      else {
        outfile << " ";
      }
    }
  }
}
