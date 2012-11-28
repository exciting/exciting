#include <fstream.h>
#include <strstream.h>
#include "lstsqr.h"
#include "getvalue.h"
#include "version.h"
#include "ridge.h"

int main(int argc, char *argv[]) {
  char *xfilename="";
  char *yfilename="";
  char *wfilename="";
  char *rfilename="";
  char *powers="";
  char *select="";
  int addone=0;
  int del_colin=0;
  int perfit=0;
  int docv=0;
  int dostd=0;
  int docopyy=0;
  int dopred=0;
  int doresid=0;
  AskStruct options[]={
    {"","Least-Square FIT " MAPS_VERSION ", by Axel van de Walle",TITLEVAL,NULL},
    {"-x","x file",STRINGVAL,&xfilename},
    {"-y","y file",STRINGVAL,&yfilename},
    {"-w","weight file",STRINGVAL,&wfilename},
    {"-r","regularization parameters file",STRINGVAL,&rfilename},
    {"-pf","Perfect Fit algorithm",BOOLVAL,&perfit},
    {"-pw","number of powers of each column to regress on",STRINGVAL,&powers},
    {"-s","select columns to regress on (one number by column, 0: ignore, 1: use)",STRINGVAL,&select},
    {"-1","add a column of ones as regressor",BOOLVAL,&addone},
    {"-colin","Ignore colinear columns in x file",BOOLVAL,&del_colin},
    {"-cv","Print crossvalidation score",BOOLVAL,&docv},
    {"-se","Print standard errors",BOOLVAL,&dostd},
    {"-ty","Print trues values of y",BOOLVAL,&docopyy},
    {"-p","Print predicted values of y",BOOLVAL,&dopred},
    {"-e","Print prediction error",BOOLVAL,&doresid}
  };
  if (!get_values(argc,argv,countof(options),options)) {
    display_help(countof(options),options);
    return 1;
  }

  {
    ostrstream cmd;
    cmd << "echo `cat " << xfilename << " | wc -l` "
	<< "`head -1 " << xfilename << " | wc -w` "
      << "`head -1 " << yfilename << " | wc -w` > tmp.tmp"
 << '\0';
    system(cmd.str());
  }
  int nb_param,nb_eqn,nb_ycol;
  {
    ifstream size("tmp.tmp");
    size >> nb_eqn >> nb_param >> nb_ycol;
  }
  unlink("tmp.tmp");
  Array2d<Real> x(nb_eqn,nb_param);
  ifstream xfile(xfilename);
  if (!xfile) {ERRORQUIT("Unable to open x file");}
  for (int eqn=0; eqn<nb_eqn; eqn++) {
    for (int param=0; param<nb_param; param++) {
      xfile >> x(eqn,param);
    }
  }
  if (strlen(powers)>0) {
    Array<int> powint(nb_param);
    istrstream powline(powers);
    int totpow=0;
    for (int i=0; i<nb_param; i++) {
      powline >> powint(i);
      if (!powline) powint(i)=1;
      totpow+=powint(i);
    }
    Array2d<Real> tmpx(nb_eqn,totpow);
    int col=0;
    for (int i=0; i<nb_param; i++) {
      for (int p=1; p<=powint(i); p++) {
	for (int j=0; j<nb_eqn; j++) {
	  tmpx(j,col)=pow(x(j,i),(Real)p);
	}
	col++;
      }
    }
    x=tmpx;
    nb_param=totpow;
  }
  if (addone) {
    Array2d<Real> tmpx(nb_eqn,nb_param+1);
    for (int j=0; j<nb_eqn; j++) {
      tmpx(j,0)=1.;
    }
    for (int i=0; i<nb_param; i++) {
      for (int j=0; j<nb_eqn; j++) {
	tmpx(j,i+1)=x(j,i);
      }
    }
    x=tmpx;
    nb_param++;
  }
  Array<int> selectcol(nb_param);
  one_array(&selectcol);
  if (strlen(select)>0) {
    zero_array(&selectcol);
    istrstream selectline(select);
    int nb_col=0;
    for (int i=0; i<nb_param; i++) {
      selectline >> selectcol(i);
      skip_delim(selectline,",; ");
      nb_col+=(selectcol(i)>0 ? 1 : 0);
    }
    Array2d<Real> tmpx(nb_eqn,nb_col);
    int i=0;
    for (int j=0; j<nb_param; j++) {
      if (selectcol(j)==1) {
	for (int k=0; k<nb_eqn; k++) {
	  tmpx(k,i)=x(k,j);
	}
	i++;
      }
    }
    x=tmpx;
    nb_param=nb_col;
  }
  ifstream yfile(yfilename);
  if (!yfile) {ERRORQUIT("Unable to open y file");}
  Array2d<Real> y_mat(nb_eqn,nb_ycol);
  for (int eqn=0; eqn<nb_eqn; eqn++) {
    for (int col=0; col<nb_ycol; col++) {
      yfile >> y_mat(eqn,col);
    }
  }
  
  Array<Real> w(nb_eqn);
  one_array(&w);
  if (strlen(wfilename)>0) {
    ifstream wfile(wfilename);
    if (!wfile) {ERRORQUIT("Unable to open weight file");}
    for (int eqn=0; eqn<nb_eqn; eqn++) {
      wfile >> w(eqn);
    }
  }

  if (perfit && strlen(rfilename)==0) {ERRORQUIT("Specify -r=filename option in addition to -pf");}
  Array<Real> r(nb_param);
  if (strlen(rfilename)>0) {
    ifstream rfile(rfilename);
    if (!rfile) {ERRORQUIT("Unable to open regularization parameters file");}
    for (int param=0; param<nb_param; param++) {
      rfile >> r(param);
    }
  }

  for (int col=0; col<nb_ycol; col++) {
    Array<Real> y;
    extract_column(&y,y_mat,col);
    if (perfit) {
      Array<Real> b;
      calc_perfect_fit(&b, x,y,r);
      int j=0;
      for (int i=0; i<selectcol.get_size(); i++) {
	if (selectcol(i)==1) {
	  cout << b(j) << endl;
	  j++;
	}
	else {
	  cout << 0. << endl;
	}
      }
      /*
      Array<Real> tmp;
      product(&tmp,x,b);
      cout << endl << tmp << endl;
      */
    }
    else if (strlen(rfilename)>0) {
      if (docv) {
        cout << calc_cv_regul(x,y,w,del_colin,r) << endl;
      }
      else if (docopyy || dopred || doresid) {
	Array<Real> yhat;
	predict_ols_regul(&yhat,x,y,w,del_colin,r);
	for (int i=0; i<yhat.get_size(); i++) {
	  if (docopyy) {cout << y(i) << " ";}
	  if (dopred)  {cout << yhat(i) << " ";}
	  if (doresid) {cout << y(i)-yhat(i) << " ";}
	  cout << endl;
	}
      }
      else if (dostd) {
	ERRORQUIT("-r with -se not implemented");
      }
      else {
        Array<Real> b;
        calc_ols_regul(&b,x,y,w,del_colin,r);
	int j=0;
        for (int i=0; i<selectcol.get_size(); i++) {
	  if (selectcol(i)==1) {
	    cout << b(j) << endl;
	    j++;
	  }
	  else {
	    cout << 0. << endl;
	  }
        }
      }
    }
    else {
      if (docv) {
        cout << calc_cv(x,y,w,del_colin) << endl;
      }
      else if (docopyy || dopred || doresid) {
	Array<Real> yhat;
	predict_ols(&yhat,x,y,w,del_colin);
	for (int i=0; i<yhat.get_size(); i++) {
	  if (docopyy) {cout << y(i) << " ";}
	  if (dopred)  {cout << yhat(i) << " ";}
	  if (doresid) {cout << y(i)-yhat(i) << " ";}
	  cout << endl;
	}
      }
      else if (dostd) {
	if (strlen(wfilename)>0) {
	  ERRORQUIT("-w with -se not implemented");
	}
	Array<Real> b;
        calc_ols(&b,x,y,w,del_colin);
	Array2d<Real> var;
	calc_ols_var(&var,x,y);
	int j=0;
        for (int i=0; i<selectcol.get_size(); i++) {
	  if (selectcol(i)==1) {
	    cout << b(j) << "\t" << sqrt(var(j,j)) << endl;
	    j++;
	  }
	  else {
	    cout << 0. << endl;
	  }
        }
	
      }
      else {
        Array<Real> b;
        calc_ols(&b,x,y,w,del_colin);
	int j=0;
        for (int i=0; i<selectcol.get_size(); i++) {
	  if (selectcol(i)==1) {
	    cout << b(j) << endl;
	    j++;
	  }
	  else {
	    cout << 0. << endl;
	  }
        }
      }
    }
    if (nb_ycol>1 && !docv) {cout << endl;}
  }
}
