#include "equil.h"

Array<Real> equildummyarray;

Equilibrator::Equilibrator(void): cur_sum(), cur_sum2(), bin_sum(), bin_sum2(), good_val(), buf_val() {}

Equilibrator::Equilibrator(Real _prec, int _which_elem, int init_granularity, int nb_bin):
  cur_sum(), cur_sum2(), bin_sum(), bin_sum2(), good_val(), buf_val() {
  init(_prec,_which_elem,init_granularity,nb_bin);
}

void Equilibrator::init(Real _prec, int _which_elem, int init_granularity, int nb_bin) {
  prec=_prec;
  which_elem=_which_elem;
  if (nb_bin%2==1) ERRORQUIT("nb_bin must be even");
  granularity=init_granularity;
  cur_corr_len=granularity/2;
  cur_bin=1;
  cur_cnt=0;
  cur_sum.resize(0); // init in new_data;
  bin_sum.resize(nb_bin+1);
  bin_sum2.resize(nb_bin+1);
  buf_corr.resize(nb_bin+1);
  buf_val.resize(granularity);
}

int Equilibrator::new_data(const Array<Real> &data) {
  if (cur_sum.get_size()==0) {
    int data_size=data.get_size();
    cur_sum.resize(data_size);
    cur_sum2.resize(data_size);
    zero_array(&cur_sum);
    zero_array(&cur_sum2);
    for (int i=0; i<bin_sum.get_size(); i++) {
      bin_sum(i).resize(data_size);
      bin_sum2(i).resize(data_size);
      zero_array(&bin_sum(i));
      zero_array(&bin_sum2(i));
    }
  }
  buf_val(cur_cnt)=data(which_elem);
  sum(&cur_sum,cur_sum,data);
  Array<Real> data2;
  product_diag(&data2,data,data);
  sum(&cur_sum2,cur_sum2,data2);
  cur_cnt++;
  if (cur_cnt==granularity) {
    bin_sum(cur_bin)=cur_sum;
    bin_sum2(cur_bin)=cur_sum2;
    
    Real corr;
    Real last_corr=0.5;
    int done=0;
    while (!done) {
      corr=0.;
      for (int i=0; i<granularity-cur_corr_len; i++) {
	corr+=buf_val(i)*buf_val(i+cur_corr_len);
      }
      corr/=(Real)(granularity-cur_corr_len);
      corr-=sqr((bin_sum(cur_bin)(which_elem)-bin_sum(cur_bin-1)(which_elem))/(Real)granularity);
      Real var=(bin_sum2(cur_bin)(which_elem)-bin_sum2(cur_bin-1)(which_elem))/(Real)granularity
	-sqr((bin_sum(cur_bin)(which_elem)-bin_sum(cur_bin-1)(which_elem))/(Real)granularity);
      if (var>0) {
	corr/=var;
      }
      else {
	corr=1.;
      }
      done=1;
      if (corr<0.) corr=0.;
      if (corr<0.25 && last_corr<=0.75 && cur_corr_len>=2) {cur_corr_len/=2; done=0;}
      if (corr>0.75 && last_corr>=0.25 && cur_corr_len<granularity/4) {cur_corr_len*=2; done=0;}
      last_corr=corr;
    }
    corr=pow(corr,1./(Real)cur_corr_len);
    buf_corr(cur_bin)=corr;
    
    int b;
    for (b=cur_bin/2; b>=1; b--) {
      Real block1=(bin_sum(cur_bin  )(which_elem)-bin_sum(cur_bin-  b)(which_elem))/(Real)(granularity*b);
      Real block2=(bin_sum(cur_bin-b)(which_elem)-bin_sum(cur_bin-2*b)(which_elem))/(Real)(granularity*b);
      if (fabs(block1-block2)<prec) {
	int bin1=cur_bin-2*b;
	Real n=(Real)(granularity*2*b);
	Real avg_corr=0.;
	for (int i=bin1+1; i<=cur_bin; i++) {
	  avg_corr+=corr;
	}
	avg_corr/=(Real)(cur_bin-bin1);
	Real var=((bin_sum2(cur_bin)(which_elem)-bin_sum2(bin1)(which_elem))/n-sqr((bin_sum(cur_bin)(which_elem)-bin_sum(bin1)(which_elem))/n))/n;
	if ((1.-avg_corr)!=0.) {
	  var*=(1.+avg_corr)/(1.-avg_corr);
	  if (var<sqr(prec) && var>0) {
	    diff(&good_val,bin_sum(cur_bin),bin_sum(bin1));
	    product(&good_val,good_val,1./n);
	    diff(&good_val2,bin_sum2(cur_bin),bin_sum2(bin1));
	    product(&good_val2,good_val2,1./n);
	    Array<Real> tmp;
	    product_diag(&tmp,good_val,good_val);
	    diff(&good_val2,good_val2,tmp);
	    equil=bin1*granularity;
	    step=(cur_bin-bin1)*granularity;
	    break;
	  }
	}
      }
    }
    
    cur_cnt=0;
    cur_bin++;
    if (cur_bin==bin_sum.get_size()) {
      int i,j;
      for (i=1,j=2; j<bin_sum.get_size(); i++,j+=2) {
	bin_sum(i)=bin_sum(j);
	bin_sum2(i)=bin_sum2(j);
	buf_corr(i)=(buf_corr(j-1)+buf_corr(j))/2.;
      }
      for (; i<bin_sum.get_size(); i++) {
	zero_array(&bin_sum(i));
	zero_array(&bin_sum2(i));
      }
      cur_bin=(bin_sum.get_size()-1)/2+1;
      granularity*=2;
      buf_val.resize(granularity);
    }
  }
  return (good_val.get_size()==0);
}

int Accumulator::new_data(const Array<Real> &data) {
  if (cur_sum.get_size()==0) {
    cur_sum.resize(data.get_size());
    zero_array(&cur_sum);
    cur_sum2.resize(data.get_size());
    zero_array(&cur_sum2);
    cur_n=0;
  }
  sum(&cur_sum,cur_sum,data);
  Array<Real> data2;
  product_diag(&data2,data,data);
  sum(&cur_sum2,cur_sum2,data2);
  cur_n++;
  return (cur_n<max_n);
}

void Accumulator::accum(void) {
  product(&good_val,cur_sum,1./(Real)cur_n);
  product(&good_val2,cur_sum2,1./(Real)cur_n);
  Array<Real> tmp;
  product_diag(&tmp,good_val,good_val);
  diff(&good_val2,good_val2,tmp);
}

GenericAccumulator *create_accum(int n_step, Real prec, int which_elem) {
  if (prec==0) {
    return new Accumulator(n_step);
  }
  else {
    return new Equilibrator(prec,which_elem);
  }
}
