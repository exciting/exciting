void read_inequalities(LinkedList<LinearInequality> *ineq_list, const Array<AutoString> &label, istream &s) {
  while (!s.eof()) {
    char dir;
    while (1) {
      LinearInequality ineq(label.get_size()+1);
      zero_array(&(ineq.v));
      Real a;
      s >> a;
      skip_delim(s,"* \t\n");
      AutoString w;
      get_string(&w,s,"\t\n +-<>");
      if (!is_in_array(label,w)) {
        cerr << "unknown atom label " << w << endl;
        ERRORQUIT("Aborting.");
      }
      ineq.v(index_in_array(label,w))+=a;
      skip_delim(s," \t\n");
      dir=s.getc();
      if (strchr("<>",dir)>=0) {break;}
      s.pushback(dir);
    }
    skip_delim(s,"=");
    s >> ineq.c;
    if (dir==">") {
      product(&(ineq.v),ineq.v,-1.);
      ineq.c=-ineq.c;
    }
    (*ineq_list) << new LinearInequality(ineq);
    skip_delim(s,", \t\n");
  }
}

Real to_real(const AutoString &s) {
  istrstream str(s);
  Real r;
  str >> r;
  return r;
}

