#!/usr/bin/awk -f

# SYNTAX:
# str.out2exciting.awk str.out

# EXPLANATION:
# Using str.out, the crystal element for an exciting 
# input.xml file is produced and written to the 
# standard output.


#-
BEGIN{n=0;
    # Initialize empty matrices BB, PP_inv (to be used in functions)
    for (i=1;i<=3;i++){
        for (j=1;j<=3;j++)  {
            BB[i,j] = 0.0; PP_inv[i,j]=0.0
        }
    }
    # Initialize empty vector vv (to be used in functions)
    for(i=1;i<=3;i++){
            vv[i]=0
    }

    # define Bohr-radius  
    au=0.529177
    
    pi=atan2(1,1)*4

    # List of chemical elements
    AT[1]="H";  AT[2]="He"; AT[3]="Li"; AT[4]="Be"; AT[5]="B";
    AT[6]="C";  AT[7]="N";  AT[8]="O";  AT[9]="F";  AT[10]="Ne";
    AT[11]="Na";AT[12]="Mg";AT[13]="Al";AT[14]="Si";AT[15]="P";
    AT[16]="S"; AT[17]="Cl";AT[18]="Ar";AT[19]="K"; AT[20]="Ca";
    AT[21]="Sc";AT[22]="Ti";AT[23]="V"; AT[24]="Cr";AT[25]="Mn";
    AT[26]="Fe";AT[27]="Co";AT[28]="Ni";AT[29]="Cu";AT[30]="Zn";
    AT[31]="Ga";AT[32]="Ge";AT[33]="As";AT[34]="Se";AT[35]="Br";
    AT[36]="Kr";AT[37]="Rb";AT[38]="Sr";AT[39]="Y"; AT[40]="Zr";
    AT[41]="Nb";AT[42]="Mo";AT[43]="Tc";AT[44]="Ru";AT[45]="Rh";
    AT[46]="Pd";AT[47]="Ag";AT[48]="Cd";AT[49]="In";AT[50]="Sn";
    AT[51]="Sb";AT[52]="Te";AT[53]="I"; AT[54]="Xe";AT[55]="Cs";
    AT[56]="Ba";AT[57]="La";AT[58]="Ce";AT[59]="Pr";AT[60]="Nd";
    AT[61]="Pm";AT[62]="Sm";AT[63]="Eu";AT[64]="Gd";AT[65]="Tb";
    AT[66]="Dy";AT[67]="Ho";AT[68]="Er";AT[69]="Tm";AT[70]="Yb";
    AT[71]="Lu";AT[72]="Hf";AT[73]="Ta";AT[74]="W"; AT[75]="Re";
    AT[76]="Os";AT[77]="Ir";AT[78]="Pt";AT[79]="Au";AT[80]="Hg";
    AT[81]="Tl";AT[82]="Pb";AT[83]="Bi";AT[84]="Po";AT[85]="At";
    AT[86]="Rn";AT[87]="Fr";AT[88]="Ra";AT[89]="Ac";AT[90]="Th";
    AT[91]="Pa";AT[92]="U"; AT[93]="Np";AT[94]="Pu";AT[95]="Am";
    AT[96]="Cm";AT[97]="Bk";AT[98]="Cf";AT[99]="Es";AT[100]="Fm";
#    print "CRYSTAL";
#    print "PRIMVEC";
}

# Read basis vecs
NR==1 {a1=$1;a2=$2;a3=$3};
NR==2 {b1=$1;b2=$2;b3=$3};
NR==3 {c1=$1;c2=$2;c3=$3

};

# Read/Determine primitive lattice vectors
NR==4 {PP[1,1]=$1*a1+$2*b1+$3*c1;PP[1,2]=$1*a2+$2*b2+$3*c2;PP[1,3]=$1*a3+$2*b3+$3*c3};
NR==5 {PP[2,1]=$1*a1+$2*b1+$3*c1;PP[2,2]=$1*a2+$2*b2+$3*c2;PP[2,3]=$1*a3+$2*b3+$3*c3};
NR==6 {PP[3,1]=$1*a1+$2*b1+$3*c1;PP[3,2]=$1*a2+$2*b2+$3*c2;PP[3,3]=$1*a3+$2*b3+$3*c3};

# Read/Determine atomic positions
NF==4 {n++;
       for(i=1;i<=100;i++){if($4~AT[i]){atom[n]=i;atom_name[n]=AT[i]}};
       x[n]=$1; 
       y[n]=$2; 
       z[n]=$3; 
       xcart[n] = x[n]*a1+y[n]*b1+z[n]*c1
       ycart[n] = x[n]*a2+y[n]*b2+z[n]*c2
       zcart[n] = x[n]*a3+y[n]*b3+z[n]*c3
};


END{
  for(i=1;i<=3;i++){
    for(j=1;j<=3;j++){
  #      printf "  %10.5f ", PP[i,j]
    }
  #    printf "  \n"
  }
  transpose(PP, BB)
  inverse3x3(BB, PP_inv)  
  #  print "PRIMCOORD";
  #  print n, "   1";
  for(i=1;i<=n;i++){
      vec[1]=xcart[i]; vec[2]=ycart[i]; vec[3]=zcart[i];    
      mat_times_vec(PP_inv, vec, vv)
      ex_coord[i,1]=vv[1];
      ex_coord[i,2]=vv[2];
      ex_coord[i,3]=vv[3];
      }
  for(i=1;i<=3;i++){
    latvec1[i]=PP[1,i]
  }
  for(i=1;i<=3;i++){
    latvec2[i]=PP[2,i]
  }
  for(i=1;i<=3;i++){
    latvec3[i]=PP[3,i]
  }  
  
  #  printf "  latvec1: %10.7f %10.7f %10.7f\n", latvec1[1], latvec1[2], latvec1[3]
  #  printf "  latvec2: %10.7f %10.7f %10.7f\n", latvec2[1], latvec2[2], latvec2[3]
  #  printf "  latvec3: %10.7f %10.7f %10.7f\n", latvec3[1], latvec3[2], latvec3[3]


# Write basis vectors
#printf "  <structure speciespath=\"__SPECIES_PATH__\">\n"
printf "    <crystal>\n"
printf "      <basevect>%12.8f %12.8f %12.8f</basevect>\n", latvec1[1]/au,latvec1[2]/au,latvec1[3]/au
printf "      <basevect>%12.8f %12.8f %12.8f</basevect>\n", latvec2[1]/au,latvec2[2]/au,latvec2[3]/au
printf "      <basevect>%12.8f %12.8f %12.8f</basevect>\n", latvec3[1]/au,latvec3[2]/au,latvec3[3]/au
printf "    </crystal>\n"


# Write atomic positionss
  for(i=1;i<=n;i++){
    aname=atom_name[i]
    aname_count[aname]++
    coord[aname, aname_count[aname], 1] = ex_coord[i,1]
    coord[aname, aname_count[aname], 2] = ex_coord[i,2]
    coord[aname, aname_count[aname], 3] = ex_coord[i,3]

  }
  for (i in aname_count) {
#    print aname_count[i], i
    aname = i
    printf "    <species speciesfile=\"%s.xml\">\n", i
    for (j=1; j<= aname_count[i]; j++){
      printf "      <atom coord=\"%12.8f  %12.8f  %12.8f\" bfcmt=\"0.0  0.0  0.0\"></atom>\n",  coord[aname,j,1], coord[aname,j,2], coord[aname,j,3] 
    }      
    printf "    </species>\n"
  }
#  printf "  </structure>\n"
}


# ==============================================================================
# Function definitions
# ==============================================================================




# ----------------------------------------------------------
# Matrix functions
# ----------------------------------------------------------


function det3x3(AA){
  det = AA[1,2]*AA[2,3]*AA[3,1]-AA[1,3]*AA[2,2]*AA[3,1]+AA[1,3]*AA[2,1]*AA[3,2] -AA[1,1]*AA[2,3]*AA[3,2]+AA[1,1]*AA[2,2]*AA[3,3]-AA[1,2]*AA[2,1]*AA[3,3]
  return det
}

function det2x2(a){
  det = AA[1,1]*AA[2,2]-AA[1,2]*AA[2,1]
  return det
}


# ----------------------------------------------------------
# Returns BB, the inverse of the input 3x3 matrix  AA

function inverse3x3(AA,            BB){
  dd=det3x3(AA)
  BB[1,1]=(AA[2,2]*AA[3,3]-AA[2,3]*AA[3,2])/dd;
  BB[1,2]=(AA[1,3]*AA[3,2]-AA[1,2]*AA[3,3])/dd;
  BB[1,3]=(AA[1,2]*AA[2,3]-AA[1,3]*AA[2,2])/dd;
  BB[2,1]=(AA[2,3]*AA[3,1]-AA[2,1]*AA[3,3])/dd;
  BB[2,2]=(AA[1,1]*AA[3,3]-AA[1,3]*AA[3,1])/dd;
  BB[2,3]=(AA[1,3]*AA[2,1]-AA[1,1]*AA[2,3])/dd;
  BB[3,1]=(AA[2,1]*AA[3,2]-AA[2,2]*AA[3,1])/dd;
  BB[3,2]=(AA[1,2]*AA[3,1]-AA[1,1]*AA[3,2])/dd;
  BB[3,3]=(AA[1,1]*AA[2,2]-AA[1,2]*AA[2,1])/dd;
}

# ----------------------------------------------------------
# Returns vv, the product of mat times vec

function mat_times_vec(mat, vec,           vv){
   vv[1]=mat[1,1]*vec[1]+mat[1,2]*vec[2]+mat[1,3]*vec[3]
   vv[2]=mat[2,1]*vec[1]+mat[2,2]*vec[2]+mat[2,3]*vec[3]
   vv[3]=mat[3,1]*vec[1]+mat[3,2]*vec[2]+mat[3,3]*vec[3]
}


# ----------------------------------------------------------
# Returns BB, the transposed of the input matrix  AA

function transpose(AA,           BB){
  for(i=1;i<=3;i++){
    for(j=1;j<=3;j++){
      BB[j,i]=AA[i,j]
    }
  }
}





