#!/usr/bin/python

import os
from pylab import *
import sys
import xml.dom.minidom as dom

import matplotlib.pyplot as plt


Arg=sys.argv
print(Arg[0])
Arg=Arg[1:]

params = {'xtick.major.size': 8,
          'ytick.major.size': 8,
          'xtick.major.width': 2,
          'ytick.major.width': 2,
          'patch.linewidth': 1.5,
          'axes.linewidth': 2.,
          'axes.formatter.limits': (-4, 6),
          'lines.linewidth': 1.0,
          'lines.markeredgewidth':1.3,
          'lines.markersize':18,
          #'font.weight':'bold',
          'font.size':14,
          'font.family':'sans',
          'axes.labelsize':17,
          'legend.fontsize':14,
          'legend.borderaxespad':1,
          'legend.borderpad':0.5}

plt.rcParams.update(params)

#Arg=[ "pattern[S,PF,ZT]","T=300,t=1.0"]
#List of files
List=os.listdir(str(os.getcwd()))
######
if(len(Arg)<2)or (len(Arg))>3:
    print("Wrong number of Input Arguments")
    sys.exit()


#Exceptions
#Check the Plotstyle Argument
def checkPlotstyle(k):
    symbols=["ZT","sigma","S","PF","kappa"]
    if (k.count("[")==0 or k.count("[")>1):
        print("Error in [ : "+k)
        return False
    elif (k.count("]")==0 or k.count("]")>1):
        print("Error in ] : "+k)
        return False
    k=k.replace("]","")
    u=k.split("[")
    for l in u[1].split(","):
        if l not in symbols:
            print("Error-Wrong Symbol input:"+l)
            return False
    return True

#Check the value Argument
def Checkvalue(k):
    if ((k.count("=")==0) or (k.count("=")>2)):
        print("Error-Wrong input:"+k)
        return False
    elif (0<k.count("T") and 0<k.count("mu")):
        print("Error- mu or T not both:  "+k)
        return False
    k=k.replace("K","")
    k=k.replace("s","")
    k=k.replace("eV","")
    u=k.split(",")
    for l in u:
        try:
            float(l.split("=")[1])
        except:
            print("Error- Check value: "+l)
            return False
    return True

for k in Arg:
    if k.count("pattern")>0 or k.count("single")>0:
        sys.exit() if not(checkPlotstyle(k))else None
    elif k.count("=")>0:
        sys.exit() if not(Checkvalue(k))else None
    else:
        print("Error-wrong input: "+k)
        sys.exit()
########


#xml read
try:
    baum=dom.parse(str(os.getcwd())+"/"+"input.xml")
    for k in baum.firstChild.childNodes:
        if k.nodeName=="properties":
            for l in k.childNodes:
                if l.nodeName=="boltzequ":
                    try:
                        if str(l.getAttribute("tsiout"))=="true":
                            Arg.append("SI")
                        elif str(l.getAttribute("tsiout"))=="false":
                            Arg.append("Au")
                    except:
                        print("Error-Check input.xml for tsiout=true or false")
                        sys.exit()
except:
    print("Error: input.xml not found")
    sys.exit()
########


#Arguments

w,u="",""#pattern or single
for k in Arg:
    if (k.count("single")>0):
        w=(k.replace("]",""))
    elif(k.count("pattern")>0):
        u=(k.replace("]",""))

tau="t=1.0"#Tau argument
for k in Arg:
    if (k.count("=")>0):
        for i in k.split(","):
            if i.count("T")>0 :
                cu=i.replace("K","")
            elif i.count("mu")>0:
                cu=i.replace("eV","")
            elif i.count("t")>0:
                tau=i.replace("s","")
        break
##########
#Translate from 12->xy...
Per_xyz=[str(a)+str(b) for b in ["x","y","z"] for a in ["x","y","z"]]
Per_123=[str(a)+str(b) for b in [1,2,3] for a in [1,2,3]]
R={Per_123[k]:Per_xyz[k] for k in range(len(Per_xyz))}
#Fermienergie
for line in open(str(os.getcwd())+"/"+"EFERMI.OUT"):
    if "SI" in Arg:
        Ef=float(line)*27.21138
    elif "Au" in Arg:
        Ef=float(line)



#List of files
List=os.listdir(str(os.getcwd()))
#Data Input function
def Read_data(file,l):
    file=str(os.getcwd())+"/"+file
    X,Y=[],[]
    value=""
    for line in open(file):
        line=line.split()
        check=False#if a X-value set then check=True else False
        value=line[0] if(cu.split("=")[0]=="T") else line[1]
        if (cu.split("=")[0]=="T") and abs(float(cu.split("=")[1])-float(line[0]))<1e-2:
            X.append(float(line[1])-Ef)
            check=True
        elif (cu.split("=")[0]=="mu") and abs(float(cu.split("=")[1])-(float(line[1])))<1e-3:
            X.append(float(line[0]))
            check=True
        if check:
            if (l=="ZT" or l=="S"):
                if l=="S":
                    Y.append(float(line[2])*1e6)
                else:
                    Y.append(float(line[2]))
            else:
                if l=="sigma" or l=="kappa":
                    Y.append(float(line[2])*float(tau.split("=")[1])*1e-2)
                else:
                    Y.append(float(line[2])*float(tau.split("=")[1]))

    if len(X)==0 and len(Y)==0:
        k=0 if (cu.split("=")[0]=="T") else 1
        print("Error- Check data in file :"+file+" or Check "+cu.split("=")[0]+"=("+cu.split("=")[1]+" - "+str(value)+")"+"<1e-3")
        sys.exit()
    return [X,Y]

#Xlabel
#Xl={"SI":{"T":"$\mathbf{\mu\ [eV]}$","mu":"$\mathbf{T [K]}$"},"Au":{"T":"$\mathbf{\mu$ $[a.u.]}$","$\mu$":"$\mathbf{T\ [a.u.]}$"}}
Xl={"SI":{"T":"$\mu\ [eV]$","mu":"$T [K]$"},"Au":{"T":"$\mu$ $[a.u.]$","$\mu$":"$T\ [a.u.]$"}}
labx=Xl["SI" if ("SI" in Arg) else "Au"][cu.split("=")[0]]
#Ylabels
if tau=="t=1.0":
    tf=False
else:
    tf=True
#K={"SEEBECK_":{True:r"$\mathbf{S\ [\mu V/K]}$",False:r"$\mathbf{S\ [\mu V/K]}$"},"ZT_":{True:r"$\mathbf{ZT}$",False:r"$\mathbf{ZT}$"},
#   "ELECTCOND_":{False:r"$\mathbf{\sigma/\tau}$ $\mathbf{[(\Omega cm\ s)^{-1}]}$",True:r"$\mathbf{\sigma$ $[(\Omega cm)^{-1}]}$"},
#   "THERMALCOND_":{False:r"$\mathbf{\kappa/\tau}$ $\mathbf{[W (cm K\ s)^{-1}]}$",True:r"$\mathbf{\kappa$ $[W (K cm\ )^{-1}]}$"},
#   "PF_":{False:r"$\mathbf{\sigma S^2/\tau}$ $\mathbf{[\mu W (cm K^2\ s)^{-1}]}$",True:r"$\mathbf{\sigma S^2$ $[\mu W/cm$ $K^2]}$"}}
K={"SEEBECK_":{True:r"${S\ [\mu V/K]$",False:r"$S\ [\mu V/K]$"},"ZT_":{True:r"$ZT$",False:r"$ZT$"},
   "ELECTCOND_":{False:r"$\sigma/\tau$ $[(\Omega cm\ s)^{-1}]$",True:r"$\sigma$ $[(\Omega cm)^{-1}]$"},
   "THERMALCOND_":{False:r"$\kappa/\tau$ $[W (cm K\ s)^{-1}]$",True:r"$\kappa$ $[W (K cm\ )^{-1}]$"},
   "PF_":{False:r"$\sigma S^2/\tau$ $[\mu W (cm K^2\ s)^{-1}]$",True:r"$\sigma S^2$ $[\mu W/cm$ $K^2]$"}}


K2={"kappa":"THERMALCOND_","sigma":"ELECTCOND_", "ZT":"ZT_","S":"SEEBECK_","PF":"PF_"}

#XYlabel font
#font = {'family': 'sans', 
#       'weight': 'bold',
#        'size': 15,
#        }
#######
#Plotroutine for pattern
u=u.split("[")
if "pattern" in u :
        #plots.reverse()
        plots=u[1].split(",")
        f,ax=plt.subplots(len(plots),sharex=True,figsize=(6.5,7))
        n=0
        for l in plots:
            plots=u[1].split(",")
            if l=="PF":#PF Factor
                for k in List:
                    if k.count("ELECTCOND_")>0:
                        for j in List:
                            if j.count("SEEBECK_")>0:
                                if j.split("_")[1]==k.split("_")[1]:
                                    Sigma=Read_data(k,"PF")[1]
                                    [X,Y]=Read_data(j,"S")
                                    try:
                                        PF=[(Sigma[j]*Y[j]**2)*1e-8 for j in range(len(X))]
                                    except:
                                        print("Error- Check data in file :"+j+"or "+k)
                                        sys.exit()
                                    ax[n].plot(X,PF,label=R[k.split("_")[1].replace(".OUT","")],linewidth=1.7,)
                                    ax[n].legend()
                                    if "Au" in Arg:
                                        ax[n].set_ylabel(K[K2[l]][tf].split("in")[0])#, fontdict=font)
                                    else:
                                        ax[n].set_ylabel((K[K2[l]][tf]))#, fontdict=font)
            else:
                for k in List:
                    if k.count(K2[l])>0:
                        [X,Y]=Read_data(k,l)
                        try:
                            if l=="sigma" or l=="kappa":
                                ax[n].plot(X,Y,label=R[k.replace(K2[l],"").replace(".OUT","")],linewidth=1.7,)
                                ax[n].set_ylim(0)
                            else:
                                ax[n].plot(X,Y,label=R[k.replace(K2[l],"").replace(".OUT","")],linewidth=1.7,)

                        except:
                            print("Error- Check data in file :"+cu+"or "+k)
                            sys.exit()

                        if "Au" in Arg:
                            ax[n].set_ylabel(K[K2[l]][tf].split("in")[0]) #, fontdict=font)
                        else:
                            ax[n].set_ylabel(K[K2[l]][tf]) #, fontdict=font)
            
            xlabel(labx) #, fontdict=font)
            ax[n].legend()#prop=font)
            
            n+=1
        f.subplots_adjust(hspace=0.08)
        f.set_figheight(9)
        f.set_figwidth(7)
        f.tight_layout()
                    
        savefig(str(os.getcwd())+"/"+"pattern_Plot.png",orientation='portrait',format='png',dpi=100)
#######
#For single plots
w=w.split("[")
if "single" in w :
    for l in w[1].split(","):
        plt.figure(figsize=(6.5,5.3))
        if l=="PF":#PF Factor
                for k in List:
                    if k.count("ELECTCOND_")>0:
                        for j in List:
                            if j.count("SEEBECK_")>0:
                                if j.split("_")[1]==k.split("_")[1]:
                                    Sigma=Read_data(k,"PF")[1]
                                    [X,Y]=Read_data(j,"S")
                                    try:
                                        PF=[(Sigma[j]*Y[j]**2)*1e-8 for j in range(len(X))]
                                    except:
                                        print("Error- Check data in file :"+k)
                                        print("or "+cu+","+tau)
                                        sys.exit()
                                    plot(X,PF,label=R[k.split("_")[1].replace(".OUT","")],linewidth=1.7)
                                    legend()
                                    if "Au" in Arg:
                                        ylabel(K[K2[l]][tf].split("in")[0])#, fontdict=font)
                                    else:
                                        ylabel((K[K2[l]][tf]))#, fontdict=font)

        else:
            for k in List:
                if k.count(K2[l])>0:
                    [X,Y]=Read_data(k,l)
                    try:
                        if l=="sigma" or l=="kappa":
                            plot(X,Y,label=R[k.replace(K2[l],"").replace(".OUT","")],linewidth=1.7,)
                            ylim(0)
                        else:
                            plot(X,Y,label=R[k.replace(K2[l],"").replace(".OUT","")],linewidth=1.7,)
                    except:
                        print("Error- Check data in file :"+k)
                        print("or "+cu+","+tau)
                        sys.exit()
                    if "Au" in Arg:
                        ylabel(K[K2[l]][tf].split("in")[0])
                    else:
                        ylabel(K[K2[l]][tf])#, fontdict=font)
        xlabel(labx)#,fontdict=font)
        legend()#prop=font)
        plt.tight_layout()

        l1=K2[l]
        l1=l1.lower().replace("_","")
        savefig(str(os.getcwd())+"/"+l1+"_Plot"+".png" ,orientation='portrait',format='png',dpi=100)
############
show()
