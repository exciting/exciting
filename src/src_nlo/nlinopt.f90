!BOP
! !ROUTINE: nlinopt
! !INTERFACE:
      subroutine nlinopt
! !USES             : Module optica
! !INPUT PARAMETERS :
! nstfv             : maximum number of bands for allocation [integer]
! banmin            : minimum band from which optical transitions start [integer]
! banmax            : maximum band to which optical transitions are allowed [integer]
! idel1             : broadening in Ry units [real]
! dw                : energy step [real]
! emesh             : number of energy mesh points [integer]
! sc                : scissors operator in Ry units (self energy) [real]
! tol               : Tolerence (to avoid singularities) [real]
! v1                : component of $\chi^{{\bf v1}v2v3}(2\omega,\omega,\omega)$ [integer]
! v2                : component of $\chi^{v1{\bf v2}v3}(2\omega,\omega,\omega)$ [integer]
! v3                : component of $\chi^{v1v2{\bf v3}}(2\omega,\omega,\omega)$ [integer]
! omega             : volume of the unit cell in a.u. (from case.scf) [real]
! nkpt              : maximum number of $k$-points [integer]
! evalfv(:,:)       : eigenvalue (band,kpoint) [real]
! wkpt(:)           : weight of the k-pt [real]
! pmat(:,:,:)       : opticalmatirix element (kpoint,frm-kpt,to-kpt,component) [complex]
! noval(:)          : number of valence bands (kpoint) [integer]
! nocond(:)         : number of comduction bands(kpoint) [integer]
! tot_ban(:)        : total number of bands (kpoint) [integer]
!
! !OUTPUT PARAMETERS:
! inter1w(:)        : 1$\omega$ interband contribution to the SHG $\chi$ in units esu and pm/V [complex]
! inter2w(:)        : 2$\omega$ interband contribution to the SCG $\chi$ in units esu and pm/V [complex]
! intra1w(:)        : 1$\omega$ intraband contribution to the SHG $\chi$ in units esu and pm/V [complex]
! intra2w(:)        : 2$\omega$ intraband contribution to the SHG $\chi$ in units esu and pm/V [complex]
! !DESCRIPTION:
! To calculates the second-order optical response of semiconductors 
!
! \begin{eqnarry}
! $\chi_{inter}^{abc}(2\omega,\omega,\omega) =
! {{e^3} \over {\hbar^2\Omega}} \sum'_{nml}\sum_{\bf k}
! \left\{
! {{2r^a_{nm}\{r^b_{ml}r^c_{ln}\}} \over{(\omega_{ln}-\omega_{ml})(\omega_{mn}-2\omega)}}
! \right. \nonumber \\
! \left.
! -{1 \over (\omega_{mn}-\omega)}
! \left[
! {{r^c_{lm}\{r^a_{mn}r^b_{nl}\}} \over{(\omega_{nl}-\omega_{mn})}}-
! {{r^b_{nl}\{r^c_{lm}r^a_{mn}\}} \over{(\omega_{lm}-\omega_{mn})}}
! \right]
! \right\}
! \end{equnarry}
!
! \begin{eqnarry}
! \chi_{intra}^{abc}(2\omega,\omega,\omega) =
! {{e^3} \over {\hbar^2\Omega}} \sum_{{\bf k}}
! \left\{
! \sum'_{nml} {1 \over {\omega^2_{mn}(\omega_{mn}-\omega)}}
! \left[
! \omega_{ln}r^b_{nl}\{ r^c_{lm}r^a_{mn}\} -\omega_{ml}r^c_{lm}\{ r^a_{mn}r^b_{nl}\}
! \right]
! \right. \nonumber \\
! \left.
! -8i\sum'_{nm} {1 \over{\omega^2_{mn}(\omega_{mn}-2\omega)}}r^a_{nm}\{ r^b_{ml}r^c_{ln}\}
! +2\sum'_{nml}{{r^a_{nm}\{ r^b_{ml}r^c_{ln}\}(\omega_{ml}-\omega_{ln})} \over {\omega^2_{mn}(\omega_{mn}-2\omega)}}
! \right\}
! \end{eqnarry}
!
! \begin{eqnarry}
! \chi_{mod}^{abc}(2\omega,\omega,\omega) = {{e^3} \over {2\hbar^2\Omega}} \sum_{{\bf k}}
! \left\{
! \sum_{nml} {1 \over {\omega^2_{mn}(\omega_{mn}-\omega)}}
! \left[
! \omega_{nl}r^a_{lm}\{r^b_{mn}r^c_{nl}\}-\omega_{lm}r^a_{nl}\{r^b_{lm}r^c_{mn}\}
! \right]
! \right. \nonumber \\
! \left.
! -i\sum_{nm} {{r^a_{nm}\{r^b_{mn}\Delta^c_{mn}\}} \over {\omega^2_{mn}(\omega_{mn}-\omega)}}
! \right\}
! \end{eqnarry}
! 
! !FILES USED:
! case.innlo : The master input file
! case.weight: contains weights of the kpts and eigenvalues.
! case.mommat: Momentum matrix elements in desired format
! case.symop : symmetry operation
! nlo.def    : names of the input/output files.
!
! !REVISION HISTORY:
! Updated and documented: 14 Nov. 2002
!
! !REMARKS:
! The factors are named accordong to the paper PRB 53 #16 10 751 1996-II, Apendix B
! I recommend paper PRB 48 #16 1993-II for the detailed formulism 
!EOP
!BOC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program: Non-linear optical suceptibility for bulk semiconductors.
! and semiconductor superlattices
!
! Author: Sangeeta Sharma
! Date of last update: Sept 23 2002
!
! !!!!!!!!!!!!! The terms B3 from PRB53 10751 have been implemented !!!!!
! !!!!!!!!!!!!! sym_mat is opimat from optic !!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!!!!!!!!!!In order to calculate in a different frame change !!!!!!!!!!
! !!!!!!!!!! opimat to imat, eg plot hex in cubic ref frame !!!!!!!!!!!!!
! !!!!!!!! Compatible with changes made to optic on 02 may 2002 !!!!!!!!!
!
! calculates the non-linear susceptibility chi(-2w,w,w) total and then 
! gives real,imaginary parts and abs value of it for plotting.
! why ksum? consult documentation for NLO
!
! using k-sum (not linear-tet!). Chi is divided into the interband
! and intraband contributions.
!   
! Terms:                 
! bxyz: equation Bx from paper PRB53, term: y, z is just the index telling which
!                                             band combination eg: 2=n <l< m
!                                                                : 3=l>m
! r(k,n,m,x): r MEs = p_nm/(i*w_nm*m)
! SYM        : symmetry tensor
! evalfv(n,s): eigenvalue for nth band at sth k-point
! wnm        : evalfv(n,k)-evalfv(m,k)
! px(n,m,x,y,z) : symmetrised r_nm wrt y and z
! py(n,m,x,y,z) : symmetrized r_nm wrt x and z
!..................................................................................................

      use modmain
      use optica
      implicit none

      real(8) wmn,wnm,wln,wnl,wml,wlm,en,em,el,e1,e2,e12,e21
      real(8) wn,wm,wl                          
      real(8) const,const_au,au2esu
      real(8) totre,totabs,e            
      real(8) emin,emax,corec1,corec2
      real(8), allocatable:: totim(:)                                      

      complex(8) b111,b121,b131,b112,b122,b132,b113,b123,b133  
      complex(8) b241,b242,b243,b221,b222,b223,b211,b212,b213,b231
      complex(8) b311,b312,b313,b331
      complex(8) b24,b21_22,b11,b12_13,b31_32                             
      complex(8) idel,w,ii
      complex(8) mat2w,mat1w1,mat1w2,mat2w_tra,mat1w3_tra 
      complex(8), allocatable :: px(:,:,:,:,:),py(:,:,:,:,:),         &    
                                 pz(:,:,:,:,:),delta(:,:,:)
      complex(8), allocatable :: inter2w(:),inter1w(:),intra2w(:),    &
                                 intra1w(:)  
      complex(8), allocatable :: pmat(:,:,:)

      integer iw, recl
      integer i,j,k,il,ik,lx,ly,lz,isym
      integer istl,istn,istm,ist1,ist2,ist
!
!.............................................................................

      allocate (px(nstfv,nstfv,3,3,3),py(nstfv,nstfv,3,3,3),     &
                pz(nstfv,nstfv,3,3,3))
      allocate (delta(3,nstfv,nstfv))
      allocate (SYM(3,3,3))
      allocate (totim(emesh))
      allocate (r(3,nstfv,nstfv,nkpt)) 
      allocate (inter2w(emesh),inter1w(emesh),intra2w(emesh),    &
                intra1w(emesh))
!..............................................................................
!bernd  const=(-9.4013754040482201E-10)*4.d0*(pi**3)*(1.d0/omega)
!sas    const=-32.d0*(1/omega)*4.6393930614E-9*(0.5d0)*(0.5d0)

! 32          : converting all energies in denominator(E^5) from Ry to au
! -ve         : sign for electron charge
! omega       : unit cell volume
! 0.5         : dividing by the sum of the weights which is 2

      const_au=-32.d0*(1/omega)*(0.5d0)

! 4.63939E-9  : au2esu : bohr*c*10^4/4pi*2*ry2ev
!               bohr: 5.29177E-11
!               c: 2.99792458   velocity of sound
!               ry2ev: 13.60569172

      au2esu=5.8300348177d-8

!    au2esu=(5.29177E-11*2.99792458*1.0E4)/(13.60569172*2)
! this const includes (e^3*hbar^3*hbar^3)/(omega*hbar^5*m_e^3)
! mass comes from converting P_mn to r_mn
! hbar^3 comes from converting all frequencies to energies in denominator
! hbar^3 comes from operator for momentum (hbar/i nabla) 
! the ME from optic doesn't include this hbar^^^
 
      const=const_au*(au2esu)                      ! in esu
!......................................................................................
      ii=(0.d0,1.d0)
      r(:,:,:,:)=0.d0
      sym(:,:,:)=0.d0
      px(:,:,:,:,:)=0.d0
      py(:,:,:,:,:)=0.d0
      pz(:,:,:,:,:)=0.d0
      delta(:,:,:)=0.d0

      mat2w=0.d0
      mat1w1=0.d0
      mat1w2=0.d0
      mat2w_tra=0.d0
      mat1w3_tra=0.d0

      inter2w(:)=0.d0
      inter1w(:)=0.d0
      intra2w(:)=0.d0
      intra1w(:)=0.d0

      idel=(0.d0,1.d0)*idel1

! .................. calculating the energy window for calculations...................     

      emin=0.d0
      emax=0.d0

      do ik=1,nkpt
        do ist=banmin,banmax
          emin=min(emin,evalfv(ist,ik))
          emax=max(emax,evalfv(ist,ik))
        end do
      end do

! ............. Making the symmetrizing tensor 'sym'..............................
      const=const*(1.d0/dble(nsymcrys))        
      do isym=1,nsymcrys
        do i=1,3
          do j=1,3
            do k=1,3
              il=lsplsymc(isym)
              sym(i,j,k)=sym(i,j,k)+ &
            &   (symlatc(i,v1,il)*symlatc(j,v2,il)*symlatc(k,v3,il))
            end do
          end do
        end do
      end do

!-----------------------------------------------------------------------------
!................ Starting the k loop ..............................................
!-----------------------------------------------------------------------------

!     Momentum matrix elements
      allocate(pmat(3,nstsv,nstsv))
      Inquire(IoLength=recl) pmat
      Open (50, File='PMAT.OUT', Action='READ', Form='UNFORMATTED', &
     &   Access='DIRECT', recl=recl)

      do ik = 1, nkpt
      
! ..... read momentum matrix elements from direct-access file ............
        read(50,rec=ik) pmat
         
        do ist1=banmin,banmax                                            
          do ist2=banmin,banmax                                          
            delta(1:3,ist1,ist2)=pmat(1:3,ist1,ist1)-pmat(1:3,ist2,ist2)
!................ Extrating the lower diagonal elements ............................
            if (ist1.lt.ist2) then
              r(1:3,ist1,ist2,ik)=pmat(1:3,ist1,ist2)
            else if (ist2.lt.ist1) then
              r(1:3,ist1,ist2,ik)=conjg(pmat(1:3,ist2,ist1))
            end if
          end do
        end do

!................... Correcting for scissors correction ............................
! check eqn # 20 to 23 of paper PRB53 10751 for this.
!
        do ist1=banmin,noval(ik)
          e1=evalfv(ist1,ik)
          do ist2=noval(ik)+1,banmax
            e2=evalfv(ist2,ik)

            e12=e1-e2
            e21=-e12
            corec1=(e12-sc)/e12
            corec2=(e21+sc)/e21

            r(1,ist1,ist2,ik)=r(1,ist1,ist2,ik)*corec1
            r(2,ist1,ist2,ik)=r(2,ist1,ist2,ik)*corec1
            r(3,ist1,ist2,ik)=r(3,ist1,ist2,ik)*corec1

            r(1,ist2,ist1,ik)=r(1,ist2,ist1,ik)*corec2
            r(2,ist2,ist1,ik)=r(2,ist2,ist1,ik)*corec2
            r(3,ist2,ist1,ik)=r(3,ist2,ist1,ik)*corec2

          end do
        end do

!................ Symmetrizing the r_nm matrix elements ...........................................
       
        do ist1=banmin,noval(ik)
          e1=evalfv(ist1,ik)
          do ist2=noval(ik)+1,banmax
            e2=evalfv(ist2,ik)
            do lx=1,3
              do ly=1,3
                do lz=1,3
                  px(ist1,ist2,lx,ly,lz)=              & ! symmetrising b c 
               &      (sym(lx,ly,lz)+sym(lx,lz,ly))*r(lx,ist1,ist2,ik)

                  py(ist2,ist1,lx,ly,lz)=              & ! symmetrising a c
               &      (sym(ly,lx,lz)+sym(ly,lz,lx))*r(lx,ist2,ist1,ik)

                  pz(ist2,ist1,lx,ly,lz)=              & ! symmetrising a b 
               &      (sym(lz,lx,ly)+sym(lz,ly,lx))*r(lx,ist2,ist1,ik)
                    
                end do ! lz
              end do ! ly
            end do ! lx
          end do
        end do

!ss.................... Start evaluating the factors ................................................
!ss The factors are named accordong to the paper PRB 53 #16 10 751 1996-II, Apendix B
!ss I recommend paper PRB 48 #16 1993-II for the detailed formulism 
!ss .................................................................................................
        b111=0.d0
        b121=0.d0
        b131=0.d0
        b112=0.d0
        b122=0.d0
        b132=0.d0
        b113=0.d0
        b123=0.d0
        b133=0.d0

        b211=0.d0
        b221=0.d0
        b212=0.d0
        b222=0.d0
        b213=0.d0
        b223=0.d0

        b231=0.d0
  
        b241=0.d0
        b242=0.d0
        b243=0.d0 
     
        b311=0.d0
        b312=0.d0
        b313=0.d0

        b331=0.d0
!ss..................start the calculation:.................................

       do istn=banmin,noval(ik)                   !ss the state n (valence)
        en=evalfv(istn,ik)
        if(abs(en).gt.1.d-5) then                 !ss to ensure band passes this k-pt

         do istm=noval(ik)+1,banmax               !ss the state m (conduction)
          em=evalfv(istm,ik)
          if(abs(em).gt.1.d-5) then

            wn=en
            wm=em+sc

            wmn=wm-wn
            wnm=-wmn

            b11=0.d0
            b12_13=0.d0
            b24=0.d0
            b31_32=0.d0
            b21_22=0.d0

            mat2w_tra=0.d0
            mat1w3_tra=0.d0

            do lx=1,3
              do ly=1,3
                do lz=1,3
                  mat2w_tra=mat2w_tra+px(istn,istm,lx,ly,lz)*r(lz,istm,istn,ik)    &
                              *delta(ly,istm,istn)

                  mat1w3_tra=mat1w3_tra+px(istn,istm,lx,ly,lz)*r(ly,istm,istn,ik)  &
                              *delta(lz,istm,istn)   

!!!!! NOTE:: lx to ly m to n in r's respectivly
! Changes are made so that this(b3) term is according to paper
! PRB48 rather than PRB53 which is incorrect!!!!

                 end do
               end do
             end do

! ............. The intraband 2w and w two band contribution....................................... 

             b231=8.d0*mat2w_tra/wmn

             b331=mat1w3_tra/(wnm)

! chi_inter^abc = sum_k sum_nm sum_l [ (r^a_nm{r^b_ml.r^c_ln}/(wln-wml))*(2/(wmn-2w)) ]............
! But remember r_nm = p_nm/(i*wnm)

!.......... l < valence band (istn)  ...The B(4) term of PRB48.....................................
           
            do istl=banmin,istn-1               !ss the valence state below n
            el=evalfv(istl,ik)
            if(abs(el).gt.1.d-7) then
              
              wl=el

              wln=wl-wn
              wml=wm-wl
              wnl=-wln
              wlm=-wml
              
              mat2w=0.d0
              mat1w1=0.d0
              mat1w2=0.d0
              mat2w_tra=0.d0

              do lx=1,3
                do ly=1,3
                  do lz=1,3
                    mat2w=mat2w+(px(istn,istm,lx,ly,lz)*r(ly,istm,istl,ik)   &
                                *r(lz,istl,istn,ik))
                               
                    mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*r(lz,istl,istm,ik) &
                                *r(ly,istn,istl,ik))
                                
                    mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*r(lz,istl,istm,ik) &
                                *r(ly,istn,istl,ik))

                   end do
                 end do
               end do
              b111=mat2w*(1.d0/(wln+wlm))*(1.d0/wlm)
              b121=mat1w1*(1.d0/(wnm+wlm))*(1.d0/wlm)
              b131=mat1w2*(1.d0/wlm)
              
              b221=0.d0
              b211=mat1w1/wml

              b241=-mat2w/wml

              b311=mat1w2/wlm

              if(dabs(wln).gt.tol) then
                b111=b111/wln
                b121=b121/wln
                b131=b131/wln

                b221=mat1w2/wln

                b241=b241+(mat2w/wln)

                b311=b311-(mat1w1/wnl)
              else
                b111=0.d0
                b121=0.d0
                b131=0.d0
    
                b221=0.d0
              end if
              

              if(dabs(wln-wnm).gt.tol) then
                b131=b131/(wln-wnm)
              else
                b131=0.d0              
              end if
              
              b11=b11-2.d0*b111        
              b12_13=b12_13+b121+b131

              b21_22=-b211+b221  

              b24=b24+2.d0*b241

              b31_32=b31_32+b311

            end if
            end do
            
!................. val(istn) < l < cond(istm)......................................................

            do istl=istn+1,istm-1
            el=evalfv(istl,ik)
            if (abs(el).gt.1.d-7) then

              mat2w=0.d0
              mat1w1=0.d0
              mat1w2=0.d0

              do lx=1,3
                do ly=1,3
                  do lz=1,3
                    mat2w=mat2w+px(istn,istm,lx,ly,lz)*r(ly,istm,istl,ik)    &
                                *r(lz,istl,istn,ik)

                    mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*r(lz,istl,istm,ik) &
                                *r(ly,istn,istl,ik))

                    mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*r(lz,istl,istm,ik) &
                                *r(ly,istn,istl,ik))

                   end do
                 end do
               end do

                if(istl.le.noval(ik)) then
                  wl=el
                  wln=wl-wn
                  wnl=-wln
                  wml=wm-wl
                  wlm=-wml
                else 
                  wl=el+sc
                  wln=wl-wn
                  wnl=-wln
                  wml=wm-wl
                  wlm=-wml
                end if

                b112=0.d0            
                b122=mat1w1*(1.d0/(wnm+wlm))
                b132=mat1w2*(1.d0/(wnm+wnl))

                b242=0.d0

                b222=0.d0
                b212=0.d0

                if (dabs(wnl).gt.tol) then
                   b112=mat2w/wln             
                   b122=b122/wnl
                   b132=b132/wnl

                   b242=mat2w/wln

                   b222=mat1w2/wln

                   b312=-1.d0*mat1w1/wnl
                else
                  b112=0.d0
                  b122=0.d0
                  b132=0.d0

                  b242=0.d0

                  b222=0.d0

                  b312=0.d0
                end if

                if (dabs(wlm).gt.tol) then
                   b112=b112/wml               
                   b122=b122/wlm
                   b132=b132/wlm

                   b242=b242-(mat2w/wml)

                   b212=mat1w1/wml

                   b312=b312+(mat1w2/wlm)  !sig +
                else
                  b112=0.d0
                  b122=0.d0
                  b132=0.d0

                  b212=0.d0
                end if                

                if (dabs(wlm-wnl).gt.tol) then
                   b112=b112/(wlm-wnl)
                else
                  b112=0.d0
                end if

                b11=b11+2.d0*b112
                b12_13=b12_13-b122+b132
                b24=b24+2.d0*b242
                b21_22=b21_22-b212+b222

                b31_32=b31_32+b312

            end if
            end do

!...................... l > cond(istm).............................................................
            do istl=istm+1,banmax
            el=evalfv(istl,ik)
            if(dabs(el).gt.1.d-7) then

              wl=el+sc

              wln=wl-wn
              wnl=-wln
              wml=wm-wl
              wlm=-wml

              mat2w=0.d0
              mat1w1=0.d0
              mat1w2=0.d0

              do lx=1,3
                do ly=1,3
                  do lz=1,3
                    mat2w=mat2w+px(istn,istm,lx,ly,lz)*r(ly,istm,istl,ik)    &
                                *r(lz,istl,istn,ik)

                    mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*r(lz,istl,istm,ik) &
                                *r(ly,istn,istl,ik)) 

                    mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*r(lz,istl,istm,ik) & 
                                *r(ly,istn,istl,ik))

                   end do
                 end do
               end do
            
               b113=mat2w*(1.d0/(wnl+wml))*(1.d0/wnl)
               b123=mat1w1*(1.d0/wnl)
               b133=mat1w2*(1.d0/wnl)*(1.d0/(wnl+wnm))

               b243=mat2w/wln
               b223=mat1w2/wln
               b213=0.d0

               b313=-1.d0*mat1w1/wnl    !sig -

               if(dabs(wml).gt.tol) then
                 b113=b113/wml
                 b123=b123/wml
                 b133=b133/wml

                 b243=b243-(mat2w/wml)
                 b213=mat1w1/wml

                 b313=b313+(mat1w2/wlm)

               else
                 b113=0.d0
                 b123=0.d0
                 b133=0.d0

                 b213=0.d0
               end if

               if(dabs(wnm-wml).gt.tol) then
                 b123=b123/(wnm-wml)
               else
                 b123=0.d0              
               end if

               b11=b11+2.d0*b113
               b12_13=b12_13+b123-b133

               b21_22=b21_22-b213+b223 

               b24=b24+2.d0*b243        

               b31_32=b31_32+b313
            end if
            end do

! ......................... l loops ends ..........................................................             
            b11=b11*ii*(1.d0/wnm)*const
            b12_13=b12_13*ii*(1.d0/wnm)*const
            b24=(b24+b231)*ii*(1.d0/(wnm**3))*const 
 
            b21_22=(b21_22)*ii*(1.d0/(wnm**3))*const   

            b31_32=(b31_32-b331)*(1.d0/(wmn**3))*const*0.5d0*ii

                 do iw=1,emesh
          	  w=(iw-1)*dw+idel
                  inter2w(iw)=inter2w(iw)+                              &
                             (wkpt(ik)*(b11/(wmn-2.d0*w)))
                  inter1w(iw)=inter1w(iw)+                              &
                             (wkpt(ik)*(b12_13/(wmn-w)))

                  intra2w(iw)=intra2w(iw)+                              &
                             (wkpt(ik)*(b24/(wmn-2.d0*w)))
                  intra1w(iw)=intra1w(iw)+                              &
                             (wkpt(ik)*((b21_22+b31_32)/(wmn-w)))       
           
                 end do
!........................................................................
          end if
          end do

        end if
        end do

!.............. End of the main k-loop ............................................................
      end do

!     PMAT.OUT
      close(50)  
      deallocate(pmat)
      
      do i=1,emesh 
        totim(i)=aimag(inter2w(i)+inter1w(i)+intra2w(i)+intra1w(i))/1.d-7
      end do
     
!............... writting the output ..............................................................
!............... File description:
!                22: Chi.Im : different contributions to Im(chi(2w,w,w))
!                23: Chi.Re : different contributions to Re(chi(2w,w,w))
!                24: Chi.Totim : Total Im(chi(2w,w,w))
!                25: Chi.Totre : Total Re(chi(2w,w,w))
!..................................................................................................
      
      write(166,*)
      write(166,*) ' Components of the non-linear susceptibility chi have &
     &  been printed into files:'
      write(166,*) '     NLOChiIm.OUT : different contributions to Im(chi(2w,w,w))'
      write(166,*) '     NLOChiRe.OUT : different contributions to Re(chi(2w,w,w))'
      write(166,*) '     NLOTotIm.OUT : Total Im(chi(2w,w,w))'
      write(166,*) '     NLOTotRe.OUT : Total Re(chi(2w,w,w))'
      write(166,*) '     NLOTotAbs.OUT : Total Abs(chi(2w,w,w))'
      write(166,*)

      open(22,file='NLOChiIm.OUT',status='UNKNOWN',form='FORMATTED')
      open(23,file='NLOChiRe.OUT',status='UNKNOWN',form='FORMATTED')
      open(24,file='NLOTotIm.OUT',status='UNKNOWN',form='FORMATTED')
      open(25,file='NLOTotRe.OUT',status='UNKNOWN',form='FORMATTED')
      open(26,file='NLOTotAbs.OUT',status='UNKNOWN',form='FORMATTED')
      
      totre=0.d0
      totabs=0.d0

      write(22,*)'# Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
      write(22,*)'# in esu'

      write(23,*)'# Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
      write(23,*)'# in esu'

      write(24,*)'# Energy      Tot-Im Chi(-2w,w,w)  Tot-Im Chi(-2w,w,w)'
      write(24,*)'# eV          *10^-7 esu        *10^-12 pm/V SI units '
      write(24,*)'#'

      write(25,*)'# Energy      Tot-Re Chi(-2w,w,w)  Tot-Re Chi(-2w,w,w)'
      write(25,*)'# eV          *10^-7 esu        *10^-12 pm/V SI units '
      write(25,*)'# '

      write(26,*)'# Energy(eV)  |TotChi(-2w,w,w)|   |Tot Chi(-2w,w,w)|'
      write(26,*)'# eV          *10^-7 esu        *10^-12 pm/V SI units '
      write(26,*)'# '

!........................... esu to SI(m/v)=(val_esu)*(4xpi)/30000........................

      do iw=1,emesh
        e=(iw-1)*dw
        e=e*13.6057
        write(22,'(f15.6,4e15.6)') e,aimag(inter2w(iw))/1.d-7,        &
                      aimag(inter1w(iw))/1.d-7,                       &
                      aimag(intra2w(iw))/1.d-7,                       &
                      aimag(intra1w(iw))/1.d-7

        write(23,'(f15.6,4e15.6)') e,dble(inter2w(iw))/1.d-7,         &
                      dble(inter1w(iw))/1.d-7,                        &
                      dble(intra2w(iw))/1.d-7,                        &
                      dble(intra1w(iw))/1.d-7 

        write(24,'(f15.6,2e15.6)') e,totim(iw),totim(iw)*4.d0*pi*(1.d0/30000.d0)*(1.d0/1.d-5)

        totre=dble(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw))/1.d-7
        write(25,'(f15.6,2e15.6)') e,totre,totre*4.d0*pi*(1.d0/30000.d0)*(1.d0/1.d-5)
        totre=0.d0
         
        totabs=abs(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw))/1.d-7
        write(26,'(f15.6,2e15.6)') e,totabs,totabs*4.d0*pi*(1.d0/30000.d0)*(1.d0/1.d-5)
        totabs=0.d0

      end do ! iw
      
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)

!      call kkinvert(totim)

!-------------------------------------------------------------------------------
!     INFO
!-------------------------------------------------------------------------------
      
      call linmsg(166,'-','Some info for debugging')
      
      write(166,*) ' '
      write(166,*) 'information about calculation just performed:'
      write(166,*) ' '
      write(166,*) 'no symmetry operations:', nsymcrys
      write(166,*) 'val',noval(1),'cond',nocond(1),'total bands',tot_ban(1),'at 1st kpt'
      write(166,*) 'energy window:',(emax-emin)*13.6,'eV'

!..................................................................................................
      return
      end subroutine

