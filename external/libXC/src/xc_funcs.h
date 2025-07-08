#define  XC_LDA_X                            1 /* Slater exchange */
#define  XC_LDA_C_WIGNER                     2 /* Wigner */
#define  XC_LDA_C_RPA                        3 /* Random Phase Approximation (RPA) */
#define  XC_LDA_C_HL                         4 /* Hedin & Lundqvist */
#define  XC_LDA_C_GL                         5 /* Gunnarson & Lundqvist */
#define  XC_LDA_C_XALPHA                     6 /* Slater's Xalpha */
#define  XC_LDA_C_VWN                        7 /* Vosko, Wilk & Nusair (VWN5) */
#define  XC_LDA_C_VWN_RPA                    8 /* Vosko, Wilk & Nusair (VWN5_RPA) */
#define  XC_LDA_C_PZ                         9 /* Perdew & Zunger */
#define  XC_LDA_C_PZ_MOD                    10 /* Perdew & Zunger (Modified) */
#define  XC_LDA_C_OB_PZ                     11 /* Ortiz & Ballone (PZ parametrization) */
#define  XC_LDA_C_PW                        12 /* Perdew & Wang */
#define  XC_LDA_C_PW_MOD                    13 /* Perdew & Wang (modified) */
#define  XC_LDA_C_OB_PW                     14 /* Ortiz & Ballone (PW parametrization) */
#define  XC_LDA_C_2D_AMGB                   15 /* AMGB (for 2D systems) */
#define  XC_LDA_C_2D_PRM                    16 /* PRM (for 2D systems) */
#define  XC_LDA_C_VBH                       17 /* von Barth & Hedin */
#define  XC_LDA_C_1D_CSC                    18 /* Casula, Sorella & Senatore */
#define  XC_LDA_X_2D                        19 /* Slater exchange */
#define  XC_LDA_XC_TETER93                  20 /* Teter 93 */
#define  XC_LDA_X_1D_SOFT                   21 /* Exchange in 1D for an soft-Coulomb interaction */
#define  XC_LDA_C_ML1                       22 /* Modified LSD (version 1) of Proynov and Salahub */
#define  XC_LDA_C_ML2                       23 /* Modified LSD (version 2) of Proynov and Salahub */
#define  XC_LDA_C_GOMBAS                    24 /* Gombas */
#define  XC_LDA_C_PW_RPA                    25 /* Perdew & Wang (fit to the RPA energy) */
#define  XC_LDA_C_1D_LOOS                   26 /* P-F Loos correlation LDA */
#define  XC_LDA_C_RC04                      27 /* Ragot-Cortona */
#define  XC_LDA_C_VWN_1                     28 /* Vosko, Wilk & Nusair (VWN1) */
#define  XC_LDA_C_VWN_2                     29 /* Vosko, Wilk & Nusair (VWN2) */
#define  XC_LDA_C_VWN_3                     30 /* Vosko, Wilk & Nusair (VWN3) */
#define  XC_LDA_C_VWN_4                     31 /* Vosko, Wilk & Nusair (VWN4) */
#define  XC_GGA_X_GAM                       32 /* Minnesota GAM exhange functional */
#define  XC_GGA_C_GAM                       33 /* Minnesota GAM correlation functional */
#define  XC_GGA_X_HCTH_A                    34 /* HCTH-A */
#define  XC_GGA_X_EV93                      35 /* Engel and Vosko */
#define  XC_HYB_MGGA_X_DLDF                 36 /* Dispersionless Density Functional */
#define  XC_MGGA_C_DLDF                     37 /* Dispersionless Density Functional */
#define  XC_GGA_X_BCGP                      38 /* Burke, Cancio, Gould, and Pittalis */
#define  XC_GGA_C_ACGGA                     39 /* acGGA, asymptotically corrected GGA correlation */
#define  XC_GGA_X_LAMBDA_OC2_N              40 /* lambda_OC2(N) version of PBE */
#define  XC_GGA_X_B86_R                     41 /* Revised Becke 86 with modified gradient correction */
#define  XC_MGGA_XC_ZLP                     42 /* Zhao, Levy & Parr, Eq. (21) */
#define  XC_LDA_XC_ZLP                      43 /* Zhao, Levy & Parr, Eq. (20) */
#define  XC_GGA_X_LAMBDA_CH_N               44 /* lambda_CH(N) version of PBE */
#define  XC_GGA_X_LAMBDA_LO_N               45 /* lambda_LO(N) version of PBE */
#define  XC_GGA_X_HJS_B88_V2                46 /* HJS screened exchange B88 corrected version */
#define  XC_GGA_C_Q2D                       47 /* Chiodo et al */
#define  XC_GGA_X_Q2D                       48 /* Chiodo et al */
#define  XC_GGA_X_PBE_MOL                   49 /* Reparametrized PBE by del Campo, Gazquez, Trickey & Vela */
#define  XC_LDA_K_TF                        50 /* Thomas-Fermi kinetic energy */
#define  XC_LDA_K_LP                        51 /* Lee and Parr Gaussian ansatz for the kinetic energy */
#define  XC_GGA_K_TFVW                      52 /* Thomas-Fermi plus von Weiszaecker correction */
#define  XC_GGA_K_REVAPBEINT                53 /* interpolated version of revAPBE */
#define  XC_GGA_K_APBEINT                   54 /* interpolated version of APBE */
#define  XC_GGA_K_REVAPBE                   55 /* revised APBE */
#define  XC_GGA_X_AK13                      56 /* Armiento & Kuemmel 2013 */
#define  XC_GGA_K_MEYER                     57 /* Meyer,  Wang, and Young */
#define  XC_GGA_X_LV_RPW86                  58 /* Berland and Hyldgaard */
#define  XC_GGA_X_PBE_TCA                   59 /* PBE revised by Tognetti et al */
#define  XC_GGA_X_PBEINT                    60 /* PBE for hybrid interfaces */
#define  XC_GGA_C_ZPBEINT                   61 /* spin-dependent gradient correction to PBEint */
#define  XC_GGA_C_PBEINT                    62 /* PBE for hybrid interfaces */
#define  XC_GGA_C_ZPBESOL                   63 /* spin-dependent gradient correction to PBEsol */
#define  XC_MGGA_XC_OTPSS_D                 64 /* oTPSS-D functional of Goerigk and Grimme */
#define  XC_GGA_XC_OPBE_D                   65 /* oPBE-D functional of Goerigk and Grimme */
#define  XC_GGA_XC_OPWLYP_D                 66 /* oPWLYP-D functional of Goerigk and Grimme */
#define  XC_GGA_XC_OBLYP_D                  67 /* oBLYP-D functional of Goerigk and Grimme */
#define  XC_GGA_X_VMT84_GE                  68 /* VMT{8,4} with constraint satisfaction with mu = mu_GE */
#define  XC_GGA_X_VMT84_PBE                 69 /* VMT{8,4} with constraint satisfaction with mu = mu_PBE */
#define  XC_GGA_X_VMT_GE                    70 /* Vela, Medel, and Trickey with mu = mu_GE */
#define  XC_GGA_X_VMT_PBE                   71 /* Vela, Medel, and Trickey with mu = mu_PBE */
#define  XC_MGGA_C_CS                       72 /* Colle and Salvetti */
#define  XC_MGGA_C_MN12_SX                  73 /* Minnesota MN12-SX correlation functional */
#define  XC_MGGA_C_MN12_L                   74 /* Minnesota MN12-L correlation functional */
#define  XC_MGGA_C_M11_L                    75 /* Minnesota M11-L correlation functional */
#define  XC_MGGA_C_M11                      76 /* Minnesota M11 correlation functional */
#define  XC_MGGA_C_M08_SO                   77 /* Minnesota M08-SO correlation functional */
#define  XC_MGGA_C_M08_HX                   78 /* Minnesota M08 correlation functional */
#define  XC_GGA_C_N12_SX                    79 /* Minnesota N12-SX correlation functional */
#define  XC_GGA_C_N12                       80 /* Minnesota N12 correlation functional */
#define  XC_HYB_GGA_X_N12_SX                81 /* Minnesota N12-SX exchange functional */
#define  XC_GGA_X_N12                       82 /* Minnesota N12 exchange functional */
#define  XC_GGA_C_REGTPSS                   83 /* regularized TPSS correlation */
#define  XC_GGA_C_OP_XALPHA                 84 /* one-parameter progressive functional (Xalpha version) */
#define  XC_GGA_C_OP_G96                    85 /* one-parameter progressive functional (G96 version) */
#define  XC_GGA_C_OP_PBE                    86 /* one-parameter progressive functional (PBE version) */
#define  XC_GGA_C_OP_B88                    87 /* one-parameter progressive functional (B88 version) */
#define  XC_GGA_C_FT97                      88 /* Filatov & Thiel correlation */
#define  XC_GGA_C_SPBE                      89 /* PBE correlation to be used with the SSB exchange */
#define  XC_GGA_X_SSB_SW                    90 /* Swart, Sola and Bickelhaupt correction to PBE */
#define  XC_GGA_X_SSB                       91 /* Swart, Sola and Bickelhaupt */
#define  XC_GGA_X_SSB_D                     92 /* Swart, Sola and Bickelhaupt dispersion */
#define  XC_GGA_XC_HCTH_407P                93 /* HCTH/407+ */
#define  XC_GGA_XC_HCTH_P76                 94 /* HCTH p=7/6 */
#define  XC_GGA_XC_HCTH_P14                 95 /* HCTH p=1/4 */
#define  XC_GGA_XC_B97_GGA1                 96 /* Becke 97 GGA-1 */
#define  XC_GGA_C_HCTH_A                    97 /* HCTH-A */
#define  XC_GGA_X_BPCCAC                    98 /* BPCCAC (GRAC for the energy) */
#define  XC_GGA_C_REVTCA                    99 /* Tognetti, Cortona, Adamo (revised) */
#define  XC_GGA_C_TCA                      100 /* Tognetti, Cortona, Adamo */
#define  XC_GGA_X_PBE                      101 /* Perdew, Burke & Ernzerhof */
#define  XC_GGA_X_PBE_R                    102 /* Revised PBE from Zhang & Yang */
#define  XC_GGA_X_B86                      103 /* Becke 86 */
#define  XC_GGA_X_B86_MGC                  105 /* Becke 86 with modified gradient correction */
#define  XC_GGA_X_B88                      106 /* Becke 88 */
#define  XC_GGA_X_G96                      107 /* Gill 96 */
#define  XC_GGA_X_PW86                     108 /* Perdew & Wang 86 */
#define  XC_GGA_X_PW91                     109 /* Perdew & Wang 91 */
#define  XC_GGA_X_OPTX                     110 /* Handy & Cohen OPTX 01 */
#define  XC_GGA_X_DK87_R1                  111 /* dePristo & Kress 87 version R1 */
#define  XC_GGA_X_DK87_R2                  112 /* dePristo & Kress 87 version R2 */
#define  XC_GGA_X_LG93                     113 /* Lacks & Gordon 93 */
#define  XC_GGA_X_FT97_A                   114 /* Filatov & Thiel 97 (version A) */
#define  XC_GGA_X_FT97_B                   115 /* Filatov & Thiel 97 (version B) */
#define  XC_GGA_X_PBE_SOL                  116 /* Perdew, Burke & Ernzerhof SOL */
#define  XC_GGA_X_RPBE                     117 /* Hammer, Hansen, and Norskov */
#define  XC_GGA_X_WC                       118 /* Wu & Cohen */
#define  XC_GGA_X_MPW91                    119 /* mPW91 of Adamo & Barone */
#define  XC_GGA_X_AM05                     120 /* Armiento & Mattsson 05 */
#define  XC_GGA_X_PBEA                     121 /* Madsen 07 */
#define  XC_GGA_X_MPBE                     122 /* Adamo & Barone modification to PBE */
#define  XC_GGA_X_XPBE                     123 /* Extended PBE by Xu & Goddard III */
#define  XC_GGA_X_2D_B86_MGC               124 /* Becke 86 with modified gradient correction for 2D */
#define  XC_GGA_X_BAYESIAN                 125 /* Bayesian best fit for the enhancement factor */
#define  XC_GGA_X_PBE_JSJR                 126 /* Reparametrized PBE by Pedroza, Silva & Capelle */
#define  XC_GGA_X_2D_B88                   127 /* Becke 88 in 2D */
#define  XC_GGA_X_2D_B86                   128 /* Becke 86 in 2D */
#define  XC_GGA_X_2D_PBE                   129 /* Perdew, Burke & Ernzerhof in 2D */
#define  XC_GGA_C_PBE                      130 /* Perdew, Burke & Ernzerhof */
#define  XC_GGA_C_LYP                      131 /* Lee, Yang & Parr */
#define  XC_GGA_C_P86                      132 /* Perdew 86 */
#define  XC_GGA_C_PBE_SOL                  133 /* Perdew, Burke & Ernzerhof SOL */
#define  XC_GGA_C_PW91                     134 /* Perdew & Wang 91 */
#define  XC_GGA_C_AM05                     135 /* Armiento & Mattsson 05 */
#define  XC_GGA_C_XPBE                     136 /* Extended PBE by Xu & Goddard III */
#define  XC_GGA_C_LM                       137 /* Langreth & Mehl */
#define  XC_GGA_C_PBE_JRGX                 138 /* Reparametrized PBE by Pedroza, Silva & Capelle */
#define  XC_GGA_X_OPTB88_VDW               139 /* opt-Becke 88 for vdW */
#define  XC_GGA_X_PBEK1_VDW                140 /* Reparametrized PBE for vdW */
#define  XC_GGA_X_OPTPBE_VDW               141 /* Reparametrized PBE for vdW */
#define  XC_GGA_X_RGE2                     142 /* Regularized PBE */
#define  XC_GGA_C_RGE2                     143 /* Regularized PBE */
#define  XC_GGA_X_RPW86                    144 /* Refitted Perdew & Wang 86 */
#define  XC_GGA_X_KT1                      145 /* Exchange part of Keal and Tozer version 1 */
#define  XC_GGA_XC_KT2                     146 /* Keal and Tozer, version 2 */
#define  XC_GGA_C_WL                       147 /* Wilson & Levy */
#define  XC_GGA_C_WI                       148 /* Wilson & Ivanov */
#define  XC_GGA_X_MB88                     149 /* Modified Becke 88 for proton transfer */
#define  XC_GGA_X_SOGGA                    150 /* Second-order generalized gradient approximation */
#define  XC_GGA_X_SOGGA11                  151 /* Second-order generalized gradient approximation 2011 */
#define  XC_GGA_C_SOGGA11                  152 /* Second-order generalized gradient approximation 2011 */
#define  XC_GGA_C_WI0                      153 /* Wilson & Ivanov initial version */
#define  XC_GGA_XC_TH1                     154 /* Tozer and Handy v. 1 */
#define  XC_GGA_XC_TH2                     155 /* Tozer and Handy v. 2 */
#define  XC_GGA_XC_TH3                     156 /* Tozer and Handy v. 3 */
#define  XC_GGA_XC_TH4                     157 /* Tozer and Handy v. 4 */
#define  XC_GGA_X_C09X                     158 /* C09x to be used with the VdW of Rutgers-Chalmers */
#define  XC_GGA_C_SOGGA11_X                159 /* To be used with HYB_GGA_X_SOGGA11_X */
#define  XC_GGA_X_LB                       160 /* van Leeuwen & Baerends */
#define  XC_GGA_XC_HCTH_93                 161 /* HCTH/93 */
#define  XC_GGA_XC_HCTH_120                162 /* HCTH/120 */
#define  XC_GGA_XC_HCTH_147                163 /* HCTH/147 */
#define  XC_GGA_XC_HCTH_407                164 /* HCTH/407 */
#define  XC_GGA_XC_EDF1                    165 /* EDF1 */
#define  XC_GGA_XC_XLYP                    166 /* XLYP */
#define  XC_GGA_XC_KT1                     167 /* Keal and Tozer, version 1 */
#define  XC_GGA_X_LSPBE                    168 /* lsPBE, a PW91-like modification of PBE exchange */
#define  XC_GGA_X_LSRPBE                   169 /* lsRPBE, a PW91-like modification of RPBE */
#define  XC_GGA_XC_B97_D                   170 /* Becke 97-D */
#define  XC_GGA_X_OPTB86B_VDW              171 /* Becke 86 reoptimized for use with vdW functional of Dion et al */
#define  XC_MGGA_C_REVM11                  172 /* Revised Minnesota M11 correlation functional */
#define  XC_GGA_XC_PBE1W                   173 /* PBE1W */
#define  XC_GGA_XC_MPWLYP1W                174 /* mPWLYP1w */
#define  XC_GGA_XC_PBELYP1W                175 /* PBELYP1W */
#define  XC_GGA_C_ACGGAP                   176 /* acGGA+, asymptotically corrected GGA correlation+ */
#define  XC_HYB_LDA_XC_LDA0                177 /* LDA hybrid exchange (LDA0) */
#define  XC_HYB_LDA_XC_CAM_LDA0            178 /* CAM version of LDA0 */
#define  XC_GGA_X_B88_6311G                179 /* Becke 88 reoptimized with the 6-311G** basis set */
#define  XC_GGA_X_NCAP                     180 /* Nearly correct asymptotic potential */
#define  XC_GGA_XC_NCAP                    181 /* NCAP exchange + P86 correlation */
#define  XC_GGA_X_LBM                      182 /* van Leeuwen & Baerends modified */
#define  XC_GGA_X_OL2                      183 /* Exchange form based on Ou-Yang and Levy v.2 */
#define  XC_GGA_X_APBE                     184 /* mu fixed from the semiclassical neutral atom */
#define  XC_GGA_K_APBE                     185 /* mu fixed from the semiclassical neutral atom */
#define  XC_GGA_C_APBE                     186 /* mu fixed from the semiclassical neutral atom */
#define  XC_GGA_K_TW1                      187 /* Tran and Wesolowski set 1 (Table II) */
#define  XC_GGA_K_TW2                      188 /* Tran and Wesolowski set 2 (Table II) */
#define  XC_GGA_K_TW3                      189 /* Tran and Wesolowski set 3 (Table II) */
#define  XC_GGA_K_TW4                      190 /* Tran and Wesolowski set 4 (Table II) */
#define  XC_GGA_X_HTBS                     191 /* Haas, Tran, Blaha, and Schwarz */
#define  XC_GGA_X_AIRY                     192 /* Constantin et al based on the Airy gas */
#define  XC_GGA_X_LAG                      193 /* Local Airy Gas */
#define  XC_GGA_XC_MOHLYP                  194 /* Functional for organometallic chemistry */
#define  XC_GGA_XC_MOHLYP2                 195 /* Functional for barrier heights */
#define  XC_GGA_XC_TH_FL                   196 /* Tozer and Handy v. FL */
#define  XC_GGA_XC_TH_FC                   197 /* Tozer and Handy v. FC */
#define  XC_GGA_XC_TH_FCFO                 198 /* Tozer and Handy v. FCFO */
#define  XC_GGA_XC_TH_FCO                  199 /* Tozer and Handy v. FCO */
#define  XC_GGA_C_OPTC                     200 /* Optimized correlation functional of Cohen and Handy */
#define  XC_MGGA_X_LTA                     201 /* Local tau approximation */
#define  XC_MGGA_X_TPSS                    202 /* Tao, Perdew, Staroverov & Scuseria */
#define  XC_MGGA_X_M06_L                   203 /* Minnesota M06-L exchange functional */
#define  XC_MGGA_X_GVT4                    204 /* GVT4 (X part of VSXC) */
#define  XC_MGGA_X_TAU_HCTH                205 /* tau-HCTH from Boese and Handy */
#define  XC_MGGA_X_BR89                    206 /* Becke-Roussel 89, gamma = 0.8 */
#define  XC_MGGA_X_BJ06                    207 /* Becke & Johnson 06 */
#define  XC_MGGA_X_TB09                    208 /* Tran & Blaha 09 */
#define  XC_MGGA_X_RPP09                   209 /* Rasanen, Pittalis & Proetto 09 */
#define  XC_MGGA_X_2D_PRHG07               210 /* Pittalis-Rasanen-Helbig-Gross 2007 */
#define  XC_MGGA_X_2D_PRHG07_PRP10         211 /* PRHG07 with Pittalis-Rasanen-Proetto 2010 correction */
#define  XC_MGGA_X_REVTPSS                 212 /* revised Tao, Perdew, Staroverov & Scuseria */
#define  XC_MGGA_X_PKZB                    213 /* Perdew, Kurth, Zupan, and Blaha */
#define  XC_MGGA_X_BR89_1                  214 /* Becke-Roussel 89, gamma = 1.0 */
#define  XC_GGA_X_ECMV92                   215 /* Engel, Chevary, Macdonald and Vosko */
#define  XC_GGA_C_PBE_VWN                  216 /* Perdew, Burke & Ernzerhof based on VWN correlation */
#define  XC_GGA_C_P86_FT                   217 /* Perdew 86 with more accurate value for ftilde */
#define  XC_GGA_K_RATIONAL_P               218 /* RATIONAL$^{p}$ by Lehtomaki and Lopez-Acevedo (by default $p=3/2$, $C_{2}=0.7687$) */
#define  XC_GGA_K_PG1                      219 /* PG1 (Pauli-Gaussian) functional by Constantin, Fabiano, and Della Sala */
#define  XC_MGGA_K_PGSL025                 220 /* PGSL025 (Pauli-Gaussian) functional by Constantin, Fabiano, and Della Sala */
#define  XC_MGGA_X_MS0                     221 /* MS exchange of Sun, Xiao, and Ruzsinszky */
#define  XC_MGGA_X_MS1                     222 /* MS1 exchange of Sun, et al */
#define  XC_MGGA_X_MS2                     223 /* MS2 exchange of Sun, et al */
#define  XC_HYB_MGGA_X_MS2H                224 /* MS2 hybrid exchange of Sun, et al */
#define  XC_MGGA_X_TH                      225 /* Tsuneda and Hirao */
#define  XC_MGGA_X_M11_L                   226 /* Minnesota M11-L exchange functional */
#define  XC_MGGA_X_MN12_L                  227 /* Minnesota MN12-L exchange functional */
#define  XC_MGGA_X_MS2_REV                 228 /* MS2 exchange of Sun, et al with revised value for c */
#define  XC_MGGA_XC_CC06                   229 /* Cancio and Chou 2006 */
#define  XC_MGGA_X_MK00                    230 /* Exchange for accurate virtual orbital energies */
#define  XC_MGGA_C_TPSS                    231 /* Tao, Perdew, Staroverov & Scuseria */
#define  XC_MGGA_C_VSXC                    232 /* VSXC (correlation part) */
#define  XC_MGGA_C_M06_L                   233 /* Minnesota M06-L correlation functional */
#define  XC_MGGA_C_M06_HF                  234 /* Minnesota M06-HF correlation functional */
#define  XC_MGGA_C_M06                     235 /* Minnesota M06 correlation functional */
#define  XC_MGGA_C_M06_2X                  236 /* Minnesota M06-2X correlation functional */
#define  XC_MGGA_C_M05                     237 /* Minnesota M05 correlation functional */
#define  XC_MGGA_C_M05_2X                  238 /* Minnesota M05-2X correlation functional */
#define  XC_MGGA_C_PKZB                    239 /* Perdew, Kurth, Zupan, and Blaha */
#define  XC_MGGA_C_BC95                    240 /* Becke correlation 95 */
#define  XC_MGGA_C_REVTPSS                 241 /* revised TPSS correlation */
#define  XC_MGGA_XC_TPSSLYP1W              242 /* TPSSLYP1W */
#define  XC_MGGA_X_MK00B                   243 /* Exchange for accurate virtual orbital energies (v. B) */
#define  XC_MGGA_X_BLOC                    244 /* functional with balanced localization */
#define  XC_MGGA_X_MODTPSS                 245 /* Modified Tao, Perdew, Staroverov & Scuseria */
#define  XC_GGA_C_PBELOC                   246 /* Semilocal dynamical correlation */
#define  XC_MGGA_C_TPSSLOC                 247 /* Semilocal dynamical correlation */
#define  XC_HYB_MGGA_X_MN12_SX             248 /* Minnesota MN12-SX hybrid exchange functional */
#define  XC_MGGA_X_MBEEF                   249 /* mBEEF exchange */
#define  XC_MGGA_X_MBEEFVDW                250 /* mBEEF-vdW exchange */
#define  XC_MGGA_C_TM                      251 /* Tao and Mo 2016 correlation */
#define  XC_GGA_C_P86VWN                   252 /* Perdew 86 based on VWN5 correlation */
#define  XC_GGA_C_P86VWN_FT                253 /* Perdew 86 based on VWN5 correlation, with more accurate value for ftilde */
#define  XC_MGGA_XC_B97M_V                 254 /* B97M-V exchange-correlation functional */
#define  XC_GGA_XC_VV10                    255 /* Vydrov and Van Voorhis */
#define  XC_MGGA_X_JK                      256 /* Jemmer-Knowles meta-GGA exchange */
#define  XC_MGGA_X_MVS                     257 /* MVS exchange of Sun, Perdew, and Ruzsinszky */
#define  XC_GGA_C_PBEFE                    258 /* PBE for formation energies */
#define  XC_LDA_XC_KSDT                    259 /* Karasiev, Sjostrom, Dufty & Trickey */
#define  XC_MGGA_X_MN15_L                  260 /* Minnesota MN15-L exchange functional */
#define  XC_MGGA_C_MN15_L                  261 /* Minnesota MN15-L correlation functional */
#define  XC_GGA_C_OP_PW91                  262 /* one-parameter progressive functional (PW91 version) */
#define  XC_MGGA_X_SCAN                    263 /* SCAN exchange of Sun, Ruzsinszky, and Perdew */
#define  XC_HYB_MGGA_X_SCAN0               264 /* SCAN hybrid exchange (SCAN0) */
#define  XC_GGA_X_PBEFE                    265 /* PBE for formation energies */
#define  XC_HYB_GGA_XC_B97_1P              266 /* version of B97 by Cohen and Handy */
#define  XC_MGGA_C_SCAN                    267 /* SCAN correlation of Sun, Ruzsinszky, and Perdew */
#define  XC_HYB_MGGA_X_MN15                268 /* Minnesota MN15 hybrid exchange functional */
#define  XC_MGGA_C_MN15                    269 /* Minnesota MN15 correlation functional */
#define  XC_GGA_X_CAP                      270 /* Correct Asymptotic Potential */
#define  XC_GGA_X_EB88                     271 /* Non-empirical (excogitated) B88 functional of Becke and Elliott */
#define  XC_GGA_C_PBE_MOL                  272 /* Reparametrized PBE by del Campo, Gazquez, Trickey & Vela */
#define  XC_HYB_GGA_XC_PBE_MOL0            273 /* PBEmol0 */
#define  XC_HYB_GGA_XC_PBE_SOL0            274 /* PBEsol0 */
#define  XC_HYB_GGA_XC_PBEB0               275 /* PBEbeta0 */
#define  XC_HYB_GGA_XC_PBE_MOLB0           276 /* PBEmolbeta0 */
#define  XC_GGA_K_ABSP3                    277 /* gamma-TFvW form by Acharya et al [$g = 1 - 1.513/N^{0.35}]$ */
#define  XC_GGA_K_ABSP4                    278 /* gamma-TFvW form by Acharya et al [$g = l = 1/(1 + 1.332/N^{1/3})$] */
#define  XC_HYB_MGGA_X_BMK                 279 /* Boese-Martin for kinetics */
#define  XC_GGA_C_BMK                      280 /* Boese-Martin correlation for kinetics */
#define  XC_GGA_C_TAU_HCTH                 281 /* correlation part of tau-hcth */
#define  XC_HYB_MGGA_X_TAU_HCTH            282 /* Hybrid version of tau-HCTH */
#define  XC_GGA_C_HYB_TAU_HCTH             283 /* correlation part of hyb-tau-hcth */
#define  XC_MGGA_X_B00                     284 /* Becke 2000 */
#define  XC_GGA_X_BEEFVDW                  285 /* BEEF-vdW exchange */
#define  XC_GGA_XC_BEEFVDW                 286 /* BEEF-vdW exchange-correlation */
#define  XC_LDA_C_CHACHIYO                 287 /* Chachiyo simple 2 parameter correlation */
#define  XC_MGGA_XC_HLE17                  288 /* high local exchange 2017 */
#define  XC_LDA_C_LP96                     289 /* Liu-Parr correlation */
#define  XC_HYB_GGA_XC_PBE50               290 /* PBE50 */
#define  XC_GGA_X_PBETRANS                 291 /* Gradient-regulated connection-based correction for the PBE exchange */
#define  XC_MGGA_C_SCAN_RVV10              292 /* SCAN + rVV10 correlation */
#define  XC_MGGA_X_REVM06_L                293 /* Minnesota revM06-L exchange functional */
#define  XC_MGGA_C_REVM06_L                294 /* Minnesota revM06-L correlation functional */
#define  XC_HYB_MGGA_X_M08_HX              295 /* Minnesota M08-HX hybrid exchange functional */
#define  XC_HYB_MGGA_X_M08_SO              296 /* Minnesota M08-SO hybrid exchange functional */
#define  XC_HYB_MGGA_X_M11                 297 /* Minnesota M11 hybrid exchange functional */
#define  XC_GGA_X_CHACHIYO                 298 /* Chachiyo exchange */
#define  XC_MGGA_X_RTPSS                   299 /* TPSS for surface adsorption */
#define  XC_MGGA_X_MS2B                    300 /* MS2beta exchange of Furness and Sun */
#define  XC_MGGA_X_MS2BS                   301 /* MS2beta* exchange of Furness and Sun */
#define  XC_MGGA_X_MVSB                    302 /* MVSbeta exchange by Furness and Sun */
#define  XC_MGGA_X_MVSBS                   303 /* MVSbeta* exchange by Furness and Sun */
#define  XC_HYB_MGGA_X_REVM11              304 /* Revised Minnesota M11 hybrid exchange functional */
#define  XC_HYB_MGGA_X_REVM06              305 /* Revised Minnesota M06 hybrid exchange functional */
#define  XC_MGGA_C_REVM06                  306 /* Revised Minnesota M06 correlation functional */
#define  XC_LDA_C_CHACHIYO_MOD             307 /* Chachiyo simple 2 parameter correlation with modified spin scaling */
#define  XC_LDA_C_KARASIEV_MOD             308 /* Karasiev reparameterization of Chachiyo */
#define  XC_GGA_C_CHACHIYO                 309 /* Chachiyo simple GGA correlation */
#define  XC_HYB_MGGA_X_M06_SX              310 /* Minnesota M06-SX short-range hybrid exchange functional */
#define  XC_MGGA_C_M06_SX                  311 /* Minnesota M06-SX correlation functional */
#define  XC_GGA_X_REVSSB_D                 312 /* Revised Swart, Sola and Bickelhaupt dispersion */
#define  XC_GGA_C_CCDF                     313 /* ccDF: coupled-cluster motivated density functional */
#define  XC_HYB_GGA_XC_HFLYP               314 /* HF + LYP correlation */
#define  XC_HYB_GGA_XC_B3P86_NWCHEM        315 /* B3P86, NWChem version */
#define  XC_GGA_X_PW91_MOD                 316 /* PW91, alternate version with more digits */
#define  XC_LDA_C_W20                      317 /* Xie, Wu, and Zhao interpolation ansatz without fitting parameters */
#define  XC_LDA_XC_CORRKSDT                318 /* Corrected KSDT by Karasiev, Dufty and Trickey */
#define  XC_MGGA_X_FT98                    319 /* Filatov and Thiel 1998 meta-GGA exchange */
#define  XC_GGA_X_PBE_MOD                  320 /* Perdew, Burke & Ernzerhof with less precise value for beta */
#define  XC_GGA_X_PBE_GAUSSIAN             321 /* Perdew, Burke & Ernzerhof with parameter values used in Gaussian */
#define  XC_GGA_C_PBE_GAUSSIAN             322 /* Perdew, Burke & Ernzerhof with parameters from Gaussian */
#define  XC_MGGA_C_TPSS_GAUSSIAN           323 /* Tao, Perdew, Staroverov & Scuseria with parameters from Gaussian */
#define  XC_GGA_X_NCAPR                    324 /* Nearly correct asymptotic potential revised */
#define  XC_HYB_GGA_XC_RELPBE0             325 /* relPBE0 a.k.a. relPBE: PBE0 refitted for actinide compounds */
#define  XC_MGGA_X_EEL                     326 /* Exact exchange-like exchange of Aschebrock et al */
#define  XC_GGA_XC_B97_3C                  327 /* Becke 97-3c by Grimme et. al. */
#define  XC_LDA_C_EPC17                    328 /* epc17(-1): electron-proton correlation 2017 */
#define  XC_LDA_C_EPC17_2                  329 /* epc17-2: electron-proton correlation 2017 for proton affinities */
#define  XC_LDA_C_EPC18_1                  330 /* epc18-1: electron-proton correlation 2018 */
#define  XC_LDA_C_EPC18_2                  331 /* epc18-2: electron-proton correlation 2018 for proton affinities */
#define  XC_GGA_X_BKL1                     338 /* Type-I band gap functional by Bhattacharjee, Koshi and Lee */
#define  XC_GGA_X_BKL2                     339 /* Type-II band gap functional by Bhattacharjee, Koshi and Lee */
#define  XC_HYB_MGGA_X_CF22D               340 /* Minnesota CF22D hybrid exchange functional */
#define  XC_MGGA_C_CF22D                   341 /* Minnesota CF22D correlation functional */
#define  XC_MGGA_X_LAK                     342 /* Lebeda-Aschebrock-Kummel meta-GGA exchange */
#define  XC_HYB_GGA_XC_OPB3LYP             386 /* opB3LYP: B3LYP reoptimized in 6-311++G(2d,2p) basis set */
#define  XC_MGGA_C_CC                      387 /* Self-interaction corrected correlation functional by Schmidt et al */
#define  XC_MGGA_C_CCALDA                  388 /* Iso-orbital corrected LDA correlation by Lebeda et al */
#define  XC_HYB_MGGA_XC_BR3P86             389 /* BR3P86 hybrid meta-GGA from Neumann and Handy */
#define  XC_HYB_GGA_XC_CASE21              390 /* CASE21: Constrained And Smoothed semi-Empirical 2021 functional */
#define  XC_MGGA_C_RREGTM                  391 /* Revised regTM correlation by Jana et al */
#define  XC_HYB_GGA_XC_PBE_2X              392 /* PBE-2X: PBE0 with 56% exact exchange */
#define  XC_HYB_GGA_XC_PBE38               393 /* PBE38: PBE0 with 3/8 = 37.5% exact exchange */
#define  XC_HYB_GGA_XC_B3LYP3              394 /* B3LYP with VWN functional 3 instead of RPA */
#define  XC_HYB_GGA_XC_CAM_O3LYP           395 /* CAM-O3LYP */
#define  XC_HYB_MGGA_XC_TPSS0              396 /* TPSS0 with 25% exact exchange */
#define  XC_MGGA_C_B94                     397 /* Becke 1994 meta-GGA correlation */
#define  XC_HYB_MGGA_XC_B94_HYB            398 /* Becke 1994 hybrid meta-GGA */
#define  XC_HYB_GGA_XC_WB97X_D3            399 /* wB97X-D3 range-separated functional */
#define  XC_HYB_GGA_XC_LC_BLYP             400 /* LC version of BLYP */
#define  XC_HYB_GGA_XC_B3PW91              401 /* The original (ACM, B3PW91) hybrid of Becke */
#define  XC_HYB_GGA_XC_B3LYP               402 /* B3LYP */
#define  XC_HYB_GGA_XC_B3P86               403 /* B3P86 */
#define  XC_HYB_GGA_XC_O3LYP               404 /* O3LYP */
#define  XC_HYB_GGA_XC_MPW1K               405 /* mPW1K */
#define  XC_HYB_GGA_XC_PBEH                406 /* PBEH (PBE0) */
#define  XC_HYB_GGA_XC_B97                 407 /* Becke 97 */
#define  XC_HYB_GGA_XC_B97_1               408 /* Becke 97-1 */
#define  XC_HYB_GGA_XC_APF                 409 /* APF hybrid functional */
#define  XC_HYB_GGA_XC_B97_2               410 /* Becke 97-2 */
#define  XC_HYB_GGA_XC_X3LYP               411 /* X3LYP */
#define  XC_HYB_GGA_XC_B1WC                412 /* B1WC */
#define  XC_HYB_GGA_XC_B97_K               413 /* Boese-Martin for Kinetics */
#define  XC_HYB_GGA_XC_B97_3               414 /* Becke 97-3 */
#define  XC_HYB_GGA_XC_MPW3PW              415 /* MPW3PW of Adamo & Barone */
#define  XC_HYB_GGA_XC_B1LYP               416 /* B1LYP */
#define  XC_HYB_GGA_XC_B1PW91              417 /* B1PW91 */
#define  XC_HYB_GGA_XC_MPW1PW              418 /* mPW1PW */
#define  XC_HYB_GGA_XC_MPW3LYP             419 /* MPW3LYP */
#define  XC_HYB_GGA_XC_SB98_1A             420 /* SB98 (1a) */
#define  XC_HYB_GGA_XC_SB98_1B             421 /* SB98 (1b) */
#define  XC_HYB_GGA_XC_SB98_1C             422 /* SB98 (1c) */
#define  XC_HYB_GGA_XC_SB98_2A             423 /* SB98 (2a) */
#define  XC_HYB_GGA_XC_SB98_2B             424 /* SB98 (2b) */
#define  XC_HYB_GGA_XC_SB98_2C             425 /* SB98 (2c) */
#define  XC_HYB_GGA_X_SOGGA11_X            426 /* Hybrid based on SOGGA11 form */
#define  XC_HYB_GGA_XC_HSE03               427 /* HSE03 */
#define  XC_HYB_GGA_XC_HSE06               428 /* HSE06 */
#define  XC_HYB_GGA_XC_HJS_PBE             429 /* HJS hybrid screened exchange PBE version */
#define  XC_HYB_GGA_XC_HJS_PBE_SOL         430 /* HJS hybrid screened exchange PBE_SOL version */
#define  XC_HYB_GGA_XC_HJS_B88             431 /* HJS hybrid screened exchange B88 version */
#define  XC_HYB_GGA_XC_HJS_B97X            432 /* HJS hybrid screened exchange B97x version */
#define  XC_HYB_GGA_XC_CAM_B3LYP           433 /* CAM version of B3LYP */
#define  XC_HYB_GGA_XC_TUNED_CAM_B3LYP     434 /* CAM version of B3LYP, tuned for excitations and properties */
#define  XC_HYB_GGA_XC_BHANDH              435 /* BHandH i.e. BHLYP */
#define  XC_HYB_GGA_XC_BHANDHLYP           436 /* BHandHLYP */
#define  XC_HYB_GGA_XC_MB3LYP_RC04         437 /* B3LYP with RC04 LDA */
#define  XC_HYB_MGGA_X_M05                 438 /* Minnesota M05 hybrid exchange functional */
#define  XC_HYB_MGGA_X_M05_2X              439 /* Minnesota M05-2X hybrid exchange functional */
#define  XC_HYB_MGGA_XC_B88B95             440 /* Mixture of B88 with BC95 (B1B95) */
#define  XC_HYB_MGGA_XC_B86B95             441 /* Mixture of B86 with BC95 */
#define  XC_HYB_MGGA_XC_PW86B95            442 /* Mixture of PW86 with BC95 */
#define  XC_HYB_MGGA_XC_BB1K               443 /* Mixture of B88 with BC95 from Zhao and Truhlar */
#define  XC_HYB_MGGA_X_M06_HF              444 /* Minnesota M06-HF hybrid exchange functional */
#define  XC_HYB_MGGA_XC_MPW1B95            445 /* Mixture of mPW91 with BC95 from Zhao and Truhlar */
#define  XC_HYB_MGGA_XC_MPWB1K             446 /* Mixture of mPW91 with BC95 for kinetics */
#define  XC_HYB_MGGA_XC_X1B95              447 /* Mixture of X with BC95 */
#define  XC_HYB_MGGA_XC_XB1K               448 /* Mixture of X with BC95 for kinetics */
#define  XC_HYB_MGGA_X_M06                 449 /* Minnesota M06 hybrid exchange functional */
#define  XC_HYB_MGGA_X_M06_2X              450 /* Minnesota M06-2X hybrid exchange functional */
#define  XC_HYB_MGGA_XC_PW6B95             451 /* Mixture of PW91 with BC95 from Zhao and Truhlar */
#define  XC_HYB_MGGA_XC_PWB6K              452 /* Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics */
#define  XC_HYB_GGA_XC_MPWLYP1M            453 /* MPW with 1 par. for metals/LYP */
#define  XC_HYB_GGA_XC_REVB3LYP            454 /* Revised B3LYP */
#define  XC_HYB_GGA_XC_CAMY_BLYP           455 /* CAMY version of BLYP */
#define  XC_HYB_GGA_XC_PBE0_13             456 /* PBE0-1/3 */
#define  XC_HYB_MGGA_XC_TPSSH              457 /* TPSSh */
#define  XC_HYB_MGGA_XC_REVTPSSH           458 /* revTPSSh */
#define  XC_HYB_GGA_XC_B3LYPS              459 /* B3LYP* */
#define  XC_HYB_GGA_XC_QTP17               460 /* Global hybrid for vertical ionization potentials */
#define  XC_HYB_GGA_XC_B3LYP_MCM1          461 /* B3LYP-MCM1 */
#define  XC_HYB_GGA_XC_B3LYP_MCM2          462 /* B3LYP-MCM2 */
#define  XC_HYB_GGA_XC_WB97                463 /* wB97 range-separated functional */
#define  XC_HYB_GGA_XC_WB97X               464 /* wB97X range-separated functional */
#define  XC_HYB_GGA_XC_LRC_WPBEH           465 /* Long-range corrected short-range hybrid PBE (LRC-wPBEh) by Rohrdanz, Martins and Herbert */
#define  XC_HYB_GGA_XC_WB97X_V             466 /* wB97X-V range-separated functional */
#define  XC_HYB_GGA_XC_LCY_PBE             467 /* LCY version of PBE */
#define  XC_HYB_GGA_XC_LCY_BLYP            468 /* LCY version of BLYP */
#define  XC_HYB_GGA_XC_LC_VV10             469 /* Vydrov and Van Voorhis */
#define  XC_HYB_GGA_XC_CAMY_B3LYP          470 /* CAMY version of B3LYP */
#define  XC_HYB_GGA_XC_WB97X_D             471 /* wB97X-D range-separated functional */
#define  XC_HYB_GGA_XC_HPBEINT             472 /* hPBEint */
#define  XC_HYB_GGA_XC_LRC_WPBE            473 /* Long-range corrected PBE (LRC-wPBE) by Rohrdanz, Martins and Herbert */
#define  XC_HYB_MGGA_X_MVSH                474 /* MVSh hybrid exchange functional */
#define  XC_HYB_GGA_XC_B3LYP5              475 /* B3LYP with VWN functional 5 instead of RPA */
#define  XC_HYB_GGA_XC_EDF2                476 /* EDF2 */
#define  XC_HYB_GGA_XC_CAP0                477 /* Correct Asymptotic Potential hybrid */
#define  XC_HYB_GGA_XC_LC_WPBE             478 /* Long-range corrected PBE (LC-wPBE) by Vydrov and Scuseria */
#define  XC_HYB_GGA_XC_HSE12               479 /* HSE12 */
#define  XC_HYB_GGA_XC_HSE12S              480 /* HSE12 (short-range version) */
#define  XC_HYB_GGA_XC_HSE_SOL             481 /* HSEsol */
#define  XC_HYB_GGA_XC_CAM_QTP_01          482 /* CAM-B3LYP retuned using ionization potentials of water */
#define  XC_HYB_GGA_XC_MPW1LYP             483 /* mPW1LYP */
#define  XC_HYB_GGA_XC_MPW1PBE             484 /* mPW1PBE */
#define  XC_HYB_GGA_XC_KMLYP               485 /* Kang-Musgrave hybrid */
#define  XC_HYB_GGA_XC_LC_WPBE_WHS         486 /* Long-range corrected PBE (LC-wPBE) by Weintraub, Henderson and Scuseria */
#define  XC_HYB_GGA_XC_LC_WPBEH_WHS        487 /* Long-range corrected short-range hybrid PBE (LC-wPBE) by Weintraub, Henderson and Scuseria */
#define  XC_HYB_GGA_XC_LC_WPBE08_WHS       488 /* Long-range corrected PBE (LC-wPBE) by Weintraub, Henderson and Scuseria */
#define  XC_HYB_GGA_XC_LC_WPBESOL_WHS      489 /* Long-range corrected PBE (LC-wPBE) by Weintraub, Henderson and Scuseria */
#define  XC_HYB_GGA_XC_CAM_QTP_00          490 /* CAM-B3LYP retuned using ionization potentials of water */
#define  XC_HYB_GGA_XC_CAM_QTP_02          491 /* CAM-B3LYP retuned using ionization potentials of water */
#define  XC_HYB_GGA_XC_LC_QTP              492 /* CAM-B3LYP retuned using ionization potentials of water */
#define  XC_MGGA_X_RSCAN                   493 /* Regularized SCAN exchange by Bartok and Yates */
#define  XC_MGGA_C_RSCAN                   494 /* Regularized SCAN correlation by Bartok and Yates */
#define  XC_GGA_X_S12G                     495 /* Swart 2012 GGA exchange */
#define  XC_HYB_GGA_X_S12H                 496 /* Swart 2012 hybrid GGA exchange */
#define  XC_MGGA_X_R2SCAN                  497 /* Re-regularized SCAN exchange by Furness et al */
#define  XC_MGGA_C_R2SCAN                  498 /* Re-regularized SCAN correlation by Furness et al */
#define  XC_HYB_GGA_XC_BLYP35              499 /* BLYP35 */
#define  XC_GGA_K_VW                       500 /* von Weiszaecker correction to Thomas-Fermi */
#define  XC_GGA_K_GE2                      501 /* Second-order gradient expansion of the kinetic energy density */
#define  XC_GGA_K_GOLDEN                   502 /* TF-lambda-vW form by Golden (l = 13/45) */
#define  XC_GGA_K_YT65                     503 /* TF-lambda-vW form by Yonei and Tomishima (l = 1/5) */
#define  XC_GGA_K_BALTIN                   504 /* TF-lambda-vW form by Baltin (l = 5/9) */
#define  XC_GGA_K_LIEB                     505 /* TF-lambda-vW form by Lieb (l = 0.185909191) */
#define  XC_GGA_K_ABSP1                    506 /* gamma-TFvW form by Acharya et al [$g = 1 - 1.412/N^{1/3}$] */
#define  XC_GGA_K_ABSP2                    507 /* gamma-TFvW form by Acharya et al [$g = 1 - 1.332/N^{1/3}$] */
#define  XC_GGA_K_GR                       508 /* gamma-TFvW form by Gazquez and Robles */
#define  XC_GGA_K_LUDENA                   509 /* gamma-TFvW form by Ludena */
#define  XC_GGA_K_GP85                     510 /* gamma-TFvW form by Ghosh and Parr */
#define  XC_GGA_K_PEARSON                  511 /* Pearson 1992 */
#define  XC_GGA_K_OL1                      512 /* Ou-Yang and Levy v.1 */
#define  XC_GGA_K_OL2                      513 /* Ou-Yang and Levy v.2 */
#define  XC_GGA_K_FR_B88                   514 /* Fuentealba & Reyes (B88 version) */
#define  XC_GGA_K_FR_PW86                  515 /* Fuentealba & Reyes (PW86 version) */
#define  XC_GGA_K_DK                       516 /* DePristo and Kress */
#define  XC_GGA_K_PERDEW                   517 /* Perdew */
#define  XC_GGA_K_VSK                      518 /* Vitos, Skriver, and Kollar */
#define  XC_GGA_K_VJKS                     519 /* Vitos, Johansson, Kollar, and Skriver */
#define  XC_GGA_K_ERNZERHOF                520 /* Ernzerhof */
#define  XC_GGA_K_LC94                     521 /* Lembarki & Chermette */
#define  XC_GGA_K_LLP                      522 /* Lee, Lee & Parr */
#define  XC_GGA_K_THAKKAR                  523 /* Thakkar 1992 */
#define  XC_GGA_X_WPBEH                    524 /* short-range part of the PBE (default w=0 gives PBEh) */
#define  XC_GGA_X_HJS_PBE                  525 /* HJS screened exchange PBE version */
#define  XC_GGA_X_HJS_PBE_SOL              526 /* HJS screened exchange PBE_SOL version */
#define  XC_GGA_X_HJS_B88                  527 /* HJS screened exchange B88 version */
#define  XC_GGA_X_HJS_B97X                 528 /* HJS screened exchange B97x version */
#define  XC_GGA_X_ITYH                     529 /* Short-range recipe for B88 functional - erf */
#define  XC_GGA_X_SFAT                     530 /* Short-range recipe for B88 functional - Yukawa */
#define  XC_HYB_MGGA_XC_WB97M_V            531 /* wB97M-V exchange-correlation functional */
#define  XC_LDA_X_REL                      532 /* Slater exchange with relativistic corrections */
#define  XC_GGA_X_SG4                      533 /* Semiclassical GGA at fourth order */
#define  XC_GGA_C_SG4                      534 /* Semiclassical GGA at fourth order */
#define  XC_GGA_X_GG99                     535 /* Gilbert and Gill 1999 */
#define  XC_LDA_XC_1D_EHWLRG_1             536 /* LDA constructed from slab-like systems of 1 electron */
#define  XC_LDA_XC_1D_EHWLRG_2             537 /* LDA constructed from slab-like systems of 2 electrons */
#define  XC_LDA_XC_1D_EHWLRG_3             538 /* LDA constructed from slab-like systems of 3 electrons */
#define  XC_GGA_X_PBEPOW                   539 /* PBE power */
#define  XC_MGGA_X_TM                      540 /* Tao and Mo 2016 exchange */
#define  XC_MGGA_X_VT84                    541 /* meta-GGA version of VT{8,4} GGA */
#define  XC_MGGA_X_SA_TPSS                 542 /* TPSS with correct surface asymptotics */
#define  XC_MGGA_K_PC07                    543 /* Perdew and Constantin 2007 */
#define  XC_GGA_X_KGG99                    544 /* Gilbert and Gill 1999 (mixed) */
#define  XC_GGA_XC_HLE16                   545 /* high local exchange 2016 */
#define  XC_LDA_X_ERF                      546 /* Short-range LDA exchange with error function kernel (erfc) */
#define  XC_LDA_XC_LP_A                    547 /* Lee-Parr reparametrization A */
#define  XC_LDA_XC_LP_B                    548 /* Lee-Parr reparametrization B */
#define  XC_LDA_X_RAE                      549 /* Rae self-energy corrected exchange */
#define  XC_LDA_K_ZLP                      550 /* Wigner including kinetic energy contribution */
#define  XC_LDA_C_MCWEENY                  551 /* McWeeny 76 */
#define  XC_LDA_C_BR78                     552 /* Brual & Rothstein 78 */
#define  XC_GGA_C_SCAN_E0                  553 /* GGA component of SCAN */
#define  XC_LDA_C_PK09                     554 /* Proynov and Kong 2009 */
#define  XC_GGA_C_GAPC                     555 /* GapC */
#define  XC_GGA_C_GAPLOC                   556 /* Gaploc */
#define  XC_GGA_C_ZVPBEINT                 557 /* another spin-dependent correction to PBEint */
#define  XC_GGA_C_ZVPBESOL                 558 /* another spin-dependent correction to PBEsol */
#define  XC_GGA_C_TM_LYP                   559 /* Takkar and McCarthy reparametrization, also known as reLYP */
#define  XC_GGA_C_TM_PBE                   560 /* Thakkar and McCarthy reparametrization */
#define  XC_GGA_C_W94                      561 /* Wilson 94 (Eq. 25) */
#define  XC_MGGA_C_KCIS                    562 /* Krieger, Chen, Iafrate, and Savin */
#define  XC_HYB_MGGA_XC_B0KCIS             563 /* Hybrid based on KCIS */
#define  XC_MGGA_XC_LP90                   564 /* Lee & Parr, Eq. (56) */
#define  XC_GGA_C_CS1                      565 /* A dynamical correlation functional */
#define  XC_HYB_MGGA_XC_MPW1KCIS           566 /* MPW1KCIS for barrier heights */
#define  XC_HYB_MGGA_XC_MPWKCIS1K          567 /* MPWKCIS1K for barrier heights */
#define  XC_HYB_MGGA_XC_PBE1KCIS           568 /* PBE1KCIS for binding energies */
#define  XC_HYB_MGGA_XC_TPSS1KCIS          569 /* TPSS1KCIS for thermochemistry and kinetics */
#define  XC_GGA_X_B88M                     570 /* Becke 88 reoptimized to be used with tau1 */
#define  XC_MGGA_C_B88                     571 /* Meta-GGA correlation by Becke */
#define  XC_HYB_GGA_XC_B5050LYP            572 /* B5050LYP */
#define  XC_LDA_C_OW_LYP                   573 /* Wigner with corresponding LYP parameters */
#define  XC_LDA_C_OW                       574 /* Optimized Wigner */
#define  XC_MGGA_X_GX                      575 /* GX functional of Loos */
#define  XC_MGGA_X_PBE_GX                  576 /* PBE-GX functional of Loos */
#define  XC_LDA_XC_GDSMFB                  577 /* Groth, Dornheim, Sjostrom, Malone, Foulkes, Bonitz */
#define  XC_LDA_C_GK72                     578 /* Gordon and Kim 1972 */
#define  XC_LDA_C_KARASIEV                 579 /* Karasiev reparameterization of Chachiyo */
#define  XC_LDA_K_LP96                     580 /* Liu-Parr kinetic */
#define  XC_MGGA_X_REVSCAN                 581 /* revised SCAN */
#define  XC_MGGA_C_REVSCAN                 582 /* revised SCAN */
#define  XC_HYB_MGGA_X_REVSCAN0            583 /* revised SCAN hybrid exchange (SCAN0) */
#define  XC_MGGA_C_SCAN_VV10               584 /* SCAN + VV10 correlation */
#define  XC_MGGA_C_REVSCAN_VV10            585 /* REVSCAN + VV10 correlation */
#define  XC_MGGA_X_BR89_EXPLICIT           586 /* Becke-Roussel 89 with an explicit inversion of x(y), gamma = 0.8 */
#define  XC_GGA_XC_KT3                     587 /* Keal and Tozer, version 3 */
#define  XC_HYB_LDA_XC_BN05                588 /* Baer and Neuhauser, gamma=1 */
#define  XC_HYB_GGA_XC_LB07                589 /* Livshits and Baer, empirical functional also used for IP tuning */
#define  XC_LDA_C_PMGB06                   590 /* Long-range LDA correlation functional */
#define  XC_GGA_K_GDS08                    591 /* Combined analytical theory with Monte Carlo sampling */
#define  XC_GGA_K_GHDS10                   592 /* As GDS08 but for an electron gas with spin */
#define  XC_GGA_K_GHDS10R                  593 /* Reparametrized GHDS10 */
#define  XC_GGA_K_TKVLN                    594 /* Trickey, Karasiev, and Vela */
#define  XC_GGA_K_PBE3                     595 /* Three parameter PBE-like expansion */
#define  XC_GGA_K_PBE4                     596 /* Four parameter PBE-like expansion */
#define  XC_GGA_K_EXP4                     597 /* Intermediate form between PBE3 and PBE4 */
#define  XC_HYB_MGGA_XC_B98                598 /* Becke 98 */
#define  XC_LDA_XC_TIH                     599 /* Neural network LDA from Tozer et al */
#define  XC_LDA_X_1D_EXPONENTIAL           600 /* Exchange in 1D for an exponentially screened interaction */
#define  XC_GGA_X_SFAT_PBE                 601 /* Short-range recipe for PBE functional - Yukawa */
#define  XC_MGGA_X_BR89_EXPLICIT_1         602 /* Becke-Roussel 89 with an explicit inversion of x(y), gamma = 1.0 */
#define  XC_MGGA_X_REGTPSS                 603 /* Regularized TPSS */
#define  XC_GGA_X_FD_LB94                  604 /* Functional derivative recovered from the stray LB94 potential */
#define  XC_GGA_X_FD_REVLB94               605 /* Revised FD_LB94 */
#define  XC_GGA_C_ZVPBELOC                 606 /* PBEloc variation with enhanced compatibility with exact exchange */
#define  XC_HYB_GGA_XC_APBE0               607 /* Hybrid based on APBE */
#define  XC_HYB_GGA_XC_HAPBE               608 /* Hybrid based in APBE and zvPBEloc */
#define  XC_MGGA_X_2D_JS17                 609 /* JS17 meta-GGA for 2D */
#define  XC_HYB_GGA_XC_RCAM_B3LYP          610 /* Similar to CAM-B3LYP, but trying to reduce the many-electron self-interaction */
#define  XC_HYB_GGA_XC_WC04                611 /* hybrid fitted to carbon NMR shifts */
#define  XC_HYB_GGA_XC_WP04                612 /* hybrid fitted to proton NMR shifts */
#define  XC_GGA_K_LKT                      613 /* Luo-Karasiev-Trickey GGA kinetic */
#define  XC_HYB_GGA_XC_CAMH_B3LYP          614 /* CAM version of B3LYP, tuned for TDDFT */
#define  XC_HYB_GGA_XC_WHPBE0              615 /* Long-range corrected short-range hybrid PBE (whPBE0) by Shao et al */
#define  XC_GGA_K_PBE2                     616 /* Three parameter PBE-like expansion */
#define  XC_MGGA_K_L04                     617 /* L0.4 by Laricchia et al */
#define  XC_MGGA_K_L06                     618 /* L0.6 by Laricchia et al */
#define  XC_GGA_K_VT84F                    619 /* VT84F by Karasiev et al */
#define  XC_GGA_K_LGAP                     620 /* LGAP by Constantin et al */
#define  XC_MGGA_K_RDA                     621 /* Reduced derivative approximation by Karasiev et al */
#define  XC_GGA_X_ITYH_OPTX                622 /* Short-range recipe for OPTX functional */
#define  XC_GGA_X_ITYH_PBE                 623 /* Short-range recipe for PBE functional */
#define  XC_GGA_C_LYPR                     624 /* Short-range LYP by Ai, Fang, and Su */
#define  XC_HYB_GGA_XC_LC_BLYP_EA          625 /* LC version of BLYP for electron affinities */
#define  XC_MGGA_X_REGTM                   626 /* Regularized Tao and Mo exchange */
#define  XC_MGGA_K_GEA2                    627 /* Second-order gradient expansion */
#define  XC_MGGA_K_GEA4                    628 /* Fourth-order gradient expansion */
#define  XC_MGGA_K_CSK1                    629 /* mGGA-rev functional by Cancio, Stewart, and Kuna (a=1) */
#define  XC_MGGA_K_CSK4                    630 /* mGGA-rev functional by Cancio, Stewart, and Kuna (a=4) */
#define  XC_MGGA_K_CSK_LOC1                631 /* mGGAloc-rev functional by Cancio, Stewart, and Kuna (a=1) */
#define  XC_MGGA_K_CSK_LOC4                632 /* mGGAloc-rev functional by Cancio, Stewart, and Kuna (a=4) */
#define  XC_GGA_K_LGAP_GE                  633 /* LGAP-GE by Constantin et al */
#define  XC_MGGA_K_PC07_OPT                634 /* Reoptimized PC07 by Mejia-Rodriguez and Trickey */
#define  XC_GGA_K_TFVW_OPT                 635 /* empirically optimized gamma-TFvW form */
#define  XC_HYB_GGA_XC_LC_BOP              636 /* LC version of B88 */
#define  XC_HYB_GGA_XC_LC_PBEOP            637 /* LC version of PBE */
#define  XC_MGGA_C_KCISK                   638 /* Krieger, Chen, and Kurth */
#define  XC_HYB_GGA_XC_LC_BLYPR            639 /* LC version of BLYP with correlation only in the short range */
#define  XC_HYB_GGA_XC_MCAM_B3LYP          640 /* Modified CAM-B3LYP by Day, Nguyen and Pachter */
#define  XC_LDA_X_YUKAWA                   641 /* Short-range LDA exchange with Yukawa attenuation */
#define  XC_MGGA_C_R2SCAN01                642 /* Re-regularized SCAN correlation with larger value for eta */
#define  XC_MGGA_C_RMGGAC                  643 /* Revised correlation energy for MGGAC exchange functional */
#define  XC_MGGA_X_MCML                    644 /* MCML exchange */
#define  XC_MGGA_X_R2SCAN01                645 /* Re-regularized SCAN exchange by Furness et al with larger value for eta */
#define  XC_HYB_GGA_X_CAM_S12G             646 /* Swart 2012 range-separated hybrid GGA exchange */
#define  XC_HYB_GGA_X_CAM_S12H             647 /* Swart 2012 range-separated hybrid GGA exchange */
#define  XC_MGGA_X_RPPSCAN                 648 /* r++SCAN: rSCAN with uniform density limit and coordinate scaling behavior */
#define  XC_MGGA_C_RPPSCAN                 649 /* r++SCAN: rSCAN with uniform density limit and coordinate scaling behavior */
#define  XC_MGGA_X_R4SCAN                  650 /* r$^{4}$SCAN, a functional that satisfies the same exact constraints that SCAN does */
#define  XC_MGGA_X_VCML                    651 /* Exchange part of VCML-rVV10 by Trepte and Voss */
#define  XC_MGGA_XC_VCML_RVV10             652 /* VCML-rVV10 by Trepte and Voss */
#define  XC_HYB_LDA_X_ERF                  653 /* Long-range corrected functional based on short-range LDA exchange (erfc) */
#define  XC_LDA_C_PW_ERF                   654 /* Short ranged correlation LDA (erfc) */
#define  XC_GGA_X_PBE_ERF_GWS              655 /* Short ranged PBE exchange (erfc) */
#define  XC_HYB_GGA_X_PBE_ERF_GWS          656 /* Short-range PBE (GWS) exchange (erfc) + long-range exact exchange */
#define  XC_GGA_C_PBE_ERF_GWS              657 /* Short ranged PBE correlation (erfc) */
#define  XC_HYB_MGGA_XC_GAS22              658 /* Google Accelerated Science 22 */
#define  XC_HYB_MGGA_XC_R2SCANH            659 /* r2SCANh: r2SCAN hybrid like TPSSh with 10% exact exchange */
#define  XC_HYB_MGGA_XC_R2SCAN0            660 /* r2SCAN0: r2SCAN hybrid like PBE0 with 25% exact exchange */
#define  XC_HYB_MGGA_XC_R2SCAN50           661 /* r2SCAN50: r2SCAN hybrid like PBE50 with 50% exact exchange */
#define  XC_HYB_GGA_XC_CAM_PBEH            681 /* CAM hybrid screened exchange PBE version */
#define  XC_HYB_GGA_XC_CAMY_PBEH           682 /* CAMY hybrid screened exchange PBE version */
#define  XC_LDA_C_UPW92                    683 /* Ruggeri, Rios, and Alavi unrestricted fit */
#define  XC_LDA_C_RPW92                    684 /* Ruggeri, Rios, and Alavi restricted fit */
#define  XC_MGGA_X_TLDA                    685 /* LDA-type exchange with tau-dependent potential */
#define  XC_MGGA_X_EDMGGA                  686 /* Tao 2001 */
#define  XC_MGGA_X_GDME_NV                 687 /* Generalized density-matrix with a=1/2 */
#define  XC_MGGA_X_RLDA                    688 /* Reparametrized local-density approximation */
#define  XC_MGGA_X_GDME_0                  689 /* Generalized density-matrix with a=0 */
#define  XC_MGGA_X_GDME_KOS                690 /* Generalized density-matrix with a=0.00638 */
#define  XC_MGGA_X_GDME_VT                 691 /* Varied-terms (VT) mGGA of Koehl, Odom, and Scuseria */
#define  XC_LDA_X_SLOC                     692 /* simple local model for Slater potential */
#define  XC_MGGA_X_REVTM                   693 /* revised Tao and Mo 2016 exchange */
#define  XC_MGGA_C_REVTM                   694 /* revised Tao and Mo 2016 exchange */
#define  XC_HYB_MGGA_XC_EDMGGAH            695 /* EDMGGA hybrid */
#define  XC_MGGA_X_MBRXC_BG                696 /* Modified Becke-Roussel for band gaps - cuspless hole */
#define  XC_MGGA_X_MBRXH_BG                697 /* Modified Becke-Roussel for band gaps - hydrogen hole */
#define  XC_MGGA_X_HLTA                    698 /* Half-and-half meta-LDAized LDA exchange by Lehtola and Marques */
#define  XC_MGGA_C_HLTAPW                  699 /* Half-and-half meta-LDAized PW correlation by Lehtola and Marques */
#define  XC_MGGA_X_SCANL                   700 /* Deorbitalized SCAN (SCAN-L) exchange */
#define  XC_MGGA_X_REVSCANL                701 /* Deorbitalized revised SCAN (revSCAN-L) exchange */
#define  XC_MGGA_C_SCANL                   702 /* Deorbitalized SCAN (SCAN-L) correlation */
#define  XC_MGGA_C_SCANL_RVV10             703 /* SCAN-L + rVV10 correlation */
#define  XC_MGGA_C_SCANL_VV10              704 /* SCAN-L + VV10 correlation */
#define  XC_HYB_MGGA_X_JS18                705 /* Jana and Samal 2018, screened range-separated TM exchange */
#define  XC_HYB_MGGA_X_PJS18               706 /* Patra, Jana and Samal 2018, screened range-separated TM exchange */
#define  XC_MGGA_X_TASK                    707 /* TASK exchange of Aschebrock and Kuemmel */
#define  XC_MGGA_X_MGGAC                   711 /* MGGAC exchange of Patra et al */
#define  XC_GGA_C_MGGAC                    712 /* beta fitted to LC20 to be used with MGGAC */
#define  XC_MGGA_X_MBR                     716 /* modified Becke-Roussel by Patra et al */
#define  XC_MGGA_X_R2SCANL                 718 /* Deorbitalized re-regularized SCAN (r2SCAN-L) exchange */
#define  XC_MGGA_C_R2SCANL                 719 /* Deorbitalized re-regularized SCAN (r2SCAN-L) correlation */
#define  XC_HYB_MGGA_XC_LC_TMLYP           720 /* Long-range corrected TM-LYP by Jana et al */
#define  XC_MGGA_X_MTASK                   724 /* modified TASK exchange */
#define  XC_GGA_X_Q1D                      734 /* Functional for quasi-1D systems */
#define  XC_MGGA_X_KTBM_0                  735 /* KTBM learned exchange - 0 */
#define  XC_MGGA_X_KTBM_1                  736 /* KTBM learned exchange - 1 */
#define  XC_MGGA_X_KTBM_2                  737 /* KTBM learned exchange - 2 */
#define  XC_MGGA_X_KTBM_3                  738 /* KTBM learned exchange - 3 */
#define  XC_MGGA_X_KTBM_4                  739 /* KTBM learned exchange - 4 */
#define  XC_MGGA_X_KTBM_5                  740 /* KTBM learned exchange - 5 */
#define  XC_MGGA_X_KTBM_6                  741 /* KTBM learned exchange - 6 */
#define  XC_MGGA_X_KTBM_7                  742 /* KTBM learned exchange - 7 */
#define  XC_MGGA_X_KTBM_8                  743 /* KTBM learned exchange - 8 */
#define  XC_MGGA_X_KTBM_9                  744 /* KTBM learned exchange - 9 */
#define  XC_MGGA_X_KTBM_10                 745 /* KTBM learned exchange - 10 */
#define  XC_MGGA_X_KTBM_11                 746 /* KTBM learned exchange - 11 */
#define  XC_MGGA_X_KTBM_12                 747 /* KTBM learned exchange - 12 */
#define  XC_MGGA_X_KTBM_13                 748 /* KTBM learned exchange - 13 */
#define  XC_MGGA_X_KTBM_14                 749 /* KTBM learned exchange - 14 */
#define  XC_MGGA_X_KTBM_15                 750 /* KTBM learned exchange - 15 */
#define  XC_MGGA_X_KTBM_16                 751 /* KTBM learned exchange - 16 */
#define  XC_MGGA_X_KTBM_17                 752 /* KTBM learned exchange - 17 */
#define  XC_MGGA_X_KTBM_18                 753 /* KTBM learned exchange - 18 */
#define  XC_MGGA_X_KTBM_19                 754 /* KTBM learned exchange - 19 */
#define  XC_MGGA_X_KTBM_20                 755 /* KTBM learned exchange - 20 */
#define  XC_MGGA_X_KTBM_21                 756 /* KTBM learned exchange - 21 */
#define  XC_MGGA_X_KTBM_22                 757 /* KTBM learned exchange - 22 */
#define  XC_MGGA_X_KTBM_23                 758 /* KTBM learned exchange - 23 */
#define  XC_MGGA_X_KTBM_24                 759 /* KTBM learned exchange - 24 */
#define  XC_MGGA_X_KTBM_GAP                760 /* KTBM learned exchange - GAP */
#define  XC_MGGA_X_MSPBEL  		   761 /* MS-PBEl, a PBE-like meta-GGA exchange */
#define  XC_MGGA_X_RMSPBEL  		   762 /* regularized MS-PBEl */
#define  XC_MGGA_X_MSRPBEL  		   763 /* MS-RPBEl, a RPBE-like meta-GGA exchange */
#define  XC_MGGA_X_RMSRPBEL  		   764 /* regularized MS-RPBEl */
#define  XC_MGGA_X_MSB86BL  		   765 /* MS-B86bl, a B86b-like meta-GGA exchange */
#define  XC_MGGA_X_RMSB86BL  		   766 /* regularized MS-B86bl */
