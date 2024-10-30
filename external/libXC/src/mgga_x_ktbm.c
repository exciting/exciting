/*
  Copyright (C) 2016 Susi Lehtola

  This Source Code Form is subject to the terms of the Mozilla Public
  License, v. 2.0. If a copy of the MPL was not distributed with this
  file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "util.h"

#define XC_MGGA_X_KTBM_0 735 /* KTBM learned 0*/
#define XC_MGGA_X_KTBM_1 736 /* KTBM learned 1*/
#define XC_MGGA_X_KTBM_2 737 /* KTBM learned 2*/
#define XC_MGGA_X_KTBM_3 738 /* KTBM learned 3*/
#define XC_MGGA_X_KTBM_4 739 /* KTBM learned 4*/
#define XC_MGGA_X_KTBM_5 740 /* KTBM learned 5*/
#define XC_MGGA_X_KTBM_6 741 /* KTBM learned 6*/
#define XC_MGGA_X_KTBM_7 742 /* KTBM learned 7*/
#define XC_MGGA_X_KTBM_8 743 /* KTBM learned 8*/
#define XC_MGGA_X_KTBM_9 744 /* KTBM learned 9*/
#define XC_MGGA_X_KTBM_10 745 /* KTBM learned 10*/
#define XC_MGGA_X_KTBM_11 746 /* KTBM learned 11*/
#define XC_MGGA_X_KTBM_12 747 /* KTBM learned 12*/
#define XC_MGGA_X_KTBM_13 748 /* KTBM learned 13*/
#define XC_MGGA_X_KTBM_14 749 /* KTBM learned 14*/
#define XC_MGGA_X_KTBM_15 750 /* KTBM learned 15*/
#define XC_MGGA_X_KTBM_16 751 /* KTBM learned 16*/
#define XC_MGGA_X_KTBM_17 752 /* KTBM learned 17*/
#define XC_MGGA_X_KTBM_18 753 /* KTBM learned 18*/
#define XC_MGGA_X_KTBM_19 754 /* KTBM learned 19*/
#define XC_MGGA_X_KTBM_20 755 /* KTBM learned 20*/
#define XC_MGGA_X_KTBM_21 756 /* KTBM learned 21*/
#define XC_MGGA_X_KTBM_22 757 /* KTBM learned 22*/
#define XC_MGGA_X_KTBM_23 758 /* KTBM learned 23*/
#define XC_MGGA_X_KTBM_24 759 /* KTBM learned 24*/
#define XC_MGGA_X_KTBM_GAP 760 /* KTBM learned for gaps*/

typedef struct{
  double ct, at, bt, a2t, b2t, xt, cb, ab, bb, a2b, b2b, xb;
} mgga_x_ktbm_params;

#define N_PAR_KTBM 12
static const char *ktbm_names[N_PAR_KTBM] = {"ct", "at", "bt", "a2t", "b2t", "xt", "cb", "ab", "bb", "a2b", "b2b", "xb"};
static const char *ktbm_desc[N_PAR_KTBM] = {"ct", "at", "bt", "a2t", "b2t", "xt", "cb", "ab", "bb", "a2b", "b2b", "xb"};

static const double par_ktbm0[N_PAR_KTBM] = { 1.0304076402855855, 0.29069155191166346, 0.09376131500450102, 0.00636399057259787, 0.2338605637430116, 0.05905616176610306, 0.912960806906729, -0.10143676314780245, 0.23304589272011897, 0.012695765733790497, 0.2120228194062504, 0.030228093216799572};
static const double par_ktbm1[N_PAR_KTBM] = { 1.0263752715156753, 0.2899411592708853, 0.09527007147746959, 0.0068639720458172964, 0.23290915428403358, 0.05947813667117197, 0.9092803457608484, -0.09990334644850286, 0.23248452710436812, 0.012457017995909028, 0.21278962441196217, 0.029719442546683742};
static const double par_ktbm2[N_PAR_KTBM] = { 1.0295936695180066, 0.29058410317853406, 0.0978868238219999, 0.007755128286066295, 0.23023677039428692, 0.06068632690160595, 0.9112142206259968, -0.09977324697775652, 0.22943895674906856, 0.011998579529446734, 0.21706408635922794, 0.028442891502834692};
static const double par_ktbm3[N_PAR_KTBM] = { 1.0379811625637196, 0.28263576038064064, 0.1086528703236551, 0.02256386879999913, 0.21951785356175765, 0.05535090451568734, 0.9062056523248283, -0.09530345147869335, 0.22208016646397852, -0.004498791216677266, 0.23786606766032575, 0.03123136106174504};
static const double par_ktbm4[N_PAR_KTBM] = { 1.0447255092727967, 0.21643762279081258, 0.11205441240705669, 0.03916064018939278, 0.218782148009531, 0.049173374763561505, 0.9111626565768056, -0.03878570821939097, 0.2223103939373648, -0.02461830223135286, 0.2420890191752141, 0.03357824594821441};
static const double par_ktbm5[N_PAR_KTBM] = { 1.0125602773047306, 0.2675776296196626, 0.09913936134558513, 0.023315761058191182, 0.24640960606519824, 0.053160334695376496, 0.9077280787511512, -0.07902833543247754, 0.23072509312940534, -0.00866483099804667, 0.2196560728349573, 0.03412753589546748};
static const double par_ktbm6[N_PAR_KTBM] = { 1.0107990187301998, 0.2831195698499447, 0.09419986422549517, 0.018516584221964772, 0.24643091238118955, 0.053222207515855054, 0.8927648294054159, -0.09162327508843797, 0.23505067511245098, -0.0038997050285408835, 0.2236142908190176, 0.03491565811591121};
static const double par_ktbm7[N_PAR_KTBM] = { 1.0247132197383149, 0.29135357108065124, 0.08982076954655123, 0.017553715624560488, 0.2389845612385746, 0.052873026336855665, 0.8796968988601903, -0.09955165062796137, 0.24024623770293027, -0.003155319989664977, 0.23357541396032006, 0.034965668999510056};
static const double par_ktbm8[N_PAR_KTBM] = { 1.0383504503761103, 0.29673875703834085, 0.10350110708250126, 0.025487794684463046, 0.22807264244165595, 0.046836991399129155, 0.8923983637114028, -0.10964831528933401, 0.22538195411613346, -0.014660499819366934, 0.2521438820727313, 0.03770293882466388};
static const double par_ktbm9[N_PAR_KTBM] = { 1.0423283694993333, 0.29520765557119005, 0.10931801428639576, 0.028253873924966876, 0.2193558860545767, 0.038107713046421104, 0.8856159764767498, -0.11135236999860085, 0.22146942269020795, -0.020769127848149577, 0.263916870673348, 0.04302665117695169};
static const double par_ktbm10[N_PAR_KTBM] = { 1.0333254850526004, 0.30912912129708214, 0.08488944601171224, 0.01761185351449011, 0.2556727204970185, 0.04186658887299335, 0.8899242792435056, -0.12160989477866269, 0.25118301876427696, -0.010401252152213812, 0.23278035355354837, 0.04218521572257066};
static const double par_ktbm11[N_PAR_KTBM] = { 1.0372855363639735, 0.31542069088069363, 0.0846516506524305, 0.017493092157697374, 0.2568451832465198, 0.04166294140424092, 0.8852633605618385, -0.12654391613327828, 0.2529563857261435, -0.011153234776641766, 0.2405626239749416, 0.04247456513199206};
static const double par_ktbm12[N_PAR_KTBM] = { 1.0328992930671033, 0.32507944455155824, 0.08867013038404557, 0.01911189167574323, 0.25267236609857674, 0.03776897033885063, 0.8678685739668615, -0.13719467942712382, 0.25010201689281614, -0.014708120975318316, 0.2562711986900481, 0.044664146022155224};
static const double par_ktbm13[N_PAR_KTBM] = { 1.0224370487563148, 0.34088085237213217, 0.09899449457477784, 0.020209081896366766, 0.24135718732652808, 0.032017852777980665, 0.8526876327199046, -0.1558957723810475, 0.23686926225137386, -0.0188712715300085, 0.27323183568634224, 0.046972503043861545};
static const double par_ktbm14[N_PAR_KTBM] = { 1.0098268364593528, 0.3712618902901585, 0.11446342331107544, 0.019271713973666525, 0.22767099588537568, 0.0186613282501236, 0.8471031407239273, -0.18372997898145776, 0.22348685551621814, -0.022330162076107284, 0.2813712594156585, 0.05447119556682036};
static const double par_ktbm15[N_PAR_KTBM] = { 1.0302976207141643, 0.3355002717573432, 0.08130545277413916, 0.02080638881254861, 0.2711168697258796, 0.026644005087231207, 0.8620530992475534, -0.15393192486729823, 0.26741916973273905, -0.02308922048967118, 0.25324767423389083, 0.05061099133357915};
static const double par_ktbm16[N_PAR_KTBM] = { 1.028944534589083, 0.34785279555635895, 0.08267108106108388, 0.01826796628619182, 0.26926390331174666, 0.025488532321196775, 0.8517262574698339, -0.16440139191706984, 0.26668246758222863, -0.020841203345624244, 0.2624707939098512, 0.05140925571333686};
static const double par_ktbm17[N_PAR_KTBM] = { 1.0292756704894197, 0.35683459177209, 0.09012873503547454, 0.02049145174816412, 0.25979318901002013, 0.020938264897749744, 0.8460317965240975, -0.17690897486541043, 0.2573678273967677, -0.02525015720470474, 0.275797970614049, 0.05316989878402506};
static const double par_ktbm18[N_PAR_KTBM] = { 1.023891416652007, 0.3695840036000602, 0.10424651498067275, 0.021389207645624835, 0.24635839037940996, 0.016303682950443012, 0.8463238571977623, -0.19156566837384478, 0.24286563818174003, -0.02916253719139352, 0.2853068266325876, 0.05439786485037524};
static const double par_ktbm19[N_PAR_KTBM] = { 1.0177509335825339, 0.3886482506196606, 0.11958299056773851, 0.01734660257009398, 0.23364313875476078, 0.013006344527878945, 0.8500063400271057, -0.19970878946385107, 0.23265545875876809, -0.02652803375774664, 0.2883152641191594, 0.05563460230305794};
static const double par_ktbm20[N_PAR_KTBM] = { 1.0165187620195661, 0.3650072557466966, 0.08954496562037528, 0.014607712335708758, 0.2632247589374224, 0.01786607619990434, 0.8287842058666036, -0.1835191673802379, 0.26577701274857213, -0.022982606699085847, 0.27472726796218805, 0.05379265178854256};
static const double par_ktbm21[N_PAR_KTBM] = { 1.0261589511390368, 0.37178251147192115, 0.09417966281307585, 0.017187315664980975, 0.25854265737744203, 0.014141308840385008, 0.8369879914727829, -0.19334218741581838, 0.261609289869405, -0.027512731460753635, 0.28028398998736687, 0.05529346995551265};
static const double par_ktbm22[N_PAR_KTBM] = { 1.0314315558545104, 0.3753257423779244, 0.09970148467691746, 0.01845516499149008, 0.25282812241932845, 0.012331220162193516, 0.842951890093288, -0.1986302049502704, 0.25728153005661014, -0.030259689597518808, 0.2837277428008583, 0.05549196058928848};
static const double par_ktbm23[N_PAR_KTBM] = { 1.0295589921206336, 0.3834598731754155, 0.1134601725160894, 0.012475919846263022, 0.2414188532038928, 0.01732362062811347, 0.8493931293985715, -0.1999873028058896, 0.24921422369400187, -0.024245285244768995, 0.2858306647480424, 0.05112243276201573};
static const double par_ktbm24[N_PAR_KTBM] = { 1.0228413140883612, 0.3920773179999707, 0.13989350865860764, -0.03442867975153089, 0.20859113727030637, 0.05898500700761443, 0.8500009175699887, -0.16852485787284124, 0.23913780777316324, 0.02878611571628562, 0.2821872346741235, 0.027565652585576854};
static const double par_ktbmgap[N_PAR_KTBM] = {1.122142370218625, 0.4354020787275538, 0.181093570244884, -0.00672327251988611, 0.1285175238165055, 0.111067440834648, 0.9723866568932589, -0.17459201010028683, 0.15992157285736205, 0.07947893639728304, 0.29944523452939364, 0.012402724366764269};

static void
mgga_x_ktbm_init(xc_func_type *p)
{
  assert(p!=NULL && p->params == NULL);
  p->params = libxc_malloc(sizeof(mgga_x_ktbm_params));
}

#include "maple2c/mgga_exc/mgga_x_ktbm.c"
#include "work_mgga.c"

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_0 = {
  XC_MGGA_X_KTBM_0,
  XC_EXCHANGE,
  "KTBM learned exchange - 0",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm0, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_1 = {
  XC_MGGA_X_KTBM_1,
  XC_EXCHANGE,
  "KTBM learned exchange - 1",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm1, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_2 = {
  XC_MGGA_X_KTBM_2,
  XC_EXCHANGE,
  "KTBM learned exchange - 2",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm2, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_3 = {
  XC_MGGA_X_KTBM_3,
  XC_EXCHANGE,
  "KTBM learned exchange - 3",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm3, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_4 = {
  XC_MGGA_X_KTBM_4,
  XC_EXCHANGE,
  "KTBM learned exchange - 4",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm4, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_5 = {
  XC_MGGA_X_KTBM_5,
  XC_EXCHANGE,
  "KTBM learned exchange - 5",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm5, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_6 = {
  XC_MGGA_X_KTBM_6,
  XC_EXCHANGE,
  "KTBM learned exchange - 6",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm6, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_7 = {
  XC_MGGA_X_KTBM_7,
  XC_EXCHANGE,
  "KTBM learned exchange - 7",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm7, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_8 = {
  XC_MGGA_X_KTBM_8,
  XC_EXCHANGE,
  "KTBM learned exchange - 8",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm8, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_9 = {
  XC_MGGA_X_KTBM_9,
  XC_EXCHANGE,
  "KTBM learned exchange - 9",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm9, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_10 = {
  XC_MGGA_X_KTBM_10,
  XC_EXCHANGE,
  "KTBM learned exchange - 10",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm10, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_11 = {
  XC_MGGA_X_KTBM_11,
  XC_EXCHANGE,
  "KTBM learned exchange - 11",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm11, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_12 = {
  XC_MGGA_X_KTBM_12,
  XC_EXCHANGE,
  "KTBM learned exchange - 12",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm12, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_13 = {
  XC_MGGA_X_KTBM_13,
  XC_EXCHANGE,
  "KTBM learned exchange - 13",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm13, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_14 = {
  XC_MGGA_X_KTBM_14,
  XC_EXCHANGE,
  "KTBM learned exchange - 14",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm14, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_15 = {
  XC_MGGA_X_KTBM_15,
  XC_EXCHANGE,
  "KTBM learned exchange - 15",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm15, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_16 = {
  XC_MGGA_X_KTBM_16,
  XC_EXCHANGE,
  "KTBM learned exchange - 16",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm16, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_17 = {
  XC_MGGA_X_KTBM_17,
  XC_EXCHANGE,
  "KTBM learned exchange - 17",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm17, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_18 = {
  XC_MGGA_X_KTBM_18,
  XC_EXCHANGE,
  "KTBM learned exchange - 18",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm18, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_19 = {
  XC_MGGA_X_KTBM_19,
  XC_EXCHANGE,
  "KTBM learned exchange - 19",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm19, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_20 = {
  XC_MGGA_X_KTBM_20,
  XC_EXCHANGE,
  "KTBM learned exchange - 20",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm20, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_21 = {
  XC_MGGA_X_KTBM_21,
  XC_EXCHANGE,
  "KTBM learned exchange - 21",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm21, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_22 = {
  XC_MGGA_X_KTBM_22,
  XC_EXCHANGE,
  "KTBM learned exchange - 22",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm22, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_23 = {
  XC_MGGA_X_KTBM_23,
  XC_EXCHANGE,
  "KTBM learned exchange - 23",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm23, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_24 = {
  XC_MGGA_X_KTBM_24,
  XC_EXCHANGE,
  "KTBM learned exchange - 24",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbm24, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

#ifdef __cplusplus
extern "C"
#endif
const xc_func_info_type xc_func_info_mgga_x_ktbm_gap = {
  XC_MGGA_X_KTBM_GAP,
  XC_EXCHANGE,
  "KTBM learned exchange - GAP",
  XC_FAMILY_MGGA,
  {&xc_ref_Kovacs2022_094110, NULL, NULL, NULL, NULL},
  XC_FLAGS_3D | XC_FLAGS_NEEDS_TAU | MAPLE2C_FLAGS,
  1e-15,
  {N_PAR_KTBM, ktbm_names, ktbm_desc, par_ktbmgap, set_ext_params_cpy},
  mgga_x_ktbm_init, NULL,
  NULL, NULL, &work_mgga
};

