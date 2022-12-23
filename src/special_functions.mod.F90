#include "cpmd_global.h"

MODULE special_functions
  USE kinds,                           ONLY: real_8

  USE, INTRINSIC :: iso_c_binding, ONLY : c_double

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: cp_erf
  PUBLIC :: cp_erfc

CONTAINS

  FUNCTION cp_erf(x) RESULT(my_res)
    REAL(real_8), INTENT(in)                 :: x
    REAL(real_8)                             :: my_res

#if defined(_HAS_EXTERNAL_C_ERF)
    ! use external c funciton
    INTERFACE
       REAL(c_double) FUNCTION ERF(x) BIND(c,name='erf')
         IMPORT :: c_double
         REAL(c_double), VALUE, INTENT(in) :: x
       END FUNCTION erf
    END INTERFACE
    my_res=ERF(REAL(x,c_double))
#else
    ! use default f08 function
    my_res=ERF(x)
#endif
  END FUNCTION cp_erf

  FUNCTION cp_erfc(x) RESULT(my_res)
    REAL(real_8), INTENT(in)                 :: x
    REAL(real_8)                             :: my_res

#if defined(_HAS_EXTERNAL_C_ERF)
    ! use external c funciton
    INTERFACE
       REAL(c_double) FUNCTION ERFC(x) BIND(c,name='erfc')
         IMPORT :: c_double
         REAL(c_double), VALUE, INTENT(in) :: x
       END FUNCTION erfc
    END INTERFACE
    my_res=ERFC(REAL(x,c_double))
#else
    ! use default f08 function
    my_res=ERFC(x)
#endif
  END FUNCTION cp_erfc

END MODULE special_functions


! *DECK DCSEVL
FUNCTION dcsevl (x, cs, n) RESULT(my_res)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  IMPLICIT NONE
  REAL(real_8)                               :: x, cs(*)
  INTEGER                                    :: n
  REAL(real_8)                               :: my_res

  INTEGER                                    :: i, ni
  LOGICAL, SAVE                              :: first = .TRUE.
  REAL(real_8)                               :: b0, b1, b2, d1mach, twox
  REAL(real_8), SAVE                         :: onepl

! ***BEGIN PROLOGUE  DCSEVL
! ***PURPOSE  Evaluate a Chebyshev series.
! ***LIBRARY   SLATEC (FNLIB)
! ***CATEGORY  C3A2
! ***TYPE      DOUBLE PRECISION (CSEVL-S, DCSEVL-D)
! ***KEYWORDS  CHEBYSHEV SERIES, FNLIB, SPECIAL FUNCTIONS
! ***AUTHOR  Fullerton, W., (LANL)
! ***DESCRIPTION
! 
! Evaluate the N-term Chebyshev series CS at X.  Adapted from
! a method presented in the paper by Broucke referenced below.
! 
! Input Arguments --
! X    value at which the series is to be evaluated.
! CS   array of N terms of a Chebyshev series.  In evaluating
! CS, only half the first coefficient is summed.
! N    number of terms in array CS.
! 
! ***REFERENCES  R. Broucke, Ten subroutines for the manipulation of
! Chebyshev series, Algorithm 446, Communications of
! the A.C.M. 16, (1973) pp. 254-256.
! L. Fox and I. B. Parker, Chebyshev Polynomials in
! Numerical Analysis, Oxford University Press, 1968,
! page 56.
! ***ROUTINES CALLED  D1MACH, XERMSG
! ***REVISION HISTORY  (YYMMDD)
! 770401  DATE WRITTEN
! 890831  Modified array declarations.  (WRB)
! 890831  REVISION DATE from Version 3.2
! 891214  Prologue converted to Version 4.0 format.  (BAB)
! 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
! 900329  Prologued revised extensively and code rewritten to allow
! X to be slightly outside interval (-1,+1).  (WRB)
! 920501  Reformatted the REFERENCES section.  (WRB)
! ***END PROLOGUE  DCSEVL
! ***FIRST EXECUTABLE STATEMENT  DCSEVL

  IF (first) onepl = 1.0_real_8 + d1mach(4)
  first = .FALSE.
  IF (n .LT. 1) STOP 100 ! CALL XERMSG ('SLATEC', 'DCSEVL',
  ! +   'NUMBER OF TERMS .LE. 0', 2, 2)
  IF (n .GT. 1000) STOP 101 ! CALL XERMSG ('SLATEC', 'DCSEVL',
  ! +   'NUMBER OF TERMS .GT. 1000', 3, 2)
  IF (ABS(x) .GT. onepl) STOP 102 ! CALL XERMSG ('SLATEC', 'DCSEVL',
  ! +   'X OUTSIDE THE INTERVAL (-1,+1)', 1, 1)
  ! 
  b1 = 0.0_real_8
  b0 = 0.0_real_8
  twox = 2.0_real_8*x
  DO i = 1,n
     b2 = b1
     b1 = b0
     ni = n + 1 - i
     b0 = twox*b1 - b2 + cs(ni)
  END DO
  ! 
  my_res = 0.5_real_8*(b0-b2)
  ! 
  RETURN
END FUNCTION dcsevl
! *DECK DE1
FUNCTION de1 (x) RESULT(my_res)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  REAL(real_8)                               :: x, my_res

  INTEGER, EXTERNAL                          :: initds
  INTEGER, SAVE                              :: ntae10, ntae11, ntae12, &
                                                ntae13, ntae14, nte11, nte12
  LOGICAL, SAVE                              :: first = .TRUE.
  REAL(real_8)                               :: d1mach, dcsevl, eta, xmaxt
  REAL(real_8), DIMENSION(25), SAVE :: e12cs = (/ &
      -0.3739021479220279511668698204827e-1_real_8,&
      +0.4272398606220957726049179176528e-1_real_8,&
      -0.130318207984970054415392055219726e+0_real_8,&
      +0.144191240246988907341095893982137e-1_real_8,&
      -0.134617078051068022116121527983553e-2_real_8,&
      +0.107310292530637799976115850970073e-3_real_8,&
      -0.742999951611943649610283062223163e-5_real_8,&
      +0.453773256907537139386383211511827e-6_real_8,&
      -0.247641721139060131846547423802912e-7_real_8,&
      +0.122076581374590953700228167846102e-8_real_8,&
      -0.548514148064092393821357398028261e-10_real_8,&
      +0.226362142130078799293688162377002e-11_real_8,&
      -0.863589727169800979404172916282240e-13_real_8,&
      +0.306291553669332997581032894881279e-14_real_8,&
      -0.101485718855944147557128906734933e-15_real_8,&
      +0.315482174034069877546855328426666e-17_real_8,&
      -0.923604240769240954484015923200000e-19_real_8,&
      +0.255504267970814002440435029333333e-20_real_8,&
      -0.669912805684566847217882453333333e-22_real_8,&
      +0.166925405435387319431987199999999e-23_real_8,&
      -0.396254925184379641856000000000000e-25_real_8,&
      +0.898135896598511332010666666666666e-27_real_8,&
      -0.194763366993016433322666666666666e-28_real_8,&
      +0.404836019024630033066666666666666e-30_real_8,&
      -0.807981567699845120000000000000000e-32_real_8/)
  REAL(real_8), DIMENSION(29), SAVE :: e11cs = (/ &
      -0.16113461655571494025720663927566180e+2_real_8,&
      +0.77940727787426802769272245891741497e+1_real_8,&
      -0.19554058188631419507127283812814491e+1_real_8,&
      +0.37337293866277945611517190865690209e+0_real_8,&
      -0.56925031910929019385263892220051166e-1_real_8,&
      +0.72110777696600918537847724812635813e-2_real_8,&
      -0.78104901449841593997715184089064148e-3_real_8,&
      +0.73880933562621681878974881366177858e-4_real_8,&
      -0.62028618758082045134358133607909712e-5_real_8,&
      +0.46816002303176735524405823868362657e-6_real_8,&
      -0.32092888533298649524072553027228719e-7_real_8,&
      +0.20151997487404533394826262213019548e-8_real_8,&
      -0.11673686816697793105356271695015419e-9_real_8,&
      +0.62762706672039943397788748379615573e-11_real_8,&
      -0.31481541672275441045246781802393600e-12_real_8,&
      +0.14799041744493474210894472251733333e-13_real_8,&
      -0.65457091583979673774263401588053333e-15_real_8,&
      +0.27336872223137291142508012748799999e-16_real_8,&
      -0.10813524349754406876721727624533333e-17_real_8,&
      +0.40628328040434303295300348586666666e-19_real_8,&
      -0.14535539358960455858914372266666666e-20_real_8,&
      +0.49632746181648636830198442666666666e-22_real_8,&
      -0.16208612696636044604866560000000000e-23_real_8,&
      +0.50721448038607422226431999999999999e-25_real_8,&
      -0.15235811133372207813973333333333333e-26_real_8,&
      +0.44001511256103618696533333333333333e-28_real_8,&
      -0.12236141945416231594666666666666666e-29_real_8,&
      +0.32809216661066001066666666666666666e-31_real_8,&
      -0.84933452268306432000000000000000000e-33_real_8/)
  REAL(real_8), DIMENSION(41), SAVE :: ae12cs = (/ &
      +0.63629589796747038767129887806803e+0_real_8,&
      -0.13081168675067634385812671121135e+0_real_8,&
      -0.84367410213053930014487662129752e-2_real_8,&
      +0.26568491531006685413029428068906e-2_real_8,&
      +0.32822721781658133778792170142517e-3_real_8,&
      -0.23783447771430248269579807851050e-4_real_8,&
      -0.11439804308100055514447076797047e-4_real_8,&
      -0.14405943433238338455239717699323e-5_real_8,&
      +0.52415956651148829963772818061664e-8_real_8,&
      +0.38407306407844323480979203059716e-7_real_8,&
      +0.85880244860267195879660515759344e-8_real_8,&
      +0.10219226625855003286339969553911e-8_real_8,&
      +0.21749132323289724542821339805992e-10_real_8,&
      -0.22090238142623144809523503811741e-10_real_8,&
      -0.63457533544928753294383622208801e-11_real_8,&
      -0.10837746566857661115340539732919e-11_real_8,&
      -0.11909822872222586730262200440277e-12_real_8,&
      -0.28438682389265590299508766008661e-14_real_8,&
      +0.25080327026686769668587195487546e-14_real_8,&
      +0.78729641528559842431597726421265e-15_real_8,&
      +0.15475066347785217148484334637329e-15_real_8,&
      +0.22575322831665075055272608197290e-16_real_8,&
      +0.22233352867266608760281380836693e-17_real_8,&
      +0.16967819563544153513464194662399e-19_real_8,&
      -0.57608316255947682105310087304533e-19_real_8,&
      -0.17591235774646878055625369408853e-19_real_8,&
      -0.36286056375103174394755328682666e-20_real_8,&
      -0.59235569797328991652558143488000e-21_real_8,&
      -0.76030380926310191114429136895999e-22_real_8,&
      -0.62547843521711763842641428479999e-23_real_8,&
      +0.25483360759307648606037606400000e-24_real_8,&
      +0.25598615731739857020168874666666e-24_real_8,&
      +0.71376239357899318800207052800000e-25_real_8,&
      +0.14703759939567568181578956800000e-25_real_8,&
      +0.25105524765386733555198634666666e-26_real_8,&
      +0.35886666387790890886583637333333e-27_real_8,&
      +0.39886035156771301763317759999999e-28_real_8,&
      +0.21763676947356220478805333333333e-29_real_8,&
      -0.46146998487618942367607466666666e-30_real_8,&
      -0.20713517877481987707153066666666e-30_real_8,&
      -0.51890378563534371596970666666666e-31_real_8/)
  REAL(real_8), DIMENSION(50), SAVE :: ae10cs = (/ &
      +0.3284394579616699087873844201881e-1_real_8,&
      -0.1669920452031362851476184343387e-1_real_8,&
      +0.2845284724361346807424899853252e-3_real_8,&
      -0.7563944358516206489487866938533e-5_real_8,&
      +0.2798971289450859157504843180879e-6_real_8,&
      -0.1357901828534531069525563926255e-7_real_8,&
      +0.8343596202040469255856102904906e-9_real_8,&
      -0.6370971727640248438275242988532e-10_real_8,&
      +0.6007247608811861235760831561584e-11_real_8,&
      -0.7022876174679773590750626150088e-12_real_8,&
      +0.1018302673703687693096652346883e-12_real_8,&
      -0.1761812903430880040406309966422e-13_real_8,&
      +0.3250828614235360694244030353877e-14_real_8,&
      -0.5071770025505818678824872259044e-15_real_8,&
      +0.1665177387043294298172486084156e-16_real_8,&
      +0.3166753890797514400677003536555e-16_real_8,&
      -0.1588403763664141515133118343538e-16_real_8,&
      +0.4175513256138018833003034618484e-17_real_8,&
      -0.2892347749707141906710714478852e-18_real_8,&
      -0.2800625903396608103506340589669e-18_real_8,&
      +0.1322938639539270903707580023781e-18_real_8,&
      -0.1804447444177301627283887833557e-19_real_8,&
      -0.7905384086522616076291644817604e-20_real_8,&
      +0.4435711366369570103946235838027e-20_real_8,&
      -0.4264103994978120868865309206555e-21_real_8,&
      -0.3920101766937117541553713162048e-21_real_8,&
      +0.1527378051343994266343752326971e-21_real_8,&
      +0.1024849527049372339310308783117e-22_real_8,&
      -0.2134907874771433576262711405882e-22_real_8,&
      +0.3239139475160028267061694700366e-23_real_8,&
      +0.2142183762299889954762643168296e-23_real_8,&
      -0.8234609419601018414700348082312e-24_real_8,&
      -0.1524652829645809479613694401140e-24_real_8,&
      +0.1378208282460639134668480364325e-24_real_8,&
      +0.2131311202833947879523224999253e-26_real_8,&
      -0.2012649651526484121817466763127e-25_real_8,&
      +0.1995535662263358016106311782673e-26_real_8,&
      +0.2798995808984003464948686520319e-26_real_8,&
      -0.5534511845389626637640819277823e-27_real_8,&
      -0.3884995396159968861682544026146e-27_real_8,&
      +0.1121304434507359382850680354679e-27_real_8,&
      +0.5566568152423740948256563833514e-28_real_8,&
      -0.2045482929810499700448533938176e-28_real_8,&
      -0.8453813992712336233411457493674e-29_real_8,&
      +0.3565758433431291562816111116287e-29_real_8,&
      +0.1383653872125634705539949098871e-29_real_8,&
      -0.6062167864451372436584533764778e-30_real_8,&
      -0.2447198043989313267437655119189e-30_real_8,&
      +0.1006850640933998348011548180480e-30_real_8,&
      +0.4623685555014869015664341461674e-31_real_8/)
  REAL(real_8), DIMENSION(50), SAVE :: ae13cs = (/ &
      -0.60577324664060345999319382737747e+0_real_8,&
      -0.11253524348366090030649768852718e+0_real_8,&
      +0.13432266247902779492487859329414e-1_real_8,&
      -0.19268451873811457249246838991303e-2_real_8,&
      +0.30911833772060318335586737475368e-3_real_8,&
      -0.53564132129618418776393559795147e-4_real_8,&
      +0.98278128802474923952491882717237e-5_real_8,&
      -0.18853689849165182826902891938910e-5_real_8,&
      +0.37494319356894735406964042190531e-6_real_8,&
      -0.76823455870552639273733465680556e-7_real_8,&
      +0.16143270567198777552956300060868e-7_real_8,&
      -0.34668022114907354566309060226027e-8_real_8,&
      +0.75875420919036277572889747054114e-9_real_8,&
      -0.16886433329881412573514526636703e-9_real_8,&
      +0.38145706749552265682804250927272e-10_real_8,&
      -0.87330266324446292706851718272334e-11_real_8,&
      +0.20236728645867960961794311064330e-11_real_8,&
      -0.47413283039555834655210340820160e-12_real_8,&
      +0.11221172048389864324731799928920e-12_real_8,&
      -0.26804225434840309912826809093395e-13_real_8,&
      +0.64578514417716530343580369067212e-14_real_8,&
      -0.15682760501666478830305702849194e-14_real_8,&
      +0.38367865399315404861821516441408e-15_real_8,&
      -0.94517173027579130478871048932556e-16_real_8,&
      +0.23434812288949573293896666439133e-16_real_8,&
      -0.58458661580214714576123194419882e-17_real_8,&
      +0.14666229867947778605873617419195e-17_real_8,&
      -0.36993923476444472706592538274474e-18_real_8,&
      +0.93790159936721242136014291817813e-19_real_8,&
      -0.23893673221937873136308224087381e-19_real_8,&
      +0.61150624629497608051934223837866e-20_real_8,&
      -0.15718585327554025507719853288106e-20_real_8,&
      +0.40572387285585397769519294491306e-21_real_8,&
      -0.10514026554738034990566367122773e-21_real_8,&
      +0.27349664930638667785806003131733e-22_real_8,&
      -0.71401604080205796099355574271999e-23_real_8,&
      +0.18705552432235079986756924211199e-23_real_8,&
      -0.49167468166870480520478020949333e-24_real_8,&
      +0.12964988119684031730916087125333e-24_real_8,&
      -0.34292515688362864461623940437333e-25_real_8,&
      +0.90972241643887034329104820906666e-26_real_8,&
      -0.24202112314316856489934847999999e-26_real_8,&
      +0.64563612934639510757670475093333e-27_real_8,&
      -0.17269132735340541122315987626666e-27_real_8,&
      +0.46308611659151500715194231466666e-28_real_8,&
      -0.12448703637214131241755170133333e-28_real_8,&
      +0.33544574090520678532907007999999e-29_real_8,&
      -0.90598868521070774437543935999999e-30_real_8,&
      +0.24524147051474238587273216000000e-30_real_8,&
      -0.66528178733552062817107967999999e-31_real_8 /)
  REAL(real_8), DIMENSION(60), SAVE :: ae11cs = (/ &
      +0.20263150647078889499401236517381e+0_real_8,&
      -0.73655140991203130439536898728034e-1_real_8,&
      +0.63909349118361915862753283840020e-2_real_8,&
      -0.60797252705247911780653153363999e-3_real_8,&
      -0.73706498620176629330681411493484e-4_real_8,&
      +0.48732857449450183453464992488076e-4_real_8,&
      -0.23837064840448290766588489460235e-5_real_8,&
      -0.30518612628561521027027332246121e-5_real_8,&
      +0.17050331572564559009688032992907e-6_real_8,&
      +0.23834204527487747258601598136403e-6_real_8,&
      +0.10781772556163166562596872364020e-7_real_8,&
      -0.17955692847399102653642691446599e-7_real_8,&
      -0.41284072341950457727912394640436e-8_real_8,&
      +0.68622148588631968618346844526664e-9_real_8,&
      +0.53130183120506356147602009675961e-9_real_8,&
      +0.78796880261490694831305022893515e-10_real_8,&
      -0.26261762329356522290341675271232e-10_real_8,&
      -0.15483687636308261963125756294100e-10_real_8,&
      -0.25818962377261390492802405122591e-11_real_8,&
      +0.59542879191591072658903529959352e-12_real_8,&
      +0.46451400387681525833784919321405e-12_real_8,&
      +0.11557855023255861496288006203731e-12_real_8,&
      -0.10475236870835799012317547189670e-14_real_8,&
      -0.11896653502709004368104489260929e-13_real_8,&
      -0.47749077490261778752643019349950e-14_real_8,&
      -0.81077649615772777976249734754135e-15_real_8,&
      +0.13435569250031554199376987998178e-15_real_8,&
      +0.14134530022913106260248873881287e-15_real_8,&
      +0.49451592573953173115520663232883e-16_real_8,&
      +0.79884048480080665648858587399367e-17_real_8,&
      -0.14008632188089809829248711935393e-17_real_8,&
      -0.14814246958417372107722804001680e-17_real_8,&
      -0.55826173646025601904010693937113e-18_real_8,&
      -0.11442074542191647264783072544598e-18_real_8,&
      +0.25371823879566853500524018479923e-20_real_8,&
      +0.13205328154805359813278863389097e-19_real_8,&
      +0.62930261081586809166287426789485e-20_real_8,&
      +0.17688270424882713734999261332548e-20_real_8,&
      +0.23266187985146045209674296887432e-21_real_8,&
      -0.67803060811125233043773831844113e-22_real_8,&
      -0.59440876959676373802874150531891e-22_real_8,&
      -0.23618214531184415968532592503466e-22_real_8,&
      -0.60214499724601478214168478744576e-23_real_8,&
      -0.65517906474348299071370444144639e-24_real_8,&
      +0.29388755297497724587042038699349e-24_real_8,&
      +0.22601606200642115173215728758510e-24_real_8,&
      +0.89534369245958628745091206873087e-25_real_8,&
      +0.24015923471098457555772067457706e-25_real_8,&
      +0.34118376888907172955666423043413e-26_real_8,&
      -0.71617071694630342052355013345279e-27_real_8,&
      -0.75620390659281725157928651980799e-27_real_8,&
      -0.33774612157467324637952920780800e-27_real_8,&
      -0.10479325703300941711526430332245e-27_real_8,&
      -0.21654550252170342240854880201386e-28_real_8,&
      -0.75297125745288269994689298432000e-30_real_8,&
      +0.19103179392798935768638084000426e-29_real_8,&
      +0.11492104966530338547790728833706e-29_real_8,&
      +0.43896970582661751514410359193600e-30_real_8,&
      +0.12320883239205686471647157725866e-30_real_8,&
      +0.22220174457553175317538581162666e-31_real_8/)
  REAL(real_8), DIMENSION(64), SAVE :: ae14cs = (/ &
      -0.1892918000753016825495679942820e+0_real_8,&
      -0.8648117855259871489968817056824e-1_real_8,&
      +0.7224101543746594747021514839184e-2_real_8,&
      -0.8097559457557386197159655610181e-3_real_8,&
      +0.1099913443266138867179251157002e-3_real_8,&
      -0.1717332998937767371495358814487e-4_real_8,&
      +0.2985627514479283322825342495003e-5_real_8,&
      -0.5659649145771930056560167267155e-6_real_8,&
      +0.1152680839714140019226583501663e-6_real_8,&
      -0.2495030440269338228842128765065e-7_real_8,&
      +0.5692324201833754367039370368140e-8_real_8,&
      -0.1359957664805600338490030939176e-8_real_8,&
      +0.3384662888760884590184512925859e-9_real_8,&
      -0.8737853904474681952350849316580e-10_real_8,&
      +0.2331588663222659718612613400470e-10_real_8,&
      -0.6411481049213785969753165196326e-11_real_8,&
      +0.1812246980204816433384359484682e-11_real_8,&
      -0.5253831761558460688819403840466e-12_real_8,&
      +0.1559218272591925698855028609825e-12_real_8,&
      -0.4729168297080398718476429369466e-13_real_8,&
      +0.1463761864393243502076199493808e-13_real_8,&
      -0.4617388988712924102232173623604e-14_real_8,&
      +0.1482710348289369323789239660371e-14_real_8,&
      -0.4841672496239229146973165734417e-15_real_8,&
      +0.1606215575700290408116571966188e-15_real_8,&
      -0.5408917538957170947895023784252e-16_real_8,&
      +0.1847470159346897881370231402310e-16_real_8,&
      -0.6395830792759094470500610425050e-17_real_8,&
      +0.2242780721699759457250233276170e-17_real_8,&
      -0.7961369173983947552744555308646e-18_real_8,&
      +0.2859308111540197459808619929272e-18_real_8,&
      -0.1038450244701137145900697137446e-18_real_8,&
      +0.3812040607097975780866841008319e-19_real_8,&
      -0.1413795417717200768717562723696e-19_real_8,&
      +0.5295367865182740958305442594815e-20_real_8,&
      -0.2002264245026825902137211131439e-20_real_8,&
      +0.7640262751275196014736848610918e-21_real_8,&
      -0.2941119006868787883311263523362e-21_real_8,&
      +0.1141823539078927193037691483586e-21_real_8,&
      -0.4469308475955298425247020718489e-22_real_8,&
      +0.1763262410571750770630491408520e-22_real_8,&
      -0.7009968187925902356351518262340e-23_real_8,&
      +0.2807573556558378922287757507515e-23_real_8,&
      -0.1132560944981086432141888891562e-23_real_8,&
      +0.4600574684375017946156764233727e-24_real_8,&
      -0.1881448598976133459864609148108e-24_real_8,&
      +0.7744916111507730845444328478037e-25_real_8,&
      -0.3208512760585368926702703826261e-25_real_8,&
      +0.1337445542910839760619930421384e-25_real_8,&
      -0.5608671881802217048894771735210e-26_real_8,&
      +0.2365839716528537483710069473279e-26_real_8,&
      -0.1003656195025305334065834526856e-26_real_8,&
      +0.4281490878094161131286642556927e-27_real_8,&
      -0.1836345261815318199691326958250e-27_real_8,&
      +0.7917798231349540000097468678144e-28_real_8,&
      -0.3431542358742220361025015775231e-28_real_8,&
      +0.1494705493897103237475066008917e-28_real_8,&
      -0.6542620279865705439739042420053e-29_real_8,&
      +0.2877581395199171114340487353685e-29_real_8,&
      -0.1271557211796024711027981200042e-29_real_8,&
      +0.5644615555648722522388044622506e-30_real_8,&
      -0.2516994994284095106080616830293e-30_real_8,&
      +0.1127259818927510206370368804181e-30_real_8,&
      -0.5069814875800460855562584719360e-31_real_8/)
  REAL(real_8), SAVE                         :: xmax

! ***BEGIN PROLOGUE  DE1
! ***PURPOSE  Compute the exponential integral E1(X).
! ***LIBRARY   SLATEC (FNLIB)
! ***CATEGORY  C5
! ***TYPE      DOUBLE PRECISION (E1-S, DE1-D)
! ***KEYWORDS  E1 FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
! SPECIAL FUNCTIONS
! ***AUTHOR  Fullerton, W., (LANL)
! ***DESCRIPTION
! 
! DE1 calculates the double precision exponential integral, E1(X), for
! positive double precision argument X and the Cauchy principal value
! for negative X.  If principal values are used everywhere, then, for
! all X,
! 
! E1(X) = -Ei(-X)
! or
! Ei(X) = -E1(-X).
! 
! 
! Series for AE10       on the interval -3.12500E-02 to  0.
! with weighted error   4.62E-32
! log weighted error  31.34
! significant figures required  29.70
! decimal places required  32.18
! 
! 
! Series for AE11       on the interval -1.25000E-01 to -3.12500E-02
! with weighted error   2.22E-32
! log weighted error  31.65
! significant figures required  30.75
! decimal places required  32.54
! 
! 
! Series for AE12       on the interval -2.50000E-01 to -1.25000E-01
! with weighted error   5.19E-32
! log weighted error  31.28
! significant figures required  30.82
! decimal places required  32.09
! 
! 
! Series for E11        on the interval -4.00000E+00 to -1.00000E+00
! with weighted error   8.49E-34
! log weighted error  33.07
! significant figures required  34.13
! decimal places required  33.80
! 
! 
! Series for E12        on the interval -1.00000E+00 to  1.00000E+00
! with weighted error   8.08E-33
! log weighted error  32.09
! approx significant figures required  30.4
! decimal places required  32.79
! 
! 
! Series for AE13       on the interval  2.50000E-01 to  1.00000E+00
! with weighted error   6.65E-32
! log weighted error  31.18
! significant figures required  30.69
! decimal places required  32.03
! 
! 
! Series for AE14       on the interval  0.          to  2.50000E-01
! with weighted error   5.07E-32
! log weighted error  31.30
! significant figures required  30.40
! decimal places required  32.20
! 
! ***REFERENCES  (NONE)
! ***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
! ***REVISION HISTORY  (YYMMDD)
! 770701  DATE WRITTEN
! 890531  Changed all specific intrinsics to generic.  (WRB)
! 891115  Modified prologue description.  (WRB)
! 891115  REVISION DATE from Version 3.2
! 891214  Prologue converted to Version 4.0 format.  (BAB)
! 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
! 920618  Removed space from variable names.  (RWC, WRB)
! ***END PROLOGUE  DE1
! ***FIRST EXECUTABLE STATEMENT  DE1

  IF (first) THEN
     eta = 0.1*REAL(d1mach(3))
     ntae10 = initds (ae10cs, 50, eta)
     ntae11 = initds (ae11cs, 60, eta)
     ntae12 = initds (ae12cs, 41, eta)
     nte11 = initds (e11cs, 29, eta)
     nte12 = initds (e12cs, 25, eta)
     ntae13 = initds (ae13cs, 50, eta)
     ntae14 = initds (ae14cs, 64, eta)
     ! 
     xmaxt = -LOG(d1mach(1))
     xmax = xmaxt - LOG(xmaxt)
  ENDIF
  first = .FALSE.
  ! 
  IF (x.GT.(-1._real_8)) go to 50
  IF (x.GT.(-32._real_8)) go to 20
  my_res = EXP(-x)/x * (1._real_8 + dcsevl (64._real_8/x+1._real_8, ae10cs, ntae10))
  RETURN
  ! 
20 IF (x.GT.(-8._real_8)) go to 30
  my_res = EXP(-x)/x * (1._real_8 + dcsevl ((64._real_8/x+5._real_8)/3._real_8, ae11cs,&
       ntae11))
  RETURN
  ! 
30 IF (x.GT.(-4._real_8)) go to 40
  my_res = EXP(-x)/x * (1._real_8 + dcsevl (16._real_8/x+3._real_8, ae12cs, ntae12))
  RETURN
  ! 
40 my_res = -LOG(-x) + dcsevl ((2._real_8*x+5._real_8)/3._real_8, e11cs, nte11)
  RETURN
  ! 
50 IF (x.GT.1.0_real_8) go to 60
  IF (x .EQ. 0._real_8) STOP 103 ! CALL XERMSG ('SLATEC', 'DE1', 'X IS 0', 2, 2)
  my_res = (-LOG(ABS(x)) - 0.6875_real_8 + x)  + dcsevl (x, e12cs, nte12)
  RETURN
  ! 
60 IF (x.GT.4.0_real_8) go to 70
  my_res = EXP(-x)/x * (1._real_8 + dcsevl ((8._real_8/x-5._real_8)/3._real_8, ae13cs,&
       ntae13))
  RETURN
  ! 
70 IF (x.GT.xmax) go to 80
  my_res = EXP(-x)/x * (1._real_8 + dcsevl (8._real_8/x-1._real_8, ae14cs, ntae14))
  RETURN
  ! 
  ! orig 80   stop 104 !CALL XERMSG ('SLATEC', 'DE1', 'X SO BIG E1 UNDERFLOWS', 1, 1)
  ! orig      DE1 = 0._real_8
80 my_res = 0._real_8
  RETURN
  ! 
END FUNCTION de1
! DECK DEI
FUNCTION dei (x) RESULT(my_res)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  REAL(real_8)                               :: x, my_res

  REAL(real_8)                               :: de1

! ***BEGIN PROLOGUE  DEI
! ***PURPOSE  Compute the exponential integral Ei(X).
! ***LIBRARY   SLATEC (FNLIB)
! ***CATEGORY  C5
! ***TYPE      DOUBLE PRECISION (EI-S, DEI-D)
! ***KEYWORDS  EI FUNCTION, EXPONENTIAL INTEGRAL, FNLIB,
! SPECIAL FUNCTIONS
! ***AUTHOR  Fullerton, W., (LANL)
! ***DESCRIPTION
! 
! DEI calculates the double precision exponential integral, Ei(X), for
! positive double precision argument X and the Cauchy principal value
! for negative X.  If principal values are used everywhere, then, for
! all X,
! 
! Ei(X) = -E1(-X)
! or
! E1(X) = -Ei(-X).
! 
! ***REFERENCES  (NONE)
! ***ROUTINES CALLED  DE1
! ***REVISION HISTORY  (YYMMDD)
! 770701  DATE WRITTEN
! 891115  Modified prologue description.  (WRB)
! 891115  REVISION DATE from Version 3.2
! 891214  Prologue converted to Version 4.0 format.  (BAB)
! ***END PROLOGUE  DEI
! ***FIRST EXECUTABLE STATEMENT  DEI

  my_res = -de1(-x)
  ! 
  RETURN
END FUNCTION dei
! *DECK INITDS
FUNCTION initds (os, nos, eta) RESULT(my_res)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  REAL(real_8)                               :: os(*)
  INTEGER                                    :: nos
  REAL(real_8)                               :: eta
  INTEGER                                    :: my_res

  INTEGER                                    :: i, ii
  REAL(real_8)                               :: err

! ***BEGIN PROLOGUE  INITDS
! ***PURPOSE  Determine the number of terms needed in an orthogonal
! polynomial series so that it meets a specified accuracy.
! ***LIBRARY   SLATEC (FNLIB)
! ***CATEGORY  C3A2
! ***TYPE      DOUBLE PRECISION (INITS-S, INITDS-D)
! ***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
! ORTHOGONAL SERIES, SPECIAL FUNCTIONS
! ***AUTHOR  Fullerton, W., (LANL)
! ***DESCRIPTION
! 
! Initialize the orthogonal series, represented by the array OS, so
! that INITDS is the number of terms needed to insure the error is no
! larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
! machine precision.
! 
! Input Arguments --
! OS     double precision array of NOS coefficients in an orthogonal
! series.
! NOS    number of coefficients in OS.
! ETA    single precision scalar containing requested accuracy of
! series.
! 
! ***REFERENCES  (NONE)
! ***ROUTINES CALLED  XERMSG
! ***REVISION HISTORY  (YYMMDD)
! 770601  DATE WRITTEN
! 890531  Changed all specific intrinsics to generic.  (WRB)
! 890831  Modified array declarations.  (WRB)
! 891115  Modified error message.  (WRB)
! 891115  REVISION DATE from Version 3.2
! 891214  Prologue converted to Version 4.0 format.  (BAB)
! 900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
! ***END PROLOGUE  INITDS
! ***FIRST EXECUTABLE STATEMENT  INITDS

  IF (nos .LT. 1) STOP 105 ! CALL XERMSG ('SLATEC', 'INITDS',
  ! +   'Number of coefficients is less than 1', 2, 1)
  ! 
  err = 0.
  DO ii = 1,nos
     i = nos + 1 - ii
     err = err + ABS(REAL(os(i)))
     IF (err.GT.eta) go to 20
  END DO
  ! 
20 IF (i .EQ. nos) STOP 106 ! CALL XERMSG ('SLATEC', 'INITDS',
  ! +   'Chebyshev series too short for specified accuracy', 1, 1)
  my_res = i
  ! 
  RETURN
END FUNCTION initds

FUNCTION d1mach(i) RESULT(my_res)
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  INTEGER                                    :: i
  REAL(real_8)                               :: my_res

  INTEGER                                    :: cray1(38), j
  INTEGER, SAVE                              :: diver(2), large(2), LOG10(2), &
                                                right(2), sc = 0, small(2)
  REAL(real_8)                               :: dmach(5)

! 
! DOUBLE-PRECISION MACHINE CONSTANTS
! D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
! D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
! D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
! D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
! D1MACH( 5) = LOG10(B)
! 

  EQUIVALENCE (dmach(1),small(1))
  EQUIVALENCE (dmach(2),large(1))
  EQUIVALENCE (dmach(3),right(1))
  EQUIVALENCE (dmach(4),diver(1))
  EQUIVALENCE (dmach(5),LOG10(1))
  COMMON /d9mach/ cray1
  ! THIS VERSION ADAPTS AUTOMATICALLY TO MOST CURRENT MACHINES.
  ! R1MACH CAN HANDLE AUTO-DOUBLE COMPILING, BUT THIS VERSION OF
  ! D1MACH DOES NOT, BECAUSE WE DO NOT HAVE QUAD CONSTANTS FOR
  ! MANY MACHINES YET.
  ! TO COMPILE ON OLDER MACHINES, ADD A C IN COLUMN 1
  ! ON THE NEXT LINE
  ! AND REMOVE THE C FROM COLUMN 1 IN ONE OF THE SECTIONS BELOW.
  ! CONSTANTS FOR EVEN OLDER MACHINES CAN BE OBTAINED BY
  ! mail netlib@research.bell-labs.com
  ! send old1mach from blas
  ! PLEASE SEND CORRECTIONS TO dmg OR ehg@bell-labs.com.
  ! 
  ! MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
  ! DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
  ! DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
  ! DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
  ! DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
  ! DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
  ! 
  ! MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
  ! 32-BIT INTEGERS.
  ! DATA SMALL(1),SMALL(2) /    8388608,           0 /
  ! DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
  ! DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
  ! DATA DIVER(1),DIVER(2) /  620756992,           0 /
  ! DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
  ! 
  ! MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
  ! DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
  ! DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
  ! DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
  ! DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
  ! DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
  ! 
  ! ON FIRST CALL, IF NO DATA UNCOMMENTED, TEST MACHINE TYPES.
  IF (sc .NE. 987) THEN
     dmach(1) = 1.e13_real_8
     IF (      small(1) .EQ. 1117925532&
          .AND. small(2) .EQ. -448790528) THEN
        ! *           *** IEEE BIG ENDIAN ***
        small(1) = 1048576
        small(2) = 0
        large(1) = 2146435071
        large(2) = -1
        right(1) = 1017118720
        right(2) = 0
        diver(1) = 1018167296
        diver(2) = 0
        LOG10(1) = 1070810131
        LOG10(2) = 1352628735
     ELSE IF ( small(2) .EQ. 1117925532&
          .AND. small(1) .EQ. -448790528) THEN
        ! *           *** IEEE LITTLE ENDIAN ***
        small(2) = 1048576
        small(1) = 0
        large(2) = 2146435071
        large(1) = -1
        right(2) = 1017118720
        right(1) = 0
        diver(2) = 1018167296
        diver(1) = 0
        LOG10(2) = 1070810131
        LOG10(1) = 1352628735
     ELSE IF ( small(1) .EQ. -2065213935&
          .AND. small(2) .EQ. 10752) THEN
        ! *** VAX WITH D_FLOATING ***
        small(1) = 128
        small(2) = 0
        large(1) = -32769
        large(2) = -1
        right(1) = 9344
        right(2) = 0
        diver(1) = 9472
        diver(2) = 0
        LOG10(1) = 546979738
        LOG10(2) = -805796613
     ELSE IF ( small(1) .EQ. 1267827943&
          .AND. small(2) .EQ. 704643072) THEN
        ! *** IBM MAINFRAME ***
        small(1) = 1048576
        small(2) = 0
        large(1) = 2147483647
        large(2) = -1
        right(1) = 856686592
        right(2) = 0
        diver(1) = 873463808
        diver(2) = 0
        LOG10(1) = 1091781651
        LOG10(2) = 1352628735
     ELSE IF ( small(1) .EQ. 1120022684&
          .AND. small(2) .EQ. -448790528) THEN
        ! *** CONVEX C-1 ***
        small(1) = 1048576
        small(2) = 0
        large(1) = 2147483647
        large(2) = -1
        right(1) = 1019215872
        right(2) = 0
        diver(1) = 1020264448
        diver(2) = 0
        LOG10(1) = 1072907283
        LOG10(2) = 1352628735
     ELSE IF ( small(1) .EQ. 815547074&
          .AND. small(2) .EQ. 58688) THEN
        ! *** VAX G-FLOATING ***
        small(1) = 16
        small(2) = 0
        large(1) = -32769
        large(2) = -1
        right(1) = 15552
        right(2) = 0
        diver(1) = 15568
        diver(2) = 0
        LOG10(1) = 1142112243
        LOG10(2) = 2046775455
     ELSE
        dmach(2) = 1.e27_real_8 + 1
        dmach(3) = 1.e27_real_8
        large(2) = large(2) - right(2)
        IF (large(2) .EQ. 64 .AND. small(2) .EQ. 0) THEN
           cray1(1) = 67291416
           DO j = 1, 20
              cray1(j+1) = cray1(j) + cray1(j)
           END DO
           cray1(22) = cray1(21) + 321322
           DO j = 22, 37
              cray1(j+1) = cray1(j) + cray1(j)
           END DO
           IF (cray1(38) .EQ. small(1)) THEN
              CALL i1mcry(small(1), j, 8285, 8388608, 0)
              small(2) = 0
              CALL i1mcry(large(1), j, 24574, 16777215, 16777215)
              CALL i1mcry(large(2), j, 0, 16777215, 16777214)
              CALL i1mcry(right(1), j, 16291, 8388608, 0)
              right(2) = 0
              CALL i1mcry(diver(1), j, 16292, 8388608, 0)
              diver(2) = 0
              CALL i1mcry(LOG10(1), j, 16383, 10100890, 8715215)
              CALL i1mcry(LOG10(2), j, 0, 16226447, 9001388)
           ELSE
              WRITE(6,9000)
              STOP 779
           ENDIF
        ELSE
           WRITE(6,9000)
           STOP 779
        ENDIF
     ENDIF
     sc = 987
  ENDIF
  ! SANITY CHECK
  IF (dmach(4) .GE. 1.0_real_8) STOP 778
  IF (i .LT. 1 .OR. i .GT. 5) THEN
     WRITE(6,*) 'D1MACH(I): I =',i,' is out of bounds.'
     STOP
  ENDIF
  my_res = dmach(i)
  RETURN
9000 FORMAT(/' Adjust D1MACH by uncommenting data statements'/&
       ' appropriate for your machine.')
  ! /* Standard C source for D1MACH -- remove the * in column 1 */
  ! #include <stdio.h>
  ! #include <float.h>
  ! #include <math.h>
  ! double d1mach_(long *i)
  ! {
  ! switch(*i){
  ! case 1: return DBL_MIN;
  ! case 2: return DBL_MAX;
  ! case 3: return DBL_EPSILON/FLT_RADIX;
  ! case 4: return DBL_EPSILON;
  ! case 5: return log10((double)FLT_RADIX);
  ! }
  ! fprintf(stderr, "invalid argument: d1mach(%ld)\n", *i);
  ! exit(1); return 0; /* some compilers demand return values */
  ! }
END FUNCTION d1mach
SUBROUTINE i1mcry(a, a1, b, c, d)
  ! *** SPECIAL COMPUTATION FOR OLD CRAY MACHINES ****
  USE kinds, ONLY: real_4, real_8, int_1, int_2, int_4, int_8
  USE error_handling, ONLY: stopgm
  IMPLICIT NONE
  INTEGER                                    :: a, a1, b, c, d

  a1 = 16777216*b + c
  a = 16777216*a1 + d
END SUBROUTINE i1mcry

