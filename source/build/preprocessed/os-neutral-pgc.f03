# 1 "pgc/os-neutral-pgc.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "pgc/os-neutral-pgc.f03" 2
# 1 "./os-config.h" 1
! Configuration file for osiris

! ----------------------------------------------------------------------------------------
! Algorithm options
! ----------------------------------------------------------------------------------------

! ----------------------------------------------------------------------------------------
! System options
! ----------------------------------------------------------------------------------------

! MPI supports MPI_IN_PLACE in global operations
!#define 1

! Use OpenMP
!#define 1

! Use SIMD optimized code
!#define SIMD
!#define SIMD_SSE
!#define SIMD_AVX
!#define SIMD_BGQ
!#define SIMD_MIC

! Use SION for checkpointing (not available on all systems)
!#define __RST_IO__ = __RST_SION__

! Use log files
!#define __USE_LOG__

! Use MPE for logging and profiling
!#define __USE_MPE__

! Compiler does not fully support the sizeof intrinsic which is a Fortran 2003 feature
!#define __NO_SIZEOF__

! Use the PAPI library for profiling
!#define __USE_PAPI__


! ----------------------------------------------------------------------------------------
! Distribution Options
! ----------------------------------------------------------------------------------------
!
! Optional modules to be removed for distribution
! These are all turned on by default unless the __DISTRO__ preprocessor macro is defined
!
! ----------------------------------------------------------------------------------------



! Include Ionization module


! Include binary collisions module


! Include particle tracking module


! Include perfectly matched layers boundary conditions for EMF


! Include spin advance module which is disabled by default
! #define __HAS_SPIN__

! Include a debug flag for cylindrical modes simulations
! #define __CYL_MODES_DEBUG__

! Compile functions for reporting B fields in the RaDiO module
! #define __HAS_RAD_BFLD__
# 2 "pgc/os-neutral-pgc.f03" 2
# 1 "./os-preprocess.fpp" 1
!
! File: os-preprocess.fpp
!
! A set of preprocessing macros for osiris, and the fpp/cpp preprocessor
!
!



! Macros for the IBM CPP/GNU preprocessor




! Assertion macro



! Debug functions
# 33 "./os-preprocess.fpp"
! SCR hostname functions







! LOG functions




! ERROR functions






! WARNING functions




! functions for restart io
# 3 "pgc/os-neutral-pgc.f03" 2

module m_neutral_pgc

# 1 "./memory/memory.h" 1
! Include file for the memory module
! Files must include this file using #include "memory/memory.h" rather than just using the module
! the #include must be placed where the module use statement would usually be i.e.
!
! module module2
!
! use module1
! #include "memory/memory.h"
!
# 19 "./memory/memory.h"
use memory
# 7 "pgc/os-neutral-pgc.f03" 2

use m_system
use m_parameters

use m_neutral
use m_emf_pgc_define
use m_species_define
use m_vdf_define

implicit none

private

! Ionization Parameters for ADK-PGC model

real(p_double), parameter, dimension(3,1) :: H_param_pgc = &
    reshape( (/ 9.839764713542517e16, 121.10599187554276, 0.4138362761594194 /), &
                                                                     (/ 3, 1 /) )
real(p_double), parameter, dimension(3,3) :: Li_param_pgc = &
    reshape( (/ -4.9924067896571425e14, 30.235882358115305, -0.10973740638330676, &
                2.5808423419428425e38, 1588.7675879635117, 5.668999676670961, &
                4.1704319938962566e63, 3272.6053071455326, 11.728084262388117 /), &
                                                                     (/ 3, 3 /) )
real(p_double), parameter, dimension(3,11) :: Na_param_pgc = &
    reshape( (/ -5.186912554033257e14, 28.135843092439796, -0.13084502302676904, &
                6.518227505635966e31, 785.3009276560339, 4.272937326824533, &
                9.609535964106987e48, 1463.8067958799568, 8.734039061391961, &
                9.797792901646653e70, 2376.642208157135, 14.254284320880808, &
                4.1065223082039873e101, 3932.207428057183, 21.552397675232996, &
                8.12908536404805e133, 5458.683593247487, 29.18963074963119, &
                6.896332546087885e170, 7270.953597628679, 37.753036848271975, &
                8.098801722845243e220, 10370.788513077201, 48.85451433990087, &
                5.0067481918939675e264, 12540.362453157853, 58.7525507454844, &
                2.922776390409296e753, 135440.4848572289, 145.75517012249176, &
                2.661653565402382e884, 161675.990841195, 170.2452126549242 /), &
                                                                     (/ 3, 11 /) )
real(p_double), parameter, dimension(3,19) :: K_param_pgc = &
    reshape( (/ -5.121839723403554e14, 21.840695763446327, -0.201209790283859, &
                3.1779008212482716e27, 429.5147527523833, 3.312210317792947, &
                9.567575141125409e39, 748.6426976844796, 6.784368069893642, &
                2.950110369365608e55, 1148.259494434296, 10.969717940248294, &
                5.4603146182687566e76, 1814.9948015862528, 16.42899624872467, &
                1.43723411241838e98, 2393.3811192216967, 21.93501952553697, &
                4.2152109180687376e122, 3078.375898636897, 28.09929721395108, &
                1.839843200321142e161, 4654.609285988508, 37.17054033523553, &
                4.21175167351257e192, 5630.216097563873, 44.753898852003836, &
                3.142739266789641e397, 27299.29363068137, 85.04541103512369, &
                2.395953238462661e465, 32486.010867344186, 99.30025962221939, &
                7.685042098480176e538, 38280.38728740257, 114.57145935818653, &
                9.13572346437244e625, 46145.25731673305, 132.24869159186383, &
                4.932759997650713e710, 53245.33418169792, 149.51015337731778, &
                8.959575268081703e800, 61006.551942052414, 167.74360738666746, &
                2.781569406962453e913, 72696.9964084871, 189.82534173755622, &
                1.105099482936575e1007, 80362.57506578414, 208.64159075238013, &
                5.14376972405305e2608, 756147.6346021439, 467.61762070575924, &
                1.264878652072477e2860, 837022.9495536234, 510.6935780431256 /), &
                                                                     (/ 3, 19 /) )
real(p_double), parameter, dimension(3,37) :: Rb_param_pgc = &
    reshape( (/ -5.023452890327682e14, 20.618115086452168, -0.21640155876083345, &
                1.4882241999234819e26, 344.2917425681459, 3.005740163009257, &
                2.5524514129168078e37, 593.8024457956027, 6.205749478880865, &
                2.672395328514122e51, 910.8314933425314, 10.080259059820088, &
                2.537154495583009e69, 1366.2067812525477, 14.85450644295686, &
                5.7632870298618255e88, 1822.905183610664, 19.945136126845977, &
                8.962772910143388e110, 2367.0638689582433, 25.65908723244124, &
                2.4776828050763747e147, 3695.978163104466, 34.34627858997096, &
                5.7170596261350214e175, 4464.774182491545, 41.34995167559996, &
                4.610252910627844e277, 11141.254362310976, 62.824674230713505, &
                6.120933222329582e325, 13380.023786262425, 73.62578122213996, &
                1.552576410011341e381, 16222.107723804585, 85.80823309957827, &
                2.069957222166411e440, 19320.673928615342, 98.68458830703847, &
                4.761550144712687e501, 22518.41783964962, 111.97557760562901, &
                1.291613032325942e577, 27163.633349147356, 127.85396570086405, &
                7.84154081578247e647, 31151.31026053306, 142.86528123287496, &
                4.981647737298132e723, 35583.11728300991, 158.78677411284053, &
                6.95295619325298e803, 40392.3300676751, 175.48835325222427, &
                1.131036049727838e886, 45304.29073373004, 192.55777587671014, &
                3.594707921245011e1046, 60590.3577781778, 223.47866879950706, &
                1.015494151985715e1133, 65784.11926656753, 241.2535927213357, &
                1.146812036114324e1226, 71712.0705336964, 260.19449040134407, &
                6.961830819703348e1330, 79137.48041160846, 281.18406536146273, &
                1.431234704298912e1431, 85717.18838892713, 301.39721776149304, &
                9.02555720497733e1524, 91129.93651861558, 320.4925708129116, &
                2.626855422562381e1681, 105773.5012103444, 350.37933638498305, &
                5.507526763227677e1786, 112491.33797235036, 371.4609115882289, &
                7.97113733320106e3157, 423579.8108890921, 599.9158956886857, &
                1.433579827158837e3357, 453880.6097729353, 635.8773219370214, &
                7.82112309366001e3568, 487695.314390501, 673.80976625396, &
                1.222742988878172e3783, 521658.1960726142, 712.1282086580444, &
                6.94229032822265e4037, 569081.0150930289, 756.7953239385334, &
                3.637052366611573e4270, 608226.0384063595, 797.998909024477, &
                9.3548459409644e4541, 660654.0185040454, 845.2152629151013, &
                9.6850122043205e4763, 694533.88497849, 884.7471581493115, &
                9.08878626796164e11472, 5.981612944959795e6, 1866.465182784449, &
                3.382208377271529e12027, 6.30781681815364e6, 1952.6135488911784 /), &
                                                                     (/ 3, 37 /) )
real(p_double), parameter, dimension(3,55) :: Cs_param_pgc = &
    reshape( (/ -4.8116271120580456e14, 18.557112620622686, -0.24343302823301816, &
                6.745823980290736e24, 269.1338887451569, 2.690031170843273, &
                7.641196102892746e34, 461.89290746438144, 5.6269281635727895, &
                7.43385038762432e46, 680.9809570005151, 9.056547071442443, &
                3.449071509036452e62, 1012.0787925143602, 13.345599084098065, &
                2.162220719347091e80, 1387.2328202415824, 18.122512156804223, &
                3.3941703436672164e100, 1822.905183610664, 23.435992147986973, &
                1.4097738742991902e132, 2789.6780909005424, 31.18241230284402, &
                7.48150453846515e157, 3399.9192896301042, 37.67312695935667, &
                2.0995777706413515e237, 7523.476473651201, 54.99510963998894, &
                3.924146985178889e272, 8589.466013234265, 63.37619699531625, &
                6.815791759097292e315, 10183.405675807735, 73.32863739819688, &
                2.353917425816354e361, 11865.308876410898, 83.7319000609827, &
                2.179849769300064e408, 13566.346585936406, 94.41710783047931, &
                4.9496422915729e464, 15949.470853392926, 106.89896092465622, &
                3.245613071319197e518, 18031.31201226235, 118.89645315164722, &
                1.959360942787137e575, 20270.16862531454, 131.45812978552215, &
                1.113832881076642e635, 22671.08475546965, 144.58183266297823, &
                9.80142339004228e695, 25080.8767000135, 157.932146508014, &
                9.76468354766127e838, 35228.46991270996, 186.35782358879885, &
                6.010003798499595e906, 38098.532906423774, 200.92927642850043, &
                1.010303132079352e981, 41509.12969772678, 216.67796415412698, &
                3.434436464486747e1054, 44727.98608012481, 232.30901902863314, &
                6.948218063643336e1147, 49902.78094909463, 251.5011738409615, &
                2.347748510406052e1229, 53727.54657541125, 268.57705344570144, &
                1.237079730944337e1357, 62509.26085738964, 293.87098878493936, &
                2.364810943186027e1446, 66961.51567140219, 312.3161390820841, &
                2.80128545315565e2097, 153407.60125356243, 427.3356626368801, &
                1.201518603944032e2234, 165115.06267560704, 453.6433112241664, &
                5.62599795203768e2377, 177864.76640179835, 481.12739257975403, &
                2.49366342798771e2529, 191860.3508854331, 509.93699437246346, &
                3.837179606771067e2681, 205726.53599189618, 538.8303859088852, &
                9.61103454819747e2840, 220726.9461185236, 568.9144060403288, &
                1.7791823998524e3016, 238582.76493692643, 601.6093276641677, &
                5.99127664978393e3181, 254325.30772045598, 632.6875103516078, &
                3.05418439952989e3355, 271447.91918698436, 665.1038083117157, &
                8.17632681824349e3526, 287868.06747088354, 697.1415165825962, &
                7.42942085558953e3848, 335631.94761776004, 753.6552619404185, &
                4.73698392525011e4019, 351516.7488036082, 785.5455135995584, &
                1.59230673952527e4200, 369193.0501716211, 819.0148940811187, &
                1.223181182973801e4378, 385975.9392898583, 852.0631418978242, &
                1.403703762392516e4637, 420481.07160218025, 898.1704323484598, &
                1.47057473751906e4830, 440049.05502615654, 933.6436579686739, &
                4.6834084858429e5100, 477107.28670525906, 981.5061327473485, &
                1.30782754582599e5297, 496862.77172424836, 1017.5176905133185, &
                2.205794683984641e8893, 1.724530643419232e6, 1575.3720061291983, &
                8.90985504470308e9236, 1.801179255811685e6, 1633.1581355005285, &
                6.79833394659372e9601, 1.887266686109058e6, 1694.1035715075022, &
                1.406208528763027e9959, 1.9685891949310298e6, 1753.924172272537, &
                1.243653169321741e10580, 2.1834824238696047e6, 1852.6617199447924, &
                2.684969394285903e10973, 2.2819391749813347e6, 1917.7369589544237, &
                2.622808659951394e11416, 2.4056715374308294e6, 1990.0982799816081, &
                1.389389857086704e11789, 2.491110838691174e6, 2052.1347972747676, &
                2.551818692651049e27536, 2.0684632930413947e7, 4234.97715760873, &
                3.752553859610269e28437, 2.146917950747324e7, 4367.292881082764 /), &
                                                                     (/ 3, 55 /) )
real(p_double), parameter, dimension(3,2) :: He_param_pgc = &
    reshape( (/ 1.632117076882469e19, 294.4427561631889, 0.9011258556818507, &
                3.9004055510620626e33, 969.4899528598119, 4.656594021244761 /), &
                                                                     (/ 3, 2 /) )
real(p_double), parameter, dimension(3,18) :: Ar_param_pgc = &
    reshape( (/ 3.109362633289577e17, 151.09515290502952, 0.5220447488727498, &
                1.9042400930425332e26, 350.7485068371556, 3.0306261316845307, &
                1.0104494709276334e38, 627.8904317895609, 6.341077123134801, &
                7.323958391028679e54, 1110.6648620999654, 10.837634114687523, &
                6.052073849168599e72, 1563.6262485606783, 15.5840869225337, &
                4.1741401703444e93, 2106.527313848287, 20.97948567122373, &
                4.74445153959699e126, 3351.395572723205, 28.935316035864872, &
                1.567630578823633e154, 4149.669785038994, 35.73711745512601, &
                8.001033588770822e324, 20980.867453320443, 69.93505335630753, &
                1.140579235626615e386, 25378.639565716032, 82.978224723549, &
                6.165182149707942e452, 30341.80815087183, 97.04311056047634, &
                4.445207613074065e532, 37189.09714848429, 113.4626249559031, &
                2.212061738976722e610, 43342.60027829773, 129.49451274890893, &
                3.308670014268435e693, 50114.66950552506, 146.50052614757402, &
                7.522862219625861e797, 60428.172493415375, 167.20864900711996, &
                8.995614498195755e884, 67214.29596899914, 184.90216358635612, &
                1.660271329912637e2310, 638827.0928776392, 417.39575905633126, &
                1.069735889002778e2546, 711184.3247804214, 458.1386017415283 /), &
                                                                     (/ 3, 18 /) )
real(p_double), parameter, dimension(3,7) :: N_param_pgc = &
    reshape( (/ 1.6416146690346304e17, 133.81840362021134, 0.4616695235333683, &
                7.775737294961277e26, 388.9532102233614, 3.1719557429798666, &
                4.0105005656852545e40, 789.263610154775, 6.922687477686023, &
                3.857442309339355e62, 1646.8803536751923, 12.498678417505777, &
                1.1217818595803872e84, 2339.056120969317, 17.966802789946712, &
                1.7150111192953165e263, 31327.10996765152, 53.05077654728761, &
                4.411388612009231e338, 41606.96603090776, 68.31554352778923 /), &
                                                                     (/ 3, 7 /) )
real(p_double), parameter, dimension(3,8) :: O_param_pgc = &
    reshape( (/ 9.948400100888674e16, 121.36818630891071, 0.414855859074839, &
                3.3463248986498183e28, 502.671929772504, 3.54431872537409, &
                2.210324165197067e43, 983.3596428856265, 7.5251617018753905, &
                3.6444245965148695e62, 1644.967563744188, 12.493450322878141, &
                1.8898698680277577e91, 2935.6993643349037, 19.45900181467824, &
                1.2585032723412369e118, 3920.2336340482675, 26.035379950568448, &
                7.666035597030094e359, 48549.74455685462, 71.97446729106068, &
                6.319381581557257e448, 62124.93982785822, 89.54321075365273 /), &
                                                                     (/ 3, 8 /) )
real(p_double), parameter, dimension(3,6) :: C_param_pgc = &
    reshape( (/ 2.4032620168838264e16, 91.25506831471063, 0.2865586086751897, &
                1.725191499410126e25, 290.80580154268296, 2.7865316802319664, &
                5.881636168670233e40, 800.3304496951282, 6.959545661477145, &
                1.1592246380143422e57, 1250.8551230851663, 11.31609146353087, &
                1.818344494801144e183, 18750.452973007043, 36.9592686040155, &
                2.998831331744547e245, 26194.935478238203, 49.921470699278885 /), &
                                                                     (/ 3, 6 /) )
real(p_double), parameter, dimension(3,10) :: Ne_param_pgc = &
    reshape( (/ 4.623217054149945e18, 241.84804541945022, 0.7804294563889254, &
                1.4147531668126254e30, 633.168471566491, 3.9077262831510717, &
                1.922911873630639e46, 1219.8477784233155, 8.160103444984605, &
                2.1659403910710488e70, 2314.006940735653, 14.119083097447419, &
                3.8364317383067267e96, 3425.814823849391, 20.539479746055953, &
                2.811827938500664e127, 4793.422790613258, 27.90971028424146, &
                1.780367456909249e170, 7206.760642772433, 37.63865345554652, &
                1.7694807639558004e208, 8928.80219852768, 46.42765273168592, &
                4.791983233424534e604, 99867.53983589928, 118.32401834381443, &
                2.524995666420888e721, 121420.82170192596, 140.50603627792864 /), &
                                                                     (/ 3, 10 /) )
real(p_double), parameter, dimension(3,54) :: Xe_param_pgc = &
    reshape( (/ 4.170127941392819e16, 102.02698640914633, 0.33531049766788823, &
                1.2105737691971357e24, 231.99853278867565, 2.5118475273557483, &
                8.872816791676921e33, 417.8541221320813, 5.409242865116984, &
                2.8741703695186756e46, 662.0655045364159, 8.962558588203676, &
                2.883936301786264e61, 962.0766211017353, 13.105347881617899, &
                6.43050410206368e78, 1315.6801866149146, 17.787915813148388, &
                2.3216961082874146e106, 2117.2663466052563, 24.686234818280607, &
                1.8333948873545727e129, 2634.842340222692, 30.575636112657435, &
                1.691508245528198e195, 5824.5412442915185, 45.27435009815517, &
                5.099590848018664e229, 6933.608855678177, 53.49169992006542, &
                5.971569229796838e269, 8370.326515246868, 62.82400631600598, &
                2.588287163481154e311, 9834.279670107457, 72.46931791259686, &
                1.790114413685948e355, 11376.057117756873, 82.55090989944684, &
                5.751909289110963e406, 13437.7561604563, 94.11467566611604, &
                2.479723519872479e457, 15341.699225482811, 105.51063912176731, &
                6.832853251705239e511, 17467.85939095684, 117.63435007646672, &
                3.739512075669454e567, 19611.207359780554, 130.0069318083264, &
                6.633971618655159e621, 21534.508673228123, 142.1071736410438, &
                4.595073955237778e759, 31066.390770677666, 169.68464169431664, &
                1.001612259921071e826, 33909.13888028995, 183.98910857015727, &
                1.72411948182608e895, 36923.541743506954, 198.83167213024134, &
                2.993712777444771e966, 40022.32543472028, 214.04731810779134, &
                3.434436464486747e1054, 44727.98608012481, 232.30901902863314, &
                7.453660342016175e1131, 48222.4204775263, 248.63462355581248, &
                8.682831721594174e1254, 56501.79428199787, 273.13932059917283, &
                8.693520681113246e1339, 60593.53933005796, 290.82737713566513, &
                5.64271119962185e1948, 139322.7363773481, 398.98926854305915, &
                3.901004636989117e2080, 150382.244194792, 424.5012061331139, &
                4.20011034918003e2218, 162308.61221681687, 451.05272995872616, &
                2.52843355074898e2365, 175591.91290455376, 479.0649530974581, &
                1.190283546280775e2511, 188444.49017318338, 506.88659699409607, &
                1.86886543144648e2667, 203022.76872499823, 536.4550368525356, &
                7.66451315885582e2835, 219748.59561165317, 568.0711290220139, &
                4.65974593071874e2995, 234574.89448713293, 598.2159116994154, &
                3.853156583949668e3163, 250741.2911613289, 629.696722173378, &
                1.293592738307522e3330, 266393.7541006866, 660.9437567572478, &
                5.70416310179541e3643, 312085.4764264122, 716.1942853686995, &
                6.73150075979773e3808, 327037.44317928125, 747.158015426826, &
                3.620174902644502e3983, 343732.8677040114, 779.696449746289, &
                2.728375485353322e4157, 359934.5632964715, 812.1021056757802, &
                9.50658460760105e4404, 391888.69843771093, 856.3971061616545, &
                1.257330119286705e4592, 410407.52673184, 890.9317826972928, &
                2.668793961338483e4856, 446018.32981189544, 937.8508450970188, &
                4.88950524334203e5047, 464880.78678410634, 973.0407338843653, &
                1.926935676259676e8492, 1.6191065828801782e6, 1509.0160735228724, &
                7.62592323391856e8826, 1.6922527096353143e6, 1565.4750473829606, &
                2.239327327624481e9183, 1.7749614762147116e6, 1625.1904715745486, &
                9.05926490739204e9532, 1.8533343231962149e6, 1683.8829567125488, &
                3.544284067886673e10124, 2.0520766046332212e6, 1778.3901411039799, &
                8.79913934144879e10507, 2.1461086926180813e6, 1842.0246866389446, &
                7.68140234303389e10940, 2.2648974578589206e6, 1912.948588485131, &
                6.09177778137702e11304, 2.3467258483170266e6, 1973.7011217758152, &
                6.44059395873296e26439, 1.951787734311023e7, 4076.8443167399187, &
                2.53185188720932e27320, 2.026995512633456e7, 4206.4787622659815 /), &
                                                                     (/ 3, 54 /) )

integer, parameter :: p_ion_max = 54

! This class overrides the ionization calculations for use with the PGC algorithm

type, extends(t_neutral) :: t_neutral_pgc

  ! no additional members

contains

  procedure :: adk_field_rates => adk_field_rates_pgc
  procedure :: read_input => read_input_neutral_pgc

end type t_neutral_pgc


public :: t_neutral_pgc

contains

!-----------------------------------------------------------------------------------------
! Get ionization rates on each cell for PGC algorithm
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_pgc( this, emf )
!-----------------------------------------------------------------------------------------
  use m_emf_define, only : t_emf
  use m_emf_pgc_define, only : t_emf_pgc

  class( t_neutral_pgc ), intent(inout) :: this
  class( t_emf ), intent(in) :: emf

  select type( emf )

  class is ( t_emf_pgc )

    ! adk_pgc_field_rates_2d only works with t_emf_pgc objects, and since emf was declared as
    ! t_emf we must use a select type construct to check that the object is indeed a
    ! t_emf_pgc instance

    select case ( p_x_dim )

      case (2)
        call adk_pgc_field_rates_2d( this, emf%env_mod, emf%omega )

      case (3)
        call adk_pgc_field_rates_3d( this, emf%env_mod, emf%omega )

      case default

        print *, "[pgc] :: ionization is not supported for 1d"
        call abort_program( p_err_invalid )

    end select

  class default

    ! This must never happen, adk_field_rates_pgc must always be called with a
    ! t_emf_pgc object
    call abort_program( p_err_invalid )

  end select


end subroutine adk_field_rates_pgc
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 2D for PGC
!-----------------------------------------------------------------------------------------
subroutine adk_pgc_field_rates_2d( this, a, omega )

  implicit none

  class( t_neutral_pgc ), intent(inout) :: this ! neutral
  type( t_vdf ), intent(in) :: a ! Laser envelope
  real( p_k_fld ), intent(in) :: omega ! Laser frequency

  integer :: i, j, l
  real(p_k_fld) :: den_center, eij

  if ( this%species1%pos_type == p_cell_low ) then

      do j = 1, a % nx_( 2 )
       do i = 1, a % nx_( 1 )

          den_center = this%multi_ion%f2(this%neut_idx,i,j)

          if (den_center > this%den_min) then
             ! Interpolate at center of the cell
             eij = 0.25_p_k_fld * omega * this%omega_p * 1.704e-12 * &
                                  ( a%f2(1,i,j ) + a%f2(1,i+1,j ) + &
                                    a%f2(1,i,j+1) + a%f2(1,i+1,j+1))
             if (eij > this%e_min) then
                ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
                do l = 1, this%multi_max
                   this%w%f2(l,i,j) = &
                      this%rate_param(1,l) * &
                      eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
                      this%omega_p
                enddo
             endif
          else
             this%w%f2(1:this%multi_max,i,j) = 0.0
          endif

       enddo
     enddo

  else

       do j = 1, a % nx_( 2 )
       do i = 1, a % nx_( 1 )

          den_center = this%multi_ion%f2(this%neut_idx,i,j)

          if (den_center > this%den_min) then
             ! Interpolate at cell corner
             eij = 0.25_p_k_fld * omega * this%omega_p * 1.704e-12 * a%f2(1,i,j)

             if (eij > this%e_min) then
                ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
                do l = 1, this%multi_max
                   this%w%f2(l,i,j) = &
                      this%rate_param(1,l) * &
                      eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
                      this%omega_p
                enddo
             endif
          else
             this%w%f2(1:this%multi_max,i,j) = 0.0
          endif

       enddo
     enddo

  endif


end subroutine adk_pgc_field_rates_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 2D for PGC
!-----------------------------------------------------------------------------------------
subroutine adk_pgc_field_rates_3d( this, a, omega )

  implicit none

  class( t_neutral_pgc ), intent(inout) :: this ! neutral
  type( t_vdf ), intent(in) :: a ! Laser envelope
  real( p_k_fld ), intent(in) :: omega ! Laser frequency

  integer :: i, j, k, l
  real(p_k_fld) :: den_center, eijk

  if ( this%species1%pos_type == p_cell_low ) then

    do k = 1, a % nx_( 3 )
      do j = 1, a % nx_( 2 )
        do i = 1, a % nx_( 1 )

          den_center = this%multi_ion%f3(this%neut_idx,i,j,k)

          if (den_center > this%den_min) then
            ! Interpolate at center of the cell
            eijk = 0.125_p_k_fld * omega * this%omega_p * 1.704e-12 * &
              ( a%f3( 1, i , j , k ) + a%f3( 1, i+1, j , k ) + &
                a%f3( 1, i , j+1, k ) + a%f3( 1, i+1, j+1, k ) + &
                a%f3( 1, i , j , k+1 ) + a%f3( 1, i+1, j , k+1 ) + &
                a%f3( 1, i , j+1, k+1 ) + a%f3( 1, i+1, j+1, k+1 ) )
            if (eijk > this%e_min) then
              ! w = r1 * eijk^-r3 * EXP(-r2/eijk) * 1/wp
              do l = 1, this%multi_max
                this%w%f3(l,i,j,k) = &
                this%rate_param(1,l) * &
                eijk**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eijk)) / &
                this%omega_p
              enddo
            endif
          else
            this%w%f3(1:this%multi_max,i,j,k) = 0.0
          endif

        enddo
      enddo
    enddo

  else

    do k = 1, a % nx_( 3 )
      do j = 1, a % nx_( 2 )
        do i = 1, a % nx_( 1 )

          den_center = this%multi_ion%f3(this%neut_idx,i,j,k)

          if (den_center > this%den_min) then
            ! Interpolate at cell corner
            eijk = omega * this%omega_p * 1.704e-12 * a%f3(1,i,j,k)

            if (eijk > this%e_min) then
              ! w = r1 * eijk^-r3 * EXP(-r2/eijk) * 1/wp
              do l = 1, this%multi_max
                this%w%f3(l,i,j,k) = &
                this%rate_param(1,l) * &
                eijk**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eijk)) / &
                this%omega_p
              enddo
            endif
          else
            this%w%f3(1:this%multi_max,i,j,k) = 0.0
          endif

        enddo
      enddo
    enddo

  endif

end subroutine adk_pgc_field_rates_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_input_neutral_pgc( this, input_file, species, def_name, mov_ions, &
                               periodic, if_move, grid, dt, sim_options )
!-----------------------------------------------------------------------------------------

  use stringutil, only : lowercase
  use m_cross
  use m_diag_neutral

  use m_input_file

  use m_grid_define

  implicit none

  class( t_neutral_pgc ), intent(inout) :: this
  class( t_species ), pointer :: species
  class( t_input_file ), intent(inout) :: input_file
  character(len = *), intent(in) :: def_name !predefined name
  type( t_options ), intent(in) :: sim_options
  logical, intent(in) :: mov_ions !if this has moving ions

  ! These are required for the associated particle species
  logical, dimension(:), intent(in) :: periodic, if_move
  class( t_grid ), intent(in) :: grid
  real(p_double), intent(in) :: dt

  character(len=p_max_spname_len) :: name ! neutral name

  ! parameters for ionization
  real(p_k_fld) :: den_min, e_min ! minimal values for profile dens and e-fld to consider ionisation
  real(p_double), dimension(3,p_ion_max) :: ion_param !custom rate parameters

  character(20) :: neutral_gas ! name of the neutral gas (selects set of hardcoded rate parameters
  logical :: if_tunnel, if_impact, inject_line
  integer :: multi_max, multi_min
  integer :: i
  character(len = p_max_spname_len) :: spname

  namelist /nl_neutral/ name, &
                        neutral_gas, ion_param, den_min, e_min, &
                        multi_max, multi_min, if_tunnel, if_impact, inject_line

  namelist /nl_neutral_mov_ions/ name, &
                        neutral_gas, ion_param, den_min, e_min, &
                        multi_max, multi_min, if_tunnel, if_impact, inject_line
  integer :: ierr


  name = "Neutral" ! neutral name

  neutral_gas = "H"
  den_min = 0.0_p_k_part
  e_min = 1.0e-6_p_k_fld

  ion_param = 0.0_p_k_fld ! custom ionization parameters this coresponds to a maximum level of 0

  if_tunnel = .true.
  if_impact = .false.

  ! multi-level by default
  multi_max = p_ion_max
  multi_min = 0

  name = trim(adjustl(def_name))

  inject_line = .true.


  if (.not. mov_ions) then

     ! read neutral without moving ions
     call get_namelist( input_file, "nl_neutral", ierr )
     if ( ierr /= 0 ) then
       if (ierr < 0) then
         print *, "Error reading neutral parameters"
       else
         print *, "Error: neutral parameters missing"
       endif
       print *, "aborting..."
       stop
     endif

     read (input_file%nml_text, nml = nl_neutral, iostat = ierr)
     if (ierr /= 0) then
       print *, "Error reading neutral parameters"
       print *, "aborting..."
       stop
     endif

     this%if_mov_ions = .false.

  else ! .not. mov_ions

     ! read neutral with moving ions
     call get_namelist( input_file, "nl_neutral_mov_ions", ierr )

     if ( ierr /= 0 ) then
       if (ierr < 0) then
         print *, "Error reading neutral moving ions parameters"
       else
         print *, "Error: neutral moving ions parameters missing"
       endif
       print *, "aborting..."
       stop
     endif

     read (input_file%nml_text, nml = nl_neutral_mov_ions, iostat = ierr)
     if (ierr /= 0) then
       print *, "Error reading neutral moving ions parameters"
       print *, "aborting..."
       stop
     endif

     this%if_mov_ions = .true.

  endif

  if (disp_out(input_file)) then
    if (mpi_node()==0) print *,"   Neutral name : ", trim(name)
  endif

  this%name = name
  this%den_min = den_min
  this%e_min = e_min

  this%if_tunnel = if_tunnel
  this%if_impact = if_impact
  this%multi_max = multi_max
  this%multi_min = multi_min

  this%inject_line = inject_line

  ! Select ionization rate calculation parameters

  select case ( lowercase(neutral_gas) )
     case ("h") ! Hydrogen
       call set_ion_parameters(this%rate_param, this%multi_max, H_param_pgc)
     case ("li") ! Lithium
       call set_ion_parameters(this%rate_param, this%multi_max, Li_param_pgc)
     case ("na") ! Sodium
       call set_ion_parameters(this%rate_param, this%multi_max, Na_param_pgc)
     case ("k") ! Potassium
       call set_ion_parameters(this%rate_param, this%multi_max, K_param_pgc)
     case ("rb") ! Rubidium
       call set_ion_parameters(this%rate_param, this%multi_max, Rb_param_pgc)
     case ("cs") ! Cesium
       call set_ion_parameters(this%rate_param, this%multi_max, Cs_param_pgc)
     case ("he") ! Helium
       call set_ion_parameters(this%rate_param, this%multi_max, He_param_pgc)
     case ("ar") ! Argon
       call set_ion_parameters(this%rate_param, this%multi_max, Ar_param_pgc)
     case ("n") ! Nitrogen
       call set_ion_parameters(this%rate_param, this%multi_max, N_param_pgc)
     case ("o") ! Oxygen
       call set_ion_parameters(this%rate_param, this%multi_max, O_param_pgc)
     case ("c") ! Carbon
       call set_ion_parameters(this%rate_param, this%multi_max, C_param_pgc)
     case ("ne") ! Neon
       call set_ion_parameters(this%rate_param, this%multi_max, Ne_param_pgc)
     case ("xe") ! Xenon
       call set_ion_parameters(this%rate_param, this%multi_max, Xe_param_pgc)
     case ("custom")
       do i=1, p_ion_max
         ! not physical, means this is multi_max + 1 (end token)
         if (ion_param(1,i) < 1.0_p_k_fld) exit
       enddo ! now i will be multi_max + 1

       if ( i == 1 ) then
          print *, "   Error reading neutral parameters"
          print *, "   You must specify the parameters for the custom gas: ion_param!"
          print *, "   aborting..."
          stop
       endif

       this%multi_max = i-1

       call set_ion_parameters(this%rate_param, this%multi_max, ion_param)

     case default

       print *, "   Error reading neutral parameters"
       print *, "   Not a valid neutral gas for PGC algorithm -> ", trim(neutral_gas)
       print *, "   Available gases:"
       print *, "    - H"
       print *, "    - Li"
       print *, "    - Na"
       print *, "    - K"
       print *, "    - Rb"
       print *, "    - Cs"
       print *, "    - He"
       print *, "    - Ar"
       print *, "    - N"
       print *, "    - O"
       print *, "    - C"
       print *, "    - Ne"
       print *, "    - Xe"
       print *, "    - Custom"
       print *, "   aborting..."
       stop

  end select



  if (this%multi_max > p_ion_max) then
     print *, "   Error reading neutral parameters"
     print *, "   The number of ionization levels selected is too high."
     print *, "   The maximum is ", p_ion_max
     print *, "   Please recompile the code with a larger p_ion_max parameter"
     print *, "   aborting..."
     stop
  endif


  if (this%multi_min > this%multi_max) then
     print *, "   Error reading neutral parameters"
     print *, "   The specified multi_min is larger than multi_max."
     print *, "   aborting..."
     stop
  endif

  this%neut_idx = this%multi_max + 2
  this%ion_idx = this%multi_max + 3

  ! if multi-level
  if ( this%multi_max > 1 ) then

     ! impact ionization must be off
     if ( this%if_impact ) then

        print *, "   Error reading neutral parameters"
        print *, "   Multi-level ionization activated, impact ionization must be off."
        print *, "   aborting..."
        stop

     endif

     ! field ionization must be on
     if (.not. this%if_tunnel) then

        print *, "   Error reading neutral parameters"
        print *, "   Multi-level ionization activated, tunnel ionization must be on."
        print *, "   aborting..."
        stop

     endif

     ! moving ions must be off
     if (this%if_mov_ions) then

        print *, "   Error reading neutral parameters"
        print *, "   Multi-level ionization activated, moving ions must be off"
        print *, "   aborting..."
        stop

     endif

  endif ! multi_max > 0

  if (this%if_impact) then
     call read_nml(this%cross_section, input_file, neutral_gas)
  endif

  ! read neutrals profile
  call this % den_neutral % read_input( input_file, grid % coordinates )

  ! read neutral diagnostics
  call read_nml( this%diag, input_file )

  ! Read associated electrons species
  if (disp_out(input_file)) then
    if (mpi_node()==0) print *,"    - reading associated electrons configuration..."
  endif
  spname = trim(adjustl(name)) // " electrons"

  call species % read_input( input_file, spname, periodic, if_move, grid, &
        dt, .false., sim_options )

  ! store pointer to electron species
  this%species1 => species

  ! If using moving ions read additional ion species
  if (this%if_mov_ions) then
    if (disp_out(input_file)) then
      if (mpi_node()==0) print *,"    - reading associated ions configuration..."
    endif
    spname = trim(adjustl(name)) // " ions"

    species => species % next
    call species % read_input( input_file, spname, periodic, if_move, grid, &
          dt, .false., sim_options )

    ! store pointer to ion species
    this%species2 => species
  endif

  contains

     subroutine set_ion_parameters(rate_param, multi_max, neut_param)
     !! Here we write the values of neut_param into rate_param
     !! truncating multi_max if necessary
        real(p_double), dimension(:,:), pointer :: rate_param
        integer, intent(out) :: multi_max
        real(p_double), dimension(:,:), intent(in) :: neut_param

        multi_max = min(multi_max, size(neut_param, 2))

        call alloc(this%rate_param, (/3,this%multi_max/),"pgc/os-neutral-pgc.f03",809)

        rate_param = neut_param(:,1:multi_max)

     end subroutine set_ion_parameters

end subroutine read_input_neutral_pgc
!-----------------------------------------------------------------------------------------


end module m_neutral_pgc
