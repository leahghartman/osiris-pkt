# 1 "os-neutral.f03"
# 1 "<built-in>" 1
# 1 "<built-in>" 3
# 467 "<built-in>" 3
# 1 "<command line>" 1
# 1 "<built-in>" 2
# 1 "os-neutral.f03" 2
!-----------------------------------------------------------------------------------------
!
! Currently disabled:
! - BSI and BSI random ionization models
! - Vacuum ionization
! - Impact ionization
!
!-----------------------------------------------------------------------------------------

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
# 11 "os-neutral.f03" 2
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
# 12 "os-neutral.f03" 2



module m_neutral

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
# 18 "os-neutral.f03" 2

use m_parameters

use m_vdf
use m_vdf_define

use m_node_conf
use m_grid_define
use m_space

use m_species
use m_species_define
use m_psource_std

use m_emf_define

use m_diag_neutral

use stringutil
use m_math
use m_random

use m_cross

use m_restart

implicit none

! restrict access to things explicitly declared public
private

! string to id restart data
character(len=*), parameter :: p_neutral_rst_id = "neutral rst data - 0x0003"

!-----------------------------------------------------------------------------------------
!
! Structure of the multi_ion vdf
!
! idx ranges comments
! +-------VDF-------+
! 1 | neutral density | [0,1] Initilaized with 1.0d0
! +-----------------+
! 2 | 1. ion. level | [0,1] Initilaized with 0.0d0
! +-----------------+
! 3 | 2. ion. level | [0,1]
! : +-----------------+
! : :::::::::::::::::::
! n-2 +-----------------+
! = m+1 | m. ion. level | [0,1]
! neut_idx +-----------------+
! = n-1 | neutral profile | > 0 Initialized with profile den_neutral
! ion_idx +-----------------+ +------VDF--------+
! = n | total charge | [0-m] => | ion-level-old |
! +-----------------+ +-----------------+
!
! - m = multi_max
! - n = multi_max + 3 where multi_max is the maximum ionization level
! - all values are normalized to the local neutral profile value (index n-1)
! - values at index neut_idx are initialized with the profile and serve for
! the normalization
!
!-----------------------------------------------------------------------------------------

! class declaration

type :: t_neutral

  ! neutral number - internal id (index in neutral array of particle object)
  integer :: neutral_id

  character(len=p_max_spname_len) :: name
  real(p_k_fld) :: omega_p, den_min, e_min
  real(p_double), dimension(:,:), pointer :: rate_param => null()

  ! neutral profile (density profile for neutral species from deck)
  type( t_psource_std ) :: den_neutral

  ! species to receive produced particles
  class( t_species ), pointer :: species1 => null() ! electrons, allways associated
  class( t_species ), pointer :: species2 => null() ! moving ions

  logical :: if_mov_ions

  ! Ionization level previous
  type( t_vdf ) :: ion_level_old

  ! Multi-level densities
  type( t_vdf ) :: multi_ion

  ! ionization rates
  type( t_vdf ) :: w

  ! diagnostic set up for neutrals
  type( t_diag_neutral ) :: diag

  ! Impact ionization
  type(t_cross) :: cross_section ! cross_section for each gas
  logical :: if_tunnel ! turn on tunnel_ionization or not
  logical :: if_analytic ! turn on semi-analytic solution
  logical :: if_impact ! turn on impact_ionization or not,
                   ! used in move and update boundary of grid density

  ! Multi level = multi_max + 2
  integer :: multi_max, ion_idx, neut_idx
  integer :: multi_min ! Ionization level after which particles are injected (0 by default)

  ! Particle injection method: line or random
  logical :: inject_line

contains

  procedure :: adk_field_rates => adk_field_rates_neutral
  procedure :: read_input => read_input_neutral
  procedure :: init => init_neutral
  procedure :: cleanup => cleanup_neutral
  procedure :: advance => advance_density
  procedure :: inject => inject_particles
  procedure :: ionize => ionize_neutral
  procedure :: move_window => move_window_neutral
  procedure :: reshape_obj => reshape_neutral
  procedure :: add_particle_load => add_particle_load_neutral
  procedure, nopass :: set_neutden_values => set_neutden_values
  procedure :: update_boundary => update_boundary_neutral
  procedure :: report => report_neutral
  procedure :: restart_write => restart_write_neutral
  procedure :: restart_read => restart_read_neutral

end type t_neutral

! parameters defining ionization models

integer, parameter :: p_adk = 1 ! Ammosove-Delone-Krainov (ADK)

! Currently disabled
!integer, parameter :: p_bsi = 2 ! Barrier Supression
!integer, parameter :: p_bsirand = 3 ! Monte Carlo Barrier Supression

! used for selectors in set_species_neutral and sp_id_neutral
integer, parameter :: p_neutral_electrons = 1
integer, parameter :: p_neutral_ions = 2

integer, parameter :: p_ion_max = 54

! parameter to controll what value we want from value-at_position
integer, parameter :: p_vap_den = 1
integer, parameter :: p_vap_ion_level_old = 2
integer, parameter :: p_vap_ion_level = 3

! Ionization Parameters for ADK model

! Data taken from:
!

! 2010 CODATA Internationally recommended values of the Fundamental Physical Constants
!

! http://physics.nist.gov/cuu/Constants/
!
! NIST Atomic Spectra Database Ionization Energies Data
! Kramida, A., Ralchenko, Yu., Reader, J., and NIST ASD Team (2014). NIST Atomic Spectra
! Database (ver. 5.2), [Online]. Available: http://physics.nist.gov/asd [2014, December].

! National Institute of Standards and Technology, Gaithersburg, MD
!
! http://physics.nist.gov/PhysRefData/ASD/ionEnergy.html

real(p_double), parameter, dimension(3,1) :: H_param = &
    reshape( (/ 8.522542995398661d19, 342.53947239007687d0, 1.0005337056631487d0 /), &
                                                                     (/ 3, 1 /) )
real(p_double), parameter, dimension(3,3) :: Li_param = &
    reshape( (/ 3.460272990838495d21, 85.51998980232813d0, 2.1770706138013733d0, &
                3.6365138642921554d20, 4493.713340713575d0, 0.6964625952167312d0, &
                2.0659396971422902d22, 9256.32561931876d0, 0.9999745128918196d0 /), &
                                                                     (/ 3, 3 /) )
real(p_double), parameter, dimension(3,11) :: Na_param = &
    reshape( (/ 4.1003766070018457d21, 79.58018178025945d0, 2.2542264609655485d0, &
                6.646980062859022d21, 2221.166444870671d0, 1.1456178592204318d0, &
                1.4450031005839113d24, 4140.270846854678d0, 1.6151368370485608d0, &
                1.1938031389836218d26, 6722.159287368323d0, 1.966696637087852d0, &
                2.363515645181666d27, 11121.962149645387d0, 2.1353951423226762d0, &
                6.533560672234857d28, 15439.488740548195d0, 2.3727930406073865d0, &
                1.2810219583986103d30, 20565.362378303853d0, 2.5763114425119764d0, &
                7.256788653171509d30, 29333.019535393767d0, 2.630951747911784d0, &
                1.2248965512782968d32, 35469.501316660346d0, 2.8341894069142315d0, &
                3.7691814856879473d24, 383083.54115896183d0, 0.9273100377897379d0, &
                1.3468138246659564d25, 457288.7579154525d0, 0.998535765107404d0 /), &
                                                                     (/ 3, 11 /) )
real(p_double), parameter, dimension(3,19) :: K_param = &
    reshape( (/ 7.216074209379994d21, 61.77481632066079d0, 2.5408885716705307d0, &
                9.406761452174264d22, 1214.8511771634942d0, 1.6236448747182823d0, &
                6.779724188614074d25, 2117.481312873944d0, 2.2701233926960898d0, &
                1.811619844825418d28, 3247.768300305309d0, 2.7807769758524734d0, &
                1.0119877391508012d30, 5133.580528079887d0, 3.0570711651756124d0, &
                9.049414272968163d31, 6769.504077462043d0, 3.4396463834481468d0, &
                5.0242638226357555d33, 8706.961891869527d0, 3.7627586361816903d0, &
                2.578271439282718d34, 13165.223159585394d0, 3.74238337717886d0, &
                7.967483576694055d35, 15924.655928532302d0, 4.007280316055681d0, &
                3.196309283214214d30, 77214.0625914301d0, 2.28713302745643d0, &
                2.3563335866942047d31, 91884.31431199581d0, 2.4121515077162674d0, &
                1.4984572622004427d32, 108273.28574947864d0, 2.524170311816701d0, &
                5.8647348952001106d32, 130518.49747304033d0, 2.5873086510013654d0, &
                3.0398926235765757d33, 150600.54746568995d0, 2.6832845094542193d0, &
                1.4304663067720226d34, 172552.5863001384d0, 2.7713790342861593d0, &
                3.4929050463473155d34, 205618.1565293413d0, 2.7944506601793733d0, &
                1.85797716283929d35, 227299.68713051578d0, 2.8991091229466286d0, &
                9.046983127409002d25, 2.1387084800213743d6, 0.9555610969933448d0, &
                1.9843676239339244d26, 2.3674584145525303d6, 0.9954563352900228d0 /), &
                                                                     (/ 3, 19 /) )
real(p_double), parameter, dimension(3,37) :: Rb_param = &
    reshape( (/ 8.142038491103064d21, 58.31683597165995d0, 2.6095364358731663d0, &
                2.5535438019000036d23, 973.8041035058765d0, 1.824374032909187d0, &
                2.7000601553750055d26, 1679.526944428912d0, 2.5327128978496343d0, &
                9.685790844349714d28, 2576.2205018430946d0, 3.0842758054317407d0, &
                1.0748558554029407d31, 3864.21631810689d0, 3.459973470196984d0, &
                1.1603503988608595d33, 5155.954467164836d0, 3.861433025510535d0, &
                7.592811710363505d34, 6695.067652968157d0, 4.198712465440701d0, &
                3.228036166231273d35, 10453.804888994271d0, 4.121312432452735d0, &
                1.2294619485882202d37, 12628.28840362558d0, 4.409748725556149d0, &
                4.594628560441955d34, 31512.226042057184d0, 3.4315574796699924d0, &
                5.099862812557592d35, 37844.42220681386d0, 3.5860783832270737d0, &
                3.66051737982471d36, 45883.04950656357d0, 3.691876466328285d0, &
                2.4532542930217856d37, 54647.11820807218d0, 3.795166356205491d0, &
                1.6954926367367167d38, 63691.70382403349d0, 3.9070049315032813d0, &
                5.741143312291227d38, 76830.35737138857d0, 3.938894193953593d0, &
                3.1772643710662465d39, 88109.210912276d0, 4.033023518460716d0, &
                1.547282384799606d40, 100644.25410629017d0, 4.115663943965691d0, &
                6.970575981163561d40, 114246.76199511335d0, 4.192469483287085d0, &
                3.2341579464837933d41, 128139.88477866953d0, 4.275232097540516d0, &
                7.882013953979101d40, 171375.41143787434d0, 4.03999269039215d0, &
                4.383609296869856d41, 186065.5873110981d0, 4.148886949420316d0, &
                2.0604060864046853d42, 202832.3654692189d0, 4.2411470329012415d0, &
                7.118253469778463d42, 223834.59618026362d0, 4.302347413111842d0, &
                3.050179866504092d43, 242444.82069622068d0, 4.387529805709947d0, &
                1.6607075650654932d44, 257754.38432565072d0, 4.4986245825104865d0, &
                1.7436005002448962d44, 299172.63990271214d0, 4.441460377264627d0, &
                7.906799692579144d44, 318173.5516199869d0, 4.535945678561607d0, &
                3.491233570113727d36, 1.1980646266135697d6, 2.690178412171373d0, &
                8.395498820952659d36, 1.283768228078111d6, 2.7349535459621346d0, &
                1.8409179480873614d37, 1.379410655833714d6, 2.772299305036374d0, &
                4.105265551141743d37, 1.4754721916179487d6, 2.811542488265903d0, &
                6.475886707043335d37, 1.6096041792672188d6, 2.822020648909448d0, &
                1.3258054862939667d38, 1.7203230250054656d6, 2.855020456297321d0, &
                2.0212930745786477d38, 1.8686117460094132d6, 2.863865259228529d0, &
                4.89932668666366d38, 1.9644384793285115d6, 2.911752011774463d0, &
                2.947497999154891d27, 1.6918556303257223d7, 0.9628968654748764d0, &
                4.4985357459233346d27, 1.784120018639596d7, 0.9820279890950028d0 /), &
                                                                     (/ 3, 37 /) )
real(p_double), parameter, dimension(3,55) :: Cs_param = &
    reshape( (/ 1.0067145666103756d22, 52.487440693139064d0, 2.7385019836913114d0, &
                7.850381786372777d23, 761.2255911152253d0, 2.0660197638382733d0, &
                1.2352289567990247d27, 1306.4304282001383d0, 2.841273587760102d0, &
                8.285899003417227d29, 1926.1050102158758d0, 3.5000370081744165d0, &
                1.3943473132068207d32, 2862.591109127988d0, 3.9290850597544393d0, &
                1.6090811519096703d34, 3923.686937109447d0, 4.324790783549418d0, &
                1.2197188038907911d36, 5155.954467164836d0, 4.671671863095624d0, &
                7.772313542644456d36, 7890.401181613261d0, 4.624790779520253d0, &
                3.3985065640224123d38, 9616.423940737586d0, 4.924077392169907d0, &
                4.8014562020006024d36, 21279.60493046487d0, 4.051203833568829d0, &
                1.3222551053571997d38, 24294.678658917313d0, 4.31624572540669d0, &
                1.6574025779488334d39, 28803.02083574891d0, 4.479631003881846d0, &
                1.9601994628561338d40, 33560.16146953232d0, 4.641372183771167d0, &
                2.3511108954259957d41, 38371.42266737041d0, 4.809982392624659d0, &
                1.4586939574159844d42, 45111.915987085304d0, 4.898074435696151d0, &
                1.3185543170683242d43, 51000.251990244644d0, 5.039189024375879d0, &
                1.0773662920890716d44, 57332.69476301885d0, 5.1711232098416176d0, &
                8.026118769310927d44, 64123.5110697902d0, 5.294812832444929d0, &
                6.132316537917955d45, 70939.4319707329d0, 5.424516464841734d0, &
                5.978655302384787d44, 99641.15986441393d0, 5.038556747870522d0, &
                4.5918747463503677d45, 107758.9238855643d0, 5.177095189338379d0, &
                2.7596043753895266d46, 117405.54836165802d0, 5.288917363302168d0, &
                1.8217938026575976d47, 126509.84906429503d0, 5.4131166262677d0, &
                5.5986254360520656d47, 141146.37923668668d0, 5.45214435668305d0, &
                3.043832153610278d48, 151964.45007995746d0, 5.557557219247648d0, &
                3.2582850917570695d48, 176802.88895687615d0, 5.484248396924973d0, &
                1.5867266730381928d49, 189395.76723911107d0, 5.580967644950393d0, &
                3.730087044273447d43, 433902.22052782355d0, 4.176983985293046d0, &
                1.0956515613976363d44, 467015.92197585426d0, 4.232029490341937d0, &
                2.962740161056216d44, 503077.5298274912d0, 4.279899983800397d0, &
                7.314590621543037d44, 542663.0206076808d0, 4.3198701538914825d0, &
                1.8934246830379147d45, 581882.5146795525d0, 4.365221097852297d0, &
                4.566171043429717d45, 624310.0815640229d0, 4.404596034427387d0, &
                9.012250277206341d45, 674813.9638445469d0, 4.425839936597809d0, &
                2.1956241197801815d46, 719340.5988659592d0, 4.467715823989638d0, &
                4.956319009281406d46, 767770.6575843783d0, 4.503108536433492d0, &
                1.2087646202679079d47, 814213.8503829134d0, 4.546320683994775d0, &
                4.2393229796198117d46, 949310.504573465d0, 4.412072205833181d0, &
                1.186145565511963d47, 994239.5071187183d0, 4.469534289314101d0, &
                3.0198512028692754d47, 1.0442356373731942d6, 4.518781954155857d0, &
                8.138731102071833d47, 1.091704816186824d6, 4.5735452197837745d0, &
                9.024016205173133d47, 1.1893000683619515d6, 4.548831754866617d0, &
                2.199734638466941d48, 1.24464668345491d6, 4.595460589784465d0, &
                2.4640844209070167d48, 1.3494631911312118d6, 4.57333408005987d0, &
                6.188697514903444d48, 1.4053401408214383d6, 4.6234319550448095d0, &
                2.0790590117445198d39, 4.877709249302955d6, 2.7966620649773297d0, &
                3.683305380261343d39, 5.094504263667927d6, 2.8233726484803303d0, &
                6.047476196981079d39, 5.337996286620714d6, 2.84442355319902d0, &
                1.0471081885558927d40, 5.568011076425189d6, 2.8697133664308336d0, &
                6.273346384614805d39, 6.175820914079348d6, 2.8146484527263547d0, &
                9.995132845792308d39, 6.454298659538151d6, 2.834157108994108d0, &
                1.3448971695667585d40, 6.804266629699228d6, 2.841129803690224d0, &
                2.4075380464334907d40, 7.0459254669033475d6, 2.8697175674767808d0, &
                1.6534330222238202d28, 5.8504976845801085d7, 0.9470580668607365d0, &
                2.1839164977677438d28, 6.072400966498238d7, 0.9586580582564923d0 /), &
                                                                     (/ 3, 55 /) )
real(p_double), parameter, dimension(3,2) :: He_param = &
    reshape( (/ 7.2207661763501d18, 832.809878216992d0, 0.48776427204592254d0, &
                2.7226733893691d21, 2742.1316798375965d0, 1.0000920088118899d0 /), &
                                                                     (/ 3, 2 /) )
real(p_double), parameter, dimension(3,18) :: Ar_param = &
    reshape( (/ 4.583155091154178d19, 427.3616288942586d0, 0.858307468844768d0, &
                2.347069522730368d23, 992.0665907024354d0, 1.8069357289301329d0, &
                1.931637188507353d26, 1775.9423286421918d0, 2.4675897958480384d0, &
                2.3004409750295724d28, 3141.43462246603d0, 2.8229627269683206d0, &
                3.468592152886133d30, 4422.602894394151d0, 3.263766733070865d0, &
                2.965675081118506d32, 5958.158993507227d0, 3.6326551045811497d0, &
                2.1249901067576396d33, 9479.178143644607d0, 3.629746649292696d0, &
                8.991061313075681d34, 11737.03857874398d0, 3.92742350307885d0, &
                7.611312705906995d29, 59342.85460567607d0, 2.229751529947494d0, &
                6.602179339057654d30, 71781.6325368281d0, 2.3680482459080237d0, &
                4.815600903407338d31, 85819.59318777094d0, 2.490706079578982d0, &
                2.055113377686043d32, 105186.65111959413d0, 2.558310025830369d0, &
                1.1945429165015238d33, 122591.38628416909d0, 2.663021333332675d0, &
                6.212435859768618d33, 141745.69057711778d0, 2.7584389081812857d0, &
                1.6060852787486324d34, 170916.68217921766d0, 2.783373249974532d0, &
                9.327897444360857d34, 190110.7378894357d0, 2.894937691775023d0, &
                6.6878740700780965d25, 1.8068758775178685d6, 0.9536895901987263d0, &
                1.5247916822996254d26, 2.0115330349032478d6, 0.9959340925415332d0 /), &
                                                                     (/ 3, 18 /) )
real(p_double), parameter, dimension(3,7) :: N_param = &
    reshape( (/ 6.441568477283914d19, 378.4956025896395d0, 0.9350660865589433d0, &
                1.4701035287257558d23, 1100.1258100528623d0, 1.711847679118431d0, &
                4.963404971976876d25, 2232.374603536868d0, 2.2130314611564854d0, &
                1.429510866064162d27, 4658.081063546513d0, 2.3525381223431663d0, &
                1.3012849307544791d29, 6615.849778853221d0, 2.7281285044063783d0, &
                2.1409091975134166d23, 88606.4475724123d0, 0.8838466900060211d0, &
                1.4212916247536375d24, 117682.27130021283d0, 0.9994495038042972d0 /), &
                                                                     (/ 3, 7 /) )
real(p_double), parameter, dimension(3,8) :: O_param = &
    reshape( (/ 8.471015773053202d19, 343.2810702373722d0, 0.9990920676510977d0, &
                4.659951179950614d22, 1421.770921017062d0, 1.4896379815554006d0, &
                1.375943467286961d25, 2781.3610873184325d0, 1.985966133301126d0, &
                1.4410401568574448d27, 4652.670876621719d0, 2.353837077475246d0, &
                2.1812887948673846d28, 8303.41171218499d0, 2.4562134926799617d0, &
                1.0859737168339012d30, 11088.09514588445d0, 2.766300924086769d0, &
                5.135265231597041d23, 137319.4144041063d0, 0.8991975448040144d0, &
                2.7649244731504574d24, 175715.8649323391d0, 0.999259077262769d0 /), &
                                                                     (/ 3, 8 /) )
real(p_double), parameter, dimension(3,6) :: C_param = &
    reshape( (/ 1.8809310466441196d20, 258.1083104918942d0, 1.1984440550739555d0, &
                5.509352426405289d23, 822.523017116882d0, 1.987881643259267d0, &
                4.572825176241475d25, 2263.6763526780164d0, 2.198152910400111d0, &
                9.835564484116034d27, 3537.9525592618193d0, 2.674447703635764d0, &
                7.536123777733801d22, 53034.28979013096d0, 0.8628040191262976d0, &
                6.587972799911127d23, 74090.46603762524d0, 0.9996157827449552d0 /), &
                                                                     (/ 3, 6 /) )
real(p_double), parameter, dimension(3,10) :: Ne_param = &
    reshape( (/ 1.2380498982127483d19, 684.0495717312216d0, 0.5886207199036224d0, &
                1.6861743227006676d22, 1790.87087951275d0, 1.3052851455522987d0, &
                4.012462791033353d24, 3450.2505445538854d0, 1.7789909006599105d0, &
                1.4251452589514637d26, 6544.999998027671d0, 1.9932260907792414d0, &
                6.670617838240815d27, 9689.667572133209d0, 2.282840577038655d0, &
                1.9391656057538128d29, 13557.847041347119d0, 2.522116807457813d0, &
                1.3881927947789378d30, 20383.797283570835d0, 2.586898525643846d0, &
                3.1075134088626144d31, 25254.466329809107d0, 2.8167466774677434d0, &
                2.1020658703689154d24, 282468.0585535283d0, 0.9200040384519763d0, &
                8.391062111505151d24, 343429.9456106982d0, 0.9988031600227598d0 /), &
                                                                     (/ 3, 10 /) )
real(p_double), parameter, dimension(3,54) :: Xe_param = &
    reshape( (/ 1.3775354155861185d20, 288.57589581574035d0, 1.1181793520578327d0, &
                1.550226815196774d24, 656.1909430408086d0, 2.2215830587335996d0, &
                2.2768338286515254d27, 1181.8699332253861d0, 2.971739666983438d0, &
                1.0222240210041819d30, 1872.6040313895708d0, 3.542491127683178d0, &
                2.1623381482112875d32, 2721.1636112083115d0, 4.013040352645607d0, &
                2.6991402989153737d34, 3721.3055273127543d0, 4.419620648906867d0, &
                2.4622862093269497d35, 5988.533565050575d0, 4.395610921298916d0, &
                1.5012400890751945d37, 7452.459544515593d0, 4.7328801021744935d0, &
                5.376912395698707d35, 16474.290444557057d0, 3.9509630414793238d0, &
                1.3066584018868181d37, 19611.207359780554d0, 4.190565038153051d0, &
                1.8497011055940987d38, 23674.858558706495d0, 4.362240665366397d0, &
                2.6681724805371856d39, 27815.543371271946d0, 4.54372243455411d0, &
                3.577640165134184d40, 32176.348524525496d0, 4.721112847931664d0, &
                2.7089693415351963d41, 38007.71401995982d0, 4.828456151144133d0, &
                2.6559087850025507d42, 43392.8782290532d0, 4.974953378510279d0, &
                2.1899377984971629d43, 49406.56731263479d0, 5.1034375243621675d0, &
                1.8601406162201934d44, 55468.870845425416d0, 5.239482352335317d0, &
                1.9313682862298425d45, 60908.78844944052d0, 5.403664925396409d0, &
                1.5426588652718157d44, 87869.02232374942d0, 4.982156226229307d0, &
                1.1590871169485781d45, 95909.52818579777d0, 5.115878165169938d0, &
                8.030370382706514d45, 104435.5470090333d0, 5.241935268399854d0, &
                5.382165292841963d46, 113200.23085498222d0, 5.365848876529456d0, &
                1.8217938026575976d47, 126509.84906429503d0, 5.4131166262677d0, &
                1.077849728298342d48, 136393.60209955153d0, 5.526234224434658d0, &
                1.145516874883936d48, 159811.20754403196d0, 5.44842538130848d0, &
                6.014132302364118d48, 171384.410225511d0, 5.5518758215736685d0, &
                1.759085923766389d43, 394064.2066635542d0, 4.1549467350722304d0, &
                5.305157185638429d43, 425345.21856075496d0, 4.211470223441139d0, &
                1.4911016576131055d44, 459078.08137395576d0, 4.2620127128282785d0, &
                3.746338544125226d44, 496648.9293453104d0, 4.3025833188760005d0, &
                1.0200735803192588d45, 533001.5075147987d0, 4.351821613265934d0, &
                2.4271623447664046d45, 574235.1060028574d0, 4.3889333565493684d0, &
                4.9683984737540496d45, 621542.8884528814d0, 4.412604825238724d0, &
                1.2484414815313788d46, 663477.9943518823d0, 4.456566977558425d0, &
                2.894901889446211d46, 709203.4692145847d0, 4.493643943279613d0, &
                7.180428844235488d46, 753475.3199613485d0, 4.53769337085137d0, &
                2.448645505910615d46, 882711.0267638017d0, 4.398978788274832d0, &
                7.090658854982765d46, 925001.57509592d0, 4.459072393688684d0, &
                1.860229435466672d47, 972223.3666808198d0, 4.510512643085584d0, &
                5.064124828052308d47, 1.0180486819614146d6, 4.5657012421984975d0, &
                5.9565042899701805d47, 1.1084286245427015d6, 4.545372106494971d0, &
                1.4905062982964533d48, 1.1608077808083333d6, 4.593864401786412d0, &
                1.660914855730416d48, 1.2615303421739575d6, 4.570386159811437d0, &
                4.2520741495048395d48, 1.3148814271135165d6, 4.621771988602169d0, &
                1.5790651376647514d39, 4.579524976873413d6, 2.7930489801003286d0, &
                2.8386327726136714d39, 4.786413465857762d6, 2.8206493017311245d0, &
                4.700763570275373d39, 5.0203491847052295d6, 2.8421056006525185d0, &
                8.200335703674566d39, 5.242021070951296d6, 2.8677440883669703d0, &
                5.20578174350226d39, 5.804149130601667d6, 2.8165062116744424d0, &
                8.391595879280377d39, 6.070112038854565d6, 2.8366647300643155d0, &
                1.135074333575459d40, 6.406097404576862d6, 2.8437495112068905d0, &
                2.0553597017159955d40, 6.637543043722892d6, 2.8730250674268727d0, &
                1.542287349895057d28, 5.520489369472207d7, 0.9483460319455781d0, &
                2.0492457824596215d28, 5.733209089671276d7, 0.9602460194755702d0 /), &
                                                                     (/ 3, 54 /) )

! General routines

interface list_algorithm
  module procedure list_algorithm_neutral
end interface

! declare things that should be public

public :: t_neutral
public :: setup, cleanup
public :: list_algorithm
public :: p_neutral_ions, p_neutral_electrons

interface alloc
  module procedure alloc_neutral
  module procedure alloc_1d_neutral
  module procedure alloc_bound_1d_neutral
end interface

interface freemem
  module procedure freemem_neutral
  module procedure freemem_1d_neutral
end interface

public :: alloc, freemem

contains

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine read_input_neutral( this, input_file, species, def_name, mov_ions, &
                               periodic, if_move, grid, dt, sim_options )
!-----------------------------------------------------------------------------------------
  use m_input_file

  implicit none

  class( t_neutral ), intent(inout) :: this
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
  logical :: if_tunnel, if_impact, if_analytic, inject_line
  integer :: multi_max, multi_min
  integer :: i
  character(len = p_max_spname_len) :: spname


  namelist /nl_neutral/ name, &
            neutral_gas, ion_param, den_min, e_min, &
            multi_max, multi_min, if_tunnel, if_impact, if_analytic, inject_line

  namelist /nl_neutral_mov_ions/ name, &
            neutral_gas, ion_param, den_min, e_min, &
            multi_max, multi_min, if_tunnel, if_impact, if_analytic, inject_line
  integer :: ierr

  neutral_gas = "H"
  den_min = 0.0_p_k_part
  e_min = 1.0e-6_p_k_fld

  ion_param = 0.0_p_k_fld ! custom ionization parameters this coresponds to a maximum level of 0

  if_tunnel = .true.
  if_impact = .false.
  if_analytic = .false.

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
  this%if_analytic = if_analytic
  this%multi_max = multi_max
  this%multi_min = multi_min

  this%inject_line = inject_line

  ! Select ionization rate calculation parameters

  select case ( lowercase(neutral_gas) )

  case ("h") ! Hydrogen
    call set_ion_parameters(this%rate_param, this%multi_max, H_param)
  case ("li") ! Lithium
    call set_ion_parameters(this%rate_param, this%multi_max, Li_param)
  case ("na") ! Sodium
    call set_ion_parameters(this%rate_param, this%multi_max, Na_param)

  case ("k") ! Potassium
    call set_ion_parameters(this%rate_param, this%multi_max, K_param)
  case ("rb") ! Rubidium
    call set_ion_parameters(this%rate_param, this%multi_max, Rb_param)
  case ("cs") ! Cesium
    call set_ion_parameters(this%rate_param, this%multi_max, Cs_param)

  case ("he") ! Helium
    call set_ion_parameters(this%rate_param, this%multi_max, He_param)

  case ("ar") ! Argon
    call set_ion_parameters(this%rate_param, this%multi_max, Ar_param)
  case ("n") ! Nitrogen
    call set_ion_parameters(this%rate_param, this%multi_max, N_param)

  case ("o") ! Oxygen
    call set_ion_parameters(this%rate_param, this%multi_max, O_param)

  case ("c") ! Carbon
    call set_ion_parameters(this%rate_param, this%multi_max, C_param)

  case ("ne") ! Neon
    call set_ion_parameters(this%rate_param, this%multi_max, Ne_param)

  case ("xe") ! Xenon
    call set_ion_parameters(this%rate_param, this%multi_max, Xe_param)
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
    print *, "   Not a valid neutral gas -> ", trim(neutral_gas)
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
    integer, intent(inout) :: multi_max
    real(p_double), dimension(:,:), intent(in) :: neut_param

    multi_max = min(multi_max, size(neut_param, 2))

    call alloc(this%rate_param, (/3,this%multi_max/),"os-neutral.f03",787)

    rate_param = neut_param(:,1:multi_max)

  end subroutine set_ion_parameters

end subroutine read_input_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine init_neutral( this, neutral_id, emf, nx_p_min, g_space, &
                          restart, restart_handle, sim_options )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_neutral ), intent(inout) :: this
  integer, intent(in) :: neutral_id
  class( t_emf),intent(in) :: emf
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  logical, intent(in) :: restart
  type( t_restart_handle ), intent(in) :: restart_handle
  type( t_options ), intent(in) :: sim_options

  integer, dimension(p_x_dim) :: nlbound, nubound
  integer, dimension(2,p_x_dim) :: lgc_num

  if ( restart ) then

    call this % restart_read( restart_handle )

  else

    ! set the internal neutral id (index in neutral array of particle object)
    this%neutral_id = neutral_id

    !set up cross section profile
    if (this%if_impact) then
      call setup( this%cross_section)
    endif

    ! check if associated species are ok
    if (.not. associated(this%species1)) then
      print *, "(*error*) species1 is not associated in setup_neutral, in os-neutral.f90"
      print *, "(*error*) unable to setup neutral"
      print *, "(*error*) bailing out..."
      call abort_program()
    endif

    if (.not. associated(this%species2) .and. this%if_mov_ions) then
      print *, "(*error*) species2 is not associated in setup_neutral, in os-neutral.f90"
      print *, "(*error*) unable to setup neutral"
      print *, "(*error*) bailing out..."
      call abort_program()
    endif

    ! Create ionization level and density vdf structures, based in vdf e
    ! copy = .false. causes vdf to be initialized with 0.
    call this%ion_level_old % new( emf%e, copy = .false., f_dim = 1)

    ! vdf with ion densities, neutral density (1st component) and

    ! total ion density (last component)
    ! copy = .false. causes vdf to be initialized with 0.

    call this % multi_ion % new( emf%e, copy = .false., f_dim = this%ion_idx)

    call this % w % new( emf%e , copy = .false. , f_dim = this%ion_idx)

    ! Fill neutral density (1st component of multi_ion) according to given neutral profile
    ! Since value at position uses guard cells we need to initialize them as well
    lgc_num = this%multi_ion%gc_num()
    nlbound = 1 - lgc_num( p_lower, 1:p_x_dim )
    nubound = emf%e%nx_(1:p_x_dim) + lgc_num( p_upper, 1:p_x_dim )
    call set_neutden_values( this%multi_ion, this%den_neutral, this%neut_idx, &
                              this%den_min, nx_p_min, g_space, nlbound, nubound)

  endif

  ! set values not present in restart info
  this%omega_p = sim_options%omega_p0

  call setup( this%diag, this%name, this%species1%interpolation )

end subroutine init_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine cleanup_neutral( this )
!-----------------------------------------------------------------------------------------

  implicit none

  class( t_neutral ), intent(inout) :: this

  call freemem(this%rate_param,"os-neutral.f03",884)
  call this % den_neutral % cleanup()

  call this%ion_level_old % cleanup()
  call this%multi_ion % cleanup()
  call this%w % cleanup()

  call cleanup( this%diag )

  if ( this%if_impact ) then
    call cleanup( this%cross_section )

  endif

end subroutine cleanup_neutral
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
! Ionizes a background source
!-----------------------------------------------------------------------------------------
subroutine ionize_neutral(this, species, emf, gdt, coordinates)
!-----------------------------------------------------------------------------------------

  implicit none

  class(t_neutral) :: this ! neutral array
  class(t_species), pointer :: species ! species array (impact ionization only)
  class( t_emf ), intent(in) :: emf ! EMF data
  real( p_double ), intent(in) :: gdt ! time step

  integer, intent(in) :: coordinates

  ! local variables

  real( p_k_fld ) :: dt
  logical :: if_impact

  ! executable statements

  dt = real( gdt, p_k_fld )

  if_impact = .false.

    ! Clear ionization rates
    call this%w %zero( )

      ! total ion density is the last component of multi_ion vdf
    call this%ion_level_old % copy( this%multi_ion, this%ion_idx)

    ! field ionization rates
    if ( this%if_tunnel ) then

        call this % adk_field_rates( emf )
      endif

      ! impact ionization rates
      ! if ( this(n)%if_impact ) call impact_rates( this(n), species, dt )

      ! Advance densities
      call advance_density( this, dt )

      ! Inject particles
      call inject_particles( this, coordinates)

end subroutine ionize_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ionization rates on each cell
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_neutral( this, emf )
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this

  class( t_emf ), intent(in) :: emf

  select case ( p_x_dim )

  case (1)
     call adk_field_rates_1d( this, emf%e_part )

  case (2)
     call adk_field_rates_2d( this, emf%e_part )

  case (3)
     call adk_field_rates_3d( this, emf%e_part )

  case default

     write(err_buf__,*) 'p_x_dim has the value:', p_x_dim;call err__("os-neutral.f03",975)
     call abort_program(p_err_invalid)

  end select

end subroutine adk_field_rates_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 1D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_1d( this, e )

  implicit none

  type( t_neutral ), intent(inout) :: this ! neutral

  type( t_vdf ), intent(in) :: e ! electrical field

  integer :: i, l
  real(p_k_fld) :: den_center, eij, e1, e2, e3

  if ( this%species1%pos_type == p_cell_low ) then

    !$omp parallel do private(den_center,e1,e2,e3,eij,l)
    do i = 1, e % nx_( 1 )

       den_center = this%multi_ion%f1(this%neut_idx,i)

        if (den_center > this%den_min) then
       ! Interpolate at center of the cell
       e1 = e%f1(1,i)
       e2 = 0.5_p_k_fld*(e%f1(2,i) + e%f1(2,i+1))

       e3 = 0.5_p_k_fld*(e%f1(3,i) + e%f1(3,i+1))

           eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12

         if (eij > this%e_min) then
        ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
        do l = 1, this%multi_max
         this%w%f1(l,i) = &
          this%rate_param(1,l) * &
          eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
          this%omega_p
        enddo
       endif
        else
           this%w%f1(1:this%multi_max,i) = 0.0d0
        endif

     enddo
     !$omp end parallel do

  else

    !$omp parallel do private(den_center,e1,e2,e3,eij,l)
    do i = 1, e% nx_( 1 )

       den_center = this%multi_ion%f1(this%neut_idx,i)

        if (den_center > this%den_min) then
       ! Interpolate at corner of the cell
       e1 = 0.5_p_k_fld*(e%f1(1,i-1) + e%f1(1,i))
       e2 = e%f1(2,i)

       e3 = e%f1(3,i)

           eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12

         if (eij > this%e_min) then
        ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
        do l = 1, this%multi_max
         this%w%f1(l,i) = &
          this%rate_param(1,l) * &
          eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
          this%omega_p
        enddo
       endif
        else
           this%w%f1(1:this%multi_max,i) = 0.0d0
        endif

     enddo
     !$omp end parallel do

  endif

end subroutine adk_field_rates_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 2D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_2d( this, e )

  implicit none

  type( t_neutral ), intent(inout) :: this ! neutral

  type( t_vdf ), intent(in) :: e ! electrical field

  integer :: i, j, l
  real(p_k_fld) :: den_center, eij, e1, e2, e3

  if ( this%species1%pos_type == p_cell_low ) then

    !$omp parallel do private(i,den_center,e1,e2,e3,eij,l)
    do j = 1, e% nx_( 2 )
     do i = 1, e% nx_( 1 )

      den_center = this%multi_ion%f2(this%neut_idx,i,j)

      if (den_center > this%den_min) then
       ! Interpolate at center of the cell
       e1 = 0.5_p_k_fld*( e%f2(1,i,j) + e%f2(1,i,j+1) )
       e2 = 0.5_p_k_fld*( e%f2(2,i,j) + e%f2(2,i+1,j) )

       e3 = 0.25_p_k_fld*( e%f2(3,i,j) + e%f2(3,i+1,j) + &
                           e%f2(3,i,j+1) + e%f2(3,i+1,j+1))

       eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12

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
       this%w%f2(1:this%multi_max,i,j) = 0.0d0
      endif

     enddo
     enddo
     !$omp end parallel do

  else

     !$omp parallel do private(i,den_center,e1,e2,e3,eij,l)
     do j = 1, e% nx_( 2 )
     do i = 1, e% nx_( 1 )

      den_center = this%multi_ion%f2(this%neut_idx,i,j)

      if (den_center > this%den_min) then
       ! Interpolate at cell corner
       e1 = 0.5_p_k_fld*(e%f2(1,i-1,j) + e%f2(1,i,j))
       e2 = 0.5_p_k_fld*(e%f2(2,i,j-1) + e%f2(2,i,j))

       e3 = e%f2(3,i,j)

       eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12

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
       this%w%f2(1:this%multi_max,i,j) = 0.0d0
      endif

     enddo
     enddo
     !$omp end parallel do

  endif

end subroutine adk_field_rates_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Get ADK ionization rates on each cell in 3D
!-----------------------------------------------------------------------------------------
subroutine adk_field_rates_3d( this, e )

  implicit none

  type( t_neutral ), intent(inout) :: this ! neutral

  type( t_vdf ), intent(in) :: e ! electrical field

  integer :: i, j, k, l
  real(p_k_fld) :: den_center, eij, e1, e2, e3

  if ( this%species1%pos_type == p_cell_low ) then

    !$omp parallel do private(j,i,den_center,e1,e2,e3,eij,l)
    do k = 1, e% nx_( 3 )
     do j = 1, e% nx_( 2 )
     do i = 1, e% nx_( 1 )

      den_center = this%multi_ion%f3(this%neut_idx,i,j,k)

      if (den_center > this%den_min) then
         ! Interpolate at center of the cell
         e1 = 0.25_p_k_fld*(e%f3(1,i,j,k) + e%f3(1,i,j+1,k) + &
                  e%f3(1,i,j,k+1) + e%f3(1,i,j+1,k+1))
         e2 = 0.25_p_k_fld*(e%f3(2,i,j,k) + e%f3(2,i+1,j,k) + &
                  e%f3(2,i,j,k+1) + e%f3(2,i+1,j,k+1))
         e3 = 0.25_p_k_fld*(e%f3(3,i,j,k) + e%f3(3,i+1,j,k) + &
                  e%f3(3,i,j+1,k) + e%f3(3,i+1,j+1,k))

         eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12

         if (eij > this%e_min) then
          ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
          do l = 1, this%multi_max
           this%w%f3(l,i,j,k) = &
            this%rate_param(1,l) * &
            eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
            this%omega_p
          enddo
         endif
      else
         this%w%f3(1:this%multi_max,i,j,k) = 0.0d0
      endif

     enddo
     enddo
     enddo
     !$omp end parallel do

  else

    !$omp parallel do private(j,i,den_center,e1,e2,e3,eij,l)
    do k = 1, e% nx_( 3 )
     do j = 1, e% nx_( 2 )
     do i = 1, e% nx_( 1 )

      den_center = this%multi_ion%f3(this%neut_idx,i,j,k)

      if (den_center > this%den_min) then
         ! Interpolate at cell corner
         e1 = 0.5_p_k_fld*(e%f3(1, i-1, j, k) + e%f3(1, i, j, k))
         e2 = 0.5_p_k_fld*(e%f3(2, i, j-1, k) + e%f3(2, i, j, k))

         e3 = 0.5_p_k_fld*(e%f3(3, i, j, k-1) + e%f3(3, i, j, k))

         eij = sqrt(e1**2 + e2**2 + e3**2)*this%omega_p*1.704e-12

         if (eij > this%e_min) then
          ! w = r1 * eij^-r3 * EXP(-r2/eij) * 1/wp
          do l = 1, this%multi_max
           this%w%f3(l,i,j,k) = &
            this%rate_param(1,l) * &
            eij**(-this%rate_param(3,l)) * exp(-this%rate_param(2,l)/(eij)) / &
            this%omega_p
          enddo
         endif
      else
         this%w%f3(1:this%multi_max,i,j,k) = 0.0d0
      endif

     enddo
     enddo
     enddo
     !$omp end parallel do

  endif

end subroutine adk_field_rates_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Advance ionization level density
!-----------------------------------------------------------------------------------------
subroutine advance_density(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  select case ( p_x_dim )

    case (1)
     if(.not. this%if_analytic) then
        call advance_density_1d(this, dt)
     else
        call advance_density_1d_analytic(this, dt)
     endif
    case (2)
     if(.not. this%if_analytic) then
        call advance_density_2d(this, dt)
     else
        call advance_density_2d_analytic(this, dt)
     endif
    case (3)
     if(.not. this%if_analytic) then
        call advance_density_3d(this, dt)
     else
        call advance_density_3d_analytic(this, dt)
     endif
    case default
     write(err_buf__,*) 'p_x_dim has the value:', p_x_dim;call err__("os-neutral.f03",1277)
     call abort_program(p_err_invalid)
   end select

end subroutine advance_density
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_1d(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  integer :: i, m, ion_idx
  real(p_k_fld) :: dens_temp, cons, inj
  logical :: shoot

  shoot = .false.

  ion_idx = this%ion_idx

  !$omp parallel do private(cons,shoot,m,inj,dens_temp)
  do i=1, this%multi_ion% nx_( 1 )

    if ( this%multi_ion%f1(this%neut_idx,i) > this%den_min ) then

  ! ionizing 0. level adding to 1. level
  cons = this%multi_ion%f1(1,i) * this%w%f1(1,i) * dt * &
       ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f1(1,i) * dt )

  ! overshoot
  if (cons > this%multi_ion%f1(1,i)) then
     shoot = .true.
     cons = this%multi_ion%f1(1,i)




  endif

  ! subtract from 0.
  this%multi_ion%f1(1,i) = this%multi_ion%f1(1,i) - cons

  do m=2, this%multi_max

     ! remember how much charge will be added to m-1. level
     inj = cons

     ! overshoot in previous level (we go to time centered densities)
     if( shoot ) then
      dens_temp = inj * 0.5_p_k_fld
      shoot = .false.
     else
     dens_temp = 0.0_p_k_fld
     endif

     ! ionizing m-1. level (index m) aadding to m. level
     cons = (this%multi_ion%f1(m,i) + dens_temp)* this%w%f1(m,i) * dt * &
        ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f1(m,i) * dt )

     ! overshoot in current level (temp might come from previous overshoot)
     if (cons > this%multi_ion%f1(m,i) + dens_temp) then
      shoot = .true.
      cons = this%multi_ion%f1(m,i) + dens_temp




     endif

     ! update density for m-1. level (index m)
     this%multi_ion%f1(m,i) = min(this%multi_ion%f1(m,i) - &
                     cons + inj, 1.0_p_k_fld)
                   !min only due to round off

  enddo ! all levels

  ! update the last level
  this%multi_ion%f1(this%multi_max+1,i) = &
          this%multi_ion%f1(this%multi_max+1,i) + cons

  ! calculate total charge
  this%multi_ion%f1(ion_idx,i) = 0.0_p_k_fld
  do m=2, this%multi_max+1
     this%multi_ion%f1(ion_idx,i) = &
             this%multi_ion%f1(ion_idx,i) + &
             (m-1)*(this%multi_ion%f1(m,i))
  enddo

  !again: just round off
  if (this%multi_ion%f1(ion_idx,i) > this%multi_max) then
     this%multi_ion%f1(ion_idx,i) = this%multi_max
  endif

  endif

  enddo
  !$omp end parallel do

end subroutine advance_density_1d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_1d_analytic(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  integer :: i, k, m, l, ion_idx
  real(p_k_fld) :: csum
  real(p_k_fld), allocatable, dimension(:) :: w_exp, coefs, multi_ion_next
  ion_idx = this%ion_idx

  ! allocate temporary buffers
  allocate(w_exp(this%multi_max + 1))
  allocate(coefs(this%multi_max + 1))
  allocate(multi_ion_next(this%multi_max + 1))

  !$omp parallel do private(m, l, w_exp,coefs,multi_ion_next)
  do i=1, this%multi_ion% nx_( 1 )

    if ( this%multi_ion%f1(this%neut_idx,i) > this%den_min ) then

     ! reset temp array level m
      multi_ion_next(1:this%multi_max+1) = 0.0_p_k_part
      l = 1

      do m= 1, this%multi_max + 1
        if(m > 1 .and. this%w%f1(m,i) == 0.0_p_k_part .and. this%w%f1(m-1,i) == 0.0_p_k_part ) then
            l = m + 1
            CYCLE
        endif

        ! calculate exp(-rate[level m] * dt) and store (only done once)
        w_exp(m) = exp(- this%w%f1(m,i) * dt )

        ! coef sum from 1 to m-1
        csum = 0.0

        do k= l, m-1
          ! update coef of exp(-rate[k] * dt)
          coefs(k) = coefs(k) * this%w%f1(m-1,i)/ &
          (this%w%f1(m,i)- this%w%f1(k,i))

          ! multiple coef by exponential exp(-rate[k] * dt) and sum to current level m
          multi_ion_next(m) = multi_ion_next(m) + w_exp(k) * coefs(k)
          csum = csum + coefs(k)
        enddo

        coefs(m) = this%multi_ion%f1(m,i) - csum
        multi_ion_next(m) = multi_ion_next(m) + w_exp(m) * coefs(m)
        this%multi_ion%f1(m,i) = multi_ion_next(m)

      enddo ! all levels


  ! calculate total charge
  this%multi_ion%f1(ion_idx,i) = 0.0_p_k_fld
  do m=2, this%multi_max+1
     this%multi_ion%f1(ion_idx,i) = &
             this%multi_ion%f1(ion_idx,i) + &
             (m-1)*(this%multi_ion%f1(m,i))
  enddo

  !again: just round off
  if (this%multi_ion%f1(ion_idx,i) > this%multi_max) then
     this%multi_ion%f1(ion_idx,i) = this%multi_max
  endif

  endif

  enddo
  !$omp end parallel do

  ! free temp buffers
  deallocate( w_exp )
  deallocate( coefs )
  deallocate( multi_ion_next )

end subroutine advance_density_1d_analytic
!-----------------------------------------------------------------------------------------




!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_2d(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  integer :: i,j, m, ion_idx
  real(p_k_fld) :: dens_temp, cons, inj
  logical :: shoot

  shoot = .false.

  ion_idx = this%ion_idx

  !$omp parallel do private(i,cons,shoot,m,inj,dens_temp)
  do j=1, this%multi_ion % nx_( 2 )
    do i=1, this%multi_ion% nx_( 1 )

        if ( this%multi_ion%f2(this%neut_idx,i,j) > this%den_min ) then

      ! 1st order method
      !cons = this(n)%multi_ion%f2(1,i,j) * w_ion(1) * dt

      ! 2nd order Runge-Kutta
      cons = this%multi_ion%f2(1,i,j) * this%w%f2(1,i,j) * dt * &
      ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f2(1,i,j) * dt )

      ! overshoot
      if (cons > this%multi_ion%f2(1,i,j)) then
        shoot = .true.
        cons = this%multi_ion%f2(1,i,j)





      endif

      ! subtract from 0.
      this%multi_ion%f2(1,i,j) = this%multi_ion%f2(1,i,j) - cons

      do m=2, this%multi_max

        ! remember how much charge will be added to m-1. level
        inj = cons

        ! overshoot in previous level (we go to time centered densities)
        if (shoot) then
          dens_temp = inj * 0.5_p_k_fld
          shoot = .false.
        else
          dens_temp = 0.0_p_k_fld
        endif

        ! ionizing m-1. level (index m) adding to m. level

        ! 1st order method
        !cons = (this(n)%multi_ion%f2(m,i,j) + dens_temp)* w_ion(m) * dt

        ! 2nd order Runge-Kutta
        cons = (this%multi_ion%f2(m,i,j) + dens_temp)* this%w%f2(m,i,j) * dt * &
        ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f2(m,i,j) * dt )

        ! overshoot in current level (temp might come from previous overshoot)
        if (cons > this%multi_ion%f2(m,i,j) + dens_temp) then
          shoot = .true.
          cons = this%multi_ion%f2(m,i,j) + dens_temp





        endif

        ! update density for m-1. level (index m)
        this%multi_ion%f2(m,i,j) = min(this%multi_ion%f2(m,i,j) - &
          cons + inj, 1.0_p_k_fld)
        !min only due to round off

      enddo ! all levels

      ! update the last level

      this%multi_ion%f2(this%multi_max+1,i,j) = &
      this%multi_ion%f2(this%multi_max+1,i,j) + cons

      ! calculate total charge
      this%multi_ion%f2(ion_idx,i,j) = 0.0_p_k_fld
      do m=2, this%multi_max+1
        this%multi_ion%f2(ion_idx,i,j) = &
        this%multi_ion%f2(ion_idx,i,j) + &
        (m-1)*(this%multi_ion%f2(m,i,j))
      enddo

      !again: just round off
      if (this%multi_ion%f2(ion_idx,i,j) > this%multi_max) then
        this%multi_ion%f2(ion_idx,i,j) = this%multi_max
      endif

      endif

    enddo
  enddo
  !$omp end parallel do

end subroutine advance_density_2d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_2d_analytic(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  integer :: i, j, k, m, l, ion_idx
  real(p_k_fld) :: csum
  real(p_k_fld), allocatable, dimension(:) :: w_exp, coefs, multi_ion_next
  ion_idx = this%ion_idx

  ! allocate temporary buffers
  allocate(w_exp(this%multi_max + 1))
  allocate(coefs(this%multi_max + 1))
  allocate(multi_ion_next(this%multi_max + 1))

  !$omp parallel do private(i,m,w_exp,coefs,multi_ion_next)
  do j=1, this%multi_ion % nx_( 2 )
    do i=1, this%multi_ion% nx_( 1 )

    if ( this%multi_ion%f2(this%neut_idx,i,j) > this%den_min ) then

     ! reset temp array level m
      multi_ion_next(1:this%multi_max+1) = 0.0_p_k_part
      l = 1

      do m= 1, this%multi_max + 1
        if(m > 1 .and. this%w%f2(m,i,j) == 0.0_p_k_part .and. this%w%f2(m-1,i,j) == 0.0_p_k_part ) then
            l = m + 1
            CYCLE
        endif

        ! calculate exp(-rate[level m] * dt) and store (only done once)
        w_exp(m) = exp(- this%w%f2(m,i,j) * dt )

        ! coef sum from 1 to m-1
        csum = 0.0

        do k= 1, m-1
          ! update coef of exp(-rate[k] * dt)
          coefs(k) = coefs(k) * this%w%f2(m-1,i,j)/ &
          (this%w%f2(m,i,j)- this%w%f2(k,i,j))

          ! multiple coef by exponential exp(-rate[k] * dt) and sum to current level m
          multi_ion_next(m) = multi_ion_next(m) + w_exp(k) * coefs(k)
          csum = csum + coefs(k)
        enddo

        coefs(m) = this%multi_ion%f2(m,i,j) - csum
        multi_ion_next(m) = multi_ion_next(m) + w_exp(m) * coefs(m)
        this%multi_ion%f2(m,i,j) = multi_ion_next(m)

      enddo ! all levels


  ! calculate total charge
  this%multi_ion%f2(ion_idx,i,j) = 0.0_p_k_fld
  do m=2, this%multi_max+1
     this%multi_ion%f2(ion_idx,i,j) = &
             this%multi_ion%f2(ion_idx,i,j) + &
             (m-1)*(this%multi_ion%f2(m,i,j))
  enddo

  !again: just round off
  if (this%multi_ion%f2(ion_idx,i,j) > this%multi_max) then
     this%multi_ion%f2(ion_idx,i,j) = this%multi_max
  endif

  endif

  enddo
  enddo
  !$omp end parallel do

  ! free temp buffers
  deallocate( w_exp )
  deallocate( coefs )
  deallocate( multi_ion_next )

end subroutine advance_density_2d_analytic
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_3d(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  integer :: i,j, k , m, ion_idx
  real(p_k_fld) :: dens_temp, cons, inj
  logical :: shoot

  shoot = .false.

  ion_idx = this%ion_idx

  !$omp parallel do private(j,i,cons,shoot,m,inj,dens_temp)
  do k=1, this%multi_ion% nx_(3)
    do j=1, this%multi_ion% nx_(2)
      do i=1, this%multi_ion% nx_(1)

        if ( this%multi_ion%f3(this%neut_idx,i,j,k) > this%den_min ) then

        ! ionizing 0. level adding to 1. level
        cons = this%multi_ion%f3(1,i,j,k) * this%w%f3(1,i,j,k) * dt * &
        ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f3(1,i,j,k) * dt )

        ! overshoot
        if (cons > this%multi_ion%f3(1,i,j,k)) then
          shoot = .true.
          cons = this%multi_ion%f3(1,i,j,k)




        endif

        ! subtract from 0.
        this%multi_ion%f3(1,i,j,k) = this%multi_ion%f3(1,i,j,k) - cons

        do m=2, this%multi_max

          ! remember how much charge will be added to m-1. level
          inj = cons

          ! overshoot in previous level (we go to time centered densities)
          if (shoot) then
            dens_temp = inj * 0.5_p_k_fld
            shoot = .false.
          else
            dens_temp = 0.0_p_k_fld
          endif

          ! ionizing m-1. level (index m) aadding to m. level
          cons = (this%multi_ion%f3(m,i,j,k) + dens_temp)* this%w%f3(m,i,j,k) * dt * &
          ( 1.0_p_k_fld + 0.5_p_k_fld * this%w%f3(m,i,j,k) * dt )

          ! overshoot in current level (temp might come from previous overshoot)
          if (cons > this%multi_ion%f3(m,i,j,k) + dens_temp) then
            shoot = .true.
            cons = this%multi_ion%f3(m,i,j,k) + dens_temp




          endif

          ! update density for m-1. level (index m)
          this%multi_ion%f3(m,i,j,k) = min(this%multi_ion%f3(m,i,j,k) - &
            cons + inj, 1.0_p_k_fld)
          !min only due to round off

        enddo ! all levels

        ! update the last level
        this%multi_ion%f3(this%multi_max+1,i,j,k) = &
        this%multi_ion%f3(this%multi_max+1,i,j,k) + cons

        ! calculate total charge
        this%multi_ion%f3(ion_idx,i,j,k) = 0.0_p_k_fld
        do m=2, this%multi_max+1
          this%multi_ion%f3(ion_idx,i,j,k) = &
          this%multi_ion%f3(ion_idx,i,j,k) + &
          (m-1)*(this%multi_ion%f3(m,i,j,k))
        enddo

        !again: just round off
        if (this%multi_ion%f3(ion_idx,i,j,k) > this%multi_max) then
          this%multi_ion%f3(ion_idx,i,j,k) = this%multi_max
        endif

        endif

      enddo

    enddo
  enddo
  !$omp end parallel do

end subroutine advance_density_3d
!-----------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine advance_density_3d_analytic(this , dt)
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this ! neutral array
  real( p_k_fld ), intent(in) :: dt

  integer :: i, j, k, l, m, n, ion_idx
  real(p_k_fld) :: csum
  real(p_k_fld), allocatable, dimension(:) :: w_exp, coefs, multi_ion_next
  ion_idx = this%ion_idx

  ! allocate temporary buffers
  allocate(w_exp(this%multi_max + 1))
  allocate(coefs(this%multi_max + 1))
  allocate(multi_ion_next(this%multi_max + 1))

  !$omp parallel do private(j,i,m,l,w_exp,coefs,multi_ion_next)
  do k=1, this%multi_ion% nx_(3)
      do j=1, this%multi_ion % nx_( 2 )
        do i=1, this%multi_ion% nx_( 1 )

        if ( this%multi_ion%f3(this%neut_idx,i,j,k) > this%den_min ) then

         ! reset temp array level m
          multi_ion_next(1:this%multi_max+1) = 0.0_p_k_part
          n = 1

          do m = 1, this%multi_max + 1
            if(m > 1 .and. this%w%f3(m,i,j,k) == 0.0_p_k_part .and. this%w%f3(m-1,i,j,k) == 0.0_p_k_part ) then
                n = m + 1
                CYCLE
            endif

            ! calculate exp(-rate[level m] * dt) and store (only done once)
            w_exp(m) = exp(- this%w%f3(m,i,j,k) * dt )

            ! coef sum from 1 to m-1
            csum = 0.0

            do l= n, m-1
              ! update coef of exp(-rate[k] * dt)
              coefs(l) = coefs(l) * this%w%f3(m-1,i,j,k)/ &
              (this%w%f3(m,i,j,k)- this%w%f3(l,i,j,k))

              ! multiple coef by exponential exp(-rate[k] * dt) and sum to current level m
              multi_ion_next(m) = multi_ion_next(m) + w_exp(l) * coefs(l)
              csum = csum + coefs(l)
            enddo

            coefs(m) = this%multi_ion%f3(m,i,j,k) - csum
            multi_ion_next(m) = multi_ion_next(m) + w_exp(m) * coefs(m)
            this%multi_ion%f3(m,i,j,k) = multi_ion_next(m)

          enddo ! all levels


      ! calculate total charge
      this%multi_ion%f3(ion_idx,i,j,k) = 0.0_p_k_fld
      do m=2, this%multi_max+1
         this%multi_ion%f3(ion_idx,i,j,k) = &
                 this%multi_ion%f3(ion_idx,i,j,k) + &
                 (m-1)*(this%multi_ion%f3(m,i,j,k))
      enddo

      !again: just round off
      if (this%multi_ion%f3(ion_idx,i,j,k) > this%multi_max) then
         this%multi_ion%f3(ion_idx,i,j,k) = this%multi_max
      endif

      endif
      enddo
    enddo
  enddo
  !$omp end parallel do

  ! free temp buffers
  deallocate( w_exp )
  deallocate( coefs )
  deallocate( multi_ion_next )

end subroutine advance_density_3d_analytic
!-----------------------------------------------------------------------------------------



!-----------------------------------------------------------------------------------------
! Inject particles from ionization
!-----------------------------------------------------------------------------------------
subroutine inject_particles( this, coordinates )
!-----------------------------------------------------------------------------------------

  class(t_neutral), intent(inout) :: this

  integer, intent(in) :: coordinates

  select case ( p_x_dim )
  case (1)
    call inject_particles_1d( this )
  case (2)
    call inject_particles_2d( this, coordinates)
  case (3)
    call inject_particles_3d( this )
  case default
    write(err_buf__,*) 'p_x_dim has the value:', p_x_dim;call err__("os-neutral.f03",1870)
    call abort_program(p_err_invalid)
  end select

end subroutine inject_particles
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine inject_particles_1d( this )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 1

  class(t_neutral), intent(inout) :: this

  integer :: n_p_c, n_p_add, i, ix, ion_idx

  real(p_k_part), dimension(rank) :: xnewpart
  integer, dimension(rank) :: ixnewpart

  real(p_k_part), dimension(p_p_dim) :: pnewpart
  real(p_k_part) :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part) :: den_center

  pnewpart = 0.0_p_k_part

  ion_idx = this%ion_idx
  n_p_c = this%species1%num_par_x(1)
  norm1 = -1.0_p_k_part / n_p_c

  ! Note that the number of particles per cell of the ion species is
  ! being ignored (which is correct )
  if (this%if_mov_ions) norm2 = +1.0_p_k_part / n_p_c

  do i=1, this%multi_ion % nx_(1)

    if ( this%multi_ion%f1(ion_idx,i) > this%multi_min ) then

      ! determine neutral density
      den_center = this%multi_ion%f1(this%neut_idx,i)

      if (den_center > this%den_min) then

        ! determine number of particles per cell
        n_p_add = int(this%multi_ion%f1(ion_idx,i)*n_p_c/this%multi_max+0.5_p_k_fld)-&
        int(this%ion_level_old%f1(1,i)*n_p_c/this%multi_max+0.5_p_k_fld)

        ! determine charge of the new particle
        qnewpart1 = real( this%multi_max*den_center*norm1, p_k_part )

        ! charge of the injected ions
        if (this%if_mov_ions) qnewpart2 = real( den_center*norm2, p_k_part )

        ! place particles in cell
        do ix = 0, n_p_add-1

          ixnewpart(1) = i

          if ( this%inject_line ) then
            xnewpart(1) = (ix + 0.5_p_k_part)/n_p_add - 0.5_p_k_part
          else
            ! using genrand_real3() makes sure that the particle is
            ! never injected in the cell boundary
            call rng % harvest_real3( xnewpart(1), ixnewpart)
            xnewpart(1) = xnewpart(1) - 0.5_p_k_part

          endif

          ! add particles to the corresponding buffers
          call this%species1%create_particle( ixnewpart, xnewpart, pnewpart, qnewpart1 )

          if (this%if_mov_ions) then

            call this%species2%create_particle( ixnewpart, xnewpart, pnewpart, qnewpart2 )
          endif

        enddo

      endif
    endif

  enddo

end subroutine inject_particles_1d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine inject_particles_2d( this, coordinates)
!-----------------------------------------------------------------------------------------
  implicit none

  integer, parameter :: rank = 2

  class(t_neutral), intent(inout) :: this

  integer, intent(in) :: coordinates

  integer, dimension(rank) :: num_par !=num_par_x
  integer :: n_p_c, n_p_add, i, j, ix, ion_idx

  real(p_k_part), dimension(rank) :: xnewpart
  integer, dimension(rank) :: ixnewpart

  real(p_k_part), dimension(p_p_dim) :: pnewpart
  real(p_k_part) :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part) :: den_center
  real(p_k_part) :: r, dr
  integer :: shift_ix2

  pnewpart = 0.0_p_k_part

  ion_idx = this%ion_idx
  num_par(1) = this%species1%num_par_x(1)
  num_par(2) = this%species1%num_par_x(2)
  n_p_c = num_par(1) * num_par(2)
  norm1 = -1.0_p_k_part / n_p_c

  if ( coordinates == p_cylindrical_b ) then
    dr = real( this%species1%dx(2), p_k_part )
    shift_ix2 = this%species1%my_nx_p(p_lower, 2) - 2
  else
    dr = 0
    shift_ix2 = 0
  endif

  ! Note that the number of particles per cell of the ion species is
  ! being ignored (which is correct )
  if (this%if_mov_ions) then
    norm2 = +1.0_p_k_part / n_p_c
  else
    norm2 = 0
  endif

  do j=1, this%multi_ion% nx_(2)

    if( coordinates == p_cylindrical_b ) then

      do i=1, this%multi_ion% nx_(1)

        if ( this%multi_ion%f2(ion_idx,i,j) > this%multi_min ) then

          ! determine neutral density
          den_center = this%multi_ion%f2(this%neut_idx,i,j)

          if (den_center > this%den_min) then

            n_p_add = int(this%multi_ion%f2(ion_idx,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
            int(this%ion_level_old%f2(1,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld)

            ! determine charge of the new particle
            qnewpart1 = real( this%multi_max*den_center*norm1, p_k_part )

            if (this%if_mov_ions) qnewpart2 = real( den_center*norm2 , p_k_part )

            ! place particles in cell
            do ix = 0, n_p_add-1

              ixnewpart(1) = i
              ixnewpart(2) = j

              if ( this%inject_line ) then
                xnewpart(1) = (ix + 0.5)/n_p_add - 0.5_p_k_part
                xnewpart(2) = 0
              else
                ! using genrand_real3() makes sure that the particle is
                ! never injected in the cell boundary
                call rng % harvest_real3( xnewpart(1), ixnewpart)
                call rng % harvest_real3( xnewpart(2), ixnewpart)

                xnewpart(1) = xnewpart(1) - 0.5_p_k_part
                xnewpart(2) = xnewpart(2) - 0.5_p_k_part
              endif

              ! get radial position normalized to cell size
              ! since g_box is adjusted based on interpolation order, this will always
              ! give the correct particle position

              ! get radial position
              r = this%species1%g_box( p_lower , p_r_dim ) + &
              ( ( ixnewpart(2) + shift_ix2 ) + xnewpart(2) ) * dr

              ! Only inject inside the box. This could be optimized for near cell positions
              ! by only injecting starting from cell 2
              if ( r > 0 ) then

                ! add particles to the corresponding buffers
                call this%species1%create_particle( ixnewpart, xnewpart, pnewpart, &
                  qnewpart1*r )

                if (this%if_mov_ions) then
                  call this%species2%create_particle( ixnewpart, xnewpart, pnewpart, &
                    qnewpart2 *r)
                endif
              endif

            enddo ! n_p_add
          endif ! den_center > den_min
        endif

      enddo ! nx1

      ! cartesian
    else

      do i=1, this%multi_ion% nx_(1)

        if ( this%multi_ion%f2(ion_idx,i,j) > this%multi_min ) then

          ! determine neutral density
          den_center = this%multi_ion%f2(this%neut_idx,i,j)

          if (den_center > this%den_min) then

            n_p_add = int(this%multi_ion%f2(ion_idx,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
            int(this%ion_level_old%f2(1,i,j)*n_p_c/this%multi_max + 0.5_p_k_fld)

            ! determine charge of the new particle
            qnewpart1 = real( this%multi_max*den_center*norm1, p_k_part )

            if (this%if_mov_ions) qnewpart2 = real( den_center*norm2 , p_k_part )

            ! place particles in cell
            do ix = 0, n_p_add-1

              ixnewpart(1) = i
              ixnewpart(2) = j

              if ( this%inject_line ) then
                xnewpart(1) = (ix + 0.5)/n_p_add - 0.5_p_k_part
                xnewpart(2) = 0
              else
                ! using genrand_real3() makes sure that the particle is
                ! never injected in the cell boundary
                call rng % harvest_real3( xnewpart(1), ixnewpart)
                call rng % harvest_real3( xnewpart(2), ixnewpart)

                xnewpart(1) = xnewpart(1) - 0.5_p_k_part
                xnewpart(2) = xnewpart(2) - 0.5_p_k_part

              endif

              ! add particles to the corresponding buffers
              call this%species1%create_particle( ixnewpart, xnewpart, pnewpart, &
                qnewpart1 )

              if (this%if_mov_ions) then
                call this%species2%create_particle( ixnewpart, xnewpart, pnewpart, &
                  qnewpart2 )
              endif
            enddo ! n_p_add
          endif ! den_center > den_min
        endif ! multi_ion > multi_min
      enddo ! nx1

    endif ! cylindrical
  enddo ! nx2

end subroutine inject_particles_2d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine inject_particles_3d( this )
!-----------------------------------------------------------------------------------------

  implicit none

  integer, parameter :: rank = 3

  class(t_neutral), intent(inout) :: this

  integer, dimension(rank) :: num_par !=num_par_x
  integer :: n_p_c, n_p_add, i, j, k, ix, ion_idx
  real(p_k_part), dimension(rank) :: xnewpart
  integer, dimension(rank) :: ixnewpart
  real(p_k_part), dimension(p_p_dim) :: pnewpart
  real(p_k_part) :: norm1, norm2, qnewpart1, qnewpart2
  real(p_k_part) :: den_center

  pnewpart = 0.0_p_k_part

  ! number of neutrals

  ion_idx = this%ion_idx

  num_par(1) = this%species1%num_par_x(1)
  num_par(2) = this%species1%num_par_x(2)
  num_par(3) = this%species1%num_par_x(3)

  n_p_c = num_par(1) * num_par(2) * num_par(3)
  norm1 = -1.0_p_k_part / n_p_c

  ! Note that the number of particles per cell of the ion species is
  ! being ignored (which is correct )
  if (this%if_mov_ions) norm2 = +1.0_p_k_part / n_p_c

  do k=1, this%multi_ion% nx_(3)
    do j=1, this%multi_ion% nx_(2)
      do i=1, this%multi_ion% nx_(1)

        if ( this%multi_ion%f3(ion_idx,i,j,k) > this%multi_min ) then

          ! determine neutral density
          den_center = this%multi_ion%f3(this%neut_idx,i,j,k)

          if (den_center > this%den_min) then

            n_p_add = int(this%multi_ion%f3(ion_idx,i,j,k)*n_p_c/this%multi_max + 0.5_p_k_fld) - &
            int(this%ion_level_old%f3(1,i,j,k)*n_p_c/this%multi_max + 0.5_p_k_fld)

            ! determine charge of the new particle
            qnewpart1 = real( this%multi_max*den_center*norm1 , p_k_part )

            if (this%if_mov_ions) qnewpart2 = real( den_center*norm2, p_k_part )

            ! place particles in cell
            do ix = 0, n_p_add-1

              ixnewpart(1) = i
              ixnewpart(2) = j
              ixnewpart(3) = k

              if ( this%inject_line ) then
                xnewpart(1) = (ix + 0.5)/n_p_add - 0.5_p_k_part
                xnewpart(2) = 0
                xnewpart(3) = 0
              else
                ! using genrand_real3() makes sure that the particle is
                ! never injected in the cell boundary
                call rng % harvest_real3( xnewpart(1), ixnewpart)
                call rng % harvest_real3( xnewpart(2), ixnewpart)
                call rng % harvest_real3( xnewpart(3), ixnewpart)
                xnewpart(1) = xnewpart(1) - 0.5_p_k_part

                xnewpart(2) = xnewpart(2) - 0.5_p_k_part

                xnewpart(3) = xnewpart(3) - 0.5_p_k_part

              endif

              ! add particles to the corresponding buffers
              call this%species1%create_particle( ixnewpart, xnewpart, pnewpart, qnewpart1 )

              if (this%if_mov_ions) then

                call this%species2%create_particle( ixnewpart, xnewpart, pnewpart, qnewpart2 )
              endif

            enddo

          endif
        endif
      enddo

    enddo
  enddo

end subroutine inject_particles_3d
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! report on ionization level diagnostic

!-----------------------------------------------------------------------------------------
subroutine report_neutral( this, g_space, grid, no_co, tstep, t )
!-----------------------------------------------------------------------------------------

  use m_time_step

  implicit none

  class( t_neutral ), intent(inout) :: this
  type( t_space ), intent(in) :: g_space
  class( t_grid ), intent(in) :: grid
  class( t_node_conf ), intent(in) :: no_co
  type( t_time_step ), intent(in) :: tstep
  real(p_double), intent(in) :: t

  call report( this%diag, this%multi_ion, this%ion_idx, this%neut_idx, g_space, &
               grid, no_co, tstep, t )

end subroutine report_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! This sets the neutral background densities on the grid
!-----------------------------------------------------------------------------------------
subroutine set_neutden_values( this, den_neutral, neut_idx, den_min, &
                               nx_p_min, g_space, lbin, ubin )

  implicit none

  ! multi_ion vdf (1st comp :: neutral density)
  type(t_vdf), intent(inout) :: this

  class(t_psource_std), intent(inout) :: den_neutral
  integer, intent(in) :: neut_idx
  real(p_k_fld), intent(in) :: den_min

  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space

  integer, dimension(:), intent(in) :: lbin
  integer, dimension(:), intent(in) :: ubin

  real(p_k_part), dimension(p_x_dim) :: xden, gxmin, ldx
  real(p_k_fld) :: den_temp
  integer :: i,j,k

  gxmin = real( xmin(g_space), p_k_part )
  ldx = real( this%dx(), p_k_part )
  select case (p_x_dim)
  case(1)

    do i = lbin(1), ubin(1)

      xden(1) = gxmin(1) + ldx(1)*( (i + nx_p_min(1) - 2) + 0.5)

      den_temp = real( den_neutral % den_value( xden ), p_k_fld )

      if (den_temp > den_min) then
        this%f1(1,i) = 1.0_p_k_fld
        this%f1(neut_idx,i) = den_temp
      endif

    enddo

  case(2)

    do j = lbin(2), ubin(2)
      xden(2) = gxmin(2) + ldx(2)*( (j + nx_p_min(2) - 2) + 0.5)
      do i = lbin(1), ubin(1)

        xden(1) = gxmin(1) + ldx(1)*( (i + nx_p_min(1) - 2) + 0.5)

        den_temp = real( den_neutral % den_value( xden ), p_k_fld )

        if (den_temp > den_min) then
          this%f2(1,i,j) = 1.0_p_k_fld
          this%f2(neut_idx,i,j) = den_temp

        endif

      enddo
    enddo

  case(3)

    do k = lbin(3), ubin(3)

      xden(3) = gxmin(3) + ldx(3)*( (k + nx_p_min(3) - 2) + 0.5)
      do j = lbin(2), ubin(2)

        xden(2) = gxmin(2) + ldx(2)*( (j + nx_p_min(2) - 2) + 0.5)
        do i = lbin(1), ubin(1)

          xden(1) = gxmin(1) + ldx(1)*( (i + nx_p_min(1) - 2) + 0.5)

          den_temp = real( den_neutral % den_value( xden ), p_k_fld )

          if ( den_temp > den_min ) then
            this%f3(1,i,j,k) = 1.0_p_k_fld
            this%f3(neut_idx,i,j,k) = den_temp
          endif

        enddo
      enddo
    enddo

  end select

end subroutine set_neutden_values
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! write object information into a restart file
!-----------------------------------------------------------------------------------------
subroutine restart_write_neutral( this, restart_handle )
!-----------------------------------------------------------------------------------------

  use m_restart

  implicit none

  class( t_neutral ), intent(in) :: this
  type( t_restart_handle ), intent(inout) :: restart_handle

  character(len=*), parameter :: err_msg = 'error writing restart data for neutral object.'
  integer :: ierr

  call restart_io_write("p_neutral_rst_id", p_neutral_rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2365)

  ! diagnostic this%diag is not saved, so that the user
  ! can change the dump factors before restarting
  call restart_io_write("this%neutral_id", this%neutral_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2370)

  call restart_io_write("this%name", this%name, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2373)

  call restart_io_write("this%if_mov_ions", this%if_mov_ions, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2376)

  call restart_io_write("this%omega_p", this%omega_p, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2379)

  call restart_io_write("this%den_min", this%den_min, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2382)

  call restart_io_write("this%if_tunnel", this%if_tunnel, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2385)

  call restart_io_write("this%if_impact", this%if_impact, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2388)

  call restart_io_write("this%multi_max", this%multi_max, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstwrt,"os-neutral.f03",2391)

  call this % multi_ion % write_checkpoint( restart_handle )
  call this % ion_level_old % write_checkpoint( restart_handle )

  call this % w % write_checkpoint( restart_handle )

  if ( this%if_impact ) then
    call restart_write( this%cross_section, restart_handle)
  endif

end subroutine restart_write_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Read checkpoint information from file
!-----------------------------------------------------------------------------------------
subroutine restart_read_neutral( this, restart_handle )
!-----------------------------------------------------------------------------------------

  use m_restart

  implicit none

  class( t_neutral ), intent(inout) :: this
  type( t_restart_handle ), intent(in) :: restart_handle

  character(len=*), parameter :: err_msg = 'error reading restart data for neutral object.'
  character(len=len(p_neutral_rst_id)) :: rst_id
  integer :: ierr

  call restart_io_read(rst_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2423)

  ! check if restart file is compatible
  if ( rst_id /= p_neutral_rst_id) then
    write(err_buf__,*) 'Corrupted restart file, or restart file, from incompatible binary (neutral)';call err__("os-neutral.f03",2427)
    call abort_program(p_err_rstrd)
  endif

  ! diagnostic this%diag is not read, so that the user
  ! can change the dump factors before restarting

  call restart_io_read(this%neutral_id, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2435)

  call restart_io_read(this%name, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2438)

  call restart_io_read(this%if_mov_ions, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2441)

  call restart_io_read(this%omega_p, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2444)

  call restart_io_read(this%den_min, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2447)

  call restart_io_read(this%if_tunnel, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2450)

  call restart_io_read(this%if_impact, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2453)

  call restart_io_read(this%multi_max, restart_handle, ierr)
  call check_error(ierr,err_msg,p_err_rstrd,"os-neutral.f03",2456)

  call this % multi_ion % read_checkpoint( restart_handle )
  call this % ion_level_old% read_checkpoint( restart_handle )

  call this % w % read_checkpoint( restart_handle )

  if ( this%if_impact ) then
    call restart_read( this%cross_section, restart_handle )
  endif

end subroutine restart_read_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Update boundaries of neutral object (communication only)
!-----------------------------------------------------------------------------------------
subroutine update_boundary_neutral( this, nx_move, no_co, send_msg, recv_msg )

  use m_vdf_comm

  implicit none

  class( t_neutral ), intent( inout ) :: this

  integer, dimension(:), intent(in) :: nx_move
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  ! update boundaries ion_level vdf
  call update_boundary( this%multi_ion, p_vdf_replace, no_co, send_msg, recv_msg, nx_move )

end subroutine update_boundary_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Move simulation window for the neutral object
!-----------------------------------------------------------------------------------------
subroutine move_window_neutral( this, nx_p_min, g_space , need_den_val )

  implicit none

  ! dummy variables

  class( t_neutral ), intent(inout) :: this
  integer, dimension(:), intent(in) :: nx_p_min
  type( t_space ), intent(in) :: g_space
  logical , intent(in) :: need_den_val

  ! local variables
  integer, dimension( 2, p_x_dim ) :: init_bnd

  integer :: i

  ! executable statements

  ! move ion_level vdf
  call move_window( this%multi_ion, g_space )
  call move_window( this%ion_level_old, g_space )

  ! initalize ionization values where required
  ! note that we do not use values from another node, we simply
  ! recalculate them (which ends up being faster that having the
  ! extra communication)
  do i = 1, p_x_dim

    if ( nx_move( g_space, i ) > 0) then

      if(need_den_val) then

        init_bnd(p_lower,1:p_x_dim) = 1 - this%multi_ion%gc_num_( p_lower, 1:p_x_dim )

        init_bnd(p_upper,1:p_x_dim) = this%multi_ion%nx_(1:p_x_dim ) + &
        this%multi_ion%gc_num_( p_upper, 1:p_x_dim )

        init_bnd(p_lower,i) = this%multi_ion%nx_(i) - nx_move( g_space, i )+ this%multi_ion%gc_num_(p_lower,i)

        call set_neutden_values( this%multi_ion, this%den_neutral, &
          this%neut_idx, this%den_min, nx_p_min, &
          g_space, init_bnd(p_lower,:), init_bnd(p_upper,:) )
      endif
    endif
  enddo

end subroutine move_window_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Reshape neutral object for dynamic load balance
!-----------------------------------------------------------------------------------------
subroutine reshape_neutral( this, old_lb, new_lb, no_co, send_msg, recv_msg )
!-----------------------------------------------------------------------------------------
  use m_vdf_comm

  implicit none

  ! dummy vars
  class(t_neutral), intent(inout) :: this
  class(t_grid), intent(in) :: old_lb, new_lb
  class( t_node_conf ), intent(in) :: no_co
  type(t_vdf_msg), dimension(2), intent(inout) :: send_msg, recv_msg

  type(t_vdf_report), pointer :: report

  call reshape_copy( this%ion_level_old, old_lb, new_lb, no_co, send_msg, recv_msg )
  call reshape_copy( this%multi_ion, old_lb, new_lb, no_co, send_msg, recv_msg )
  call reshape_nocopy( this%w, new_lb )

  ! Reshape tavg_data
  report => this%diag%reports
  do
    if ( .not. associated(report) ) exit
    if ( report%tavg_data%x_dim_ > 0 ) then
      call reshape_copy( report%tavg_data, old_lb, new_lb, no_co, send_msg, recv_msg )
    endif
    report => report%next
  enddo

end subroutine reshape_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Add ionization calculation load
!-----------------------------------------------------------------------------------------
subroutine add_particle_load_neutral( this, grid )
!-----------------------------------------------------------------------------------------

    implicit none

    ! dummy vars
    class(t_neutral), intent(in) :: this
    class(t_grid), intent(inout) :: grid

    integer, dimension( p_max_dim ) :: dim_off
    real( p_double ), dimension(:), pointer :: int_load
    integer :: i, j, k, i1, i2, i3

    ! base_load is the basic load for calculating adk rates, doing injection, etc.
    ! multi_fac is the additional load for each ionization level (mostly in advance_density)
    ! note that the real load will depend on the nearby field strength, but that is ignored
    ! these parameters were calculated based on a few sample simulations by K. Miller
    real(p_double), parameter :: base_load = 0.35, multi_fac = 0.15
    real(p_double) :: load

    ! Get the offset of the int_load array for each direction
    dim_off(1) = 0
    do i = 2, grid%x_dim
        dim_off(i) = dim_off(i-1) + grid%g_nx(i-1)
    enddo

    ! Add the correction for global cell positions to dim_off
    do i = 1, grid%x_dim
        dim_off(i) = dim_off(i) + ( this%species1%my_nx_p( p_lower, i ) - 1 )
    enddo

    ! just for clarity
    int_load => grid%int_load

    ! Approximate load for each cell
    load = base_load + this%multi_max * multi_fac

    select case( grid%x_dim )
    case(1)

        do i = 1, this%multi_ion%nx_(1)
            if ( this%multi_ion%f1(this%neut_idx,i) > this%den_min ) then
                i1 = dim_off(1) + i
                int_load(i1) = int_load(i1) + load
            endif
        enddo

    case(2)

        do j = 1, this%multi_ion%nx_(2)
            do i = 1, this%multi_ion%nx_(1)
                if ( this%multi_ion%f2(this%neut_idx,i,j) > this%den_min ) then
                    i1 = dim_off(1) + i
                    i2 = dim_off(2) + j
                    int_load(i1) = int_load(i1) + load
                    int_load(i2) = int_load(i2) + load
                endif
            enddo
        enddo

    case(3)

        do k = 1, this%multi_ion%nx_(3)
            do j = 1, this%multi_ion%nx_(2)
                do i = 1, this%multi_ion%nx_(1)
                    if ( this%multi_ion%f3(this%neut_idx,i,j,k) > this%den_min ) then
                        i1 = dim_off(1) + i
                        i2 = dim_off(2) + j
                        i3 = dim_off(3) + k
                        int_load(i1) = int_load(i1) + load
                        int_load(i2) = int_load(i2) + load
                        int_load(i3) = int_load(i3) + load
                    endif
                enddo
            enddo
        enddo

    end select

end subroutine add_particle_load_neutral
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
! Printout the algorithm used by the pusher
!-----------------------------------------------------------------------------------------
subroutine list_algorithm_neutral( this )

  implicit none
  class( t_neutral ), intent(in) :: this

  print *, ' '
  print *, trim(this%name),' :'
  print *, "  - Ionization model: ADK"

  print *, "   Maximum ionization level used: ", this%multi_max

end subroutine list_algorithm_neutral
!-----------------------------------------------------------------------------------------

! ****************************************************************************************
! Impact Ionization code
! - This is currently offline
! ****************************************************************************************
# 3122 "os-neutral.f03"
!---------------------------------------------------------------------------------------------------
! Generate memory allocation / deallocation routines





# 1 "./memory/mem-template.h" 1
! Generated automatically by memory.py on 2019-03-06 18:36:58
# 24 "./memory/mem-template.h"
subroutine alloc_neutral( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_neutral ), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize







  ! Nullify pointer
  p => null()

  ! Allocate memory



  allocate( p, stat = ierr )


  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         "t_neutral" // ', ' // &
           trim(strfileline(file,line))
    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif






  ! log allocation
  addr = loc( p )
  bsize = sizeof( p )
  call log_allocation( addr, "t_neutral" , &
                        bsize, &
                        file, line )

end subroutine alloc_neutral



subroutine alloc_1d_neutral( p, n, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_neutral ), dimension(:), pointer :: p
  integer, dimension(:), intent(in) :: n

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: i, ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_dims = -1001
# 111 "./memory/mem-template.h"
  ! Nullify pointer
  p => null()

  ! Verify requested size
  do i = 1, 1
    if ( n(i) <= 0 ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid dimensions requested, ' //&
             "t_neutral" // '(' // &
             trim(strrange(1,n)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_dims
        return
      else
        call exit(p_err_invalid_dims)
      endif

    endif
  enddo

  ! Allocate memory
# 143 "./memory/mem-template.h"
  allocate( p( n(1) ), stat = ierr )


  ! Check allocation
  if ( ierr /= 0 ) then
  write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
         "t_neutral" // '(' // &
         trim(strrange(1,n)) // '), ' // &
           trim(strfileline(file,line))
  if ( present(stat) ) then
    stat = ierr
    return
  else
    call exit(ierr)
  endif
  endif






  ! log allocation
  addr = loc(p)
  bsize = sizeof( p )
  call log_allocation( addr, "t_neutral" // '('// trim(strrange(1,n)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_1d_neutral
# 270 "./memory/mem-template.h"
subroutine alloc_bound_1d_neutral( p, lb, ub, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_neutral ), dimension(:), pointer :: p
  integer, dimension(:), intent(in) :: lb, ub

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: i,ierr
  integer( p_int64 ) :: addr, bsize

  ! Error constants
  integer, parameter :: p_err_invalid_bounds = -1002
# 298 "./memory/mem-template.h"
  ! Nullify pointer
  p => null()

  ! Verify requested boundaries
  do i = 1, 1
    if ( ub(i) < lb(i) ) then
      write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
             'invalid boundaries requested, ' // &
             "t_neutral" // '(' // &
             trim(strrange(1,lb,ub)) // '), ' // &
             trim(strfileline(file,line))

      if ( present(stat) ) then
        stat = p_err_invalid_bounds
        return
      else
        call exit(p_err_invalid_bounds)
      endif

    endif
  enddo

  ! Allocate memory
# 329 "./memory/mem-template.h"
  allocate( p( lb(1):ub(1) ), stat = ierr )


  ! Check allocation
  if ( ierr /= 0 ) then
    write( p_stderr, '(A)' ) '(*error*) Unable to allocate memory, ' // &
           "t_neutral" // '(' // &
           trim(strrange(1,lb,ub)) // '), ' // &
           trim(strfileline(file,line))

    if ( present(stat) ) then
      stat = ierr
      return
    else
      call exit(ierr)
    endif
  endif







  ! log allocation
  addr = loc(p)
  bsize = sizeof(p)
  call log_allocation( addr, "t_neutral" // '('// trim(strrange(1,lb,ub)) // ')', &
                        bsize, &
                        file, line )

end subroutine alloc_bound_1d_neutral



subroutine freemem_neutral( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_neutral ), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
   call log_deallocation( addr, "t_neutral", &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_neutral" // ', ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_neutral



subroutine freemem_1d_neutral( p, file, line, stat )

  use iso_c_binding

  implicit none

  type( t_neutral ), dimension(:), pointer :: p

  character(len=*), intent(in) :: file
  integer, intent(in) :: line
  integer, intent(out), optional :: stat

  integer :: ierr
  integer( p_int64 ) :: addr, bsize

  ! check if p is associated, otherwise return silently
  if ( associated(p) ) then
   ! get memory size to deallocate
   bsize = sizeof(p)
   addr = loc(p)

   ! Decrease global memory counter
     call log_deallocation( addr, "t_neutral" // '(:)', &
                        bsize, file, line )

   ! deallocate memory





   deallocate( p, stat = ierr )


   ! Check allocation
   if ( ierr /= 0 ) then
     write( p_stderr, * ) '(*error*) Unable to deallocate memory, ' // &
            "t_neutral" // '(' // &
        trim(strrange(1,lbound(p),ubound(p))) // '), ' // &
        trim(strfileline(file,line))
     if ( present(stat) ) then
     stat = ierr
     return
     else
     call exit(ierr)
     endif
   endif

   ! Nullify pointer
   p => null()

  endif

end subroutine freemem_1d_neutral
# 3130 "os-neutral.f03" 2

!---------------------------------------------------------------------------------------------------

end module m_neutral
