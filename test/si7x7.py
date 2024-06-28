from ase import Atoms
from ase.dft.kpoints import monkhorst_pack
from ase.constraints import FixAtoms
from fireballpy import Fireball
import time

atoms = Atoms(symbols=['Si','Si','H','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'H','Si','Si','Si','H','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','H','Si',
                       'Si','Si','H','Si','Si','H',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','H','H','Si','Si',
                       'Si','Si','H','Si','Si','Si',
                       'H','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','H','Si','H',
                       'Si','H','H','Si','Si','Si',
                       'Si','Si','Si','Si','H','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','H','H','Si',
                       'Si','Si','Si','Si','H','H',
                       'H','H','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','H','H','H','H','H',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','H','Si','Si',
                       'Si','Si','H','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','H','H',
                       'Si','Si','Si','H','H','H',
                       'H','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','H','Si','H',
                       'H','Si','Si','Si','Si','Si',
                       'H','H','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','H','H','Si','Si',
                       'Si','H','H','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','H',
                       'Si','H','Si','Si','H','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','Si','Si',
                       'Si','Si','Si','Si','H','Si',
                       'Si','H','Si','Si','Si','Si',
                       'Si','Si','H','Si','Si','Si'],
        positions=[(    -11.647951,     19.969287,      2.445701),
                  (    -11.520889,     19.954767,      0.000000),
                  (    -11.520880,     22.171962,     -2.185591),
                  (    -11.520880,     22.171962,     -0.783900),
                  (    -11.487294,     17.707865,      3.224180),
                  (    -11.208755,     17.108104,      5.555609),
                  (    -11.169513,     24.855104,      5.541515),
                  (     -9.997961,     14.994170,      5.600854),
                  (     -9.941084,     26.966564,      5.590516),
                  (     -9.611841,     14.393798,      3.281327),
                  (     -9.448434,     18.713650,      6.311693),
                  (     -9.600739,     18.846167,     -2.185591),
                  (     -9.600739,     18.846167,     -0.783900),
                  (     -9.600739,     23.280563,      0.000000),
                  (     -9.600737,     16.628969,      0.000000),
                  (     -9.600735,     25.497748,     -2.185591),
                  (     -9.600735,     25.497748,     -0.783900),
                  (     -9.664839,     16.538436,      2.320885),
                  (     -9.612825,     20.971526,      3.209639),
                  (     -9.637164,     23.186626,      2.349647),
                  (     -9.412095,     23.320723,      6.369627),
                  (     -8.960226,     21.037849,      5.543046),
                  (     -7.776423,     15.642481,      6.335392),
                  (     -7.685365,     19.963186,      2.318801),
                  (     -7.671195,     26.583300,      2.351363),
                  (     -7.682878,     17.674347,      2.833572),
                  (     -7.703153,     13.215638,      2.445256),
                  (     -7.680595,     15.520372,     -2.185591),
                  (     -7.680595,     15.520372,     -0.783900),
                  (     -7.680594,     13.303176,      0.000000),
                  (     -7.680594,     26.606352,      0.000000),
                  (     -7.680588,     22.171966,     -2.185591),
                  (     -7.680588,     22.171966,     -0.783900),
                  (     -7.680588,     19.954767,      0.000000),
                  (     -7.680583,     28.823546,     -2.185591),
                  (     -7.680583,     28.823546,     -0.783900),
                  (     -7.784849,     26.335422,      6.462397),
                  (     -7.673847,     24.311948,      2.886947),
                  (     -7.661257,     17.670432,      5.142476),
                  (     -7.644122,     10.976570,      3.294407),
                  (     -7.616939,     24.372108,      5.213560),
                  (     -7.318979,     10.398198,      5.628949),
                  (     -7.317891,     31.589406,      5.620383),
                  (     -6.516556,     21.054337,      5.568565),
                  (     -6.098780,      8.298049,      5.609314),
                  (     -6.080422,     33.680083,      5.613546),
                  (     -5.864944,     14.227297,      5.732363),
                  (     -5.799030,     20.985571,      3.253390),
                  (     -7.680658,     17.735130,      7.658666),
                  (     -5.813638,      9.880331,      2.318426),
                  (     -5.706119,     23.184560,      2.350696),
                  (     -5.957231,     18.812778,      6.304960),
                  (     -5.799191,     29.846879,      2.357079),
                  (     -5.694391,     16.526392,      2.382782),
                  (     -5.770556,      7.742365,      3.289187),
                  (     -5.760443,     18.846170,     -2.185591),
                  (     -5.760443,     25.497751,     -2.185591),
                  (     -5.760443,     18.846170,     -0.783900),
                  (     -5.760443,     25.497751,     -0.783900),
                  (     -5.760443,     16.628972,      0.000000),
                  (     -5.760443,     29.932150,      0.000000),
                  (     -5.760442,     12.194575,     -2.185591),
                  (     -5.760442,     12.194575,     -0.783900),
                  (     -5.760442,      9.977375,      0.000000),
                  (     -5.760442,     23.280556,      0.000000),
                  (     -5.760435,     32.149343,     -2.185591),
                  (     -5.760435,     32.149343,     -0.783900),
                  (     -5.755852,     27.647417,      3.279643),
                  (     -5.778495,     14.315329,      3.290808),
                  (     -5.917542,     27.730232,      5.734441),
                  (     -5.971504,     23.236321,      6.431085),
                  (     -7.721485,     24.276167,      7.898124),
                  (     -5.566898,     11.894747,      6.454932),
                  (     -5.603494,     30.056118,      6.435933),
                  (     -3.902759,     15.500149,      5.865428),
                  (     -3.945304,      9.097583,      6.458573),
                  (     -3.833549,     33.237194,      2.358042),
                  (     -3.873641,      6.564390,      2.431709),
                  (     -3.852568,     11.016282,      2.889278),
                  (     -3.814557,     19.921365,      2.442820),
                  (     -3.840300,      8.868781,     -2.185591),
                  (     -3.840300,      8.868781,     -0.783900),
                  (     -3.840298,     22.171959,     -2.185591),
                  (     -3.840298,     22.171959,     -0.783900),
                  (     -3.840293,     15.520371,     -2.185591),
                  (     -3.840293,     28.823550,     -2.185591),
                  (     -3.840293,     15.520371,     -0.783900),
                  (     -3.840293,     28.823550,     -0.783900),
                  (     -3.840293,      6.651589,      0.000000),
                  (     -3.840293,     19.954771,      0.000000),
                  (     -3.840293,     26.606351,      0.000000),
                  (     -3.840293,     33.257947,      0.000000),
                  (     -3.840291,     13.303172,      0.000000),
                  (     -3.840284,     35.475141,     -2.185591),
                  (     -3.840284,     35.475141,     -0.783900),
                  (     -3.916223,     32.925643,      6.450167),
                  (     -3.910788,     26.484514,      5.955012),
                  (     -3.837704,     26.538716,      2.408424),
                  (     -3.834534,     13.293680,      2.372180),
                  (     -3.821475,     24.309173,      3.264768),
                  (     -3.817248,     17.697527,      3.271529),
                  (     -3.815974,     30.968734,      2.870905),
                  (     -3.839610,     11.041588,      8.084635),
                  (     -3.805886,      4.317441,      3.260621),
                  (     -3.794620,     17.841567,      5.712035),
                  (     -3.850305,     30.950441,      7.989633),
                  (     -3.788145,     11.050788,      5.170735),
                  (     -3.821667,     24.156891,      5.717844),
                  (     -3.800704,     30.919319,      5.197587),
                  (     -3.476160,      3.695662,      5.579888),
                  (     -3.396596,     38.260562,      5.573082),
                  (     -2.246271,      1.594062,      5.523670),
                  (     -2.189887,     40.385553,      5.535936),
                  (     -2.111390,     12.012415,      6.465541),
                  (     -2.088368,     29.942485,      6.449295),
                  (     -2.053561,     34.350712,      5.720613),
                  (     -2.080665,      7.685565,      5.704270),
                  (     -1.972495,      3.204727,      2.313848),
                  (     -1.869505,      1.023724,      3.221292),
                  (     -1.954679,     36.506388,      2.351395),
                  (     -1.912688,     34.305908,      3.272177),
                  (     -1.943173,     23.174340,      2.351842),
                  (     -1.920149,     12.194577,     -2.185591),
                  (     -1.920149,     25.497755,     -2.185591),
                  (     -1.920149,     12.194577,     -0.783900),
                  (     -1.920149,     25.497755,     -0.783900),
                  (     -1.920149,     23.280557,      0.000000),
                  (     -1.920148,     36.583736,      0.000000),
                  (     -1.920147,      9.977380,      0.000000),
                  (     -1.920141,      5.542989,     -2.185591),
                  (     -1.920141,     18.846170,     -2.185591),
                  (     -1.920141,     32.149349,     -2.185591),
                  (     -1.920141,     38.800927,     -2.185591),
                  (     -1.920141,      5.542989,     -0.783900),
                  (     -1.920141,     18.846170,     -0.783900),
                  (     -1.920141,     32.149349,     -0.783900),
                  (     -1.920141,     38.800927,     -0.783900),
                  (     -1.920141,      3.325793,      0.000000),
                  (     -1.920141,     16.628969,      0.000000),
                  (     -1.920141,     29.932148,      0.000000),
                  (     -1.942283,     16.535726,      2.377934),
                  (     -1.883963,     14.329064,      5.731525),
                  (     -1.952237,      7.690628,      3.267860),
                  (     -1.914333,     27.632733,      3.274582),
                  (     -1.899197,     14.330690,      3.289902),
                  (     -1.875435,      9.899727,      2.361647),
                  (     -1.850260,     29.837074,      2.356548),
                  (     -1.848878,     20.986854,      3.292212),
                  (     -1.874438,     27.630081,      5.733343),
                  (     -1.782532,     36.678198,      6.414673),
                  (     -1.825417,      5.324047,      6.295498),
                  (     -1.667327,     18.803415,      6.461889),
                  (     -1.649273,     23.278855,      6.474871),
                  (     -1.202863,     21.045780,      5.614173),
                  (      0.018802,     39.623781,      6.383702),
                  (      0.014066,     19.941030,      2.324688),
                  (      0.000000,      2.217197,     -2.185591),
                  (      0.000000,     15.520376,     -2.185591),
                  (      0.000000,     22.171956,     -2.185591),
                  (      0.000000,     28.823554,     -2.185591),
                  (      0.000000,     35.475135,     -2.185591),
                  (      0.000000,      2.217197,     -0.783900),
                  (      0.000000,     15.520376,     -0.783900),
                  (      0.000000,     22.171956,     -0.783900),
                  (      0.000000,     28.823554,     -0.783900),
                  (      0.000000,     35.475135,     -0.783900),
                  (      0.000000,      0.000000,      0.000000),
                  (      0.000000,     13.303174,      0.000000),
                  (      0.000000,     26.606355,      0.000000),
                  (      0.000000,     39.909533,      0.000000),
                  (      0.000002,      8.868779,     -2.185591),
                  (      0.000002,      8.868779,     -0.783900),
                  (      0.000002,      6.651581,      0.000000),
                  (      0.000002,     19.954760,      0.000000),
                  (      0.000006,     33.257947,      0.000000),
                  (      0.000009,     42.126727,     -2.185591),
                  (      0.000009,     42.126727,     -0.783900),
                  (      0.014146,     13.239869,      2.396719),
                  (      0.011723,     24.315705,      2.890102),
                  (      0.011080,     -2.284746,      3.666583),
                  (      0.016024,     17.680005,      2.935077),
                  (      0.014451,     26.582032,      2.362753),
                  (      0.022253,     33.259957,      5.983881),
                  (      0.025723,     39.905759,      2.348739),
                  (      0.010312,     33.205469,      2.411856),
                  (      0.003604,      6.660704,      2.381307),
                  (      0.014192,     37.634551,      2.872691),
                  (      0.001096,      4.355310,      2.811383),
                  (      0.011321,     30.991987,      3.291043),
                  (      0.005761,     11.026828,      3.275030),
                  (      0.021002,     -0.221245,      2.446056),
                  (      0.017082,     15.676256,      6.423733),
                  (      0.013096,     17.684352,      5.227605),
                  (      0.033054,     37.551605,      5.197095),
                  (      0.004573,      8.741139,      5.870956),
                  (      0.014986,     17.704830,      7.962494),
                  (      0.009723,     11.080391,      5.723563),
                  (      0.012421,      4.436251,      5.135342),
                  (      0.007763,      2.341448,      6.270571),
                  (      0.017680,     37.635466,      7.854242),
                  (      0.022402,     30.929293,      5.746208),
                  (      0.014210,     24.356362,      5.215722),
                  (      0.016154,     24.299911,      8.083577),
                  (      0.012811,     26.304345,      6.508129),
                  (     -0.015632,      4.358036,      7.635386),
                  (      1.226515,     21.040787,      5.614905),
                  (      1.695816,     18.801853,      6.459240),
                  (      1.669461,     23.271543,      6.476849),
                  (      1.832190,     36.682828,      6.429818),
                  (      1.885952,     20.981076,      3.276854),
                  (      1.854182,      5.313657,      6.328939),
                  (      1.878465,     29.840786,      2.366546),
                  (      1.886254,      9.893233,      2.375223),
                  (      1.920144,      5.542986,     -2.185591),
                  (      1.920144,     18.846164,     -2.185591),
                  (      1.920144,      5.542986,     -0.783900),
                  (      1.920144,     18.846164,     -0.783900),
                  (      1.920149,      3.325794,      0.000000),
                  (      1.920150,     12.194575,     -2.185591),
                  (      1.920150,     25.497755,     -2.185591),
                  (      1.920150,     32.149352,     -2.185591),
                  (      1.920150,     38.800933,     -2.185591),
                  (      1.920150,     12.194575,     -0.783900),
                  (      1.920150,     25.497755,     -0.783900),
                  (      1.920150,     32.149352,     -0.783900),
                  (      1.920150,     38.800933,     -0.783900),
                  (      1.920150,     29.932154,      0.000000),
                  (      1.920150,     36.583733,      0.000000),
                  (      1.920152,      9.977377,      0.000000),
                  (      1.920152,     16.628972,      0.000000),
                  (      1.920152,     23.280553,      0.000000),
                  (      1.927518,     14.333273,      3.282423),
                  (      1.899699,      1.046469,      3.220428),
                  (      1.964232,     23.172097,      2.349666),
                  (      1.937619,     34.301563,      3.281571),
                  (      1.936860,     27.636949,      3.280473),
                  (      1.963071,      7.683497,      3.272746),
                  (      1.983659,     36.496439,      2.352004),
                  (      1.976165,     16.534980,      2.377345),
                  (      1.923667,     14.337669,      5.729294),
                  (      1.983455,      3.215158,      2.312851),
                  (      1.902077,     27.632369,      5.738588),
                  (      2.110466,     34.355443,      5.736637),
                  (      2.092871,      7.675109,      5.711157),
                  (      2.134182,     29.941362,      6.446964),
                  (      2.129376,     12.028309,      6.474015),
                  (      3.822944,      4.314657,      3.264921),
                  (      3.853225,     19.924914,      2.438693),
                  (      3.840293,     22.171959,     -2.185591),
                  (      3.840293,     22.171959,     -0.783900),
                  (      3.840294,      8.868783,     -2.185591),
                  (      3.840294,     35.475138,     -2.185591),
                  (      3.840294,      8.868783,     -0.783900),
                  (      3.840294,     35.475138,     -0.783900),
                  (      3.840294,      6.651582,      0.000000),
                  (      3.840294,     33.257941,      0.000000),
                  (      3.840294,     19.954761,      0.000000),
                  (      3.840301,     15.520374,     -2.185591),
                  (      3.840301,     28.823553,     -2.185591),
                  (      3.840301,     15.520374,     -0.783900),
                  (      3.840301,     28.823553,     -0.783900),
                  (      3.840301,     26.606352,      0.000000),
                  (      3.840303,     13.303176,      0.000000),
                  (      3.846944,     30.969950,      2.870324),
                  (      3.862941,     13.281719,      2.375074),
                  (      3.888755,      6.567273,      2.443350),
                  (      3.852564,     17.703008,      3.273537),
                  (      3.833400,     30.938576,      5.193938),
                  (      3.850982,     26.538716,      2.404132),
                  (      3.847611,     24.309519,      3.263666),
                  (      3.860706,     11.013199,      2.921276),
                  (      3.861079,     33.240583,      2.355124),
                  (      3.841637,     11.066603,      5.209299),
                  (      3.828440,     17.848655,      5.714159),
                  (      3.837532,     24.144850,      5.713784),
                  (      3.938892,     15.507543,      5.878523),
                  (      3.944350,     26.497157,      5.942478),
                  (      3.870550,     30.936416,      8.022840),
                  (      3.958301,     32.920937,      6.478004),
                  (      3.879093,     11.016476,      8.035417),
                  (      3.964815,      9.070092,      6.446088),
                  (      5.616981,     30.064613,      6.442306),
                  (      5.635843,     11.909520,      6.467952),
                  (      5.760443,     25.497758,     -2.185591),
                  (      5.760443,     12.194580,     -2.185591),
                  (      5.760443,     25.497758,     -0.783900),
                  (      5.760443,     12.194580,     -0.783900),
                  (      5.760443,      9.977380,      0.000000),
                  (      5.760446,     18.846160,     -2.185591),
                  (      5.760446,     32.149340,     -2.185591),
                  (      5.760446,     18.846160,     -0.783900),
                  (      5.760446,     32.149340,     -0.783900),
                  (      5.760446,     23.280559,      0.000000),
                  (      5.760448,     16.628962,      0.000000),
                  (      5.760453,     29.932151,      0.000000),
                  (      5.726544,     23.183198,      2.343344),
                  (      5.715014,     16.525397,      2.375410),
                  (      5.819143,     29.848728,      2.348579),
                  (      5.798988,     14.314132,      3.290223),
                  (      5.772197,     27.648527,      3.258058),
                  (      5.793891,      7.733766,      3.286894),
                  (      5.814696,      9.878135,      2.322134),
                  (      5.828900,     20.993655,      3.258970),
                  (      5.934053,     27.744796,      5.711900),
                  (      5.910345,     14.240372,      5.728980),
                  (      5.996868,     18.812432,      6.310628),
                  (      5.985248,     23.227827,      6.426345),
                  (      6.546255,     21.054991,      5.571336),
                  (      7.664266,     10.969578,      3.281873),
                  (      7.706665,     19.962554,      2.316180),
                  (      7.680589,     15.520369,     -2.185591),
                  (      7.680589,     15.520369,     -0.783900),
                  (      7.680595,     28.823557,     -2.185591),
                  (      7.680595,     28.823557,     -0.783900),
                  (      7.680595,     13.303179,      0.000000),
                  (      7.680596,     22.171959,     -2.185591),
                  (      7.680596,     22.171959,     -0.783900),
                  (      7.680596,     19.954760,      0.000000),
                  (      7.680596,     26.606356,      0.000000),
                  (      7.670384,     17.659721,      5.132284),
                  (      7.693526,     24.308851,      2.876243),
                  (      7.694543,     26.580576,      2.343653),
                  (      7.726969,     13.211192,      2.443420),
                  (      7.702892,     17.675619,      2.807052),
                  (      7.632005,     24.377783,      5.206293),
                  (      7.800004,     26.350380,      6.463902),
                  (      7.748516,     17.711213,      7.662213),
                  (      7.725043,     24.290178,      7.880434),
                  (      7.831095,     15.634598,      6.322861),
                  (      8.987030,     21.056746,      5.534080),
                  (      9.440374,     23.327680,      6.365109),
                  (      9.487631,     18.712938,      6.270648),
                  (      9.627314,     20.971142,      3.215069),
                  (      9.600739,     18.846164,     -2.185591),
                  (      9.600739,     18.846164,     -0.783900),
                  (      9.600740,     16.628965,      0.000000),
                  (      9.600747,     25.497757,     -2.185591),
                  (      9.600747,     25.497757,     -0.783900),
                  (      9.600747,     23.280557,      0.000000),
                  (      9.665440,     23.188195,      2.344954),
                  (      9.639561,     14.400465,      3.255675),
                  (      9.690419,     16.549828,      2.313899),
                  (     11.535206,     17.709230,      3.203472),
                  (     11.520889,     22.171963,     -2.185591),
                  (     11.520889,     22.171963,     -0.783900),
                  (     11.520889,     19.954764,      0.000000),
                  (     11.668673,     19.982752,      2.441314)])
atoms.set_cell([(    -13.44079900,     23.28040000,      0.00000000),
                  (     13.44079900,     23.28040000,      0.00000000),
                  (      0.00000000,      0.00000000,     99.00000000)])

FRAGMENTS = FixAtoms(indices=[1,2,3,11,12,13,14,15,16,27,28,29,30,31,32,33,34,35,55,56,57,58,59,60,61,62,63,64,65,66,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,283,284,285,286,287,288,289,290,291,292,293,294,310,311,312,313,314,315,316,317,318,333,334,335,336,337,338,343,344,345])
atoms.set_constraint(FRAGMENTS)


start_time = time.time()

kpoints = monkhorst_pack([4, 4, 1])
C='Mulliken_dipole_preserving'
atoms.calc = Fireball(charges_method=C,kpts=kpoints)
#atoms.calc = Fireball(charges=C)
atoms.get_potential_energy()
print('ETOT =',atoms.calc.results['energy'],'eV')
end_time = time.time()
print(f"Tiempo de ejecución: {end_time - start_time} segundos")