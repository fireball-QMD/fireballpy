import unittest
import json
from numpy.testing import assert_allclose
from importlib.util import find_spec
from ase.build import molecule
from fireballpy import Fireball


class Test(unittest.TestCase):
    TOL = 1e-8

    def setUp(self):
        self.assertAllClose = assert_allclose  # Better syntax
        self.atoms = molecule('CH4')
        self.data = json.loads(RAW_DATA)
        self.isdftd3 = find_spec('dftd3') is not None

    def test_lowdin_improved(self):
        c = 'lowdin'
        d = 'improved'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_mulliken_improved(self):
        c = 'mulliken'
        d = 'improved'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_mulliken_dipole_improved(self):
        c = 'mulliken_dipole'
        d = 'improved'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_npa_improved(self):
        c = 'npa'
        d = 'improved'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_lowdin_legacy(self):
        c = 'lowdin'
        d = 'legacy'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_mulliken_legacy(self):
        c = 'mulliken'
        d = 'legacy'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_mulliken_dipole_legacy(self):
        c = 'mulliken_dipole'
        d = 'legacy'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_npa_legacy(self):
        c = 'npa'
        d = 'legacy'
        data = self.data[c][d]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d, correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, dipole_method=d)
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)


RAW_DATA = """
{
    "lowdin": {
        "improved": {
            "none": {
                "energy": -208.38005830684278,
                "charges": [
                    -0.4450714274312948,
                    0.1112678533461885,
                    0.11126785334618816,
                    0.1112678603694588,
                    0.11126786036945924
                ],
                "forces": [
                    [
                        -7.419119030065145e-16,
                        9.528858120669776e-16,
                        -4.383524904941055e-07
                    ],
                    [
                        -1.5560034019433733,
                        -1.5560034019433635,
                        -1.5560034005338483
                    ],
                    [
                        1.5560034019433702,
                        1.5560034019433746,
                        -1.5560034005338528
                    ],
                    [
                        -1.5560036183006147,
                        1.5560036183006085,
                        1.5560036197100908
                    ],
                    [
                        1.556003618300621,
                        -1.5560036183006183,
                        1.5560036197101081
                    ]
                ]
            },
            "dftd3": {
                "energy": -208.56448015852982,
                "charges": [
                    -0.4450714274312948,
                    0.1112678533461885,
                    0.11126785334618816,
                    0.1112678603694588,
                    0.11126786036945924
                ],
                "forces": [
                    [
                        -7.417376782663482e-16,
                        9.527115873268115e-16,
                        -4.3835249049393126e-07
                    ],
                    [
                        -1.5631785328050307,
                        -1.563178532805021,
                        -1.5631785313955058
                    ],
                    [
                        1.5631785328050276,
                        1.563178532805032,
                        -1.5631785313955102
                    ],
                    [
                        -1.5631787491622722,
                        1.563178749162266,
                        1.5631787505717483
                    ],
                    [
                        1.5631787491622784,
                        -1.5631787491622757,
                        1.5631787505717656
                    ]
                ]
            }
        },
        "legacy": {
            "none": {
                "energy": -208.38005830684327,
                "charges": [
                    -0.44507142743128947,
                    0.11126785334618661,
                    0.11126785334618672,
                    0.11126786036945813,
                    0.1112678603694589
                ],
                "forces": [
                    [
                        -4.676566848369638e-06,
                        2.518922759784538e-13,
                        -9.791486003996148e-06
                    ],
                    [
                        -1.5560004314808533,
                        -1.5560006999621616,
                        -1.5560001615900876
                    ],
                    [
                        1.5560027697641572,
                        1.556002501282849,
                        -1.5560019629107877
                    ],
                    [
                        -1.5560042504796092,
                        1.5560045189609073,
                        1.5560050573329733
                    ],
                    [
                        1.556006588763145,
                        -1.5560063202818433,
                        1.5560068586538978
                    ]
                ]
            },
            "dftd3": {
                "energy": -208.5644801585303,
                "charges": [
                    -0.44507142743128947,
                    0.11126785334618661,
                    0.11126785334618672,
                    0.11126786036945813,
                    0.1112678603694589
                ],
                "forces": [
                    [
                        -4.676566848369463e-06,
                        2.5189218886608373e-13,
                        -9.791486003995887e-06
                    ],
                    [
                        -1.5631755623425108,
                        -1.563175830823819,
                        -1.563175292451745
                    ],
                    [
                        1.5631779006258146,
                        1.5631776321445063,
                        -1.5631770937724452
                    ],
                    [
                        -1.5631793813412667,
                        1.5631796498225647,
                        1.5631801881946308
                    ],
                    [
                        1.5631817196248023,
                        -1.5631814511435007,
                        1.5631819895155552
                    ]
                ]
            }
        }
    },
    "mulliken": {
        "improved": {
            "none": {
                "energy": -207.13130160816561,
                "charges": [
                    -0.5181323325955076,
                    0.12953307793008373,
                    0.12953307793008317,
                    0.12953308836767008,
                    0.1295330883676702
                ],
                "forces": [
                    [
                        1.933292875832063e-14,
                        -1.0308823447509072e-14,
                        -5.401885519841289e-07
                    ],
                    [
                        -1.618871128712557,
                        -1.6188711287125561,
                        -1.6188711354458944
                    ],
                    [
                        1.6188711287125446,
                        1.6188711287125437,
                        -1.6188711354458896
                    ],
                    [
                        -1.618871412273578,
                        1.6188714122735717,
                        1.6188714055401747
                    ],
                    [
                        1.6188714122735681,
                        -1.6188714122735637,
                        1.6188714055401614
                    ]
                ]
            },
            "dftd3": {
                "energy": -207.31572345985265,
                "charges": [
                    -0.5181323325955076,
                    0.12953307793008373,
                    0.12953307793008317,
                    0.12953308836767008,
                    0.1295330883676702
                ],
                "forces": [
                    [
                        1.9333190095430878e-14,
                        -1.0308823447509072e-14,
                        -5.401885519838676e-07
                    ],
                    [
                        -1.6260462595742144,
                        -1.6260462595742136,
                        -1.626046266307552
                    ],
                    [
                        1.626046259574202,
                        1.6260462595742011,
                        -1.626046266307547
                    ],
                    [
                        -1.6260465431352353,
                        1.6260465431352291,
                        1.6260465364018322
                    ],
                    [
                        1.6260465431352256,
                        -1.6260465431352211,
                        1.6260465364018188
                    ]
                ]
            }
        },
        "legacy": {
            "none": {
                "energy": -207.1313016081658,
                "charges": [
                    -0.5181323325955027,
                    0.12953307793008106,
                    0.12953307793008106,
                    0.12953308836767152,
                    0.12953308836766997
                ],
                "forces": [
                    [
                        -5.558287046496967e-06,
                        2.574590663971391e-13,
                        -1.1656762555104588e-05
                    ],
                    [
                        -1.618867598198805,
                        -1.618867917299604,
                        -1.6188672858314552
                    ],
                    [
                        1.618870377342216,
                        1.618870058241408,
                        -1.6188694267732764
                    ],
                    [
                        -1.6188721636436934,
                        1.6188724827444778,
                        1.6188731142126151
                    ],
                    [
                        1.6188749427873383,
                        -1.6188746236865463,
                        1.6188752551546686
                    ]
                ]
            },
            "dftd3": {
                "energy": -207.31572345985285,
                "charges": [
                    -0.5181323325955027,
                    0.12953307793008106,
                    0.12953307793008106,
                    0.12953308836767152,
                    0.12953308836766997
                ],
                "forces": [
                    [
                        -5.558287046496706e-06,
                        2.5745897928476903e-13,
                        -1.1656762555104327e-05
                    ],
                    [
                        -1.6260427290604624,
                        -1.6260430481612616,
                        -1.6260424166931127
                    ],
                    [
                        1.6260455082038734,
                        1.6260451891030654,
                        -1.6260445576349338
                    ],
                    [
                        -1.6260472945053508,
                        1.6260476136061353,
                        1.6260482450742726
                    ],
                    [
                        1.6260500736489958,
                        -1.6260497545482038,
                        1.626050386016326
                    ]
                ]
            }
        }
    },
    "mulliken_dipole": {
        "improved": {
            "none": {
                "energy": -217.0911637299205,
                "charges": [
                    -0.004199778908609009,
                    0.0010499401956590715,
                    0.0010499401956575172,
                    0.001049949258646432,
                    0.0010499492586454329
                ],
                "forces": [
                    [
                        6.306806702222358e-14,
                        9.202700568687676e-17,
                        -4.958156209086564e-07
                    ],
                    [
                        -0.6599332537683553,
                        -0.659933253768366,
                        -0.6599332537057876
                    ],
                    [
                        0.6599332537683269,
                        0.6599332537683198,
                        -0.659933253705745
                    ],
                    [
                        -0.6599335015510586,
                        0.6599335015510683,
                        0.6599335016135925
                    ],
                    [
                        0.659933501551023,
                        -0.6599335015510239,
                        0.6599335016135561
                    ]
                ]
            },
            "dftd3": {
                "energy": -217.27558558160754,
                "charges": [
                    -0.004199778908609009,
                    0.0010499401956590715,
                    0.0010499401956575172,
                    0.001049949258646432,
                    0.0010499492586454329
                ],
                "forces": [
                    [
                        6.306824124696375e-14,
                        9.193989331679363e-17,
                        -4.958156209083951e-07
                    ],
                    [
                        -0.6671083846300128,
                        -0.6671083846300234,
                        -0.667108384567445
                    ],
                    [
                        0.6671083846299843,
                        0.6671083846299772,
                        -0.6671083845674024
                    ],
                    [
                        -0.667108632412716,
                        0.6671086324127258,
                        0.66710863247525
                    ],
                    [
                        0.6671086324126805,
                        -0.6671086324126814,
                        0.6671086324752136
                    ]
                ]
            }
        },
        "legacy": {
            "none": {
                "energy": -217.09116372992088,
                "charges": [
                    -0.0041997789086081205,
                    0.0010499401956579613,
                    0.0010499401956574062,
                    0.00104994925864621,
                    0.001049949258646432
                ],
                "forces": [
                    [
                        -3.768235144279904e-08,
                        6.481426007256542e-15,
                        -5.711802338125906e-07
                    ],
                    [
                        -0.659933229833221,
                        -0.6599332319965465,
                        -0.6599332276073615
                    ],
                    [
                        0.6599332486743665,
                        0.6599332465110357,
                        -0.6599332421218289
                    ],
                    [
                        -0.6599335066449075,
                        0.6599335088082636,
                        0.6599335131974486
                    ],
                    [
                        0.6599335254861085,
                        -0.659933523322751,
                        0.659933527711972
                    ]
                ]
            },
            "dftd3": {
                "energy": -217.2755855816079,
                "charges": [
                    -0.0041997789086081205,
                    0.0010499401956579613,
                    0.0010499401956574062,
                    0.00104994925864621,
                    0.001049949258646432
                ],
                "forces": [
                    [
                        -3.768235144262482e-08,
                        6.481338894886459e-15,
                        -5.711802338123293e-07
                    ],
                    [
                        -0.6671083606948784,
                        -0.667108362858204,
                        -0.667108358469019
                    ],
                    [
                        0.667108379536024,
                        0.6671083773726931,
                        -0.6671083729834864
                    ],
                    [
                        -0.6671086375065649,
                        0.6671086396699211,
                        0.6671086440591061
                    ],
                    [
                        0.667108656347766,
                        -0.6671086541844085,
                        0.6671086585736294
                    ]
                ]
            }
        }
    },
    "npa": {
        "improved": {
            "none": {
                "energy": -208.3800583068429,
                "charges": [
                    -0.4450714274312926,
                    0.11126785334618794,
                    0.11126785334618738,
                    0.11126786036945824,
                    0.11126786036945846
                ],
                "forces": [
                    [
                        -4.949066601412505e-16,
                        2.465295678950685e-14,
                        -4.3835250313196697e-07
                    ],
                    [
                        -1.5560034019433617,
                        -1.556003401943359,
                        -1.5560034005338341
                    ],
                    [
                        1.5560034019433537,
                        1.5560034019433466,
                        -1.556003400533827
                    ],
                    [
                        -1.556003618300605,
                        1.5560036183006023,
                        1.5560036197100748
                    ],
                    [
                        1.5560036183006156,
                        -1.5560036183006107,
                        1.55600361971009
                    ]
                ]
            },
            "dftd3": {
                "energy": -208.56448015852993,
                "charges": [
                    -0.4450714274312926,
                    0.11126785334618794,
                    0.11126785334618738,
                    0.11126786036945824,
                    0.11126786036945846
                ],
                "forces": [
                    [
                        -4.946453230310011e-16,
                        2.465295678950685e-14,
                        -4.383525031317056e-07
                    ],
                    [
                        -1.5631785328050192,
                        -1.5631785328050165,
                        -1.5631785313954916
                    ],
                    [
                        1.5631785328050112,
                        1.563178532805004,
                        -1.5631785313954845
                    ],
                    [
                        -1.5631787491622624,
                        1.5631787491622597,
                        1.5631787505717323
                    ],
                    [
                        1.563178749162273,
                        -1.5631787491622682,
                        1.5631787505717474
                    ]
                ]
            }
        },
        "legacy": {
            "none": {
                "energy": -208.38005830684313,
                "charges": [
                    -0.44507142743128947,
                    0.11126785334618716,
                    0.11126785334618716,
                    0.11126786036945846,
                    0.11126786036945868
                ],
                "forces": [
                    [
                        -4.676566846249537e-06,
                        2.476009457830938e-13,
                        -9.791486015616593e-06
                    ],
                    [
                        -1.5560004314808857,
                        -1.5560006999621967,
                        -1.5560001615901102
                    ],
                    [
                        1.5560027697641958,
                        1.5560025012828804,
                        -1.5560019629108073
                    ],
                    [
                        -1.5560042504796354,
                        1.5560045189609357,
                        1.5560050573330049
                    ],
                    [
                        1.5560065887631742,
                        -1.5560063202818752,
                        1.556006858653927
                    ]
                ]
            },
            "dftd3": {
                "energy": -208.56448015853016,
                "charges": [
                    -0.44507142743128947,
                    0.11126785334618716,
                    0.11126785334618716,
                    0.11126786036945846,
                    0.11126786036945868
                ],
                "forces": [
                    [
                        -4.6765668462492755e-06,
                        2.4760085867072373e-13,
                        -9.791486015616332e-06
                    ],
                    [
                        -1.5631755623425432,
                        -1.5631758308238541,
                        -1.5631752924517677
                    ],
                    [
                        1.5631779006258533,
                        1.5631776321445379,
                        -1.5631770937724647
                    ],
                    [
                        -1.5631793813412929,
                        1.5631796498225932,
                        1.5631801881946623
                    ],
                    [
                        1.5631817196248317,
                        -1.5631814511435327,
                        1.5631819895155845
                    ]
                ]
            }
        }
    }
}
"""

if __name__ == '__main__':
    unittest.main()
