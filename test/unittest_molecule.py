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
    "npa": {
        "improved": {
            "dftd3": {
                "energy": -208.56448015852973,
                "forces": [
                    [
                        3.188257808363033e-14,
                        3.4936269309991296e-14,
                        -4.169105326112646e-07
                    ],
                    [
                        -1.5631785301247925,
                        -1.563178530124796,
                        -1.563178536756041
                    ],
                    [
                        1.5631785301247627,
                        1.5631785301247645,
                        -1.563178536756013
                    ],
                    [
                        -1.563178751842529,
                        1.5631787518425246,
                        1.5631787452112829
                    ],
                    [
                        1.5631787518425326,
                        -1.563178751842529,
                        1.5631787452112909
                    ]
                ],
                "charges": [
                    -0.445071427431293,
                    0.1112678533461886,
                    0.11126785334618772,
                    0.11126786036945857,
                    0.1112678603694588
                ]
            },
            "none": {
                "energy": -208.3800583068427,
                "forces": [
                    [
                        3.1882316746520076e-14,
                        3.4936269309991296e-14,
                        -4.16910532611526e-07
                    ],
                    [
                        -1.556003399263135,
                        -1.5560033992631386,
                        -1.5560034058943835
                    ],
                    [
                        1.5560033992631053,
                        1.556003399263107,
                        -1.5560034058943555
                    ],
                    [
                        -1.5560036209808716,
                        1.5560036209808672,
                        1.5560036143496254
                    ],
                    [
                        1.5560036209808752,
                        -1.5560036209808716,
                        1.5560036143496334
                    ]
                ],
                "charges": [
                    -0.445071427431293,
                    0.1112678533461886,
                    0.11126785334618772,
                    0.11126786036945857,
                    0.1112678603694588
                ]
            }
        },
        "legacy": {
            "dftd3": {
                "energy": -208.56448015853022,
                "forces": [
                    [
                        -4.676566806159515e-06,
                        2.5564544485834995e-13,
                        -9.77004403865443e-06
                    ],
                    [
                        -1.5631755596622778,
                        -1.5631758281435821,
                        -1.563175297812266
                    ],
                    [
                        1.5631778979455493,
                        1.5631776294642552,
                        -1.5631770991329403
                    ],
                    [
                        -1.563179384021539,
                        1.5631796525028399,
                        1.563180182834152
                    ],
                    [
                        1.5631817223050788,
                        -1.5631814538237727,
                        1.5631819841550796
                    ]
                ],
                "charges": [
                    -0.4450714274312908,
                    0.11126785334618683,
                    0.1112678533461865,
                    0.11126786036945857,
                    0.11126786036945868
                ]
            },
            "none": {
                "energy": -208.38005830684318,
                "forces": [
                    [
                        -4.676566806159777e-06,
                        2.5564553197072e-13,
                        -9.77004403865469e-06
                    ],
                    [
                        -1.5560004288006204,
                        -1.5560006972819247,
                        -1.5560001669506085
                    ],
                    [
                        1.5560027670838918,
                        1.5560024986025978,
                        -1.556001968271283
                    ],
                    [
                        -1.5560042531598817,
                        1.5560045216411824,
                        1.5560050519724946
                    ],
                    [
                        1.5560065914434213,
                        -1.5560063229621153,
                        1.5560068532934221
                    ]
                ],
                "charges": [
                    -0.4450714274312908,
                    0.11126785334618683,
                    0.1112678533461865,
                    0.11126786036945857,
                    0.11126786036945868
                ]
            }
        }
    },
    "mulliken": {
        "improved": {
            "dftd3": {
                "energy": -207.31572345985268,
                "forces": [
                    [
                        1.4876722568635236e-14,
                        8.706296086916617e-15,
                        -5.172476321129665e-07
                    ],
                    [
                        -1.626046256706599,
                        -1.6260462567065963,
                        -1.6260462720428162
                    ],
                    [
                        1.6260462567065754,
                        1.6260462567065794,
                        -1.6260462720428015
                    ],
                    [
                        -1.6260465460028282,
                        1.62604654600283,
                        1.6260465306666194
                    ],
                    [
                        1.626046546002829,
                        -1.626046546002822,
                        1.6260465306666296
                    ]
                ],
                "charges": [
                    -0.5181323325955081,
                    0.12953307793008484,
                    0.12953307793008417,
                    0.12953308836766997,
                    0.12953308836766997
                ]
            },
            "none": {
                "energy": -207.13130160816564,
                "forces": [
                    [
                        1.487654834389507e-14,
                        8.706470311656783e-15,
                        -5.172476321131408e-07
                    ],
                    [
                        -1.6188711258449415,
                        -1.6188711258449389,
                        -1.6188711411811587
                    ],
                    [
                        1.618871125844918,
                        1.618871125844922,
                        -1.618871141181144
                    ],
                    [
                        -1.6188714151411707,
                        1.6188714151411725,
                        1.618871399804962
                    ],
                    [
                        1.6188714151411716,
                        -1.6188714151411645,
                        1.6188713998049722
                    ]
                ],
                "charges": [
                    -0.5181323325955081,
                    0.12953307793008484,
                    0.12953307793008417,
                    0.12953308836766997,
                    0.12953308836766997
                ]
            }
        },
        "legacy": {
            "dftd3": {
                "energy": -207.31572345985302,
                "forces": [
                    [
                        -5.558287094779126e-06,
                        3.2694108638284385e-13,
                        -1.163382167674267e-05
                    ],
                    [
                        -1.626042726192815,
                        -1.6260430452936172,
                        -1.626042422428333
                    ],
                    [
                        1.626045505336209,
                        1.6260451862354073,
                        -1.6260445633701315
                    ],
                    [
                        -1.6260472973729156,
                        1.6260476164737099,
                        1.6260482393390057
                    ],
                    [
                        1.6260500765166352,
                        -1.626049757415834,
                        1.6260503802811244
                    ]
                ],
                "charges": [
                    -0.5181323325955001,
                    0.12953307793008018,
                    0.12953307793007984,
                    0.12953308836766975,
                    0.12953308836767063
                ]
            },
            "none": {
                "energy": -207.13130160816598,
                "forces": [
                    [
                        -5.558287094779301e-06,
                        3.269411734952139e-13,
                        -1.163382167674293e-05
                    ],
                    [
                        -1.6188675953311575,
                        -1.6188679144319598,
                        -1.6188672915666755
                    ],
                    [
                        1.6188703744745516,
                        1.6188700553737498,
                        -1.618869432508474
                    ],
                    [
                        -1.6188721665112582,
                        1.6188724856120524,
                        1.6188731084773482
                    ],
                    [
                        1.6188749456549778,
                        -1.6188746265541765,
                        1.618875249419467
                    ]
                ],
                "charges": [
                    -0.5181323325955001,
                    0.12953307793008018,
                    0.12953307793007984,
                    0.12953308836766975,
                    0.12953308836767063
                ]
            }
        }
    },
    "mulliken_dipole": {
        "improved": {
            "dftd3": {
                "energy": -217.27558558160763,
                "forces": [
                    [
                        3.5321016295830704e-14,
                        5.633392911827495e-15,
                        -4.783248062767548e-07
                    ],
                    [
                        -0.6671083824436492,
                        -0.6671083824436439,
                        -0.6671083889401401
                    ],
                    [
                        0.6671083824436146,
                        0.6671083824436332,
                        -0.6671083889401273
                    ],
                    [
                        -0.6671086345990498,
                        0.6671086345990436,
                        0.6671086281025429
                    ],
                    [
                        0.6671086345990391,
                        -0.6671086345990285,
                        0.6671086281025322
                    ]
                ],
                "charges": [
                    -0.004199778908609453,
                    0.0010499401956591825,
                    0.0010499401956584054,
                    0.001049949258645766,
                    0.001049949258645766
                ]
            },
            "none": {
                "energy": -217.0911637299206,
                "forces": [
                    [
                        3.5320842071090537e-14,
                        5.633480024197578e-15,
                        -4.783248062770161e-07
                    ],
                    [
                        -0.6599332515819918,
                        -0.6599332515819865,
                        -0.6599332580784827
                    ],
                    [
                        0.6599332515819571,
                        0.6599332515819758,
                        -0.6599332580784698
                    ],
                    [
                        -0.6599335037373923,
                        0.6599335037373861,
                        0.6599334972408855
                    ],
                    [
                        0.6599335037373817,
                        -0.659933503737371,
                        0.6599334972408748
                    ]
                ],
                "charges": [
                    -0.004199778908609453,
                    0.0010499401956591825,
                    0.0010499401956584054,
                    0.001049949258645766,
                    0.001049949258645766
                ]
            }
        },
        "legacy": {
            "dftd3": {
                "energy": -217.27558558160797,
                "forces": [
                    [
                        -3.76823643420357e-08,
                        -2.1761631730976832e-14,
                        -5.536893155117133e-07
                    ],
                    [
                        -0.6671083585085167,
                        -0.6671083606718413,
                        -0.6671083628417349
                    ],
                    [
                        0.6671083773496866,
                        0.6671083751863625,
                        -0.6671083773562201
                    ],
                    [
                        -0.667108639692924,
                        0.6671086418562748,
                        0.667108639686375
                    ],
                    [
                        0.667108658534125,
                        -0.6671086563707644,
                        0.6671086542008997
                    ]
                ],
                "charges": [
                    -0.0041997789086081205,
                    0.0010499401956580723,
                    0.0010499401956587384,
                    0.001049949258645988,
                    0.001049949258645877
                ]
            },
            "none": {
                "energy": -217.09116372992094,
                "forces": [
                    [
                        -3.7682364342297035e-08,
                        -2.1761544618606748e-14,
                        -5.536893155119746e-07
                    ],
                    [
                        -0.6599332276468592,
                        -0.6599332298101839,
                        -0.6599332319800775
                    ],
                    [
                        0.6599332464880292,
                        0.659933244324705,
                        -0.6599332464945626
                    ],
                    [
                        -0.6599335088312666,
                        0.6599335109946174,
                        0.6599335088247176
                    ],
                    [
                        0.6599335276724676,
                        -0.659933525509107,
                        0.6599335233392423
                    ]
                ],
                "charges": [
                    -0.0041997789086081205,
                    0.0010499401956580723,
                    0.0010499401956587384,
                    0.001049949258645988,
                    0.001049949258645877
                ]
            }
        }
    },
    "lowdin": {
        "improved": {
            "dftd3": {
                "energy": -208.56448015852968,
                "forces": [
                    [
                        4.473689387034572e-16,
                        4.274963547894674e-15,
                        -4.1691051130273486e-07
                    ],
                    [
                        -1.563178530124765,
                        -1.563178530124749,
                        -1.5631785367560005
                    ],
                    [
                        1.5631785301247485,
                        1.5631785301247516,
                        -1.563178536755997
                    ],
                    [
                        -1.5631787518424949,
                        1.5631787518424882,
                        1.5631787452112482
                    ],
                    [
                        1.5631787518425064,
                        -1.5631787518424993,
                        1.5631787452112567
                    ]
                ],
                "charges": [
                    -0.44507142743129346,
                    0.11126785334618794,
                    0.11126785334618794,
                    0.11126786036945846,
                    0.1112678603694588
                ]
            },
            "none": {
                "energy": -208.38005830684264,
                "forces": [
                    [
                        4.471076015932079e-16,
                        4.275050660264757e-15,
                        -4.169105113029962e-07
                    ],
                    [
                        -1.5560033992631075,
                        -1.5560033992630915,
                        -1.556003405894343
                    ],
                    [
                        1.556003399263091,
                        1.5560033992630942,
                        -1.5560034058943395
                    ],
                    [
                        -1.5560036209808374,
                        1.5560036209808308,
                        1.5560036143495908
                    ],
                    [
                        1.556003620980849,
                        -1.5560036209808419,
                        1.5560036143495992
                    ]
                ],
                "charges": [
                    -0.44507142743129346,
                    0.11126785334618794,
                    0.11126785334618794,
                    0.11126786036945846,
                    0.1112678603694588
                ]
            }
        },
        "legacy": {
            "dftd3": {
                "energy": -208.56448015852988,
                "forces": [
                    [
                        -4.676566830995935e-06,
                        2.513748414019776e-13,
                        -9.770044031376343e-06
                    ],
                    [
                        -1.563175559662272,
                        -1.5631758281435768,
                        -1.563175297812246
                    ],
                    [
                        1.5631778979455604,
                        1.563177629464245,
                        -1.5631770991329288
                    ],
                    [
                        -1.5631793840215187,
                        1.5631796525028154,
                        1.563180182834144
                    ],
                    [
                        1.5631817223050561,
                        -1.563181453823738,
                        1.5631819841550556
                    ]
                ],
                "charges": [
                    -0.44507142743129213,
                    0.11126785334618716,
                    0.11126785334618705,
                    0.11126786036945868,
                    0.11126786036945946
                ]
            },
            "none": {
                "energy": -208.38005830684284,
                "forces": [
                    [
                        -4.676566830996197e-06,
                        2.5137492851434765e-13,
                        -9.770044031376604e-06
                    ],
                    [
                        -1.5560004288006146,
                        -1.5560006972819194,
                        -1.5560001669505885
                    ],
                    [
                        1.556002767083903,
                        1.5560024986025875,
                        -1.5560019682712714
                    ],
                    [
                        -1.5560042531598612,
                        1.556004521641158,
                        1.5560050519724866
                    ],
                    [
                        1.5560065914433987,
                        -1.5560063229620806,
                        1.5560068532933982
                    ]
                ],
                "charges": [
                    -0.44507142743129213,
                    0.11126785334618716,
                    0.11126785334618705,
                    0.11126786036945868,
                    0.11126786036945946
                ]
            }
        }
    }
}
"""

if __name__ == '__main__':
    unittest.main()
