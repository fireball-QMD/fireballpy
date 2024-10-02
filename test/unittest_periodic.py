import unittest
import json
from numpy.testing import assert_allclose
from importlib.util import find_spec
from ase.build import diamond111
from fireballpy import Fireball


class Test(unittest.TestCase):
    TOL = 1e-8

    def setUp(self):
        self.assertAllClose = assert_allclose  # Better syntax
        self.atoms = diamond111('C', (2,2,1))
        self.data = json.loads(RAW_DATA)
        self.isdftd3 = find_spec('dftd3') is not None

    def test_lowdin(self):
        c = 'lowdin'
        data = self.data[c]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2), correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2))
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_mulliken(self):
        c = 'mulliken'
        data = self.data[c]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2), correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2))
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_mulliken_dipole(self):
        c = 'mulliken_dipole'
        data = self.data[c]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2), correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2))
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)

    def test_npa(self):
        c = 'npa'
        data = self.data[c]
        self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2), correction='off')
        self.assertAlmostEqual(self.atoms.get_potential_energy(), data['none']['energy'], delta=self.TOL)
        self.assertAllClose(self.atoms.get_charges(), data['none']['charges'], atol=self.TOL)
        self.assertAllClose(self.atoms.get_forces(), data['none']['forces'], atol=self.TOL)
        if self.isdftd3:
            self.atoms.calc = Fireball(fdata='biology', charges_method=c, kpts=(2,2,2))
            self.assertAlmostEqual(self.atoms.get_potential_energy(), data['dftd3']['energy'], delta=self.TOL)
            self.assertAllClose(self.atoms.get_charges(), data['dftd3']['charges'], atol=self.TOL)
            self.assertAllClose(self.atoms.get_forces(), data['dftd3']['forces'], atol=self.TOL)


RAW_DATA = """
{
    "npa": {
        "dftd3": {
            "energy": -5614.527608118373,
            "forces": [
                [
                    2.8774511221890733,
                    1.6745455031138277,
                    0.025308230797985123
                ],
                [
                    -1.0766765619471095,
                    3.571654563778994,
                    0.001971622708770125
                ],
                [
                    2.2145286024755704,
                    -2.6798804316717915,
                    -0.0020035586074445416
                ],
                [
                    -4.015303162717506,
                    -2.566319635221016,
                    -0.025276294900120178
                ]
            ],
            "charges": [
                -0.014397901375988287,
                0.006478995747229543,
                0.011388726764153922,
                -0.003469821135395179
            ]
        },
        "none": {
            "energy": -5614.176180153041,
            "forces": [
                [
                    2.8153948246387683,
                    1.6387172830182477,
                    0.025308230797985123
                ],
                [
                    -1.0410854031131076,
                    3.5100088683782484,
                    0.001971622708770125
                ],
                [
                    2.1789374436415683,
                    -2.618234736271046,
                    -0.0020035586074445416
                ],
                [
                    -3.953246865167202,
                    -2.530491415125436,
                    -0.025276294900120178
                ]
            ],
            "charges": [
                -0.014397901375988287,
                0.006478995747229543,
                0.011388726764153922,
                -0.003469821135395179
            ]
        }
    },
    "mulliken": {
        "dftd3": {
            "energy": -5628.870602030968,
            "forces": [
                [
                    2.9433716304102218,
                    1.7116635296572338,
                    0.027073642216087088
                ],
                [
                    -1.1156801723673373,
                    3.6496720668945173,
                    0.002514456209727692
                ],
                [
                    2.2612731938916766,
                    -2.7515211354409885,
                    -0.0027620695999456615
                ],
                [
                    -4.088964651934621,
                    -2.6098144611108913,
                    -0.026826028826241466
                ]
            ],
            "charges": [
                -0.014292676711373797,
                0.006347903561197299,
                0.011280545893961769,
                -0.0033357727437852702
            ]
        },
        "none": {
            "energy": -5628.519174065636,
            "forces": [
                [
                    2.881315332859917,
                    1.6758353095616538,
                    0.027073642216087088
                ],
                [
                    -1.0800890135333354,
                    3.588026371493772,
                    0.002514456209727692
                ],
                [
                    2.2256820350576745,
                    -2.689875440040243,
                    -0.0027620695999456615
                ],
                [
                    -4.026908354384316,
                    -2.5739862410153114,
                    -0.026826028826241466
                ]
            ],
            "charges": [
                -0.014292676711373797,
                0.006347903561197299,
                0.011280545893961769,
                -0.0033357727437852702
            ]
        }
    },
    "mulliken_dipole": {
        "dftd3": {
            "energy": -5628.995438829362,
            "forces": [
                [
                    2.9452514413679864,
                    1.7123059829535547,
                    0.027025937033575823
                ],
                [
                    -1.1181970790574864,
                    3.649758427969985,
                    0.002564177293457509
                ],
                [
                    2.26136165194304,
                    -2.753168288063421,
                    -0.0027220562120513422
                ],
                [
                    -4.088416014253496,
                    -2.6088961228600795,
                    -0.02686805811437569
                ]
            ],
            "charges": [
                -0.014365854525192567,
                0.006257766461340797,
                0.011331254286017778,
                -0.003223166222165119
            ]
        },
        "none": {
            "energy": -5628.64401086403,
            "forces": [
                [
                    2.8831951438176815,
                    1.6764777628579748,
                    0.027025937033575823
                ],
                [
                    -1.0826059202234846,
                    3.5881127325692397,
                    0.002564177293457509
                ],
                [
                    2.225770493109038,
                    -2.691522592662676,
                    -0.0027220562120513422
                ],
                [
                    -4.026359716703191,
                    -2.5730679027644996,
                    -0.02686805811437569
                ]
            ],
            "charges": [
                -0.014365854525192567,
                0.006257766461340797,
                0.011331254286017778,
                -0.003223166222165119
            ]
        }
    },
    "lowdin": {
        "dftd3": {
            "energy": -5614.527608118368,
            "forces": [
                [
                    2.8774511221892016,
                    1.6745455031139467,
                    0.025308230798295306
                ],
                [
                    -1.0766765619471765,
                    3.571654563778921,
                    0.0019716227090920066
                ],
                [
                    2.2145286024754753,
                    -2.6798804316718687,
                    -0.002003558607164016
                ],
                [
                    -4.01530316271743,
                    -2.5663196352209616,
                    -0.025276294900164573
                ]
            ],
            "charges": [
                -0.014397901375988287,
                0.0064789957472299875,
                0.011388726764157031,
                -0.0034698211353965114
            ]
        },
        "none": {
            "energy": -5614.1761801530365,
            "forces": [
                [
                    2.8153948246388967,
                    1.6387172830183667,
                    0.025308230798295306
                ],
                [
                    -1.0410854031131747,
                    3.5100088683781756,
                    0.0019716227090920066
                ],
                [
                    2.1789374436414732,
                    -2.6182347362711234,
                    -0.002003558607164016
                ],
                [
                    -3.953246865167125,
                    -2.5304914151253817,
                    -0.025276294900164573
                ]
            ],
            "charges": [
                -0.014397901375988287,
                0.0064789957472299875,
                0.011388726764157031,
                -0.0034698211353965114
            ]
        }
    }
}
"""

if __name__ == '__main__':
    unittest.main()
