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
    "lowdin": {
        "none": {
            "energy": -5614.176180153027,
            "charges": [
                -0.014397901375988731,
                0.006478995747228211,
                0.011388726764155255,
                -0.0034698211353969555
            ],
            "forces": [
                [
                    2.8153948246389526,
                    1.6387172830185168,
                    0.02530823079822178
                ],
                [
                    -1.04108540311322,
                    3.5100088683780872,
                    0.00197162270872156
                ],
                [
                    2.178937443641406,
                    -2.6182347362713396,
                    -0.0020035586072675304
                ],
                [
                    -3.9532468651670376,
                    -2.5304914151252422,
                    -0.025276294900924604
                ]
            ]
        },
        "dftd3": {
            "energy": -5614.527608118359,
            "charges": [
                -0.014397901375988731,
                0.006478995747228211,
                0.011388726764155255,
                -0.0034698211353969555
            ],
            "forces": [
                [
                    2.8774511221892576,
                    1.6745455031140968,
                    0.02530823079822178
                ],
                [
                    -1.0766765619472218,
                    3.5716545637788326,
                    0.00197162270872156
                ],
                [
                    2.2145286024754083,
                    -2.679880431672085,
                    -0.0020035586072675304
                ],
                [
                    -4.015303162717342,
                    -2.566319635220822,
                    -0.025276294900924604
                ]
            ]
        }
    },
    "mulliken": {
        "none": {
            "energy": -5628.519174065628,
            "charges": [
                -0.014292676711373353,
                0.0063479035611950785,
                0.011280545893967986,
                -0.0033357727437879348
            ],
            "forces": [
                [
                    2.881315332860082,
                    1.6758353095617102,
                    0.027073642216281293
                ],
                [
                    -1.0800890135335441,
                    3.5880263714938985,
                    0.0025144562092419553
                ],
                [
                    2.225682035057557,
                    -2.68987544004053,
                    -0.0027620696000611733
                ],
                [
                    -4.026908354384132,
                    -2.5739862410152012,
                    -0.026826028826162654
                ]
            ]
        },
        "dftd3": {
            "energy": -5628.87060203096,
            "charges": [
                -0.014292676711373353,
                0.0063479035611950785,
                0.011280545893967986,
                -0.0033357727437879348
            ],
            "forces": [
                [
                    2.943371630410387,
                    1.7116635296572902,
                    0.027073642216281293
                ],
                [
                    -1.115680172367546,
                    3.649672066894644,
                    0.0025144562092419553
                ],
                [
                    2.261273193891559,
                    -2.7515211354412754,
                    -0.0027620696000611733
                ],
                [
                    -4.088964651934437,
                    -2.609814461110781,
                    -0.026826028826162654
                ]
            ]
        }
    },
    "mulliken_dipole": {
        "none": {
            "energy": -5628.644010864022,
            "charges": [
                -0.014365854525190791,
                0.006257766461347902,
                0.011331254286008452,
                -0.003223166222164675
            ],
            "forces": [
                [
                    2.883195143817764,
                    1.6764777628580654,
                    0.027025937033443567
                ],
                [
                    -1.082605920223751,
                    3.588112732569195,
                    0.0025641772941653734
                ],
                [
                    2.22577049310909,
                    -2.6915225926629094,
                    -0.0027220562122520775
                ],
                [
                    -4.026359716703076,
                    -2.5730679027643824,
                    -0.02686805811476456
                ]
            ]
        },
        "dftd3": {
            "energy": -5628.995438829354,
            "charges": [
                -0.014365854525190791,
                0.006257766461347902,
                0.011331254286008452,
                -0.003223166222164675
            ],
            "forces": [
                [
                    2.945251441368069,
                    1.7123059829536453,
                    0.027025937033443567
                ],
                [
                    -1.118197079057753,
                    3.64975842796994,
                    0.0025641772941653734
                ],
                [
                    2.261361651943092,
                    -2.7531682880636548,
                    -0.0027220562122520775
                ],
                [
                    -4.08841601425338,
                    -2.6088961228599623,
                    -0.02686805811476456
                ]
            ]
        }
    },
    "npa": {
        "none": {
            "energy": -5614.17618015304,
            "charges": [
                -0.014397901375987843,
                0.006478995747228655,
                0.011388726764152146,
                -0.003469821135394291
            ],
            "forces": [
                [
                    2.815394824638861,
                    1.6387172830185266,
                    0.02530823079793236
                ],
                [
                    -1.0410854031131938,
                    3.5100088683782484,
                    0.001971622709333619
                ],
                [
                    2.1789374436414413,
                    -2.6182347362714333,
                    -0.0020035586067469885
                ],
                [
                    -3.9532468651670998,
                    -2.5304914151253755,
                    -0.025276294900588497
                ]
            ]
        },
        "dftd3": {
            "energy": -5614.527608118372,
            "charges": [
                -0.014397901375987843,
                0.006478995747228655,
                0.011388726764152146,
                -0.003469821135394291
            ],
            "forces": [
                [
                    2.877451122189166,
                    1.6745455031141065,
                    0.02530823079793236
                ],
                [
                    -1.0766765619471956,
                    3.571654563778994,
                    0.001971622709333619
                ],
                [
                    2.2145286024754434,
                    -2.6798804316721787,
                    -0.0020035586067469885
                ],
                [
                    -4.015303162717404,
                    -2.5663196352209554,
                    -0.025276294900588497
                ]
            ]
        }
    }
}
"""

if __name__ == '__main__':
    unittest.main()
