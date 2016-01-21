#!/usr/bin/env python
# $Id: healalm.py,v 1.1 2008/07/04 22:12:03 oxon Exp $
# Author: Akira Okumura 2008/06/26

import unittest

import ROOT

class TestHealAlm(unittest.TestCase):
    """
    Unit test class for THealAlm.
    """
    def setUp(self):
        """
        Set up function. libRootHealPix must be loaded.
        """
        ROOT.std.complex("double") # hack
        ROOT.std.complex("float")
        ROOT.gSystem.Load("libRootHealPix")

    def testGetSet(self):
        """
        Check basic getter/setter.
        """
        lmax = 5
        mmax = 3
        alm = ROOT.THealAlm('Double_t')(lmax, mmax)
        p = alm(1, 1)
        for l in range(lmax + 1):
            for m in range(mmax + 1):
                alm.SetAlm(l, m, ROOT.std.complex('double')(l, m))

        for l in range(lmax + 1):
            for m in range(mmax + 1):
                if (l <= lmax) and (m <= min(l, mmax)):
                    rhs = ROOT.std.complex('double')(l, m)
                else:
                    rhs = ROOT.std.complex('double')(0, 0)
                self.assertEqual(alm(l, m).real(), rhs.real())
                self.assertEqual(alm(l, m).imag(), rhs.imag())

        alm = ROOT.THealAlm('Double_t')(-1, 3)
        self.assertEqual(alm.GetLmax(), 0)
        self.assertEqual(alm.GetMmax(), 0)

        alm = ROOT.THealAlm('Double_t')(3, 5)
        self.assertEqual(alm.GetLmax(), 3)
        self.assertEqual(alm.GetMmax(), 3)

        alm = ROOT.THealAlm('Double_t')(3, -1)
        self.assertEqual(alm.GetLmax(), 3)
        self.assertEqual(alm.GetMmax(), 0)

    def testArithmeticOperations(self):
        """
        Check basic operations
        """

        lmax = 5
        mmax = 3

        alm1 = ROOT.THealAlm('Double_t')(lmax, mmax)
        alm2 = ROOT.THealAlm('Double_t')(lmax, mmax)

        for l in range(lmax + 1):
            for m in range(mmax + l):
                alm1.SetAlm(l, m, ROOT.std.complex('Double_t')(l*l, m*m))
                alm2.SetAlm(l, m, ROOT.std.complex('Double_t')(l, m))

        alm = alm1 + alm2
        for l in range(lmax + 1):
            for m in range(mmax + l):
                if (l <= lmax) and (m <= min(l, mmax)):
                    self.assertEqual(alm(l, m).real(), l*l + l)
                    self.assertEqual(alm(l, m).imag(), m*m + m)
                else:
                    self.assertEqual(alm(l, m).real(), 0)
                    self.assertEqual(alm(l, m).imag(), 0)

        alm = alm1 - alm2
        for l in range(lmax + 1):
            for m in range(mmax + l):
                if (l <= lmax) and (m <= min(l, mmax)):
                    self.assertEqual(alm(l, m).real(), l*l - l)
                    self.assertEqual(alm(l, m).imag(), m*m - m)
                else:
                    self.assertEqual(alm(l, m).real(), 0)
                    self.assertEqual(alm(l, m).imag(), 0)

        alm = alm1 * alm2
        for l in range(lmax + 1):
            for m in range(mmax + l):
                if (l <= lmax) and (m <= min(l, mmax)):
                    self.assertEqual(alm(l, m).real(), l*l*l - m*m*m)
                    self.assertEqual(alm(l, m).imag(), l*l*m + l*m*m)
                else:
                    self.assertEqual(alm(l, m).real(), 0)
                    self.assertEqual(alm(l, m).imag(), 0)

        alm = alm1 / alm2
        for l in range(lmax + 1):
            for m in range(mmax + l):
                if (l <= lmax) and (m <= min(l, mmax)) and (l != 0 or m != 0):
                    rhs = complex(l*l, m*m)/complex(l, m)
                else:
                    rhs = complex(0, 0)
                self.assertAlmostEqual(alm(l, m).real(), rhs.real)
                self.assertAlmostEqual(alm(l, m).imag(), rhs.imag)

        alm = 3.*alm1
        for l in range(lmax + 1):
            for m in range(mmax + l):
                if (l <= lmax) and (m <= min(l, mmax)):
                    self.assertEqual(alm(l, m).real(), 3*l*l)
                    self.assertEqual(alm(l, m).imag(), 3*m*m)
                else:
                    self.assertEqual(alm(l, m).real(), 0)
                    self.assertEqual(alm(l, m).imag(), 0)

        alm = alm1*3.
        for l in range(lmax + 1):
            for m in range(mmax + l):
                if (l <= lmax) and (m <= min(l, mmax)):
                    self.assertEqual(alm(l, m).real(), 3*l*l)
                    self.assertEqual(alm(l, m).imag(), 3*m*m)
                else:
                    self.assertEqual(alm(l, m).real(), 0)
                    self.assertEqual(alm(l, m).imag(), 0)

if __name__=="__main__":
   suite = unittest.TestLoader().loadTestsFromTestCase(TestHealAlm)
   unittest.TextTestRunner(verbosity=2).run(suite)
