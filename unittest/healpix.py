#!/usr/bin/env python
# $Id: healpix.py,v 1.5 2008/07/03 07:56:55 oxon Exp $
# Author: Akira Okumura 2008/06/25

import unittest

import ROOT

class TestHealPix(unittest.TestCase):
    """
    Unit test class for THealPix familly.
    """
    def setUp(self):
        """
        Set up function. libRootHealPix must be loaded.
        """
        ROOT.gSystem.Load("libRootHealPix")

    def testGetSet(self):
        """
        Check basic getter/setter.
        """
        order = 1
        nside = 2**order
        name = ("hpf", "hpd")
        title = ("THealPixF", "THealPixD")
        hp = (ROOT.THealPixF(name[0], title[0], order, False),\
              ROOT.THealPixD(name[1], title[1], order, True))
        for i in range(2):
            self.assertEqual(hp[i].GetName(), name[i])
            self.assertEqual(hp[i].GetTitle(), title[i])
            self.assertEqual(hp[i].GetOrder(), order)
            self.assertEqual(hp[i].GetNside(), nside)
            self.assertEqual(hp[i].GetNpix(), 12*(nside**2))
            self.assertEqual(hp[i].GetDirectory(), ROOT.gROOT)
            if i==0:
                self.assertEqual(hp[i].IsNested(), False)
                self.assertEqual(hp[i].GetSchemeString(), "RING")
            else:
                self.assertEqual(hp[i].IsNested(), True)
                self.assertEqual(hp[i].GetSchemeString(), "NESTED")

        hp[0].SetUnit("MeV")
        self.assertEqual(hp[0].GetUnit(), "MeV")

    def testBinContent(self):
        """
        Check if bin contents are set correctly.
        """
        hpd = ROOT.THealPixD("hpd", "title", 0);
        for i in range(12):
            hpd.SetBinContent(i, i)

        # Check bin contents
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), i)

        self.assertEqual(hpd.GetEntries(), 12)

    def testRebin(self):
        """
        Compare bin contents between before/after of rebinning.
        """
        ring3 = ROOT.THealPixD("ring3", "", 3, ROOT.kFALSE) # RING
        nest3 = ROOT.THealPixD("nest3", "", 3, ROOT.kTRUE)  # NESTED
        ring3.SetDegree()
        nest3.SetDegree()
        for theta in range(180):
            for phi in range(360):
                ring3.Fill(theta, phi)
                nest3.Fill(theta, phi)
        
        ring4 = ring3.Rebin(4, "ring4") # higher resolution
        nest4 = nest3.Rebin(4, "nest4") # higher resolution
        ring2 = ring3.Rebin(2, "ring2") # lower resolution
        nest2 = nest3.Rebin(2, "nest2") # lower resolution

        for i in range(ring3.GetNpix()):
            # Check RING 3->4 rebinning
            inest = ring3.Ring2Nest(i)
            before = ring3.GetBinContent(i)
            for j in range(4):
                after = ring4.GetBinContent(ring4.Nest2Ring(inest*4 + j))
                self.assertEqual(before, after*4.)

        for i in range(nest3.GetNpix()):
            # Check NESTED 3->4 rebinning
            before = nest3.GetBinContent(i)
            for j in range(4):
                after = nest4.GetBinContent(i*4 + j)
                self.assertEqual(before, after*4.)

        for i in range(ring2.GetNpix()):
            # Check RING 3->2 rebinning
            inest = ring2.Ring2Nest(i)
            after = ring2.GetBinContent(i)
            before = 0.
            for j in range(4):
                before += ring3.GetBinContent(ring3.Nest2Ring(inest*4 + j))
            self.assertEqual(before, after)

        for i in range(nest2.GetNpix()):
            # Check NESTED 3->2 rebinning
            after = nest2.GetBinContent(i)
            before = 0.
            for j in range(4):
                before += nest3.GetBinContent(i*4 + j)
            self.assertEqual(before, after)

        ring4.Rebin(3)

        for i in range(ring4.GetNpix()):
            self.assertEqual(ring4.GetBinContent(i), ring3.GetBinContent(i))

        nest4.Rebin(3)

        for i in range(ring2.GetNpix()):
            self.assertEqual(nest4.GetBinContent(i), nest3.GetBinContent(i))

            
    def testWeight(self):
        ring3 = ROOT.THealPixD("ring3", "title", 3, ROOT.kFALSE);
        nest3 = ROOT.THealPixD("nest3", "title", 3, ROOT.kTRUE);
        for i in range(ring3.GetNpix()):
            ring3.SetBinContent(i, i)
            nest3.SetBinContent(i, i)

        ring3.Sumw2()
        nest3.Sumw2()
        # Check bin errors
        for i in range(ring3.GetNpix()):
            self.assertEqual(ring3.GetBinError(i), i**0.5)
            self.assertEqual(nest3.GetBinError(i), i**0.5)

        ring4 = ring3.Rebin(4, "ring4")
        for i in range(ring3.GetNpix()):
            e3 = ring3.GetBinError(i)
            for j in range(4):
                e4 = ring4.GetBinError(ring4.Nest2Ring(ring3.Ring2Nest(i)*4 + j))
                self.assertAlmostEqual(e3, e4*2.)

        nest4 = nest3.Rebin(4, "nest4")
        for i in range(nest3.GetNpix()):
            e3 = nest3.GetBinError(i)
            for j in range(4):
                e4 = nest4.GetBinError(i*4 + j)
                self.assertAlmostEqual(e3, e4*2.)

        ring2 = ring3.Rebin(2, "ring2")
        for i in range(ring2.GetNpix()):
            e2 = ring2.GetBinError(i)
            e3 = 0.
            for j in range(4):
                e3 += ring3.GetBinError(ring3.Nest2Ring(ring2.Ring2Nest(i)*4 + j))**2
            self.assertAlmostEqual(e3, e2**2)

        nest2 = nest3.Rebin(2, "nest2")
        for i in range(nest2.GetNpix()):
            e2 = nest2.GetBinError(i)
            e3 = 0.
            for j in range(4):
                e3 += nest3.GetBinError(i*4 + j)**2
            self.assertAlmostEqual(e3, e2**2)

    def testArithmeticOperations(self):
        hpd1 = ROOT.THealPixD("hpd1", "", 0)
        hpd2 = ROOT.THealPixD("hpd2", "", 0)
        for i in range(12):
            hpd1.SetBinContent(i, i)
            hpd2.SetBinContent(i, i*i)
        
        hpd = hpd1 + hpd2
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), i + i*i)

        hpd = hpd1 - hpd2
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), i - i*i)

        hpd = hpd1 * hpd2
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), i*i*i)

        hpd = hpd1 / hpd2
        for i in range(12):
            if i==0:
                self.assertEqual(hpd.GetBinContent(i), 0)
            else:
                self.assertEqual(hpd.GetBinContent(i), 1./i)

        hpd = 3.*hpd1
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), 3.*i)

        hpd = hpd1*3.
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), i*3.)
            
if __name__=="__main__":
   suite = unittest.TestLoader().loadTestsFromTestCase(TestHealPix)
   unittest.TextTestRunner(verbosity=2).run(suite)
