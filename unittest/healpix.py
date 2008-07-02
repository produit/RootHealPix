#!/usr/bin/env python
# $Id: healpix.py,v 1.4 2008/07/02 22:40:59 oxon Exp $
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
        hpd1 = ROOT.THealPixD("hpd1", "", 3, ROOT.kFALSE) # RING
        hpd2 = ROOT.THealPixD("hpd2", "", 3, ROOT.kTRUE)  # NESTED
        hpd1.SetDegree()
        hpd2.SetDegree()
        for theta in range(180):
            for phi in range(360):
                hpd1.Fill(theta, phi)
                hpd2.Fill(theta, phi)
        
        hpd3 = hpd1.Rebin(4, "hpd3") # higher resolution
        hpd4 = hpd2.Rebin(4, "hpd4") # higher resolution
        hpd5 = hpd1.Rebin(2, "hpd5") # lower resolution
        hpd6 = hpd2.Rebin(2, "hpd6") # lower resolution

        for i in range(hpd1.GetNpix()):
            # Check RING 3->4 rebinning
            inest = hpd1.Ring2Nest(i)
            before = hpd1.GetBinContent(i)
            for j in range(4):
                after = hpd3.GetBinContent(hpd3.Nest2Ring(inest*4 + j))
                self.assertEqual(before, after*4.)

            # Check NESTED 3->4 rebinning
            before = hpd2.GetBinContent(i)
            for j in range(4):
                after = hpd4.GetBinContent(i*4 + j)
                self.assertEqual(before, after*4.)
        
        for i in range(hpd5.GetNpix()):
            # Check RING 3->2 rebinning
            inest = hpd5.Ring2Nest(i)
            after = hpd5.GetBinContent(i)
            before = 0
            for j in range(4):
                before += hpd1.GetBinContent(hpd1.Nest2Ring(inest*4 + j))
            self.assertEqual(before, after)

            # Check NESTED 3->2 rebinning
            after = hpd6.GetBinContent(i)
            before = 0
            for j in range(4):
                before += hpd2.GetBinContent(i*4 + j)
            self.assertEqual(before, after)

        hpd3.Rebin(3)
        hpd4.Rebin(3)

        for i in range(hpd3.GetNpix()):
            self.assertEqual(hpd3.GetBinContent(i), hpd4.GetBinContent(i))

        hpd5.Rebin(3)
        hpd6.Rebin(3)

        for i in range(hpd5.GetNpix()):
            self.assertEqual(hpd5.GetBinContent(i), hpd6.GetBinContent(i))

            
    def testWeight(self):
        hpd1 = ROOT.THealPixD("hpd1", "title", 0);
        for i in range(12):
            hpd1.SetBinContent(i, i)

        hpd1.Sumw2()
        # Check bin errors
        for i in range(12):
            self.assertEqual(hpd1.GetBinError(i), i**0.5)

        hpd2 = hpd1.Rebin(1, "hpd2")
        for i in range(hpd1.GetNpix()):
            e1 = hpd1.GetBinError(i)
            for j in range(4):
                e2 = hpd2.GetBinError(hpd2.Nest2Ring(hpd1.Ring2Nest(i)*4 + j))
                self.assertEqual(e1, e2*2.)

        hpd3 = hpd2.Rebin(0, "hpd3")
        for i in range(hpd3.GetNpix()):
            e3 = hpd3.GetBinError(i)
            for j in range(4):
                e2 = hpd2.GetBinError(hpd2.Nest2Ring(hpd3.Ring2Nest(i)*4 + j))
                self.assertEqual(e3, e2*2.)

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
