#!/usr/bin/env python
# $Id: healpix.py,v 1.1 2008/06/26 00:52:18 oxon Exp $
# Author: Akira Okumura 2008/06/25

import unittest
from test import test_support

import ROOT

class TestHealPix(unittest.TestCase):
    def setUp(self):
        ROOT.gSystem.Load("libRootHealPix")

    def testGetSet(self):
        order = 1
        nside = 2**order
        name = ("hpf", "hpd")
        title = ("THealPixF", "THealPixD")
        hp = (ROOT.THealPixF(name[0], title[0], order),\
              ROOT.THealPixD(name[1], title[1], order))
        for i in range(2):
            self.assertEqual(hp[i].GetName(), name[i])
            self.assertEqual(hp[i].GetTitle(), title[i])
            self.assertEqual(hp[i].GetOrder(), order)
            self.assertEqual(hp[i].GetNside(), nside)
            self.assertEqual(hp[i].GetNpix(), 12*(nside**2))
            self.assertEqual(hp[i].IsNested(), False)
            self.assertEqual(hp[i].GetDirectory(), ROOT.gROOT)

    def testBin(self):
        hpd = ROOT.THealPixD("hpd", "title", 0);
        for i in range(12):
            hpd.SetBinContent(i, i)

        # Check bin contents
        for i in range(12):
            self.assertEqual(hpd.GetBinContent(i), i)

        self.assertEqual(hpd.GetEntries(), 12)

    def testWeight(self):
        hpd = ROOT.THealPixD("hpd", "title", 0);
        for i in range(12):
            hpd.SetBinContent(i, i)

        hpd.Sumw2()
        # Check bin errors
        for i in range(12):
            self.assertEqual(hpd.GetBinError(i), i**0.5)

    def testOperations(self):
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
            
def main():
    test_support.run_unittest(TestHealPix)

if __name__=="__main__":
    main()
