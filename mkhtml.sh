#!/bin/sh
root -l <<EOF
gSystem->Load("libRootHealPix")
THtml html
html.SetProductName("RootHealPix")
html.MakeAll()
EOF