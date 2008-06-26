#!/bin/sh
root -l <<EOF
gSystem->Load("libRootHealPix")
THtml html
html.MakeAll()
EOF