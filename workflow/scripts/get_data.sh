#!/usr/bin/bash

cat data/model_proteomes.tsv | while read row; do
     pdb=$(echo "$row" | cut -f5)
     fa=$(echo "$row" | cut -f6)
     
     if [ -f data/fastas/$(basename $fa) ]; then
          echo "already done"
          else
          wget $fa -P data/fastas/
     fi

     if [ -f data/pdbs/$(basename $pdb) ]; then
          echo "already done"
          else
          wget $pdb -P data/pdbs/
     fi
done
