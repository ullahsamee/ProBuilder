#########################################################################
# File Name: example_monomer_scaffold_generation.sh
# Author: longxing
# mail: caolongxing@westlake.edu.cn
# Created Time: 2025年07月27日 星期日 20时24分00秒
#########################################################################
#!/bin/bash

# -frag_path: the fragment library. We used 7mer for fragment insertion.
# -rpx_db: rpx score database, the spatial resolution is 2.0 for cart and 26.0 for angle.
# -gzip: dump the pdb in gzip format to save dist space
# -dump_polyG: dump backbone structure as poly glycine
# -use_ss: use the secondary structure information for rpx scoring
# -hallucinate: randomly generate secondary structure configuration based on the number of helix
# -nstruct: the number of structures to genrate
# -total_score_cutoff: the total score cutoff for output, the score is the (rpx score + clash score) / protein_length
# -sidechain_neighbor_cutoff: sidechain neighbor cutoff for output
# -num_helix: the number of helix for output

# the fragment database and rpx score database are too large to be deposited 

../build/scaffold_generator -len 88 -num_repeats 1 -symmetry C1 -frag_path ../database/awesome_7HL_frags_20220427.bin -rpx_db ../database/rpx_cart2.0_ang26.0_ss_ALLAA.phmap -rpx_cart_resl 2.0 -rpx_ang_resl 26.0 -gzip -dump_polyG -use_ss -hallucinate -nstruct 100 -total_score_cutoff -15 -sidechain_neighbor_cutoff 3.0 -num_helix 6
