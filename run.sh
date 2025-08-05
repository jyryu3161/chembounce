python chembounce.py \
    -o ./output_test \
    -i "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl" \
    -n 100 -t 0.5 \
    --cand_max_n__rplc 10 

python chembounce.py \
    -o ./output_test2 \
    -i "CC1=C2C(=C(N(C1=O)C)NC3=C(C=C(C=C3)I)F)C(=O)N(C(=O)N2C4=CC=CC(=C4)NC(=O)C)C5CC5" \
    -n 100 -t 0.5 \
    --cand_max_n__rplc 10 
    --wo_lipinski

    