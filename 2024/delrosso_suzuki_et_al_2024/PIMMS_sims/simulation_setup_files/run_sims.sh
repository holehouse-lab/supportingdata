#!/bin/zsh

for valency in Nb; do
    cd $valency

    ### INTERACTION STRENGTH ###
    for AM in "AM_-1" "AM_-110" "AM_-160" "AM_-210"; do
        cd $AM
        AM_int=$(echo $AM | cut -d '_' -f 2)

        ### INTERACTION STRENGTH ###
        for BM in "BM_-1" "BM_-110" "BM_-160" "BM_-210"; do
            BM_int=$(echo $BM | cut -d '_' -f 2)

            if [[ $BM_int -ge $AM_int ]]; then
                cd $BM

                for AN in "AN_-1" "AN_-110" "AN_-160" "AN_-210"; do
                    cd $AN
                    AN_int=$(echo $AN | cut -d '_' -f 2)

                    for BN in "BN_-1" "BN_-110" "BN_-160" "BN_-210"; do
                        cd $BN
                        BN_int=$(echo $BN | cut -d '_' -f 2)

                        for coa_conc in 5 10 20 40 80 160 250 500 1000 2000 4000; do
                            cd $coa_conc          

                            for rep in 1 2 3; do
                                cd $rep

                                lsf_run_pimms_beta "doubleTADxdoubleCoA${coa_conc}_AM${AM_int}_BM${BM_int}_AN${AN_int}_BN${BN_int}_${rep}" low
                                
                                cd ..
                            done
                            cd ..
                        done
                        cd ..
                    done
                    cd ..
                done
                cd ..
            fi
        done
        cd ..
    done 
    cd ..
done