#!/bin/zsh

for valency in Nb; do
    mkdir $valency
    cd $valency
    echo $valency

    # Set Ab chains
    if [[ $valency == "Ab" ]]; then
        Abchains="100"
    else
        Abchains="0"
    fi

    TAD_seq="FGGGMGGGGNGG"
    CoA_seq="ASSSSSSSSSSSSSSSSB"

    ### INTERACTION STRENGTH ###
    for AM in "AM_-1" "AM_-110" "AM_-160" "AM_-210"; do
        mkdir $AM
        cd $AM
        AM_int=$(echo $AM | cut -d '_' -f 2)

        ### INTERACTION STRENGTH ###
        for BM in "BM_-1" "BM_-110" "BM_-160" "BM_-210"; do
            BM_int=$(echo $BM | cut -d '_' -f 2)

            if [[ $BM_int -ge $AM_int ]]; then
                mkdir $BM
                cd $BM
                echo "${AM}  ${BM}"

                for AN in "AN_-1" "AN_-110" "AN_-160" "AN_-210"; do
                    mkdir $AN
                    cd $AN
                    AN_int=$(echo $AN | cut -d '_' -f 2)

                    for BN in "BN_-1" "BN_-110" "BN_-160" "BN_-210"; do
                        mkdir $BN
                        cd $BN
                        BN_int=$(echo $BN | cut -d '_' -f 2)

                        for coa_conc in 5 10 20 40 80 160 250 500 1000 2000 4000; do
                            mkdir $coa_conc
                            cd $coa_conc          

                            for rep in 1 2 3; do
                                mkdir $rep
                                cd $rep

                                # Copy files to directory
                                cp /work/degriffith/pimms/fordyce_collab/multichain_sims/fixed_surface_set4_affinity_titration/doubleTAD_doubleCoA/full_keyfile.kf ./KEYFILE.kf
                                cp /work/degriffith/pimms/fordyce_collab/multichain_sims/fixed_surface_set4_affinity_titration/doubleTAD_doubleCoA/doubleTAD_doubleCoA_params.prm ./full_params.prm

                                cp "/work/degriffith/pimms/fordyce_collab/multichain_sims/fixed_surface_set3_TAD_CoA_arch/setup/freeze_files/full/${valency}_double_full_freezefile.pimms" ./freezefile.pimms
                                cp "/work/degriffith/pimms/fordyce_collab/multichain_sims/fixed_surface_set3_TAD_CoA_arch/equil/${valency}/double/restart.pimms" ./equil_restart.pimms

                                # Edit fields in keyfile
                                sed -i "s/<~Ab_chains~>/${Abchains}/g" KEYFILE.kf
                                sed -i "s/<~CoA_chains~>/${coa_conc}/g" KEYFILE.kf
                                sed -i "s/<~TAD_seq~>/${TAD_seq}/g" KEYFILE.kf
                                sed -i "s/<~CoA_seq~>/${CoA_seq}/g" KEYFILE.kf
                                sed -i "s/<~AM_affinity~>/${AM_int}/g" full_params.prm
                                sed -i "s/<~BM_affinity~>/${BM_int}/g" full_params.prm
                                sed -i "s/<~AN_affinity~>/${AN_int}/g" full_params.prm
                                sed -i "s/<~BN_affinity~>/${BN_int}/g" full_params.prm

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