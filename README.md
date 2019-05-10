# MCDC

DCDB: contains drug pair ID.

DCDB_DrugCom:  contains (DCDB, drug1, drug2), drug1 and drug2 can be found in DrugBank.

OMIM_ID: disease ID in OMIM dataset.




data: contain the (drug pair, disease).

drug similarity: drugSimSmap, drugSimSmile, drugSimTargetSW.

disease similarity: HPOSim, SimMimMiner, DoSim.




instanceDrugKernel: constructs the pairwise kernel.

MultipleViewComplete: multi-view learning to complete the kernel.

combineKernel: linear combination of multiple kernels.

kronrls: kernel-based algorithm.

auc: measures AUC and AUPR.

Demo:  run main.m

