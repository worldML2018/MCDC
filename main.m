% % the main function for the drug combination project
 clear all;
 clc;
%function [] = main(alphapara,betapara,ranks,filename,variablename)
load DCDB;
load OMIM_ID;

load drugSimSmile;
load drugSimSmap;
load drugSimTargetGO;
load drugSimTargetSW;


load DoSim;
load HPOSim;
load SimMimMiner;



[dcdb, omim, label] =  textread('data','%q %f %f');

dataMatrix = zeros(length(DCDB), length(OMIM_ID));

for i = 1:length(dcdb)
    index1 = find(strcmp(DCDB(:,1),dcdb(i,1)));
    index2 = find(OMIM_ID == omim(i,1));
    dataMatrix(index1, index2) = 1;
end
        

drugKernel{1} = drugSimSmile;
drugKernel{2} = drugSimSmap;
drugKernel{3} = drugSimTargetGO;
drugKernel{4} = drugSimTargetSW;


diseaseKernel{1} = DoSim;
diseaseKernel{2} = HPOSim;
diseaseKernel{3} = SimMimMiner;



% complete kernel
[drugKernel, objs1] = MultipleViewComplete(drugKernel, 0.01);
[diseaseKernel, objs2] = MultipleViewComplete(diseaseKernel, 0.01);

% combine kernel


 drugKernel = combineKernel(drugKernel,0.3,0.3,0.2,0.2);
 diseaseKernel = combineKernel(diseaseKernel,0.2,0.2,0.3);

% construct pairwise kernel
 pairKernel = zeros(length(DCDB),length(DCDB));
 for i = 1:length(DCDB)
     for j = 1:i
            pairKernel(i,j) = instanceDrugKernel(i, j,drugKernel);
            pairKernel(j,i) = pairKernel(i,j);
     end
 end
 



 



 
 % constrain training data and test data
 folds = 5; % 50% as test data
 
 positiveId = find(dataMatrix);
 crossval_id = crossvalind('Kfold',positiveId(:),folds);
 
 AUC = zeros(folds,1);
 AUPR = zeros(folds,1);
 lambda = 1;

 for  fold = 1:folds
     
     
     yy = dataMatrix;
     PtrainID = positiveId(find(crossval_id~=fold));
	 PtestID  = positiveId(find(crossval_id==fold));

     % sample equal amount of negative sample
     negativeID = find(yy==0);
	 num = numel(negativeID);
     Nidx = randperm(num);
     NtestID = negativeID(Nidx(1:length(PtestID)));
     
     % construct test datasize
     testID = [PtestID;PtestID];
     testlabel = [ones(length(PtestID),1); zeros(length(PtestID),1)];
     
     yy(testID) = 0;
     y2 = kronrls( pairKernel, diseaseKernel, yy, lambda);
     predict = y2(testID);
      
     [testroc, testpr, rocx, rocy, prx, pry] = auc( yy(:), y2(:));
     
     
      AUC(fold,1) = testroc;
      AUPR(fold,1) = testpr;
      
       fprintf('--------------- FOLD %d:  AUC %d, AUPR %d \n', fold, testroc, testpr)
     
    

end
    
     
    


