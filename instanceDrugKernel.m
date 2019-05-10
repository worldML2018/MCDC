function Sdrug = instanceDrugKernel(dcdb1,dcdb2,drugKernel);
    load DCDB_DrugCom.mat;
    drug1 = DCDB_DrugCom(dcdb1,:);
    drug2 = DCDB_DrugCom(dcdb2,:);
    Sdrug = drugKernel(drug1(1),drug2(1))+drugKernel(drug1(2),drug2(1))+ ...
        drugKernel(drug1(1),drug2(2))+drugKernel(drug1(2),drug2(2));
end
    
