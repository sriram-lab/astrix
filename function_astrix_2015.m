function [TRN,correl,rmsdval,betas,predictors] =  function_astrix_2015(data,probeids,tfprobes,N,z,RMSDTHRESHOLD)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function can be run in two modes; when z = 1, inference is run for the entire
% data set.
%when z = 0 the script performs N fold cross validation on the
% data with N specified by the user.
%Data - expression data - gene expression values for all genes ( both TFs and non-TFs)
% probeids - ids for the data. could be numbers from 1 to length of data
% tfprobes - the ids of the TFs in data
% output
% correl - this structure has the train and test set correlations
% rmsdval has rmsd values (Root mean square deviation)
% test_expression - is a 3d structure (N x no. of to-be-predicted variables
% x size of testset expression)
% betas has the coefficients of regression - dimension -(N x no. of to-be-predicted variables
% x 10); 10 is the max. no of predictors allowed for each variable
%predictors - same dimension as betas; it gives the position of the chosen
%predictors for that to-be-predicted variable
% predictors has the position of the top predictor (TF) for each gene in
% the data set. the top tfs are sorted based on
% MI with the target gene
% TRN - THIS STRUCTURE HAS THE TRN.   it has the list of tfs and the
% corresponding targets, and the RMSD value and the predictive correlation.
%
% RMSDTHRESHOLD - this is the threshold used to determine the genes
% selected to be in the final TRN. genes that can be predicted accurately
% below the RMSDTHRESHOLD will be chosen to be part of the final TRN.
% default value = 0.5

% NOTE: change the system call to JAVA based on your OS
% to DOS OR UNIX version.  also change the address to java based on its
% location in your computer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('z','var') || isempty(z)
    z = 1;
end

if (z == 1)
    N = 1;
    disp('running inference for the whole data')
end

if ~exist('RMSDTHRESHOLD','var') || isempty(RMSDTHRESHOLD)
    RMSDTHRESHOLD = 0.5;
end

TRN = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
translabel = probeids;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalizing and standardizing data
data = quantilenorm(data);
data = zscore(data,[],2);
% temp = (data - repmat(ignoreNaN(double(data),@mean,2),[1,size(data,2)]));
% kk = repmat(ignoreNaN(double(data),@var,2),[1,size(data,2)]);
% temp = temp./kk;
% data = single(temp);
% clear temp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = size(data,2);
idx = randsample(n,n);
fold= 0:ceil(n/N):n;
fold(end+1)=n;
m = size(data,1);

[ix, idsec] = ismember(tfprobes,probeids);
prot = data(idsec,:);  %TF expression
data1 = data;
data_probes = translabel;



for times = 1:N,
    
    if z == 1
        trainl = idx((fold(times)+1):fold(times+1));
        traintrans = data1(:,trainl);
        trainprot = prot(:,trainl);
        traintrans = knnimpute(traintrans);
        trainprot = knnimpute(trainprot);
        
    else
        testid = idx((fold(times)+1):fold(times+1));
        testtrans = data1(:,testid);
        testprot = prot(:,testid);
        trainl = idx(~ismember(idx,testid));
        traintrans = data1(:,trainl);
        trainprot = prot(:,trainl);
        
        traintrans = knnimpute(traintrans);
        testtrans = knnimpute(testtrans);
        trainprot = knnimpute(trainprot);
        testprot = knnimpute(testprot);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % splitting data for cluster
    li = 0:ceil(length(traintrans)/1):length(traintrans);
    li(end+1) = length(traintrans);
    split = 1;
    indi = li(split)+ 1:li(split+1);
    st = clock;
    sd = [date,'___',num2str(st(4))];
    
    % preallocate
    %memory
    if (times == 1)
        if (z ~= 1)
            
            predictorpos = single(zeros(N,length(indi),10));
            res12 = repmat(NaN,[length(indi),N]);
            res_train12 = repmat(NaN,[length(indi),N]);
            corrs12 = repmat(NaN,[length(indi),N]);
            corrs12_test = repmat(NaN,[length(indi),N]);
            test_expression = repmat(NaN,[N,length(indi),length(testprot(1,:))]);
            b_opt12 = repmat(NaN,[N,length(indi),10]);
        else
            predictorpos = single(zeros(length(indi),10));
            res_train12 = repmat(NaN,[length(indi),N]);
            corrs12 = repmat(NaN,[length(indi),N]);
            b_opt12 = repmat(NaN,[length(indi),10]);
        end
    end
    trainprot1 = trainprot;
    for i = 1:10%length(indi)
        trainprot = trainprot1;
        disp('GENE NO:');
        disp(i)
        clear spredrna smiscore miscore predrna X y
        y = traintrans(indi(i),:);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if any(ismember(tfprobes,data_probes(indi(i)))) % if its a tf remove it from the predictor list
            p1 = find(ismember(tfprobes,data_probes(indi(i))));
            trainprot(p1,:) = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % formatting input to ARACNE
        
        filenamess = ['nonsecgene',sd,num2str(split),'.txt'];
        fid = fopen(filenamess,'w');
        for j = 1:length(y),
            header{j} = ['sample',num2str(j)];
        end
        colnames = ['ID',header];
        
        
        fprintf(fid,'%s\t',colnames{:});
        towrite = [y;trainprot];
        towrite = [[1:size(towrite,1)]',towrite];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% calling aracne
        dlmwrite(filenamess,towrite,'delimiter','\t','-append','roffset',1);
        %arac_call = [ 'unix(''/usr/java/jdk1.6.0_13/bin/java -jar -Xms1024m -Xmx1024m ARACNE-java.jar -i nonsecgene',sd,num2str(split),'.txt -o singleout' ,sd,num2str(split),'.txt -h 1 -p 1e-6 -e 0.1'')'];
        if ispc
            arac_call = [ 'dos(''C:\Java\jre1.8.0_121\bin\java -jar -Xms1024m -Xmx1024m ARACNE-java.jar -i nonsecgene',sd,num2str(split),'.txt -o singleout' sd,num2str(split),'.txt -h 1 -p 1e-6 -e 0.1 '')'];
        elseif ismac
            arac_call = ['system(''unset DYLD_FRAMEWORK_PATH DYLD_LIBRARY_PATH; java -jar  -Xms1024m -Xmx1024m ARACNE-java.jar -i nonsecgene',sd,'.txt -o singleout' sd,'.txt -h 1 -p 1e-6 -e 0.1'')'];
        end
        eval(arac_call)
        %%%%% reading aracne output %%%%%%
        fclose(fid);
        filenamess1 = ['singleout',sd,num2str(split),'.txt'];
        fid = fopen(filenamess1);
        aracneoutfile = textscan(fid,'%n','delimiter','\t','commentstyle','>');
        fclose(fid);
        tempoutfile = aracneoutfile{:};
        tempoutfile = tempoutfile(2:end);
        temp = tempoutfile(tempoutfile < 2);
        if ~isempty(temp)
            miscore = temp;
            predrna = tempoutfile(tempoutfile >= 2);
            
            [smiscore,yy] = sort(miscore,1,'descend');
            
            for j = 1:size(miscore,2)
                spredrna(:,j) = predrna(yy(:,j),j);
            end
            
            boud = min(length(spredrna),10);
            idd = spredrna(1:boud);
            idd = idd - 1;
            X = trainprot(idd,:);
            
            
            if (z ~= 1)
                [n p] = size(X);
                ytest = testtrans(indi(i),:);
                Xtest = testprot(idd,:);
                predictorpos(times,i,1:length(idd)) = idd;
            else
                predictorpos(i,1:length(idd)) = idd;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %% LARS
            try [s_opt1(i,:), b, res_mean1, res_std1] = crossvalidate(@lars, 10, 1000, X', y');
                temp = b*X;
                res_train12(i,times) = sqrt(sum((y - temp).^2)/length(y));
                corrs12(i,times)  = corr2(temp(:),y(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if (z ~= 1)
                    res12(i,times) = sqrt(sum((ytest - b*Xtest).^2)/length(ytest));
                    temp1 = b*Xtest;
                    corrs12_test(i,times) = corr2(temp1(:),ytest(:));
                    test_expression(times,i,1:length(temp1)) = temp1(:);
                    b_opt12(times,i,1:length(b)) = b; clear b;
                else
                    b_opt12(i,1:length(b)) = b; clear b;
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            catch MM
                res12(i,times) = NaN;
                res_train12(i,times) = NaN;
                corrs12(i,times) = NaN;
                corrs12_test(i,times) = NaN;
                
            end
        else
            res12(i,times) = NaN;
            res_train12(i,times) = NaN;
            corrs12(i,times) = NaN;
            corrs12_test(i,times) = NaN;
            
            
        end
        
        
    end
    
end

if (z ~= 1)
    %tempname = ['save NI_alltfs_astrix_',sd,'_',num2str(split),' res12 res_train12 b_opt12 predictorpos corrs12 corrs12_test data_probes'];
    correl = struct('train_corr',corrs12,'test_corr',corrs12_test);
    rmsdval = struct('train_rmsd',res_train12,'test_rmsd',res12);
else
    correl = struct('train_corr',corrs12);
    rmsdval = struct('train_rmsd',res_train12);
end

betas =  b_opt12;
predictors = predictorpos;
if z == 1
    rmsdx = (rmsdval.train_rmsd < RMSDTHRESHOLD);sum(rmsdx)  % 
    c = 1;
    for i = 1:length(data_probes)
        if (rmsdx(i))
            lx = (abs(betas(i,:)) > 0.1);
            p1 = predictors(i,:); p1 = p1(lx);
            for j = 1:length(p1)
                TRN.targets(c,1) = data_probes(i);
                TRN.tfs(c,1) = tfprobes(p1(j));
                TRN.rmsd(c,1) = rmsdval.train_rmsd(i);
                TRN.corr(c,1) = correl.train_corr(i);
                c = c + 1;
            end
        end
    end
    

end

end
