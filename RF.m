%run plethora of tests
clc
close all

%compile everything
if strcmpi(computer,'PCWIN') |strcmpi(computer,'PCWIN64')
   compile_windows
else
   compile_linux
end

total_train_time=0;
total_test_time=0;

%diabetes

 
%modify so that training data is NxD and labels are Nx1, where N=#of
%examples, D=# of features

X = data;
Y = time;

[N D] =size(X);
%randomly split into 400 examples for training and 42 for testing
% rand('seed',2)
randvector = randperm(N);
trn_rate=0.75;%训练集的百分比，就是用70%的数据用于训练，其他的用于验证
X_trn = X(randvector(1:ceil(N*trn_rate)),:);
Y_trn = Y(randvector(1:ceil(N*trn_rate)));
X_tst = X(randvector(ceil(N*trn_rate)+1:end),:);
Y_tst = Y(randvector(ceil(N*trn_rate)+1:end));
   clear extra_options
    extra_options.importance = 1; %(0 = (Default) Don't, 1=calculate)
   
    model = regRF_train(X_trn,Y_trn, 100, 4, extra_options);
    Y_hat = regRF_predict(X_tst,model);
    fprintf('\nexample 8: MSE rate %f\n',   sum((Y_hat-Y_tst).^2));
    

    %model will have 3 variables for importance importanceSD and localImp
    %importance = a matrix with nclass + 2 (for classification) or two (for regression) columns.
    %           For classification, the first nclass columns are the class-specific measures
    %           computed as mean decrease in accuracy. The nclass + 1st column is the
    %           mean decrease in accuracy over all classes. The last column is the mean decrease
    %           in Gini index. For Regression, the first column is the mean decrease in
    %           accuracy and the second the mean decrease in MSE. If importance=FALSE,
    %           the last measure is still returned as a vector.
    figure('Name','Importance Plots')
    subplot(3,1,1);
    bar(model.importance(:,end-1));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Accuracy');
    
    subplot(3,1,2);
    bar(model.importance(:,end));xlabel('feature');ylabel('magnitude');
    title('Mean decrease in Gini index');
    
    
    %importanceSD = The ?standard errors? of the permutation-based importance measure. For classification,
    %           a D by nclass + 1 matrix corresponding to the first nclass + 1
    %           columns of the importance matrix. For regression, a length p vector.
    model.importanceSD
    subplot(3,1,3);
    bar(model.importanceSD);xlabel('feature');ylabel('magnitude');
    title('Std. errors of importance measure');
R2 = 1 - norm(Y_tst-Y_hat)^2/norm(Y_tst - mean(Y_tst))^2
%%
  

    n=0;
     bardata=nan(fenzu,size(data,2));
    for i=1:fenzu
    m=n+1;
    n=m+size(model.importance(find(idx==i)),1)-1;  
   bardata(i,m:n)=sort(model.importance(find(idx==i),end-1),'descend');
  end
    
    
      figure('Name','Importance Plots')
%     subplot(3,1,1);
    hold on
   
    
    for i=1:fenzu
    bar(bardata(i,:));
    xlabel('feature');
    ylabel('magnitude');
    colormap summer
    end
    title('Mean decrease in Accuracy');
        figure
    violin(bardata','facecolor','b')
    
    %%
       n=0;
     bardata2=nan(fenzu,size(data,2));
    for i=1:fenzu
    m=n+1;
    n=m+size(model.importance(find(idx==i)),1)-1;  
   bardata2(i,m:n)=sort(model.importance(find(idx==i),end),'descend');
  end
    
    
      figure('Name','Importance Plots')
%     subplot(3,1,1);
    hold on
   
    
    for i=1:fenzu
    bar(bardata2(i,:));
    xlabel('feature');
    ylabel('magnitude');
    colormap summer
    end
   title('Mean decrease in Gini index');;
    
        figure
    violin(bardata2','facecolor','b')
%%
    
%       figure('Name','Importance Plots')
% %     subplot(3,1,1);
%     hold on
%    
%     
%     for i=1:fenzu
%     hist(model.importance(find(idx==i),end-1));
%     xlabel('feature');
%     ylabel('magnitude');
%     colormap summer
%     end
%     title('Mean decrease in Accuracy');
    %% 
