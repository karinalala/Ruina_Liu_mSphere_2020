    clear
    clc
    load('E:\user\Desktop\data\data_origin.mat')
    filepath=uigetdir('C:\','choose destination folder ');
 
    clear data_filter1
    j = 1;
    for i = 1:size(data_origin,2)
        if size(find(data_origin(:,i)>5),1)>=24 | sum(data_origin(:,i))>=20
            selleteddata1(:,j)=data_origin(:,i);
            selletedotu1(:,j)=otu_origin(:,i);
            selletedtaxonomy1(:,j)=taxonomy(:,i);
            j=j+1;
        end
    end
    clear i j m n

    n0=arrayfun(@(x) x==0,selleteddata1);
    [t,b,c]=unique(time);
    for i=1:size(t,1)
       t0(i,:)=sum(n0(find(c==i),:))/24>0.6;
    end

    selleteddata3=selleteddata1(:,find(sum(t0)~=size(t,1)));
    selletedotu3=selletedotu1(:,find(sum(t0)~=size(t,1)));
    selletedtaxonomy3=selletedtaxonomy1(:,find(sum(t0)~=size(t,1)));
    
    [m,n]=size(selleteddata3);
    mx=mean(selleteddata3);
    stdx=std(selleteddata3);
    data=(selleteddata3-mx(ones(m,1),:))./stdx(ones(m,1),:);

    otu=selletedotu3;
    clear i m n mx stdx  j n0 p  selleteddata1   selleteddata2   selletedotu1 selletedotu2;
    tax=selletedtaxonomy3;

       for i=1:size(data,2)
        p(:,i) = kruskalwallis(data(:,i),time,'off');
    end
    KWresult=p<0.001;

    
    for i=1:size(t,1)
        med(i,:)= median(data(find(c==i),:));
    end

  
    for i=1:size(t,1)
        mea(i,:)= mean(data(find(c==i),:));
    end
   
    for m=1:size(data,2)
        for i=1:size(t,1)-1
            for j=i+1:size(t,1)
                [~,h]=ranksum(data(find(c==i),m),data(find(c==j),m));%U-test
                location(i,j,m)=h*sign(med(j,m)-med(i,m));
            end
        end
    end

    location2=zeros(size(t,1),size(data,2));
    for m=1:size(data,2)
        for i=1:size(t,1)-1
            j=i+1;
            h=0;
            while h~=1 & j<=10
                [~,h]=ranksum(data(find(c==i),m),data(find(c==j),m));
                j=j+1;
            end
            j=j-1;
            location2(j,m)=h*sign(med(j,m)-med(i,m));
        end
    end
    clear h  i  j  location  c  b m a ans  i
   
    for i=1:size(location2,1)
        qushi(i,:)=sum(location2(1:i,:),1); 
    end
    qushi2=qushi;
    

    for i=1:size(qushi2,2)
      if sum(abs(location2(:,i)))>0
     qushi2(end,i)=qushi2(find((qushi2(1:end,i))~=0,1,'last'),i);
      end
    a=find((location2(2:end-1,i))==0);
    b=find((location2(2:end-1,i))~=0);
    qushi2(a+1,i)=interp1(t([1;b+1;end]),qushi2([1;b+1;end],i),t(a+1),'Linear');
    clear a b
    end

    clear  i  j
    
    x=qushi2;l=[-6,4];
    
    for fenzu=7
     treeCluster=linkage(squareform(pdist(qushi','cityblock')),'ward');
    close all
    

figure
[H,~,O]=dendrogram(treeCluster,0,'Orientation', 'left','ColorThreshold',1600);

h= gca;
set(h,'YTickLabel', []);
set(h,'YTick', []);
set(H,'Linewidth',2);

% 
% for i=1:1810
% O2(i)=idx(O(i));
% end
 h = figure();	
    hold on
    SIZE = get(0,'ScreenSize');			
    set(h, 'position', SIZE);	
    for i=1:fenzu
        m=find(idx==i & KWresult'==1);
        n=find(idx==i & KWresult'==0);
        hold on
        subplot(ceil(fenzu/4),4,i)
        set(gca,'YLim',l)
          set(gca,'XLim',[0 360])
          set(gca,'XLim',[0 360])
            set(gca,'fontsize',[14])
            set(gca,'fontsize',[14])
            set(gca,'linewidth',1.5)
            set(gca,'linewidth',1.5)
            xlabel('PMI (h)','FontSize',16)
            ylabel('Relative trend','FontSize',16)
          box on
                      hold on
            plot(t,x(:,m),'Color',[213 62 79]/213,'linewidth',1.5) 
            if n~=-1
            plot(t,x(:,n),'Color',[69 134 198]/255,'linewidth',1.5) 
            end
        hold on
        scatter(t,mean(x(:,m),2),'fill','sizedata',[100],'markeredgeColor',[50 50 50]/255,...
     'markerfaceColor',[50 50 50]/255)
        plot(t,mean(x(:,m),2),'Color',[50 50 50]/255,'LineWidth',3)
        title(['Class ' num2str(i) ' (n=' num2str(size(find(idx==i),1)) ')'])
    end
    clear h  i j a  b  c m n 
    %pause
    end
  
PLS=pls(data,time,20,'center');
PLS.VIP;
PLS.SR;
PLS.tploadings;

% uninformative variable elimination(UVE)
  UVE=mcuvepls(data,time,15,'center');
  UVE.RI;
% Competitive Adaptive Reweighted Sampling£¨CARS£©	
CARS=carspls(data,time,20,10,'center',50); 
CARS.vsel;
CARSvesl=zeros(1810,1);
CARSvesl(CARS.vsel)=1;

% Random Frog
Frog=randomfrog_pls(data,time,10,'center',2000,2);
Frog.probability;
% spearman 
spearman=corr(data,time,'type','spearman');
   


pathout=[filepath '\' num2str(fenzu) 'cluster.xlsx']
xlswrite(pathout,{'otu','annotation','class','KW','VIP','SR','tploadings','UVE','CARS','Random Frog','Spearman'},1,'A1');
xlswrite(pathout,otu',1,'A2');
xlswrite(pathout,tax',1,'B2');
xlswrite(pathout,idx,1,'C2');
xlswrite(pathout,KWresult',1,'D2');


xlswrite(pathout,PLS.VIP',1,'E2');
xlswrite(pathout,PLS.SR',1,'F2');
xlswrite(pathout,PLS.tploadings,1,'G2');
xlswrite(pathout,UVE.RI',1,'H2');
xlswrite(pathout,CARSvesl,1,'I2');
xlswrite(pathout,Frog.probability',1,'J2');
xlswrite(pathout,spearman,1,'K2');

[a1,~,c1]=unique(name1810(:,3));
[a2,~,c2]=unique(name1810(:,7));
[a3,~,c3]=unique(time);
OTUSUM3=[];

for    j=1:size(a3,1)
       for    i=1:size(a1,1)
            OTUSUM3(i,j)=sum(sum(selleteddata3(find(c3==j),find(c1==i))));
       end 
       OTUSUM3(:,j)=OTUSUM3(:,j)/sum(OTUSUM3(:,j));
end
  pathout3=[filepath '\' num2str(fenzu) 'classes_p.xlsx']
xlswrite(pathout3,OTUSUM3,1,'B2');    
  xlswrite(pathout3,a1,1,'A2'); 
    xlswrite(pathout3,t',1,'B1'); 
    
   OTUSUM7=[];

 for j=1:size(a3,1)
       for i=1:size(a2,1)     
            OTUSUM7(i,j)=sum(sum(selleteddata3(find(c3==j),find(c2==i))));
       end 
            OTUSUM7(:,j)=OTUSUM7(:,j)/sum(OTUSUM7(:,j));
end
    
pathout7=[filepath '\' num2str(fenzu) 'classes_g.xlsx']
xlswrite(pathout7,OTUSUM7,1,'B2');    
  xlswrite(pathout7,a2,1,'A2'); 
    xlswrite(pathout7,t',1,'B1'); 
    
for k=1:fenzu

[a1,~,c1]=unique(name1810(find(idx==k),3));
[a2,~,c2]=unique(name1810(find(idx==k),7));
[a3,~,c3]=unique(time);
OTUSUM3=[];

for    j=1:size(a3,1)
       for    i=1:size(a1,1)
            OTUSUM3(i,j)=sum(sum(selleteddata3(find(c3==j),find(c1==i))));
       end 
       OTUSUM3(:,j)=OTUSUM3(:,j)/sum(OTUSUM3(:,j));
end
  pathout3=[filepath '\' num2str(fenzu) 'cluster_class_' num2str(k) '_p.xlsx']
xlswrite(pathout3,OTUSUM3,1,'B2');    
  xlswrite(pathout3,a1,1,'A2'); 
    xlswrite(pathout3,t',1,'B1'); 
    xlswrite(pathout3,{['class_' num2str(k)]},1,'A1');  
   OTUSUM7=[];

 for j=1:size(a3,1)
       for i=1:size(a2,1)     
            OTUSUM7(i,j)=sum(sum(selleteddata3(find(c3==j),find(c2==i))));
       end 
            OTUSUM7(:,j)=OTUSUM7(:,j)/sum(OTUSUM7(:,j));
end
    
pathout7=[filepath '\' num2str(fenzu) 'cluster_class_' num2str(k) '_g.xlsx']
xlswrite(pathout7,OTUSUM7,1,'B2');    
  xlswrite(pathout7,a2,1,'A2'); 
    xlswrite(pathout7,t',1,'B1'); 
      xlswrite(pathout7,{['class_' num2str(k)]},1,'A1');  
    
    
end
    