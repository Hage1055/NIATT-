clc
clear all
close all
% %%
% %we begin by the 2D case
% %construct a rectangular gridding
% N=5000 %number of points or steps on the X and Y directions
% 
% tx=-N/2:1:N/2 ;%gridding
% ty=tx;
% %generate random poiints with coordiantes in the specified intervall
% C=-N/2+randi([0,N],3,N); 
% Px=C(1,:);
% Py=C(2,:);
% Pz=C(3,:);
% scatter(Px,Py,3,'filled');
% hold on
%%
%%
%get the data from the las file
% fid=fopen('20171215_flight_5_edit.las'); % reading las file 
% if fid == -1
%     error('Error opening file')
% end
% fseek(fid, 96, 'bof');
% OffsetToPointData = fread(fid,1,'uint16');
% fseek(fid, 131, 'bof');
% XScaleFactor = fread(fid,1,'double');
% YScaleFactor = fread(fid,1,'double');
% ZScaleFactor = fread(fid,1,'double');
% XOffset = fread(fid,1,'double');
% YOffset = fread(fid,1,'double');
% ZOffset = fread(fid,1,'double');
% 
% % The number of bytes from the beginning of the file to the first point record
% % data field is used to access the attributes of the point data
% %
% c = OffsetToPointData;
% 
% % Read in the X coordinates of the points
% %
% % Reads in the X coordinates of the points making use of the 
% % XScaleFactor and XOffset values in the header.
% fseek(fid, c, 'bof');
% X1=fread(fid,inf,'int16',24);
% X=X1*XScaleFactor+XOffset;
% Xt=unique(X);
% 
% % Read in the Y coordinates of the points
% fseek(fid, c+4, 'bof');
% Y1=fread(fid,inf,'int16',24);
% Y=Y1*YScaleFactor+YOffset;
% Yt=unique(Y);
% 
% % Read in the Z coordinates of the points
% fseek(fid, c+8, 'bof');
% Z1=fread(fid,inf,'int16',24);
% Z=Z1*ZScaleFactor+ZOffset;
% Zt=unique(Z);
% 
% figure;
% pcshow([X(1:size(Z)),Y(1:size(Z)),Z(:)]);
% title('Point Cloud Visualization');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% hold on;
%%
%open the file and load the data (2)
c = lasdata([cd,'\LIDAR_Data\LAS_files\20171215_flight_6_edit.las']);
c
X = c.x;
Y = c.y;
Z = c.z;
[X,Is0]=sort(X);
Y=Y(Is0);
Z=Z(Is0);
figure;
scatter3(X,Y,Z,'.');
c = lasdata([cd,'\vehicle_pointcloud.las']);
c
Xr = c.x;
Yr = c.y;
Zr = c.z;
[Xr,Is]=sort(Xr)
Yr=Yr(Is)
Zr=Zr(Is)
figure;
scatter3(Xr,Yr,Zr,'.');
%figure;
%scatter(X,Y);
%extracting a region
% a = 0.5;
% X = X( (X > min(Xr)-a & X < max(Xr)+a ) & (Y > min(Yr)-a & Y < max(Yr)+a) );
% Y = Y( (X > min(Xr)-a & X < max(Xr)+a ) & (Y > min(Yr)-a & Y < max(Yr)+a) );
% Z = Z( (X > min(Xr)-a & X < max(Xr)+a ) & (Z > min(Zr)-a & Z < max(Zr)+a) );
%%
%get the vehicle point cloud
% fid=fopen('vehicle_pointcloud.las'); % reading las file 
% if fid == -1
%     error('Error opening file')
% end
% fseek(fid, 96, 'bof');
% OffsetToPointData = fread(fid,1,'uint16');
% fseek(fid, 131, 'bof');
% XScaleFactor = fread(fid,1,'double');
% YScaleFactor = fread(fid,1,'double');
% ZScaleFactor = fread(fid,1,'double');
% XOffset = fread(fid,1,'double');
% YOffset = fread(fid,1,'double');
% ZOffset = fread(fid,1,'double');
% 
% % The number of bytes from the beginning of the file to the first point record
% % data field is used to access the attributes of the point data
% %
% c = OffsetToPointData;
% 
% % Read in the X coordinates of the points
% %
% % Reads in the X coordinates of the points making use of the 
% % XScaleFactor and XOffset values in the header.
% fseek(fid, c, 'bof');
% X1r=fread(fid,inf,'int16',24);
% Xr=X1r*XScaleFactor+XOffset;
% Xtr=unique(Xr);
% 
% % Read in the Y coordinates of the points
% fseek(fid, c+4, 'bof');
% Y1r=fread(fid,inf,'int16',24);
% Yr=Y1r*YScaleFactor+YOffset;
% Ytr=unique(Yr);
% 
% % Read in the Z coordinates of the points
% fseek(fid, c+8, 'bof');
% Z1r=fread(fid,inf,'int16',24);
% Zr=Z1r*ZScaleFactor+ZOffset;
% Ztr=unique(Zr);
% 
% %figure;
% pcshow([Xr(1:size(Xr)-1),Yr(1:size(Yr)-1),Zr(:)],'red');
% title('Point Cloud Visualization');
% xlabel('Xr');
% ylabel('Yr');
% zlabel('Zr');
% 
% s=randi([0 2400],1,1)  %window size
%%
%gridding
%delta_r=min(abs(Xtr(1)-Xtr(2)),abs(Ytr(1)-Ytr(2)));
%compute and display some statistics
disp('compute and display some statistics')
disp('average X= ');
Xavg=sum(X)/length(X)
disp('var of X')
Vx=var(X)
disp('min absolute value of X');
min(abs(X))
disp('max absolute value of X')
max(abs(X))
disp('average Y= ');
Yavg=sum(Y)/length(Y)
disp('var of Y')
Vy=var(Y)
disp('min absolute value of Y');
min(abs(Y))
disp('max absolute value of Y')
max(abs(Y))
%for the filter
disp('statistics for the filter')
disp('average Xr= ');
Xavg=sum(Xr)/length(Xr)
disp('var of Xr')
Vxr=var(Xr)
disp('min absolute value of Xr');
min(abs(Xr))
disp('max absolute value of Xr')
max(abs(Xr))
disp('average Yr= ');
Yravg=sum(Yr)/length(Yr)
disp('var of Yr')
Vye=var(Yr)
disp('min absolute value of Yr');
min(abs(Yr))
disp('max absolute value of Yr')
max(abs(Yr))


%%
s=1000;
% the reference coordinates intervall
% S=randi([-N/2+s N/2-s],1,2);
% Del_x=[S(1) S(1)+s]; %x-compoment of the window
% Del_y=[S(2) S(2)+s]; %y-component of the window
% a=1;
% X=X(X>=min(Xr)-a & X<=max(Xr)+a);
% 
% I1=find(X>=min(Xr)-a & X<=max(Xr)+a); %find the indeces of the points belonging to the reference cloud
% Y=Y(I1);
% I2=find(Y>min(Yr)-a & Y<max(Yr)+a);
% I=intersect(I1,I2);
% X=X(I);
% Y=Y(I);
% %Pz(I)=Pz(I)+500; %to make it a significant cloud of a particular distribution of heights
% Z=Z(I); 
%Constructing the matrix
% Xs=unique(Px);   %Xs contains a sorted version of Px without duplicates
% Ys=unique(Py);   %Ys contains a sorted version of Py without duplicates
% Xmin1=Xs(1);     %minimum of Px
% Xmax=max(Px);
% Xmin2=Xs(2);     %second minimum of Px
% Ymin1=Ys(1);     %minimum of Py
% Ymax=max(Py);
% Ymin2=Ys(2);     %second minimum of Py
% Zmin=min(Pz);
% delta=min(Xmin2-Xmin1,Ymin2-Ymin1);
% Max=max(Xmax,Ymax);
% M=0.5*Zmin*ones(Max,Max);
% for i=1:max(Xmax,Ymax)
%      M(floor(N/2+Px(i))+1,floor(N/2+Py(i))+1)=Pz(i); %note that Px(i) and Py(i) are different numbers
% end
% scatter(Xr,Yr,3,'green')
% title(['size: ',num2str(s)])
% % [x,y]=meshgrid(M);
% % surf(x,y,M(x,y))
% %construction/extracting the reference matrix for the convolution
% Xrs=unique(Xr);   %Xs contains a sorted version of Px without duplicates
% Yrs=unique(Yr);   %Ys contains a sorted version of Py without duplicates
% Xrmin=Xrs(1)     %minimum of Px
% Xrmax=max(Xr)
% Yrmin=Yrs(1)     %minimum of Py
% Yrmax=max(Yr)
% Maxr=max(Xrmax,Yrmax);
% Mr=zeros(length(Zr),length(Zr));
% for i=1:length(Zr)
%      Mr(floor(N/2+Xr(i))+1,floor(N/2+Yr(i))+1)=Zr(i); %note that Px(i) and Py(i) are different numbers
% end
% figure;
% [v1,v2]=meshgrid(M);
% surf(v1,v2,M)
%% 
%Sliding the window and convolving matrices:Mr and M
rsl=100; %resolution
X=round(rsl*X);
Y=round(rsl*Y);
Z=round(rsl*Z);
Xr=round(rsl*Xr);
Yr=round(rsl*Yr);
Zr=round(rsl*Zr);
a=100;
Xsub=X(X>=min(Xr)-a & X<=max(Xr)+a);
I1=find(X>=min(Xr)-a & X<=max(Xr)+a); %find the indeces of the points belonging to the reference cloud
Ysub=Y(I1);
I2=find(Y>min(Yr)-a & Y<max(Yr)+a);
I=intersect(I1,I2);
Xsub=X(I);
Ysub=Y(I);
Zsub=Z(I); 


%%
%building the matrices
%L=max(max(X)-min(X)+1,max(Y)-min(Y)+1)
L=max(max(Xsub)-min(Xsub)+1,max(Ysub)-min(Ysub)+1)
disp('size of the general matrix')
M=0.5*min(Z)*ones(L,L);

disp('size of the filter matrix')
Lr=max(max(Xr)-min(Xr),max(Yr)-min(Yr))


%just to see if we can enhance the gridding efficiency
minimum=Xr(2)-Xr(1);
for i=1:length(Xr)-1
    if (Xr(i+1)-Xr(i)< minimum && Xr(i+1)-Xr(1)<abs(Yr(i+1)-Yr(i)))
        minimum = Xr(i+1)-Xr(i);
    elseif (abs(Yr(i+1)-Yr(i))<Xr(i+1)-Xr(i) && abs(Yr(i+1)-Yr(i))< minimum)
        minimum = abs(Yr(i+1)-Yr(1));
    end
end
display('minimum distance');
minimum

%building the filter matrix
Xroff=0;
if min(Xr)<0
    Xroff=abs(min(Xr))+1
elseif min(Xr)>0
    Xroff=-min(Xr)+1
end
Yroff=0;
if min(Yr)<0
    Yroff=abs(min(Yr))+1
elseif min(Yr)>0
    Yroff=-min(Yr)+1
end



Mr=zeros(Lr,Lr);
for i=1:length(Zr)
     Mr(Xr(i)+Xroff,Yr(i)+Yroff)=Zr(i); %note that Px(i) and Py(i) are different numbers
end

%building the tile matrix
% Xoff=0;
% if min(X)<0
%     Xoff=abs(min(X))+1
% elseif min(X)>0
%     Xoff=-min(X)+1
% end
% Yoff=0;
% if min(Y)<0
%     Yoff=abs(min(Y))+1
% elseif min(Y)>0
%     Yoff=-min(Y)+1
% end
Xoff=0;
if min(Xsub)<0
    Xoff=abs(min(Xsub))+1
elseif min(Xsub)>0
    Xoff=-min(Xsub)+1
end
Yoff=0;
if min(Ysub)<0
    Yoff=abs(min(Ysub))+1
elseif min(Ysub)>0
    Yoff=-min(Ysub)+1
end

for i=1:length(Zsub)
     M(Xsub(i)+Xoff,Ysub(i)+Yoff)=Zsub(i); %note that X(i) and Y(i) are different numbers
end
%%
%dist=55555*ones(L-Lr+1,L-Lr+1);
t=0
dist=5555555;
summation=0;
x0=0;
y0=0;
%convolving M with the filter matrix Mr;
for n=1:L-Lr
    for m=1:L-Lr
        fprintf('Window %d, %d \n', n,m)
        for i=1:Lr
            summation=0;
            for j=1:Lr
                summation=summation+(1/(Lr^2))*(Mr(i,j)-M(i,j))^2;
                %dist(n,m)=dist(n,m)+(1/(Lr^2))*(Mr(i,j)-M(i,j))^2; %for a
                %matrix setting
            end
            if i==Lr-1
                    %j
                    %disp('distance so far: ')
               summation
            end
        end
        if summation<dist
            dist=summation
            x0=n
            y0=m
        end
    end      
end
% [DistMin,IdxRow]=min(dist);
% [DistMin2,IdxCln]=min(DistMin);
disp('least distance measure:');
dist
%DistMin2
disp('region');
Pos=[n-Xoff,m-Yoff;n-Xoff+Lr m-Yoff+Lr]
disp('error');
(abs(Xr(1)-Pos(1,1))/Xr(1))+(abs(Yr(1)-Pos(1,2))/Yr(1))
% [IdxRow-Xoff,IdxCln-Yoff]
% [IdxRow-Xoff+Lr,IdxCln-Yoff+Lr]
%scatter(10*max(Selec(1,:)),10*max(Selec(2,:)),100,'red','x');
