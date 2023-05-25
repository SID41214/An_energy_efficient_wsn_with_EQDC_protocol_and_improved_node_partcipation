close all;
clear;
clc;
xm=100;
ym=100;
x=0; % added for better display results of the plot
y=0; % added for better display results of the plot
p=0.1; % Probability of node transmitting data
Eo=0.5; % initial energy of node
ETX=50*0.000000001; %
ERX=50*0.000000001;
Efs=10e-12;
Emp=0.0013e-12;
EDA=5*0.000000001;
rmax=2000;
do=sqrt(Efs/Emp);
Et=0;
% Number of Nodes in the field %
global n;
n=100;
global mn;
mn=10;
%loop variable to determine the termination of the network
operatingNodes=n+mn;
% Number of Dead Nodes in the beggining %
dead_nodes=0;
% Coordinates of the Sink (location is predetermined in this simulation) %
sinkx=50;
sinky=120;
 
%energy details
 
% Initial Energy of a Node (in Joules) % 
% E0=0.05; % units in Joules
% ME0=0.1;
% deadNodeThreshold=.000025;
 
LeaderNodeThreshold=.03;
global SN;
global MN;
global totalNodes;
totalNodes=n+mn+1;
 
 
 
            %%% Creation of the Wireless Sensor Network %%%
            
pts=[n,2];
pts=[rand([n,1])*xm,rand([n,1])*xm];
disp(pts);
 
%%% randomly generating mobile node points
mnpts=[mn,2];
mnpts=[rand([mn,1])*xm,rand([mn,1])*xm];
%disp(pts);
 
 
 
 
 
% X-axis coordinates of sensor node
 
 
for i=1:mn
    
    MN(i).x=mnpts(i,1);
    MN(i).y=mnpts(i,2);
    MN(i).radius=50; %assuming 20meters
    MN(i).sel=0;
    
end
 
 
%% Adding base station to the sensor node list
 
 
 %adding sink point to pts array
 
%aading sink point to pts array
sink=[sinkx,sinky];
disp('sink is..');
disp(sink);
pts=[sink;pts];
% pts(end+1,1)=sinkx;
% pts(end,2)=sinky;
disp(pts);
 
 
 
 
 
 
 
 
 
 
 
% Plotting the WSN %
for i=1:n+1
    
    
      SN(i).id=i; % sensor's ID number
     SN(i).radius=20; %assuming 20 meters
   SN(i).x=pts(i,1);
   SN(i).y=pts(i,2);
%    SN(i).chainID=0;
%  
%     SN(i).E=E0;     % nodes energy levels (initially set to be equal to "Eo"
%         
%     SN(i).order=0;
%     SN(i).sel=0;    % states if the node has already operated for this round or not (if 0 then no, if 1 then yes) 
%     SN(i).rop=0;    % number of rounds node was operational
%     SN(i).tel=0;    % states how many times the node was elected as a Cluster Head
%     order(i)=0;
%     SN(i).weight=0;
 
 
       hold on;
    figure(1)
    plot(x,y,xm,ym,SN(i).x,SN(i).y,'ob',sinkx,sinky,'*r');
     
    title 'Wireless Sensor Network';
    xlabel '(m)';
    ylabel '(m)';
    
end
 
%trying to create coordinates array and y coordinates array seperately to
%give coordinates in the grpah
xCoord=(n+1+mn);
yCoord=(n+1+mn);
 
for i=1:n+1
   xCoord(i)=SN(i).x;
   yCoord(i)=SN(i).y;
end
%xCoord(end+1)=sinkx;
%xCoord=[sinkx,xCoord];
%yCoord(end+1)=sinky;
%yCoord=[sinky,yCoord];
% for i=1:mn
%     
%      plot(x,y,xm,ym,MN(i).x,MN(i).y,'+r');
%      
%      
%     
% end
 
 
 
 
 
 
 
   
 
 
 % Calculates Euclidaean Distance Between Each Node and the Sink (Base Station) %
 for i=1:n+1
 
     
    dts(i)=sqrt((sinkx-xCoord(i))^2 + (sinky-yCoord(i))^2);
%    SN(i).Esink=Eelec*k + Eamp*k*(SN(i).dts)^2;
    T(i)=dts(i);
 end
% disp('t of one is..');
 %disp(T(1));
 %sdisp('t of 21 is..')
 %disp(T(21));
 %euclidean distance from each node to the sink node now is stored in T
 %array T   can access by T(i) in for loop
 
 
 
%checking operatingNodes loop
 
 
%Delaunay triangulating all the sensor nodes
 
 
 
DT = delaunayTriangulation(pts); 
figure(2)
 
e=DT.edges;
plot(sinkx, sinky, 'r*', 'markersize', 5);
hold on
 
pl=triplot(DT,'-ob');
hold on
plot(pts(1),pts(2),'-r','LineWidth',2);
 
%nodes to be displayed
 
 
 
%plot(xy(:, 1), xy(:, 2), 'ro', MarkerSize=12, MarkerFaceColor='r');
title 'Delaunay Triangulation of the network with Sink and Static sensor nodes';
 
 
 
 
DT = delaunayTriangulation(pts); 
figure(3)
hold on
for i=1:mn
    
     plot(x,y,xm,ym,MN(i).x,MN(i).y,'+r');
     
     
    
end
pl=triplot(DT,'-ob');
title 'Current positions of mobile sensor nodes'
 
%identifying gaps
 
%identifying edge lengths of DT
 
edgelengths = sqrt(sum((DT.Points(e(:,1),:) - DT.Points(e(:,2),:)).^2,2));
%disp('edge lengths');
disp(edgelengths);
%disp(edgelengths);
noOfEdges=size(edgelengths);
 
%disp(pts);
%disp((DT.Points(e(:,1),:)));
%disp((DT.Points(e(:,2),:)));
 
 
mnx=(mn);
mny=(mn);
pending_edge=(0);
%checking if it is greater than 2*radius of static node
 global k;
 k=0;
 in=1;
 selectedEdge=[];
for i=1:noOfEdges
 
  edgeX1=DT.Points(e(i,1),1);
 
 edgeY1=DT.Points(e(i,1),2);
 edgeX2=DT.Points(e(i,2),1);
 edgeY2=DT.Points(e(i,2),2);
    
    
    if k<=mn
   
if edgelengths(i)> 2*SN(1).radius
  % disp('found a hole');
   %k=k+1;
   
     %check if the found hole is patchable or not   
   if 2*MN(1).radius>abs(edgelengths(i)-2*SN(1).radius)
       
     %  disp('patchable hole found');
       selectedEdge=[selectedEdge,i];
       
       % edgeX1=DT.Points(e(i,1),1);
 
 %edgeY1=DT.Points(e(i,1),2);
 %edgeX2=DT.Points(e(i,2),1);
 %edgeY2=DT.Points(e(i,2),2);
 if k~=mn
  k=k+1;
 end
 
 
 
 mnx(k)=(edgeX1+edgeX2)/2;
 mny(k)=(edgeY1+edgeY2)/2;
% newPoint=[1,2];
%newPoint(1,1)=mnx(k);
%newPoint(1,2)=mny(k);
 
 %DT.Points(end+1,:)=newPoint;
 
 %if not patchable
          
   end
   
   
   
    
   
end
%disp(k);
 
 
 
 
 
 
 
 
 
    end 
   
   MNx=(mn);
   MNy=(mn);
   
     
   
 
 
end
disp(k);
 
if k<mn
    for i=1:length(edgelengths)
        
  edgeX1=DT.Points(e(i,1),1);
 
 edgeY1=DT.Points(e(i,1),2);
 edgeX2=DT.Points(e(i,2),1);
 edgeY2=DT.Points(e(i,2),2);
        if edgelengths(i)>30 && k<mn && ismember(i,selectedEdge)==0
            k=k+1;
            selectedEdge=[selectedEdge,i];
             mnx(k)=(edgeX1+edgeX2)/2;
               mny(k)=(edgeY1+edgeY2)/2;
        end
        
        
        
        
        
    end
   
end
 
%disp('k now is....');
%disp(k);
%disp('selected edges are.....');
%disp(selectedEdge);
 
 
 
 
 %if k<mn
      
    
     
  %   for i=k+1:mn
   %       disp('getting in.......');
         
          
    %      y = randsample(n,1);
         % z = randsample(n,1);
          
     %     edgeX1=DT.Points(e(y,1),1);
 
% edgeY1=DT.Points(e(y,1),2);
 %edgeX2=DT.Points(e(y,2),1);
 %edgeY2=DT.Points(e(y,2),2);
          
          
  %       mnx(i)=(edgeX1+edgeX2)/2;
   %     mny(i)=(edgeY1+edgeY2)/2;
    %     disp(mnx(i));
     %    disp(mny(i));
         
          
    % end
     % k=mn;
   
% end
 
 
 
 
 
 
 
 
 
%replot with mobile nodes moved to the positions according to the
%conditions
 
 
DT = delaunayTriangulation(pts); 
figure(5)
 
%disp(pts);
 
pl=triplot(DT,'-ob');
 title 'Movable sensor nodes placed at the mid points of the longest edges';
hold on;
 
for i=1:k
    
     plot(mnx(i),mny(i),'+r');
     pts(end+1,1)=mnx(i);
pts(end,2)=mny(i);

    
end
 disp(pts);
 
%graph is to be created here
 
 
 
 
 
 
 
 
 
 
newDT = delaunayTriangulation(pts);
%calling edges to get this delaunay graph
e=newDT.edges;
%creating grpah object to get graph of DT and minimum spanning tree with
%all mobile nodes placed finally
 
%xCood and yCoord list to be updated with mobile nodes coordinates
 
for i=1:mn
    
    xCoord(n+1+i)=mnx(i);
    yCoord(n+1+i)=mny(i);
    
end
 
%disp('size of points array');
%disp(size(pts));
 
 
for i=1:mn
%     
%     %disp('n+i+1   is..........................................................');
%     %disp(n+i+1);
%     SN(n+1+i).id=n+1+i; % sensor's ID number
%     SN(n+1+i).radius=40; %assuming 20 meters
   SN(n+1+i).x=pts(n+1+i,1);
   SN(n+1+i).y=pts(n+1+i,2);
%  
%     SN(n+1+i).E=ME0;     % nodes energy levels (initially set to be equal to "Eo"
%     SN(n+1+i).cond=1;   % States the current condition of the node. when the node is operational its value is =1 and when dead =0
%     SN(n+1+i).dts=0;    % nodes distance from the sink
%     SN(n+1+i).role=0;   % node acts as normal if the value is '0', if elected as a cluster head it  gets the value '1' (initially all nodes are normal)
%     SN(n+1+i).pos=0;
%     SN(n+1+i).closest=0;
%     SN(n+1+i).prev=0;
%     SN(n+1+i).dis=0;    % distance between two nodes headin towards to the cluster head from position 1
%     SN(n+1+i).dis2=0;   % distance between two nodes headin towards to the cluster head from position 2
%     SN(n+1+i).order=0;
%     SN(n+1+i).sel=0;    % states if the node has already operated for this round or not (if 0 then no, if 1 then yes) 
%     SN(n+1+i).rop=0;    % number of rounds node was operational
%     SN(n+1+i).tel=0;    % states how many times the node was elected as a Cluster Head
%     order(n+1+i)=0;
%      SN(i).chainID=0;
%      SN(i).weight=0;
%  
%  
%     
 end
 
 
newDT = delaunayTriangulation(pts); 
%calling edges to get this delaunay graph
e=newDT.edges;
%creating grpah object to get graph of DT and minimum spanning tree with
%all mobile nodes placed finally
 
%xCood and yCoord list to be updated with mobile nodes coordinates
 
for i=1:mn
    
    xCoord(n+1+i)=mnx(i);
    yCoord(n+1+i)=mny(i);
    
end
 
figure(6)
pl=triplot(newDT,'-ob');
title 'Retriangulating the network with mobile nodes'
hold on;
 

figure(7)
hold on
for i=1:n+1+mn
    plot(SN(i).x, SN(i).y, 'ob');
end


% % Concatenate all points into a single array
% all_x = [];
% all_y = [];
% 
% for i=1:n+1
%     all_x = [all_x; SN(i).x];
%     all_y = [all_y; SN(i).y];
% end
% 
% for i=1:mn
%     all_x = [all_x; xCoord(n+1+i)];
%     all_y = [all_y; yCoord(n+1+i)];
% end
% 
% all_x = [all_x; sinkx];
% all_y = [all_y; sinky];
% 
% figure(8)
% % Plot all points together
% plot(all_x, all_y, '*k');
% 
% % Add legend
% legend('SN points', 'MN points', 'Sink point', 'All points');

disp(pts);
% Remove sink coordinates from pts array
pts = pts(pts(:,1) ~= sink(1) | pts(:,2) ~= sink(2), :);
disp(pts);
n=n+mn
sink=[sinkx,sinky];
disp('sink is..');
disp(sink);
% pts=[sink;pts];
pts(end+1,1)=sinkx;
pts(end,2)=sinky;
disp(pts);
 
disp(n);
for i=1:n+1
   S(i).xd=pts(i,1);
   S(i).yd=pts(i,2);
      hold on;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%                           LEACH                               %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
for h=1:1
    S(n+1).xd=sinkx;
    S(n+1).yd=sinky;
    Et=0;
    %% composing sensors
    for i=1:1:n
        XR(i)=S(i).xd;
        YR(i)=S(i).yd;
        distance=sqrt( (S(i).xd-(S(n+1).xd) )^2 + (S(i).yd-(S(n+1).yd) )^2 );
        S(i).distance=distance;
        S(i).G=0;
        %initially there are no cluster heads only nodes
        S(i).type='N';
        S(i).E=Eo;
        Et=Et+S(i).E;
        %figure(h*10)
        %  plot(S(i).xd,S(i).yd,'bo');
        %  text(S(i).xd+1,S(i).yd-0.5,num2str(i));
        %  hold on;
    end

%plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
%text(S(n+1).xd+1,S(n+1).yd-0.5,num2str(n+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
countCHs=0;  %variable, counts the cluster head
cluster=1;  %cluster is initialized as 1
flag_first_dead=0; %flag tells the first node dead
flag_half_dead=0;  %flag tells the 10th node dead
flag_all_dead=0;  %flag tells all nodes dead
first_dead=0;
half_dead=0;
all_dead=0;
allive=n;
%counter for bit transmitted to Bases Station and to Cluster Heads
packets_TO_BS=0;
packets_TO_CH=0;
packets_TO_BS_per_round=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r=0:1:rmax
    figure(20)
    clf
    r
    title('Normal nodes:Blue||CH:Red||Dead:Empty circle||Round no:'+string(r));
    hold on

     line([50,50],[0,120],'color','green')
%      line([0,100],[100,100],'color','green')
    line([100,100],[100,0],'color','green')
    line([0,100],[0,0],'color','green')
    line([0, 100], [60, 60], 'color', 'green');

    line([0,0],[0,100],'color','green')
    line([0, 0], [0, 120], 'color', 'green');
%     line([0, 120], [100, 120], 'color', 'green');
%     line([0, 120], [100, 120], 'color', 'green');
line([0, 100], [120, 120], 'color', 'green');
line([100, 100], [0, 120], 'color', 'green');





    packets_TO_BS_per_round=0;
    %Operations for epochs
    if(mod(r, round(1/p) )==0)
        for i=1:1:n
            S(i).G=0;
            S(i).cl=0;
        end
    end
   
    %Number of dead nodes
    dead=0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     for i=1:1:n
        %checking if there is a dead node
        if (S(i).E<=0)
            plot(S(i).xd,S(i).yd,'o');            
            dead=dead+1;
            if (dead==1)
              if(flag_first_dead==0)
                 first_dead=r;
                 flag_first_dead=1;
              end
            end
            if(dead==0.5*n)
              if(flag_half_dead==0)
                  half_dead=r;
                  flag_half_dead=1;
              end
            end
            if(dead==n)
              if(flag_all_dead==0)
                  all_dead=r;
                  flag_all_dead=1;
              end
            end
        end
        if S(i).E>0
            S(i).type='N';
        end
    end
    STATISTICS.DEAD(h,r+1)=dead;
    STATISTICS.ALLIVE(h,r+1)=allive-dead;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %head finding
    ch1=0;ch2=0;ch3=0;ch4=0;ch=0; 
    countCHs=0;
    cluster=1;
    for i=1:1:n
        if(S(i).E>0)
            temp_rand=rand;
            if ( (S(i).G)<=0)
                
                %Election of Cluster Heads for normal nodes
                if ( temp_rand <= ( p/ ( 1 - p * mod(r,round(1/p)) )) )
                    
                    %countCHs=countCHs+1;
                    %packets_TO_BS=packets_TO_BS+1;
                    %packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                    %PACKETS_TO_BS(r+1)=packets_TO_BS;
                    
                    if(S(i).xd>0 && S(i).xd<=50 && S(i).yd>0 && S(i).yd<=60 && ch1==0)
                        ch1=i;%msgbox('ch1')
                        ch=i;
                        %line([S(i).xd,SN(ch1).x],[S(i).yd,SN(ch1).y])
                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS;
                    end
                    
                    if(S(i).xd>0 && S(i).xd<=50 && S(i).yd>60 && S(i).yd<=120 && ch2==0)
                        %line([S(i).xd,SN(ch2).x],[S(i).yd,SN(ch2).y])
                        ch2=i;%msgbox('ch2')
                        ch=i;
                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS;
                    end
                    
                    if(S(i).xd>50 &&S(i).xd<=100 && S(i).yd>60 &&S(i).yd<=120 && ch3==0)
%                         line([S(i).xd,SN(ch4).x],[S(i).yd,SN(ch4).y])
                        ch3=i;%msgbox('ch3')
                        ch=i;
                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS;
                    end

                    if(S(i).xd>50 && S(i).xd<=100 && S(i).yd>0 && S(i).yd<=60&& ch4==0)
                        %line([S(i).xd,SN(ch3).x],[S(i).yd,SN(ch3).y])
                        ch4=i;%msgbox('ch4')
                        ch=i;
                        countCHs=countCHs+1;
                        packets_TO_BS=packets_TO_BS+1;
                        packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                        PACKETS_TO_BS(r+1)=packets_TO_BS;
                    end
                    if (ch>0)
                        S(ch).type='C';
                        S(ch).G=round(1/p)-1;
                        C(cluster).xd=S(ch).xd;
                        C(cluster).yd=S(ch).yd;
                        %plot(S(ch).xd,S(ch).yd,'k*');

                        distance=sqrt( (S(ch).xd-(S(n+1).xd) )^2 + (S(ch).yd-(S(n+1).yd) )^2 );
                        C(cluster).distance=distance;
                        C(cluster).id=ch;
                        X(cluster)=S(ch).xd;
                        Y(cluster)=S(ch).yd;
                        cluster=cluster+1;

                        %Calculation of Energy dissipated
                        distance;
                        if (distance>do)
                            S(ch).E=S(ch).E- ( (ETX+EDA)*(4000) + Emp*4000*( distance*distance*distance*distance ));
                        end
                        if (distance<=do)
                            S(ch).E=S(ch).E- ( (ETX+EDA)*(4000)  + Efs*4000*( distance * distance ));
                        end
                    end
                    ch=0;
                    if (ch1>0 && ch2>0 && ch3>0 &&ch4>0)
                      % msgbox('break');
                      %x=[ch1,ch2,ch3,ch4, cluster]
                   
                       break;
                    end
                    %x=[ch1,ch2,ch3,ch4, cluster]
                end
    
            end     
        end
    end
   % Distance matrix initialization
distMatrix = inf(cluster-1, cluster-1);

% Calculate distances between cluster heads
for i = 1:cluster-1
    for j = i+1:cluster-1
        distMatrix(i, j) = sqrt((C(i).xd - C(j).xd)^2 + (C(i).yd - C(j).yd)^2);
        distMatrix(j, i) = distMatrix(i, j);
    end
end

% Minimum spanning tree algorithm
visited = zeros(1, cluster-1);
visited(1) = 1;
routeMatrix = zeros(cluster-1, cluster-1);

while sum(visited) < cluster-1
    minDist = inf;
    minIdx = 0;
    
    for i = 1:cluster-1
        if visited(i)
            % Find the cluster head with the shortest distance to the visited set
            for j = 1:cluster-1
                if ~visited(j) && distMatrix(i, j) < minDist
                    minDist = distMatrix(i, j);
                    minIdx = j;
                end
            end
        end
    end
    
    % Connect the cluster head to the visited set
    routeMatrix(minIdx, find(visited == 1, 1)) = 1;
    visited(minIdx) = 1;
end

% Plotting the route
for i = 1:cluster-1
    for j = 1:cluster-1
        if routeMatrix(i, j) == 1
             line([C(i).xd, C(j).xd], [C(i).yd, C(j).yd], 'Color', 'r', 'LineStyle', '--');
        end
    end
end


    %head connection
%     chs=[ch1,ch2,ch3,ch4];
%     for i=1:3
%         if chs(i)==0
%             continue;
%         end
%         for j=i+1:4
%             if(chs(j)==0)
%                 continue;
%             end
%              plot([S(chs(i)).xd,S(chs(j)).xd],[S(chs(i)).yd,S(chs(j)).yd],'r--');break;
%         end
%     end
    STATISTICS.COUNTCHS(h,r+1)=countCHs;
    % or STATISTICS.COUNTCHS(h,r+1)=clster-1;
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Election of Associated Cluster Head for Normal Nodes
    plot(S(n+1).xd,S(n+1).yd,'o', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
   for i=1:1:n
       if( S(i).type=='N' && S(i).E>0 )
           if(cluster-1>=1)
                min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
                min_dis_cluster=0;
                for c=1:1:cluster-1
                    temp=min(min_dis,sqrt( (S(i).xd-C(c).xd)^2 + (S(i).yd-C(c).yd)^2 ) );
                   if ( temp<min_dis )
                       min_dis=temp;
                       min_dis_cluster=c;
                   end
                end
               %Calculating the culsterheads%
               if(min_dis_cluster~=0)    
                    min_dis;
                    if (min_dis>do)
                        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
                    end
                    if (min_dis<=do)
                        S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
                    end

                    S(C(min_dis_cluster).id).E = S(C(min_dis_cluster).id).E- ( (ERX + EDA)*4000 ); 

                    if(S(i).type=='N')
                        plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'black');      
                        line([S(i).xd,S(C(min_dis_cluster).id).xd],[S(i).yd,S(C(min_dis_cluster).id).yd]);
                    elseif(S(i).type=='C')
                        plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
                       

                    else
                        plot(S(i).xd,S(i).yd,'bo');
                        line([S(i).xd,S(n+1).xd],[S(i).yd,S(n+1).yd]);
                    end
                    packets_TO_CH=packets_TO_CH+1;
               else
                    min_dis;
                    if (min_dis>do)
                        S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
                    end
                    if (min_dis<=do)
                        S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
                    end
                    packets_TO_BS=packets_TO_BS+1;
                    packets_TO_BS_per_round=packets_TO_BS_per_round+1;
                    PACKETS_TO_BS(r+1)=packets_TO_BS;
                     % Plot the line to the sink
                plot([S(i).xd, S(n+1).xd], [S(i).yd, S(n+1).yd], 'g--');
               end
                S(i).min_dis=min_dis;
               S(i).min_dis_cluster=min_dis_cluster;
        else
            min_dis=sqrt( (S(i).xd-S(n+1).xd)^2 + (S(i).yd-S(n+1).yd)^2 );
            if (min_dis>do)
                S(i).E=S(i).E- ( ETX*(4000) + Emp*4000*( min_dis * min_dis * min_dis * min_dis)); 
            end
            if (min_dis<=do)
                S(i).E=S(i).E- ( ETX*(4000) + Efs*4000*( min_dis * min_dis)); 
            end
            packets_TO_BS=packets_TO_BS+1;
            packets_TO_BS_per_round=packets_TO_BS_per_round+1;
            
       end
  end
       
  if(S(i).type=='C'&&S(i).E>0)
        plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
  elseif(S(i).type~='N'&&S(i).E>0)
        plot(S(i).xd,S(i).yd,'o', 'MarkerSize', 5, 'MarkerFaceColor', 'black');
  elseif S(i).type == 'N' && S(i).E > 0 && S(i).min_dis_cluster == 0
        plot(S(i).xd, S(i).yd, 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'black');
        line([S(i).xd,S(n+1).xd],[S(i).yd,S(n+1).yd],'Color','black');
    
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Connection of Cluster head and Sink
    cl_list=[0 0 0 0];
    cl_sequence=[0 0 0 0];
    cl_number=0;
    for ii=1:n
        if (S(ii).type == 'C' && S(ii).E>0)
            cl_number=cl_number+1;
            cl_list(cl_number)=ii;            
        end
    end
    min_dis_cluster=0;
    min_dis=999999999;
    for c=1:1:cluster-1
        %min_dis=sqrt( (S(cl_list(c)).xd-S(n+1).xd)^2 + (S(cl_list(c)).yd-S(n+1).yd)^2 );   
        if cl_list(c)==0
            break;        
        end
        temp=min(min_dis,sqrt( (S(cl_list(c)).xd-S(n+1).xd)^2 + (S(cl_list(c)).yd-S(n+1).yd)^2 ));
           
           if ( temp<min_dis )
               min_dis=temp;
               min_dis_cluster=c;
           end
           
    end
    cl_sequence(1)=min_dis_cluster;
    if min_dis_cluster==0 
        continue
    end
    plot([S(n+1).xd,S(cl_list(min_dis_cluster)).xd],[S(n+1).yd,S(cl_list(min_dis_cluster)).yd],'r--');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
STATISTICS.PACKETS_TO_CH(h,r+1)=packets_TO_CH;
STATISTICS.PACKETS_TO_BS(h,r+1)=packets_TO_BS;
STATISTICS.PACKETS_TO_BS_PER_ROUND(h,r+1)=packets_TO_BS_per_round;
STATISTICS.THROUGHPUT(h,r+1)=STATISTICS.PACKETS_TO_BS(h,r+1)+STATISTICS.PACKETS_TO_CH(h,r+1);

 En=0;
for i=1:n
    if S(i).E<=0
        continue;
    end
    En=En+S(i).E;
end
ENERGY(r+1)=En;
STATISTICS.ENERGY(h,r+1)=En;


end
first_dead_LEACH(h)=first_dead
half_dead_LEACH(h)=half_dead
all_dead_LEACH(h)=all_dead


% cluster head display-------

end 
for r=0:rmax
    STATISTICS.DEAD(h+1,r+1)=sum(STATISTICS.DEAD(:,r+1))/h;
    STATISTICS.ALLIVE(h+1,r+1)=sum(STATISTICS.ALLIVE(:,r+1))/h;
    STATISTICS.PACKETS_TO_CH(h+1,r+1)=sum(STATISTICS.PACKETS_TO_CH(:,r+1))/h;
    STATISTICS.PACKETS_TO_BS(h+1,r+1)=sum(STATISTICS.PACKETS_TO_BS(:,r+1))/h;
    STATISTICS.PACKETS_TO_BS_PER_ROUND(h+1,r+1)=sum(STATISTICS.PACKETS_TO_BS_PER_ROUND(:,r+1))/h;
    STATISTICS.THROUGHPUT(h+1,r+1)=sum(STATISTICS.THROUGHPUT(:,r+1))/h;
    STATISTICS.COUNTCHS(h+1,r+1)=sum(STATISTICS.COUNTCHS(:,r+1))/h;
    STATISTICS.ENERGY(h+1,r+1)=sum(STATISTICS.ENERGY(:,r+1))/h;
end

first_dead=sum(first_dead_LEACH)/h;
half_dead=sum(half_dead_LEACH)/h;
all_dead=sum(all_dead_LEACH)/h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0:rmax;
figure(21)
plot(r,STATISTICS.DEAD(h+1,r+1));
title('Dead Nodes')
xlabel('No of Rounds')
ylabel('Number of Dead Nodes')
figure(22)
plot(r,STATISTICS.ALLIVE(h+1,r+1));
title('Live Nodes')
xlabel('No of Rounds')
ylabel('Number of live Nodes')
figure(23)
plot(r,STATISTICS.PACKETS_TO_BS(h+1,r+1));
title('pkts to BS')
xlabel('No of Rounds')
ylabel('Number of packets')
figure(24)
plot(r,STATISTICS.PACKETS_TO_BS_PER_ROUND(h+1,r+1));
title('pkts to BS per round')
xlabel('No of Rounds')
ylabel('Number of packets')
figure(25)
plot(r,STATISTICS.PACKETS_TO_CH(h+1,r+1));
title('pkts to CH')
xlabel('No of Rounds')
ylabel('Number of packets')
figure(26)
plot(r,STATISTICS.THROUGHPUT(h+1,r+1));
title('THROUGHPUT')
xlabel('No of Rounds')
ylabel('No of bits')
figure(27)
plot(r,STATISTICS.COUNTCHS(h+1,r+1));
title('COUNTCHS')
xlabel('No of Rounds')
ylabel('Number of cluster head ')
figure(28)
plot(r,STATISTICS.ENERGY(h+1,r+1));
title('Average Residual Energy') 
xlabel('No of Rounds')
ylabel('Energy in joules')