function solver_heterogeneous_creeping_ctm(model,test)


if nargin<1 model =1; test = 2; end  % by default, creeping model in the creeping simulation
%-------------------------------------------------------

clc;  % clear the command window
close all % close all figures

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Examples and comparisons with the n-populations model.
% For the input parameters
%     model = 1: the creeping model
%     model = 2: the n-populations model
%     test = 1: the overtaking simulation
%     test = 2: the creeping simulation
% For other parameters in the code
%     rho1, rho2 are densities of first and second vehicle class
%     rm1, rm2 are effective jam densities of rho1 and rho2
%     vm1, vm2 are the maximum traffic veloicty of rho1 and rho2
%------------------------------------------------
% This solver is based on the finite volume Godunov method, and have
% applied the sending and receiving of vehicles on each vehicle class 
% to find the numerical flux.
% In this simulation, first vehicle class rho1 is assumed to be small
% vehicles, while the second vehicle class is larger.
%--------------------------------------------------
% Aug 29 2013
% Shimao Fan
% University of Illinois at Urbana-Champaign
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\



% parameter of maximum occupied space and maximim velocity
switch model  % gives better results for larger rm1!!
    case 1
       model_name = 'Creeping model';
       rm1 = 1.8; rm2 = 1.0;
       vm1 = 1.8; vm2 = 1.8;
    case 2
       model_name = '$n$-populations model';
       rm1 = 1.8; rm2 = 1.8;
       vm1 = 1.8; vm2 = 1.0;
end
% ----------------------------------

tfinal = 150;                % final time
len = 50;                    % length of study area
x = linspace(0,len,1000);    % number of grids 1000
dx = x(2)-x(1);              % size space step
% time step
k = (x(2)-x(1))/vm1;         % by CFL condition
lambda = k/dx;
dt = k;
t = 0:dt:tfinal;
M = length(t);               % length of time vector
N = length(x);               % length of space vector
%==========================================
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% boundary: traffic light, red
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
switch test
    case 1   % overtaking

     U(1,1:ceil(1*N/len)) = 0;
     U(1,ceil(1*N/len)+1:ceil(10*N/len)) = 0.9;
     U(1,ceil(10*N/len)+1:N) = 0;

     U(2,1:ceil(11*N/len)) = 0;
     U(2,ceil(11*N/len)+1:ceil(20*N/len)) = 0.9;
     U(2,ceil(20*N/len)+1:N) = 0;
     
     % boundary condition 
     d1l = U(1,1)+0*t;    % upstream
     d2l = U(2,1)+0*t;    % downstream

     d1r = U(1,end)+0*t;
     d2r = U(2,end)+0*t;
    
    case 2 % creeping
     U(1,1:ceil(1*N/50)) = 0;
     U(1,ceil(1*N/50)+1:ceil(19*N/50)) = .7;
     U(1,ceil(19*N/50)+1:N) = 0;

     U(2,1:ceil(20*N/50)) = 0;
     U(2,ceil(20*N/50)+1:N) = .7;
     
     % boundary condition
     d1l = U(1,1)+0*t;   % upstream 
     d2l = U(2,1)+0*t;   % downstream
      
     d1r = 0.9+0*t;
     d2r = 0.9+0*t; 
end
%----------------------------------------
U_0 = U;                % initial condition
inter = 250;            % printing figures every 250 steps


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%        The main algorithm
%/////////////////////////////////////


for n = 0:M-1
    % the boundary condition

    lbc = [d1l(n+1);d2l(n+1)];
    rbc = [d1r(n+1);d2r(n+1)];
    % insert the boundary condition

    Up1=[U(:,2:end),rbc];    % right boundary
    Um1=[lbc,U(:,1:end-1)];  % left boundary 
    
    
    
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % phase transition happens
    %////////////////////////////
    if model ==1
      rl = Um1(1,:) + Um1(2,:);   % total density on left
      r = U(1,:) + U(2,:);        % total density in middle
      rr = Up1(1,:) + Up1(2,:);   % total density on right 

      case1 = rl<=rm2;        % left state in the domain D1
      case2 = r<=rm2;         % middle state in the domain D1
      case3 = rr<=rm2;        % right state in the domain D1

      id_l = find(case1-case2); % find all phase changes ul and u
      id_r = find(case2-case3); % find all phase changes u and ur
    
      idl = case1(id_l)-case2(id_l);% -1 means left large density
      idr = case2(id_r)-case3(id_r);% +1 mean left state small
    
% 
      for i = 1:length(idr)-1
        if idr(i)<0
           U(1,id_r(i)) = rm2-U(2,id_r(i));   % insert a intermediate state
        else
           Up1(1,id_r(i)) = rm2-Up1(2,id_r(i));   % insert a intermediate state
        end
      end
    end
    
    %============================================
    % updating with godunov methods/heterogeneous extension of Cell
    % Transmission Model
     U = U+lambda*(Num_Flux(Um1,U,vm1,vm2,rm1,rm2)- ...
        Num_Flux(U,Up1,vm1,vm2,rm1,rm2));    
    %============================================
   
   

    % \\\\\\\\\\\\\\\
    % plotting part
    %/////////////////
    if mod(n,inter) ==0
    subplot(2,1,1),
    plot(x(2:1:end),U_0(1,2:1:end),'--','color',[0.8,0,0],'linewidth',2), hold on
    plot(x(2:end),U(1,2:end),'-','color',[.6,.6,.6],'linewidth',6)
    title(sprintf('%s',model_name),'fontsize',30,'interpreter', 'latex')
    axis([0 len 0 1.9])

    h = legend('Initial state','Density of small vehicles');
    set(h,'Location','NorthWest','interpreter', 'latex',...
        'position',[0.11 0.78 0.5 0.14])    
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',30)    
    set(gca,'xtick',[])     
    set(gca,'position',[.10 .52 .885 .41])
    ylabel('$\rho_1$','interpreter', 'latex')
    
    %****************
    subplot(2,1,2),
    %****************
    plot(x(2:1:end),U_0(2,2:1:end),'--','color',[0.8,0,0],'linewidth',2), hold on
    plot(x(2:end),U(2,2:end),'-.','color',[0,0,.1],'linewidth',6)
    h = legend('Initial state','Density of large vehicles');
    set(h,'Location','NorthWest','interpreter', 'latex',...
        'position',[0.11 0.35 0.5 0.14])
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',30)
    axis([0 len 0 1.9])

    ylabel('$\rho_2$','interpreter', 'latex')
    xlabel('x','interpreter', 'latex')
   
    
    set(gca,'position',[.10 .09 .885 .41])
    res = 800;
    set(gcf,'paperpositionmode','auto')
    set(gcf,'position',[10  50 res res*.8])
     
%      print figure
%      filename_save = sprintf('fig_multi_class_%1.0f_%s_%1.0f',test,model_name(1:3),ceil(n/inter));
%      print(gcf,'-dpng',filename_save,'-r290') 
%      print(gcf,'-depsc',filename_save,'-r290','-painters')
%      fixPSlinestyle([filename_save,'.eps'])
    end
    drawnow
end


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% subfunctionals
%///////////////////////////////////////////////////////
function y = new_flux(Ul,Ur,vm1,vm2,rm1,rm2)
    rl = Ul(1,:) + Ul(2,:); % total density on left
    rr = Ur(1,:) + Ur(2,:); % total density on right
%     n = length()
    inst1 = (rl==0 | rr==0); % case either initial data are zero
    inst2 = 1-inst1;   % other cases
    
    % update the numerical flux
    flux_hll = HLL(Ul,Ur,vm1,vm2,rm1,rm2);      % case inst2
    flux_sr = Num_Flux(Ul,Ur,vm1,vm2,rm1,rm2);  % sending receiving 
    indx1 = find(isnan(flux_hll(1,:)));
    indx2 = find(isnan(flux_hll(2,:)));
    flux_hll(1,indx1) = flux_sr(1,indx1);
    flux_hll(2,indx2) = flux_sr(2,indx2);    
%     y = flux_hll;
    y = (ones(2,1)*inst1).*flux_sr + (ones(2,1)*inst2).*flux_hll;
    
function F = HLL(Ul,Ur,vm1,vm2,rm1,rm2) % the HLL Riemann solver
    % where v contains boundary condition of velocity
    % eigenvalues:
    eigns1 = lambda(Ul,vm1,vm2,rm1,rm2);  % left bound
    eigns2 = lambda(Ur,vm1,vm2,rm1,rm2);  % right bound
    % another way to define c1 and c2
    c1 = min(eigns1(1,:),eigns2(1,:));
    c2 = max(eigns1(2,:),eigns2(2,:));

    case1 = c1>=0;
    case2 = c1<0 & c2>=0;
    case3 = c2<0;
    % left flux
    Fl = flux(Ul,vm1,vm2,rm1,rm2);  % plug into left boundary velocity information
    % right flux
    Fr = flux(Ur,vm1,vm2,rm1,rm2);
    %-----------------------------------------------
    % intermediate flux
    F_hll(1,:) = (c2.*Fl(1,:)-c1.*Fr(1,:)+(Ur(1,:)-Ul(1,:)).*c1.*c2)./...
            (c2-c1);
    F_hll(2,:) = (c2.*Fl(2,:)-c1.*Fr(2,:)+(Ur(2,:)-Ul(2,:)).*c1.*c2)./...
            (c2-c1);
    F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
    F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);
%--------------------------------------------------------------------------   
function y = lambda(U,vm1,vm2,rm1,rm2)  % U = (density, w), calculate eigenvalues
       r = U(1,:)+U(2,:);    % sum
       u1 = vel(r,vm1,rm1);
       u2 = max(vel(r,vm2,rm2),0);
       dv1 = DiffVel(r,vm1,rm1);
       dv2 = DiffVel(r,vm2,rm2);
       indx = find(r>rm2);
       dv2(indx) = 0;
       % characteristic speed
       c1 = U(1,:).*dv1;
       c2 = U(2,:).*dv2;
       q1 = c1+u1; q2 = u2+c2;
       delta = sqrt((q1-q2).^2+4*c1.*c2);
       
       y(1,:) = .5*(q1+q2-delta);   % slower chareateritic field
       y(2,:) = .5*(q1+q2+delta);   % faster characteristic field
%--------------------------------------------------------------
function  F = Num_Flux(Ul,Ur,vm1,vm2,rm1,rm2) % numerical flux
   % find sending function/first class
   s_1 = sending_1(Ul,vm1,rm1);
   % receiving of flux/first class
   r_1 = receiving_1(Ul,Ur,vm1,rm1);
   
   % find sending function
   s_2 = sending_2(Ul,vm2,rm2);
   % receiving of flux
   r_2 = receiving_2(Ul,Ur,vm2,rm2);
   % conservation of mass
   F(1,:) = min(s_1,r_1);   % take the minimum, flow of first class
   % conservation of momentum
   F(2,:) = min(s_2,r_2);   % take the minimum, flow of the second class
   
%-------------------------------------------------   
% define the receiving function for 1st vehicle class
function y = receiving_1(Ul,Ur,vm1,rm1)

      rhoc = (rm1-Ur(2,:))/2;  % critical density
      case1 = Ur(1,:)<=rhoc;  % maximum possible
      case2 = Ur(1,:)>rhoc;   % maximum
      u0(1,:) = rhoc;
      u0(2,:) = Ur(2,:);
      q_max = traffic_flux_1(u0,vm1,rm1); % maximum flow rate
%       q_max = max(traffic_flux_1(Ur,vm1,rm1)); % maximum flow rate      
      q = traffic_flux_1(Ur,vm1,rm1);
      y = case2.*q + case1.*q_max;
%--------------------------------------------------  
% define the receiving function for 2nd vehicle class
function y = receiving_2(Ul,Ur,vm2,rm2)
      % when r<rm2, non-creeping phase
      r = Ur(1,:) + Ur(2,:); % total density

      rhoc = max((rm2-Ur(1,:))/2,0);  % critical density
      case1 = Ur(2,:)<=rhoc;  % maximum possible
      case2 = Ur(2,:)>rhoc;   % maximum
      u0(2,:) = rhoc;
      u0(1,:) = Ur(1,:);
      q_max = traffic_flux_2(u0,vm2,rm2); % maximum flow rate
      q = traffic_flux_2(Ur,vm2,rm2);
      y = case2.*q + case1.*q_max;    
%----------------------------------------------
% define sending function for 1st vehicle class
function y = sending_1(U,vm1,rm1)% sending of creeping vehicles/for two phases
      
      rhoc = (rm1-U(2,:))/2;  % a vector of rhoc/3-params fd
      case1 = U(1,:)<=rhoc;  % maximum possible
      case2 = U(1,:)>rhoc;   % maximum
      u0(1,:) = rhoc;
      u0(2,:) = U(2,:);
      q_max = traffic_flux_1(u0,vm1,rm1); % maximum flow rate
%       q_max = max(traffic_flux_1(U,vm1,rm1));
      q = traffic_flux_1(U,vm1,rm1);
      y = case1.*q + case2.*q_max;
%------------------------------------------------      
% define sending function for 2nd vehicle class
function y = sending_2(U,vm2,rm2)% U = [\rho, w] two state variables

      rhoc = max((rm2-U(1,:))/2,0);  % critical density

      case1 = U(2,:)<=rhoc;    % maximum possible
      case2 = U(2,:)>rhoc;     % maximum
      u0(2,:) = rhoc;
      u0(1,:) = U(1,:);
      q_max = traffic_flux_2(u0,vm2,rm2); % maximum flow rate
      q = traffic_flux_2(U,vm2,rm2);
      y = case1.*q + case2.*q_max;  
%-----------------------------------------------
function y = flux(U,vm1,vm2,rm1,rm2)  % define flux function        
      u1 = vel(U(1,:)+U(2,:),vm1,rm1);
      u2 = vel(U(1,:)+U(2,:),vm2,rm2);
      y(1,:) = U(1,:).*u1;
      y(2,:) = U(2,:).*u2;
%-------------------------------------------------------
function y = traffic_flux_1(U,vm1,rm1)  % traffic flux of 1st class
         rho = U(1,:)+U(2,:);
         v = vel(rho,vm1,rm1);
         y = U(1,:).*v; 
%-------------------------------------------------------         
function y = traffic_flux_2(U,vm2,rm2) % traffic flux of 2nd class
         rho = U(1,:)+U(2,:);
         v = vel(rho,vm2,rm2);
         y = U(2,:).*v; 
         y = max(y,0*rho);
%-------------------------------------------------------
function y = DiffVel(r,vm,rm)% define derivative of velocity        
        y = 0*r-vm./rm;
%-------------------------------------------------------
function y = vel(rho,vm,rhom)  % define velocity functions 
     y = vm*(1-(rho/rhom));    % which is the quadratic form     
%-------------------------------------------------------
