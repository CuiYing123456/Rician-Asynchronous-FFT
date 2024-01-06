clear all;
clc;
close all;

N =100  ;   % number of devices
M =32;     % number of anttenas
p = 0.1;    % device actice probability
Pt = 1;        % transmit power
n_var = 0.1;   %  noise power
 
% CFOrange =  2*pi*1;
  
  NumGrid = 2^7;      %    discretization size $Q$
 
Omega = 2^6;         %   integer CFO

CCFOrange =  2*(Omega.*( 2*pi/  NumGrid ));  % maximum CFO

LL =8;          % pilot length
Loop =5;        % number of experiments

TT_g = 2;%2-1;   % maximum STO $D$

         
PSU_NEW_tg = zeros( length(LL),length( TT_g) );         % Estimated detection error probability
 

 
CT_Prop_tg = zeros( length(LL),length( TT_g) );
 
PData_TR_tg= zeros( length(LL),length( TT_g) );
PData_NEW_tg= zeros( length(LL),length( TT_g) );


AandD = zeros(Loop,N, length( TT_g)  );
AandDest_Prop = zeros(Loop,N, length( TT_g)  );

for loop_tg = 1: length(TT_g)
    
    T_g = TT_g( loop_tg);          
    
    CFOrange = CCFOrange( 1);
    
    X = zeros(Loop,N,M);%靽⊿��拚
    A = zeros(Loop,N);%�???瘣餌嚙??
    B =  zeros(N,  T_g + 1, Loop );
    Tau =  zeros(Loop,N);
    H  = zeros(N,M,Loop);
    H_bar  = zeros(N,M,Loop);
     Pind = zeros(Loop,N);
    
     
     
   %%  this part is for generating $Loop$ i.i.d. 
    for ii = 1:Loop
        ii;
        h = normrnd(0,sqrt(1/2),[N M]) + 1i * normrnd(0,sqrt(1/2),[N M]);  %    Rayleigh fading
        
        kappa_all = 1e-1;                 % Rician factor of all devices
        Kappa = kappa_all *ones(N,1);    %  vector consisting of Rician factors 
        gamma_sq = sqrt(1./(Kappa + 1));
%         h_hat = 1*exp(j*  2*pi*rand(N,M)   );
              h_hat = 1*exp(1j*  diag(2*pi*rand(N,1 ))* ones(N,1)* [0:M-1]   );  % normalized LoS channel matrix
        h_bar = diag(sqrt(Kappa))*h_hat;               % LoS channel matrix
        
        h = diag( gamma_sq ) *(h_bar +    normrnd(0,sqrt(1/2),[N M]) + 1i * normrnd(0,sqrt(1/2),[N M]));  % channel matrix
        
        
        
        a = randsrc(N,1,[[0 1];[1-p p]]);     %  activity states $\mathbf a$, the active probability is $p$
        T_delay =      floor(  (T_g + 1)* rand( 1,N)) ;   % 
          P_index =   T_delay +1;
          Pind(ii,:) = P_index;
        
        Beta = zeros( T_g + 1, N );
        for  jj = 1: N
            Beta (   T_delay(jj )   +  1,jj  ) =1;
        end
        
        %     a1 = diag(a);
        %     x = a.*h;
        Tau(ii,:) =  T_delay.*a'  ;
        B(:, :,ii   ) = Beta';
        H(:,:,ii) = h;
        H_bar(:,:,ii) = h_bar;
        
        A(ii,:) = a';
    end
    
       AandD(:,:, loop_tg) =   Tau;
    
    
    
    
    
    
    
    %% LSC �典�?
    
    A1 = sum(A,2);
    num1 = size(A,1);
    num2 = size(A,2);


    if Omega == NumGrid/2
         num_EP = NumGrid*(T_g+1);
    else
        num_EP  = (2*Omega+1)*(T_g+1);
    end
   
    X_est = zeros(num2,     NumGrid,  num1);
    
    X_est_lasso = zeros( num2,     NumGrid,  num1 );
    
    X_est_AMP = zeros( num2,     NumGrid,  num1 );
    
    X_est_Liuliang = zeros(  num1, num2);
    
    
  
    X_est_NEW = zeros(num1,num2);

    X_est_NEWb = zeros(num1,num2);

    X_est_NEWc = zeros(num1,num2);
    
    
    A_est_NEW  = zeros(num1,num2);
    
    X_est1 = zeros(num1,num2,M);
    
    
    
    MN = zeros(length(LL),1);
    t = zeros(length(LL),1);
    
    
    %% define the metric and recorder
    PSU_TR= zeros(length(LL),1);
    MSE_TR = zeros(length(LL),1);
    
    PSU_NEW= ones(length(LL),1);
    MSE_NEW= zeros(length(LL),1);
    
    PSU_ANK = zeros(length(LL),1);
    MSE_ANK = zeros(length(LL),1);
    
    
    PSU_LASSO= zeros(length(LL),1);
    MSE_LASSO = zeros(length(LL),1);
    
    
    PSU_Liuliang= zeros(length(LL),1);
    MSE_Liuliang = zeros(length(LL),1);
   
    PSU_AMP= zeros(length(LL),1);
     
    
    for l = 1:length(LL)
        tto_start = clock;
        
        L = LL( l);
        CT_Prop  = 0;
       
        
        
        for ii = 1:num1
            %                 tic
            ii;
            %        T_delay  =   Tau(ii,:);
            
            h  = H(:,:,ii);
            h_bar = H_bar(:,:,ii);
            a =  A(ii,:)';
            Beta   = B(:, :,ii   )';
           
            
             [ P_array, P_matrix,P_array_CFO, P_matrix_CFO, P_orig,  MatrixtriuPPH, Magic,CFO]   = Fx_P_Gaussian_Pilotpool_CFO ( T_g , L, N, CFOrange ) ;  %generate effective pilots
            
       
            
            L_tilde  = size(P_matrix ,1);  % length of effective pilot $L_{i}$
            
            
            
            Z = normrnd(0,sqrt(n_var/2),[ L_tilde ,M])+1i * normrnd(0,sqrt(n_var/2),[  L_tilde ,M]);  %AWGN
            
            P_index = Pind( ii, :);
            P_effct = rand( L_tilde,N ); 
            active_device_index = find(a);
            if isempty(active_device_index)  
                Y = Z;
            else
            for i_act = 1:length(active_device_index)
               P_effct(:, active_device_index(i_act)) = P_array_CFO(:, P_index( active_device_index(i_act)),active_device_index(i_act));
            end
             Y =  P_effct  *   diag(  a) * h    + Z;   %  received signal $Y_i$,  in which large-scale fading powers, $g_n, n\in \mathcal N$, are set to 1.
            end
          
          
        %%  4 Proposed BCD algorithms: Prop-MLE-L and Prop-MLE-S. 
          %   If TT_g>0, CCFOrange=0, the case is Asynchronous case-t
          %    if TT_g=0, CCFOrange>0, the case is Asynchronous case-f 
           %    if TT_g>0, CCFOrange>0, the case is Asynchronous case-(t,f) 
          
                 Sigma_bar = Y * Y' ./ M;
   
          [gamma, Indpilot ,w, time_CD_prop] = Prop_FT_Ri_Fourcases_LcfoLsto  (P_matrix,  P_orig, Sigma_bar, Y, n_var,  Kappa, ones(N,1),h_bar,Omega,  NumGrid);  % Prop-MLE-L
          
%           [gammb, Indpilot ,w, time_CD_prop] = Prop_FT_Ri_Fourcases_LcfoSsto (P_matrix,  Sigma_bar, Y, n_var,  Kappa, ones(N,1),h_bar,N ,Omega,NumGrid); % Omega, NumGrid/2
 
          [gammc, Indpilot ,w, time_CD_prop] = Prop_FT_Ri_Fourcases_ScfoSsto (P_matrix,  Sigma_bar, Y, n_var,  Kappa, ones(N,1),h_bar,N ,Omega,NumGrid);  % Prop-MLE-S

 
           
            CT_Prop =   CT_Prop  +  time_CD_prop;
%           CT_gradFFT =     CT_gradFFT  + time_CD_prop_gradFFT;
            %                 toc
            
             X_est_NEW ( ii, :) = gamma';
%             X_est_NEWb ( ii, :) = gammb';
             X_est_NEWc ( ii, :) = gammc';

            Indpilot_NEW (ii,:) = Indpilot';


            
        end
          
        
         %%  2. Hard decision
        
         
        
        %%  2.4  Proposed
         
         [PSU_NEW,  PdataU_NEW, ~] = Hard_decision_opt_thresh(  X_est_NEW , A  , AandD(:,:, loop_tg) , Indpilot_NEW );
%         [PSU_NEWb,  PdataU_NEW, ~] = Hard_decision_opt_thresh(  X_est_NEWb , A  , AandD(:,:, loop_tg) , Indpilot_NEW );
        [PSU_NEWc,  PdataU_NEW, ~] = Hard_decision_opt_thresh(  X_est_NEWc , A  , AandD(:,:, loop_tg) , Indpilot_NEW );
         
    end
%       PSU_AMP_tg(1,loop_tg) = PSU_AMP;
%     PSU_LASSO_tg(1,loop_tg) =PSU_LASSO;
%     PSU_TR_tg(1,loop_tg) = PSU_TR;
    PSU_NEW_tg(1,loop_tg) = PSU_NEW;
    t_end= clock;
%     t_T_g = etime(t_end, t_start)
    
end

% close all
% figure(1)
% plot(TT_g, (PSU_TR_tg)','b-+')
% hold on
% plot(TT_g, (PSU_NEW_tg)','r-*')
% xlabel('Delay/Length of pilot ')
% ylabel('Error probability')



% figure(2)
% plot(MN,(MSE_TR),'b-+')
% % hold on
% % plot(MN, (PSU_NEW),'r-*')
% xlabel('撖潮��踹漲瘥�??嚙賜�瑟')
% ylabel('mse')

% save('LASSO\10db\1000\MSE010','MSE')
% save('LASSO\10db\1000\PSU010','PSU')
% save('LASSO\10db\1000\t010','t')
