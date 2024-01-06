function [y, idx ,w, time ] = Prop_FT_Ri_Fourcases_ScfoSsto(A, sampCov ,  Y, sigma2,  Kappa,  lsfc, h_bar, N , Omega, Loop_search)

% This function is the proposed Prop-MLE-S. In this version:
%1. We carefully designed a process to quickly compute the function 
%   $\boldsymbol  \Xi_i (\mathbf X , \mathbf p )$.





[L, Nall] = size(A);
Q = Nall/N;

[~,M] = size(h_bar);

L_tilde  = L+0;

gamma = zeros(N,1);

omega = zeros(N,1);
Indpilot = ones(N,1);
Indw  = ones(N,1);

n_var = sigma2;
inv_Sigma= eye(L) /n_var;


Sigma_bar = sampCov;

Y_tilde = Y;




 if  Omega == Loop_search/2
        T_gsamlple  = 1: 1: Loop_search;
         else
          T_gsamlple = [ 1:Omega+1, (Loop_search +1- Omega):1:(Loop_search)  ];
  end
  L_T_gsamlple = length(T_gsamlple);


%  1.
%                K_power = nextpow2(L_tilde);

% Loop_search = 2^8;%  2^  K_power;
%  2.
% numpilot = Q;
% Loop_search =   L_tilde +  (  numpilot  -mod( L_tilde, numpilot ));




%%

for  n =1:N
    B(:, (n-1)*Q+1 : (n-1)*Q+Q) = sqrt( lsfc(n)./(1 + Kappa(n)  )) * A(:, (n-1)*Q+1 : (n-1)*Q+Q)  ;
end


Y_tilde = Y;


%%   consturct the tridline
tridline = [];
for i = 1 :L
    tridline = [tridline , (i-1)*L+i  :  i*L];
end
tridline =  tridline';

%%   consturct triuline

L = L_tilde;
lineL = 1:L^2;
LineL = reshape(lineL, L,L  );
LineL = triu(LineL);
LineL_for_triu = LineL;



%%   consturct the triuline
triuline = [];
for i = 1 :L
    triuline = [triuline , diag(LineL_for_triu,i-1)'];
end
triuline =  triuline';



triulinePPH = zeros(L*(L+1)/2, Q  ,N);
for nn = 1: N
    for t =1 : Q
        P  =    B(  : ,  Q*(nn-1) +t );
        conjP =  conj(P );
        
        PPH  = conjP*conjP';
        triulinePPH(:, t,nn) = PPH( triuline );  %clounm vector
        
    end
end




n_var = sigma2;
inv_Sigma= eye(L) /n_var;
Sigma_bar = sampCov;
SigbarSig =  inv_Sigma *  Sigma_bar * inv_Sigma ;

SigbarSig_back_triuline  = SigbarSig (triuline) ;
inv_Sigma_back_triuline  =     inv_Sigma(triuline) ;



%%




time = 0;
loopbcd = 0;


totalincrement  = 0;
Delta_totalincrement = 10;

inv_sigmaPPH = zeros(L,L );
SigbarSigPPH = zeros(L,L );


epsilon = 1e-8;

invinv_Sigma = inv(inv_Sigma);
flike0 = log( det( invinv_Sigma )) + 1/M * trace( inv_Sigma* (Y_tilde  *  Y_tilde'));

Indupdate_sigbarsig= 0;
val = 0;
while   loopbcd < 1000
    loopbcd =   loopbcd + 1;
    %     g = grad(gamma, A, sampCov, sigma2);
    
    %     if norm(max(gamma - g, 0) - gamma) < epsilon
    %         break
    %     end
    
    if  abs(Delta_totalincrement)  <  abs(totalincrement) * epsilon
        %       if abs(Delta_totalincrement)  <  abs(LastDelta) * epsilon
        break
    end
    
    TQ =  dftmtx( Loop_search )';
    TQ =  TQ(1:L,:);
    wt  =  zeros(Q,1);
    MinInct  =  zeros(Q,1);
    dt  = zeros(Q,1);
    indw  = zeros(Q,1);
    
    c1sc  = zeros( 1,L_T_gsamlple );c2sc  = zeros(1, L_T_gsamlple);c3sc  = zeros( 1,L_T_gsamlple);
    
    tic;
    
    last_tot_incre= totalincrement;
    %        LastDelta = Delta_totalincrement;
    
    
    for j = 1:N
        
        
        
        Indwold = Indw (j);
        
        IndPold = Indpilot(j);
        omegaold = omega(j);
        older_gamma  = gamma(j);
        pary = B(:,  (j-1)*Q + IndPold );
        %       p_iterate =    exp( 1j* omega(j)*(0: L -1)'     ) .* pary;
        para_k_n = Kappa(j) ;
        hbar = h_bar(j,:);
        p_iterate =   TQ(:,Indwold ) .* pary;
        
        %%    compute back, corresponds to steps 5, 6
        if  gamma(j) ==0
            inv_Sigma_back = inv_Sigma;
            
            
            Y_tilde_back =   Y_tilde;
            
            
            
            
            
        else
%             Indupdate_sigbarsig= Indupdate_sigbarsig+1;
            vc1 = inv_Sigma* p_iterate;
            
            c1im =  (p_iterate' *vc1 );
            c1  = real(c1im);
            c1c1 = ( vc1*vc1');
            
            factorc1c1 = (0 -older_gamma  )/( 1+  (0- older_gamma  )   *c1)  ;
            
            inv_Sigma_back =  inv_Sigma -   factorc1c1 * c1c1 ;
            
            %               update  Y_tilde_back
            Y_tilde_back = Y_tilde + ( older_gamma * p_iterate)*  hbar;
            
            
            
            
            
        end
        
              
%         if  Indupdate_sigbarsig~=0
%             AA= inv_Sigma_back   *Y_tilde_back;
%             SigbarSig_back  =1/M* (AA*AA');
%             
%             SigbarSig_back_triuline  =    SigbarSig_back(triuline) ;
%             
%             
%             inv_Sigma_back_triuline  =     inv_Sigma_back(triuline) ;
%         else
%             
%         end
        
        
        
        
        %    fft
        
        
%         
%         Yh =    Y_tilde_back*  hbar';
%         SigYh =  inv_Sigma_back*Yh;
        for tt = 1:Q
            
            
%             QQtriuline = triulinePPH(:,tt,j);
%             pary =   B(:,  (j-1)*Q + tt);
%             
%             inv_Sigma_backdot = inv_Sigma_back_triuline.*QQtriuline;
%             SigbarSig_backdot =    SigbarSig_back_triuline.*QQtriuline;
%             
%             inv_sigmaPPH(tridline) =  inv_Sigma_backdot;
%             SigbarSigPPH (tridline)=  SigbarSig_backdot;
%             
%             
%             VMC1_up=   sum(  inv_sigmaPPH  );
%             VMC2_up =   sum( SigbarSigPPH   );
%             
%             VMC1_up(1) = VMC1_up(1)/2;
%             VMC2_up(1) = VMC2_up(1)/2;
%             
%             VMC3_up = conj( pary).* SigYh ;
%             VMC3_up=   VMC3_up';
%             
%             
%             C1_sample =  2*real ( Loop_search* ifft( VMC1_up  ,Loop_search    ));
%             C2_sample = 2*real ( Loop_search* ifft( VMC2_up  ,Loop_search    ));
%             C3_sample = 2/M*( real (  Loop_search* ifft( VMC3_up  ,Loop_search    )));
%             
%             
%             %         T_gsamlple  = 1: Loop_search/( numpilot) : Loop_search;
% %             T_gsamlple  = 1: 1: Loop_search;
%             C1_samp0=  C1_sample(T_gsamlple);
%             C2_samp0=  C2_sample(T_gsamlple);
%             C3_samp0=  C3_sample(T_gsamlple);

            %%   SEPARATELY COMPUTE,  corresponds to steps 7-11
             pary =   B(:,  (j-1)*Q + tt);
         for ww = 1:length(T_gsamlple)
          www =  T_gsamlple(ww);
          p_ary_tt_ww =    TQ(:, www) .* pary;
          vcsc = inv_Sigma_back*  p_ary_tt_ww ;
            
            c1imsc =  ( p_ary_tt_ww' *vcsc );
            c1sc(ww)  = real(c1imsc);

            YSP = Y_tilde_back'*vcsc;
            c2sc(ww) =  (YSP'* YSP)/M;

            c3sc(ww) =  2*real(hbar* YSP)/M;

         end
            C1_samp =  c1sc;
            C2_samp =  c2sc;
            C3_samp =  c3sc;
%          trace(  (   c1sc  - C1_samp0)*(   c1sc  - C1_samp0)' + ... 
%               (   c2sc  - C2_samp0)*(   c2sc  - C2_samp0)' + ...
%                (   c3sc  - C3_samp0)*(   c3sc  - C3_samp0)')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            d_zeroderev = -1/2/  para_k_n*ones(1, L_T_gsamlple) -  1./C1_samp + 1/2/  para_k_n*  ...
                sqrt( max ( 1 + 4*para_k_n *(para_k_n + C2_samp + C3_samp)./(C1_samp).^2 ,  zeros(   1, L_T_gsamlple) )   ) ;
            d_samp =  min( max (   zeros(   1, L_T_gsamlple),    d_zeroderev         ) ,ones(1, L_T_gsamlple)); %d_samp shuold be the same size of C1_samp
            %           d_samp =    max (   zeros(   1, Loop_search), ( C2_samp -  C1_samp)./  C1_samp.^2         )  ;
            
            fml_samp = log(1 + d_samp.*  ( C1_samp) ) +   ( para_k_n* d_samp.^2 .*C1_samp -  d_samp.*  ( C2_samp + C3_samp  )   )./  (1 +d_samp.*   C1_samp );
            
            
            
            
            [ Minincre , ind] =  min (    fml_samp);
            Indpi = ind(1);
            dd =    d_samp(Indpi );
            
            wt(tt) =  ( T_gsamlple(Indpi)-1) *2*pi/Loop_search ;
            MinInct(tt) = Minincre;
            dt(tt) = dd;
            indw(tt)= T_gsamlple(Indpi);

          


            
        end
        
        [~,idt] = min(    MinInct) ;
        idt = idt(1);
        gamma(j) = dt(idt);
        omega(j) =   wt(idt);
        Indpit = idt;
        Indpilot(j) = idt;
        Indw(j) =   indw(idt);
        dd = gamma(j);
        
        %%    Variable update,  corresponds to steps 12, 13
%         Indupdate_sigbarsig =   0;
        if  dd==0
            inv_Sigma =inv_Sigma_back ;
           

            
            Y_tilde  = Y_tilde_back  ;
            
            
        elseif Indpit==IndPold &&   older_gamma~=0   &&      omegaold ==  omega(j)
            
            factorc1c1 = (dd - older_gamma  )/( 1+ (dd - older_gamma)   *c1   );
            inv_Sigma =  inv_Sigma -  factorc1c1 *  c1c1;
            
            
            
            Y_tilde  = Y_tilde_back  - dd *  p_iterate*  hbar;
            
            
%             Indupdate_sigbarsig= Indupdate_sigbarsig+1;
            
        else
            
            
            p_ary =    TQ(:, Indw(j)) .* B(:,  (j-1)*Q +          Indpit );
            
            
            vc1 =inv_Sigma_back *  p_ary;
            c1= p_ary'*  vc1  ;% no simplify
            c1   =real( c1);
            
            c1c1 =(vc1*vc1');
            factorc1c1 =(dd )/( 1+dd*c1);
            
            inv_Sigma =inv_Sigma_back - factorc1c1*  c1c1;
            Y_tilde = Y_tilde_back - dd * p_ary * hbar;
            
%             Indupdate_sigbarsig= Indupdate_sigbarsig+1;
            
            
            
            
            
            
            
        end
        
        
        
        
    end
    
    
    invinv_Sigma = inv(inv_Sigma);
    totalincrement = real( flike0 - log( det( invinv_Sigma )) - 1/M * trace( inv_Sigma* (Y_tilde  *  Y_tilde')));
    
    
    %             val = [val, totalincrement];
    
    
    
    
    Delta_totalincrement = totalincrement - last_tot_incre ;
    
    
    
    
    
    time  = time  + toc;
end

fprintf('Prop_FT_Ri_ScfoSsto:A total of %d iterations were run, ',  loopbcd ) ;
fprintf('costs  %d  seconds \n',  time ) ;

% close all
% plot(1:length(val), val)
y = gamma;
idx =  Indpilot;
w = omega;
end