function [y, idx ,w, time ] = Prop_FT_Ri_Fourcases_LcfoLsto(A, A_orig, sampCov ,  Y, sigma2,  Kappa,  lsfc, h_bar,  Omega, Loop_search)

% This function is the proposed Prop-MLE-L. In this version:
%1. We carefully designed a process to quickly compute the function $\boldsymbol  \Xi_i (\mathbf X
% , \mathbf p )$. 
%2. We simply use the basic matrix-matrix multiplication (O(ML^2)), rather than
% the propsed matrix-vector multiplication and matrix-matrix addition in our paper (O(L^2)), to compute the updation of
%  the matrix
  %${\bf\Sigma}_{\text{t},n }^{-1}   \widetilde {\mathbf Y}_{\text{t} ,n}    \widetilde {\mathbf Y}_{\text{t},n }^H      {\bf\Sigma}_{\text{t} ,n}^{-1}$
 % This is because we find that in the simulation setup of our paper, the basic
 % matrix-matrix multiplication (O(ML^2)) costs almost the same time as 
 %the propsed matrix-vector multiplication and matrix-matrix addition in our paper (O(L^2)), 
 %and the basic matrix-matrix multiplication is easier to debug.

 
% input
% lsfc   large-scale fading power
% h_bar  normalized LoS
 
 
[L, Nall] = size(A);
[L0, N] = size(A_orig);
Q = Nall/N;

[~,M] = size(h_bar);

STO = L-L0;

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

 
for  n =1:N
    A_orig(:,n) = sqrt( lsfc(n)./(1 + Kappa(n)  )) * A_orig(:,n) ;
A(:, (n-1)*Q+1 : (n-1)*Q+Q) = sqrt( lsfc(n)./(1 + Kappa(n)  )) * A(:, (n-1)*Q+1 : (n-1)*Q+Q)  ;
end





%%   consturct the tridline
tridline = [];
for i = 1 :L
    tridline = [tridline , (i-1)*L+i  :  i*L];
end
 tridline =  tridline';  % column vector
 
 %%   consturct the tridline_flipcir
tridline2 = [];
for i = 1 :L
    tridline2 = [tridline2,(i-1)*L+1  :  i*L-i+1];
end
 tridline2 =  tridline2';  % column vector
 
 tridline_flipcir = [];
for i = 1 :L0
    tridline_flipcir = [tridline_flipcir,  (i-1)*L+1,   i* L :-1:(i-1)*L+ STO+i+1 ];
end
 tridline_flipcir =  tridline_flipcir';  % column vector
 
 
%%   consturct the magic matrix
lineL = 1:L0^2;
LineL = reshape(lineL, L0,L0  );
LineL = triu(LineL);
LineL_for_triu = LineL;
 triulineL0 = [];
for i = 1 :L0
   triulineL0 = [ triulineL0 , diag(LineL_for_triu,i-1)'];
end
 triulineL0 =   triulineL0'; % column vector

 
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
 
 
 
 
%%  construct the pilot
 PPHfft_N = zeros(L,L0,N);
 FPsiP = zeros(L,N);
for n = 1: N
     PPHcc  = zeros(L,L ); 
     p_orig = conj(  A_orig(:,n)) ;
     PPH    =         (p_orig*        p_orig') ;
     PPH_lineL0  =   PPH( triulineL0);
     PPHcc( tridline_flipcir)= PPH_lineL0 ;
      PPHcc(:,1) = 1/2*PPHcc(:,1);
     PPHfft_N(:,:,n) =  fft(PPHcc(:,1:L0));
      
%      PsiP  =  [ conj(A_orig(:,n)); zeros(STO,1)];
 PsiP  =  [  conj(A_orig(:,n)); zeros(STO,1)];
     PsiP(2:end ) = flipud( PsiP(2:end));
  FPsiP(:,n) =   fft ( PsiP);
     
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
        P  =   A(  : ,  Q*(nn-1) +t );
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

inv_Sigma_backcc = zeros(L,L );
inv_Sigma_backcc(tridline2)=     inv_Sigma(triuline) ;
fftinv_Sigma_backcc =fft(  inv_Sigma_backcc(:, 1:L0));

SigbarSig_backcc = zeros(L,L );
SigbarSig_backcc(tridline2) =  SigbarSig(triuline) ;
fftSigbarSig_backcc =fft( SigbarSig_backcc(:, 1:L0)); 
 


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
    
    
    tic;
    
    last_tot_incre= totalincrement;
    %        LastDelta = Delta_totalincrement;
    
    
    for j = 1:N
        
        
        
        Indwold = Indw (j);
        
        IndPold = Indpilot(j);
        omegaold = omega(j);
        older_gamma  = gamma(j);
        
        para_k_n = Kappa(j) ;
        hbar = h_bar(j,:);
        p_iterate =   TQ(:,Indwold ) .* A(:,  (j-1)*Q + IndPold );
        
        %%    compute back, corresponds to setps 5-7, but step 7 is replaced with a 
                 %process with similar computation time
        if  gamma(j) ==0
            inv_Sigma_back = inv_Sigma;
            
            
            Y_tilde_back =   Y_tilde;
            
            
            
            
            
        else
            Indupdate_sigbarsig= Indupdate_sigbarsig+1;
            vc1 = inv_Sigma* p_iterate;
            
            c1im =  (p_iterate' *vc1 );
            c1  = real(c1im);
            c1c1 = ( vc1*vc1');
            
            factorc1c1 = (0 -older_gamma  )/( 1+  (0- older_gamma  )   *c1)  ;
            
            inv_Sigma_back =  inv_Sigma -   factorc1c1 * c1c1 ;
            
            %               update  Y_tilde_back
            Y_tilde_back = Y_tilde + ( older_gamma * p_iterate)*  hbar;
            
            
            
            
            
        end
        
        
      if  Indupdate_sigbarsig~=0
         AA= inv_Sigma_back   *Y_tilde_back;
         SigbarSig_back  =1/M* (AA*AA');
        SigbarSig_backcc(tridline2) =  SigbarSig_back(triuline) ;
        fftSigbarSig_backcc =fft( SigbarSig_backcc(:, 1:L0)); 
        
        
      inv_Sigma_backcc(tridline2)=     inv_Sigma_back(triuline) ;
        fftinv_Sigma_backcc =fft(  inv_Sigma_backcc(:, 1:L0));
    else

    end
        
        
        
        
        %%    fft,   corresponds to steps 8, 9
         PPHfft = PPHfft_N (:,:,j); 
        
        Yh =    Y_tilde_back*  hbar';
        SigYh =  inv_Sigma_back*Yh;
        VMC3_up = conj( A(:,Q*(j-1)+1 :  Q*(j-1)+Q  )).* SigYh ;
        VMC3_up=   conj(VMC3_up);

  
        
      ccInv_SigPPH   = ifft( fftinv_Sigma_backcc.*  PPHfft)     ;
      ccSigbarSigPPH   = ifft( fftSigbarSig_backcc.* PPHfft )     ;
      
      ccInv_SigPPHo   =   ccInv_SigPPH( 1:Q,:).'   ;

      ccSigbarSigPPHo  = ccSigbarSigPPH( 1:Q,:) .'    ;     
   

        
     C1_samp=   2*real ( Loop_search* ifft(ccInv_SigPPHo  ,Loop_search    )); 
     C2_samp=   2*real ( Loop_search* ifft( ccSigbarSigPPHo  ,Loop_search    ));   
     C3_samp= 2/M*( real (  Loop_search* ifft( VMC3_up  ,Loop_search    ))); 

     C1_samp =  C1_samp (T_gsamlple,:);
      C2_samp =  C2_samp (T_gsamlple,:);
       C3_samp = C3_samp (T_gsamlple,:);

        
    d_zeroderev = -1/2/  para_k_n*ones(   L_T_gsamlple,Q) -  1./C1_samp + 1/2/  para_k_n*  ...
                sqrt( max ( 1 + 4*para_k_n *(para_k_n + C2_samp + C3_samp)./(C1_samp).^2 ,  zeros(     L_T_gsamlple,Q) )   ) ;
            d_samp =  min( max (   zeros(    L_T_gsamlple,Q),    d_zeroderev         ) ,ones(  L_T_gsamlple,Q)); %d_samp shuold be the same size of C1_samp
            %           d_samp =    max (   zeros(   1, Loop_search), ( C2_samp -  C1_samp)./  C1_samp.^2         )  ;
            
            fml_samp = log(1 + d_samp.*  ( C1_samp) ) +   ( para_k_n* d_samp.^2 .*C1_samp -  d_samp.*  ( C2_samp + C3_samp  )   )./  (1 +d_samp.*   C1_samp );
            
      
      [minf] = min(min(fml_samp));
      [indcfo ,indsto] = find( fml_samp ==minf);
       dd  =   d_samp( indcfo(1), indsto(1) );
       
         Minincre  = minf ;
         Indpiw = T_gsamlple(indcfo(1));
         Indpit = indsto(1);
         
        gamma(j) = dd;
        omega(j) =  ( Indpiw-1) *2*pi/Loop_search ; 
        Indpilot(j) = Indpit;
        Indw(j) =    Indpiw;
        dd = gamma(j);
      
   
        
        %%    Variable update, corresponds to steps 10-12, but step 12 is replaced with a 
                 %process with similar computation time
        Indupdate_sigbarsig =   0;
        if  dd==0
            inv_Sigma =inv_Sigma_back ;
           

            
            Y_tilde  = Y_tilde_back  ;
            
            
        elseif Indpit==IndPold &&   older_gamma~=0   &&      omegaold ==  omega(j)
            
            factorc1c1 = (dd - older_gamma  )/( 1+ (dd - older_gamma)   *c1   );
            inv_Sigma =  inv_Sigma -  factorc1c1 *  c1c1;
            
            
            
            Y_tilde  = Y_tilde_back  - dd *  p_iterate*  hbar;
            
            
            Indupdate_sigbarsig= Indupdate_sigbarsig+1;
            
        else
            
            
            p_ary =    TQ(:, Indw(j)) .* A(:,  (j-1)*Q +          Indpit );
            
            
            vc1 =inv_Sigma_back *  p_ary;
            c1= p_ary'*  vc1  ;% no simplify
            c1   =real( c1);
            
            c1c1 =(vc1*vc1');
            factorc1c1 =(dd )/( 1+dd*c1);
            
            inv_Sigma =inv_Sigma_back - factorc1c1*  c1c1;
            Y_tilde = Y_tilde_back - dd * p_ary * hbar;
            
            Indupdate_sigbarsig= Indupdate_sigbarsig+1;
            
            
            
            
            
            
            
        end
        
        
        
        
    end
    
    
    invinv_Sigma = inv(inv_Sigma);
    totalincrement = real( flike0 - log( det( invinv_Sigma )) - 1/M * trace( inv_Sigma* (Y_tilde  *  Y_tilde')));
    
    
    %             val = [val, totalincrement];
    
    
    
    
    Delta_totalincrement = totalincrement - last_tot_incre ;
    
    
    
    
    
    time  = time  + toc;
end

fprintf(' Prop_FT_Ri_Fourcases_LcfoLsto:A total of %d iterations were run, ',  loopbcd ) ;
fprintf('costs  %d  seconds \n',  time ) ;

% close all
% plot(1:length(val), val)
y = gamma;
idx =  Indpilot;
w = omega;
end