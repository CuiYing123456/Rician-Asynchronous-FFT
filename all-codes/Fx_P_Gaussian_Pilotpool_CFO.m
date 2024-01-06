function [ P_array, P_matrix,P_array_CFO, P_matrix_CFO, P_orig, MatrixtriuPPH, LineL,CFO] = Fx_P_Gaussian_Pilotpool_CFO( T_g , L, N, CFOrange )
% This function generates the effective pilots with uniformly distributed
% CFOs and STOs


%input:
%T_g   maximum STO $D$
% N    number of devices
%L      pilot length
% CFOrange \in [0,2*pi]

%output:

%  P_orig   preassigned pilots
%P_array   3-dimensional-array effective pilots which only suffer STO              
%P_matrix    matrix consisting of effective pilots which only suffer STO
% P_array_CFO  3-dimensional-array effective pilots which suffer both STO and CFO 
%P_matrix_CFO matrix consisting of effective pilots which suffer both STO and CFO 
% MatrixtriuPPH, LineL : these two quantaties are generated for quickly
                          % computing $ \boldsymbol  \Xi_i (\mathbf X , \mathbf p ) $ in the proposed algortihm and are not used in generating the received signal, and the
                          % details of quickly computing $ \boldsymbol  \Xi_i (\mathbf X , \mathbf p ) $ are omitted here.
% CFO         generated CFO
 
 

L_tilde = L + T_g ;
CFO =   CFOrange* rand(N,1)-  CFOrange /2 ;   %  2*pi/ L_tilde*  ceil( L_tilde*rand(N,1)) ;  
% CFO =  zeros(N,1); %2*pi/ 2^8 *ceil( 2^8*rand(N,1)-1 ) ;

% P =normrnd(0,sqrt(1/2),[L N]) + 1i * normrnd(0,sqrt(1/2),[L N]); %未归一化


% for nn = 1:N
%     P(:,nn)= sqrt(L)*  P (:,nn)/( norm(P (:,nn))   );
% end

P_array = zeros(    L_tilde  ,T_g +1,N);
P_matrix= zeros(    L_tilde  ,N*(T_g +1) );

P_array_CFO = zeros(    L_tilde  ,T_g +1,N);
P_matrix_CFO = zeros(    L_tilde  ,N*(T_g +1) );

CFO_Matrix = zeros(    L_tilde,N      );


P =normrnd(0,sqrt(1/2),[L N]) + 1i * normrnd(0,sqrt(1/2),[L N]); %未归一化


for nn = 1:N
    P(:,nn)= sqrt(L)*  P (:,nn)/( norm(P (:,nn))   );
    
    for l = 1: L_tilde
        CFO_Matrix(l,nn)  =  exp(1j*(l-1)*CFO(nn)   );
    end
    
end
P_orig =P;




for   n  =1 :N
    
    for  tau = 1: T_g+1
        P_array(  tau: tau+L-1  ,tau   ,  n  ) = P(:,n);
        P_array_CFO(  :  ,tau   ,  n  ) = diag( CFO_Matrix(:,n) ) *   P_array(  :  ,tau   ,  n  );
        
    end
    
    for j = (n-1)*(T_g+1) +  1:  (n-1)*(T_g+1) +  T_g+1
        
        P_matrix( j - (n-1)*(T_g+1)     : j - (n-1)*(T_g+1) +L-1  ,   j         ) = P(:,n);
        P_matrix_CFO(:,j) =  diag( CFO_Matrix(:,n) ) *P_matrix(:,j);
    end
end



%%   consturct the magic matrix

L = L_tilde;  
lineL = 1:L^2;
LineL = reshape(lineL, L,L  );
LineL = triu(LineL);


for i = 1 :L
    LineL(:,i) = circshift(  LineL(:,i), L-i );
end
LineL = flip( LineL);
L2tril = tril( (2)* ones(L) );  %    due to this step, an extra step which sets the element of PPH indexed by '2' shuold be set to zero
L2tril = L2tril - 2*eye(L);
LineL =   LineL + L2tril;
LineL = LineL.';


%%  Compute PP^H

MatrixtriuPPH = zeros(L,L,  1+T_g,  N );
MatrixPPH = zeros(L,L,  1+T_g, N);
LinePPH = zeros(L^2,  1+T_g  ,N);






for nn = 1: N
    for t =1 : T_g +1
        P  =  P_array(  : , t  ,  nn  );
        conjP =  conj(P );
        
        PPH  = conjP*conjP';
        
        MatrixPPH(:,:,t,nn) =     PPH;
        
        PPH(2,1) = 0;       % this is for the time saving in the inner product of invsigma and SigbarSig
        Linepph =  reshape(  PPH, L^2,1 );
        LinePPH(:,t,nn) = Linepph;
        
        MatrixtriuPPH(:,:,t,nn) =     Linepph( LineL);
        
    end
end


