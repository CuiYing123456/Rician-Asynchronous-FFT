function [PSU,  PdataU , AandDest] = Hard_decision_opt_thresh( X_est,  A  , AandDreal  ,Indpilot_NEW)

Num = size( X_est );






if     nargin == 4    % the constraints on one-data transimmsion are not released
    %     strcmp(type, 'R')
    
  Q  = size( Indpilot_NEW,2);
    N = Num(2);
   
        
    Loop = Num(1);
  
    
    A_est = zeros(    Loop,  N );
    Delay_est = zeros(  Loop,  N);
    
    B_est = zeros(N,  Q,    Loop);
    
    th1 = 0.01:0.01:1.5;
    th2 = [1e-5,1e-4,1e-3,1e-2,1e-1];
    th = [th2, th1];
    
    PS = zeros(length(th),1);
    for thres = 1:length(th)
        for ii=1:  Loop
            for jj = 1: N
                if(abs(X_est(ii,jj)) > th(thres))
                    A_est(ii,jj)=1;
                else
                    A_est(ii,jj)=0;
                end
            end
        end
        
        PS(thres) = sum(sum(1-abs(A_est - A)))/  Loop/ N;
        
    end
    [aaaa, idx] = max(PS);
    PSU = 1-aaaa;
    optimal_thres = th(idx);
    %     optimal_thres = 0.3 ;
    A_est = X_est>optimal_thres;
    
    %          ErDatNum =0;
    AandDest = zeros(  Loop, N);
    for ii=1:   Loop
        for jj = 1: N
            
            if abs(X_est(ii,jj)) > optimal_thres
                A_est(ii,jj)=1;
                AandDest(ii,jj ) = Indpilot_NEW(ii,jj);
                
            else
                A_est(ii,jj)=0;
                AandDest(ii,jj ) =0;
            end
            
            
            
            
        end
        %              ErDatNum =  ErDatNum +  sum(sum(abs((  B_est(:,:,ii) - DB(:, :,ii   )))));
    end
    
    PdataU  = 1- length( find(  AandDest ==  AandDreal  ) )/  Loop/ N   ;
    %         AandDest_Prop(:,:,loop_tg) =  AandDest;
    
    
    
    
    
elseif  nargin == 3    % the constraints on one-data transimmsion are released
    %%  2.1 Yu Wei   ML
    
    Q = Num(2);
    N = Num(1);
  if length(Num ) ==2
        Loop =1;
    else
        
    Loop = Num(3);
    end
    
    
    A_est = zeros(N, Q ,  Loop);
    Av_est = zeros(Loop, N);
    Delay_est = zeros(Loop,  N);
    
    
    th1 = 0.01:0.01:1.5;
    th2 = [1e-5,1e-4,1e-3,1e-2,1e-1];
    th = [th2, th1];
    
    PS = zeros(length(th),1);
    for thres = 1:length(th)
        for ii=1:Loop
            for jj = 1:N
                [ maxXest , indexXest  ] = max (  abs( X_est(jj,:,ii) ));
                if   ( maxXest > th(thres))
                    A_est(jj,indexXest ,ii)=1;
                    Av_est(ii,jj) = 1;
                    Delay_est(ii,jj) =  indexXest-1;
                else
                    A_est(jj,indexXest ,ii)=0;
                    Av_est(ii,jj) = 0;
                    Delay_est(ii,jj) =  indexXest-1;
                end
            end
            
        end
        
        
        PS(thres) = sum(sum(1- abs( Av_est - A)))/N/Loop;
        
    end
    [ aaaa, idx] = max(PS);
    PSU  = 1-aaaa;
    optimal_thres = th(idx);
    
    
    %         ErDatNum =0;
    AandDest = zeros(Loop,N);
    for ii=1:Loop
        for jj = 1:N
            [ maxXest , indexXest  ] = max (  abs( X_est(jj,:,ii) ));
            if    (maxXest > optimal_thres)
                A_est(jj,indexXest ,ii)=1;
                Av_est(ii,jj) = 1;
                Delay_est(ii,jj) =  indexXest;
                AandDest(ii,jj ) =  indexXest;
                
            else
                A_est(jj,indexXest ,ii)=0;
                Av_est(ii,jj) = 0;
                Delay_est(ii,jj) =  indexXest;
                AandDest(ii,jj ) =  0;
                
            end
        end
        %             ErDatNum =  ErDatNum +  sum(sum(abs((  A_est(:,:,ii) - DB(:, :,ii   )))));
        
    end
    
    
    PdataU = 1-length( find(  AandDest == AandDreal   ) )/Loop/N     ;
    
    %       AandDest_ML(:,:,loop_tg) =  AandDest;
    
    
  elseif  nargin == 2    % the constraints on one-data transimmsion are released  
      
     %%  2.1 Yu Wei   ML    for case-f
    
    Q = Num(2);
    N = Num(1);
    
    if length(Num ) ==2
        Loop =1;
    else
        
    Loop = Num(3);
    end
    
    AandDest = zeros(Loop,N);  
    PdataU  = 1;
    
    
    A_est = zeros(N, Q ,  Loop);
    Av_est = zeros(Loop, N);
    Delay_est = zeros(Loop,  N);
    
    
    th1 = 0.01:0.01:1.5;
    th2 = [1e-5,1e-4,1e-3,1e-2,1e-1];
    th = [th2, th1];
    
    PS = zeros(length(th),1);
    for thres = 1:length(th)
        for ii=1:Loop
            for jj = 1:N
                [ maxXest , indexXest  ] = max (  abs( X_est(jj,:,ii) ));
                if   ( maxXest >  th(thres) )
                    A_est(jj,indexXest ,ii)=1;
                    Av_est(ii,jj) = 1;
                    Delay_est(ii,jj) =  indexXest-1;
                else
                    A_est(jj,indexXest ,ii)=0;
                    Av_est(ii,jj) = 0;
                    Delay_est(ii,jj) =  indexXest-1;
                end
            end
            
        end
        
        
        PS(thres) = sum(sum(1- abs( Av_est - A)))/N/Loop;
        
    end
    [ aaaa, idx] = max(PS);
    PSU  = 1-aaaa;
    optimal_thres = th(idx);
    
    
else
    error('error');
    
end