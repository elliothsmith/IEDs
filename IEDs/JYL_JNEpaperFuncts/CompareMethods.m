function [ DD, DS, DH, pD, pS, pH, dD, dS, DetectRate ] = CompareMethods( Comp, Ref, Per  )
% [ dD, dS ] = CompareMethods( Comp, Ref, Per  )
% 
% DD = Direction error
% DS = Speed error
% DH = XOR, whether they agree or not
% pD = paired circular test (cm-test), see whether the two has the same
%      median direction or not
% pS = signrank test
% pH = McNemar test
% dD, dS, dH = depends on the percentage
%
% If you want to use old comparison method, use 'SpatialError.m'
%
% Jyun-you Liou 2016/09/14

[M,N] = size(Comp);
Speed = @(x) sqrt(sum(x.^2,2));

if nargin<3;Per = 0.5;end

for m = M:-1:1
    for n = N:-1:1
        Sel = Ref.Significance; 
        % Sel = Ref.Significance & Comp(m,n).Significance;       
        Vcomp = Comp(m,n).V(Sel,:);
        Vref = Ref.V(Sel,:);
        DS{m,n} = Speed(Vcomp) - Speed(Vref);
        pS(m,n) = signrank(Speed(Vcomp),Speed(Vref));
        dS(m,n,:) = quantile(DS{m,n},Per);
        
        DD{m,n} = acos( sum(Vcomp .* Vref,2) ./ Speed(Vcomp) ./ Speed(Vref));
        Dcomp = Comp(m,n).direction(Sel,:);
        Dref = Ref.direction(Sel,:);        
        pD(m,n) = circ_cmtest(Dcomp, Dref);
        dD(m,n,:) = quantile(DD{m,n},Per);        
        
        DH{m,n} = ~xor(Comp(m,n).Significance, Ref.Significance);
        pH(m,n) = testcholdout(Comp(m,n).Significance,Ref.Significance,true(numel(Ref.Significance),1));
        DetectRate(m,n) = sum(Comp(m,n).Significance)/numel(Comp(m,n).Significance);
    end
end


end
