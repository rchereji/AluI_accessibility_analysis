function Ratios = CorrectRatios_sacCer3(Ratios, Sites_with_problems_Right, Sites_with_problems_Left, Sites_with_problems_Both_Sides)
% Function to correct the problem of missing short DNA fragments. If 2
% cleavage sites are close and the resulting DNA fragment is shorter than
% 50 bp, then this fragment will not be sequenced (or adapters will
% compromise the read). To correct this problem, if the left side of an
% AluI site has this problem, the right side estimation is used, and
% vice-versa. If both sides have this problem (close AluI sites on both
% sides) then this site is discarded, as the estimation of cut ratio is not
% reliable.

noChr = numel(Ratios);
for chr = 1:noChr
    ind = strfind(Sites_with_problems_Left{chr}, [1,0]);
    Ratios{chr}(ind) = Ratios{chr}(ind+1);
    
    ind = strfind(Sites_with_problems_Right{chr}, [0,1]) + 1;
    Ratios{chr}(ind) = Ratios{chr}(ind-1);
    
    ind = strfind(Sites_with_problems_Both_Sides{chr}, [1,1]);
    Ratios{chr}(ind) = NaN;
    Ratios{chr}(ind+1) = NaN;
end
