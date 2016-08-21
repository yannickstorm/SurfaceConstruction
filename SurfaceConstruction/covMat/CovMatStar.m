function covMatStar = CovMatStar(Prior,X1,X2)
switch Prior.kernel
    case 'TP'
    covMatStar = CovMatStarTP(Prior,X1,X2);
    
    case 'SE'
    covMatStar = CovMatStarSE(Prior,X1,X2);
    
    otherwise
        error('The kernel defined in Prior does not correspond to any of the possible ones.')
end

end

