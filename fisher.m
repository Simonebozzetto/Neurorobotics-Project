function [Fisher,P,Pk,P_2D] = fisher(PSD,Ck,CFbK,cue_labels)
    P = PSD(CFbK == 781, :, :);  
    Pk = Ck(CFbK==781);

    n_wind = size(P,1);
    n_frequencies = size(P,2);
    n_channels = size(P,3);
    
    P_2D = reshape(P, [n_wind, n_frequencies * n_channels]);
    
    P_771 = P_2D(Pk==cue_labels(2),:);
    P_773 = P_2D(Pk==cue_labels(1),:);
    
    avg_P_771 = mean(P_771,1);
    avg_P_773 = mean(P_773,1);
    
    var_P_771 = var(P_771,0,1);
    var_P_773 = var(P_773,0,1);
    
    Fisher_score = ((avg_P_771 - avg_P_773).^2)./(var_P_771 + var_P_773);

    Fisher = reshape(Fisher_score,[n_frequencies,n_channels]);
end

