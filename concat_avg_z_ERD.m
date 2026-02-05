function z_cal_ERD = concat_avg_z_ERD(cal_ERD)
    n_cal_runs = size(cal_ERD,2);    

    length = zeros(n_cal_runs,1);
    for i = 1:n_cal_runs
        length(i) = size(cal_ERD{i},1);
    end
    clear i
    min_length = min(length);
    clear length

    % concatenation over files
    cal_ERD_concat = [];
    for i = 1:n_cal_runs
        temp = cal_ERD{i}(1:min_length,:,:);
    
        cal_ERD_concat = cat(3, cal_ERD_concat, temp);
    end
    clear i
    clear temp
    
    % averaging over trials
    avg_cal_ERD  = mean(cal_ERD_concat,3);
    
    % z-score 
    %z_cal_ERD  = zscore(avg_cal_ERD);
    z_cal_ERD = avg_cal_ERD;
end

