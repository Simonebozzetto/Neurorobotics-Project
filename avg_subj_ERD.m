function subj_avg_cal_ERD = avg_subj_ERD(cel_z_cal_ERD,n_subjects)
    % averaging among subjects for each ERD 
    % before averaging among subject, we have again to check the minimum length 
    % of the samples, and cut the longer ones
    cut_cel_z_cal_ERD = cut_samples_concatenate(cel_z_cal_ERD,0);
    
    
    % averaging among subjects
    subj_avg_cal_ERD = zeros(size(cut_cel_z_cal_ERD{1}));

    for i = 1:n_subjects
        subj_avg_cal_ERD = subj_avg_cal_ERD + cut_cel_z_cal_ERD{i};
    end
    clear i
    subj_avg_cal_ERD = subj_avg_cal_ERD/n_subjects;
end

