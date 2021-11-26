function [FL1, FL2, FL3, FL1_std, FL2_std, FL3_std] = ft_std()
    start_path = pwd;
    uiwait(msgbox('Please navigate to the .h5 file to plot.',...
                    'Loading files...','help'));
    [baseFileName, folder] = uigetfile('*.h5');
    cd(folder);
    F1 = h5read(baseFileName, '/NEO/ParticleData/Xe1_FluorPeak');
    F2 = F1(2, :);
    F1 = F1(1,:);
    F3 = h5read(baseFileName, '/NEO/ParticleData/Xe2_FluorPeak');
    F3 = F3(2,:);
    
    FL1 = mean(F1);
    FL2 = mean(F2);
    FL3 = mean(F3);
    
    FL1_std = std(F1);
    FL2_std = std(F2);
    FL3_std = std(F3);
    
    
    cd (start_path);
end