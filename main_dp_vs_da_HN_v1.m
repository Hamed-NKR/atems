clc
close all
warning('off')

%% initialize dpp vs. da figure %%

figure;
h = gcf;
h.Position = [100, 100, 500, 500];
set(h, 'color', 'white');

% plot universal correlation
r0 = (2e4 / 2e1) ^ (1 / (1e4 - 1));
da0 = 2e1 * ones(1e4, 1) .* r0 .^ (((1:1e4)-1)');
dpp0 = 17.8 * (da0 / 100) .^ (0.35);
plt_0 = plot(da0, dpp0, 'Color', [0.4940 0.1840 0.5560],...
    'LineStyle', '-.', 'LineWidth', 2);
hold on

%% low agglomeration level data %%

fd_im_lal = 'D:\Hamed\CND\PhD\TEM\PFA_new_ET+NIT\26AUG24\PFA_ET+NIT_HAL_20AUG24_R1\19AUG-LAL-Start';

fd_pp_lal = 'D:\Hamed\CND\PhD\TEM\PFA_new_ET+NIT\ImageJ\ET+NIT_19AUG24_LAL_R1\CSV';

f_info_lal = dir([fd_im_lal, '\*.tif']);

n_im_lal = length(f_info_lal);

if ~exist('i_im_lal', 'var') || isempty(i_im_lal) 
    i_im_lal = sort(randi(n_im_lal, 10, 1));
end

if ~exist('imgs_binary_lal', 'var')

    [Imgs_lal, imgs_lal, pixsizes_lal] = tools.load_imgs(fd_im_lal, i_im_lal);
    fname_lal = {Imgs_lal.fname};
    
    imgs_binary_lal = agg.seg_slider(imgs_lal);
    
    Aggs_lal = agg.analyze_binary(imgs_binary_lal,...
        pixsizes_lal, imgs_lal, fname_lal);
    
    Aggs_lal = pp.pcm(Aggs_lal);
    Aggs_lal = pp.hough_kook2(Aggs_lal);

end

n_agg_lal = length(Aggs_lal);

pp_manu_lal = cell(1,n_agg_lal);
fname_pp_lal = cell(1,n_agg_lal);
dpp_lal = cell(1,n_agg_lal);
dbarpp_manu_lal = zeros(n_agg_lal,1);
sigmapp_manu_lal = zeros(n_agg_lal,1);
npp_manu_lal = zeros(n_agg_lal,1);

for i = 1 : n_agg_lal

    fname_pp_lal{i} = char(strcat(fd_pp_lal, {'\'}, Aggs_lal(i).fname(1:end-4), '.csv'));

    if isfile(fname_pp_lal{i})

        % opts_pp = detectImportOptions(fname_pp{i});
        % opts_pp = setvartype(opts_pp, 'char');
        pp_manu_lal{i} = readtable(fname_pp_lal{i});
        
        dpp_lal{i} = sqrt(4 * pp_manu_lal{i}.Area(2:end) / pi);

        dbarpp_manu_lal(i) = geomean(dpp_lal{i});
        sigmapp_manu_lal(i) = morph.geostd(dpp_lal{i});
        npp_manu_lal = length(dpp_lal{i});
        
        plt_lal = scatter(Aggs_lal(i).da, dbarpp_manu_lal(i), 15,...
            hex2rgb('#006989'), '^');
        
    end

end

set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 12,...
    'TickLength', [0.02 0.02], 'XScale', 'log', 'YScale', 'log')
xlim([20,2000])
ylim([10,60])
xlabel('$d_\mathrm{a}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
ylabel('$d_\mathrm{pp}$ [nm]', 'interpreter', 'latex', 'FontSize', 14)
legend(cat(2, plt_0, plt_lal),...
    cat(2, {'Olfert and Rogak (2019)'}, {'Low agglomeration'}),...
    'interpreter', 'latex', 'FontSize', 14, 'location', 'northwest')
