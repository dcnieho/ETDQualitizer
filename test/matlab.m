p = fileparts(mfilename("fullpath"));
addpath(fullfile(p,'..','matlab'))

screen = ScreenConfiguration(528.0, 296.9997253417969, 1920, 1080, 650);

% get all tsv files in data folder
fold = fullfile(p,'data');
files= dir(fold);
files= files(~[files.isdir]);
files= files(arrayfun(@(x) endsWith(x.name,'.tsv'),files));

all_dq = {};
for f=1:length(files)
    fprintf('----------\n%s\n',files(f).name);
    gaze = readtable(fullfile(fold,files(f).name),'FileType','text','Delimiter','\t');

    dq = compute_data_quality_from_validation(gaze, 'pixels', screen, false, true)     % include_data_loss for testing, this is probably *not* what you want
    all_dq{end+1,1} = addvars(dq, repmat(string(files(f).name),size(dq,1),1), 'Before', 1, 'NewVariableNames',{'file'});

    eyes = {'left','right'};
    for e=1:length(eyes)
        if ~any(strcmp([eyes{e} '_x'],gaze.Properties.VariableNames))
            continue
        end
        % and RMS S2S calculated in two ways over the whole datafile
        dq_calc = DataQuality(gaze.([eyes{e} '_x']),gaze.([eyes{e} '_y']),gaze.timestamp/1000,'pixels',screen); % timestamps are in ms in the file
    
        fp = split(files(f).name,'_');
        fs = str2double(fp{end}(1:end-length('Hz.tsv')));
        window_len = round(.2*fs);  % 200 ms
    
        fprintf('RMS S2S using median (%s eye): %.4f deg\n', eyes{e}, dq_calc.precision_RMS_S2S(@(x) median(x, 'omitnan')))
        fprintf('RMS S2S using moving window (%s eye): %.4f deg\n', eyes{e}, dq_calc.precision_using_moving_window(window_len,"RMS_S2S"))
    
        % data loss and effective frequency
        fprintf('Data loss (%s eye): %.1f%%\n', eyes{e}, dq_calc.data_loss_percentage())
        fprintf('Data loss nominal frequency (%s eye): %.1f%%\n', eyes{e}, dq_calc.data_loss_percentage_nominal(fs))
        fprintf('Effective frequency (%s eye): %.1f Hz\n', eyes{e}, dq_calc.effective_frequency())
    end
end
all_dq = vertcat(all_dq{:});
