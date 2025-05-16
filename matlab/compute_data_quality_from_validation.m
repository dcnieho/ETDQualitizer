function dq = compute_data_quality_from_validation(gaze, unit, screen, advanced, include_data_loss)
arguments
    gaze                           {mustBeA(gaze,'table')}
    unit                (1,:) char {mustBeMember(unit,{'pixels','degrees'})}
    screen                         {mustBeScalarOrEmpty,mustBeA(screen,{'ScreenConfiguration','double'})} = []
    advanced            (1,1)      {mustBeA(advanced         ,'logical')} = false
    include_data_loss   (1,1)      {mustBeA(include_data_loss,'logical')} = false
end

% get all targets
targets         = unique(gaze.target_id);
targets(targets==-1) = [];
target_locations= nan(length(targets),2);
for t=1:length(targets)
    qTarget = gaze.target_id==targets(t);
    iTarget = find(qTarget,1);
    target_locations(t,:) = [gaze.tar_x(iTarget) gaze.tar_y(iTarget)];
end

% ensure we have target locations in degrees
if unit=='pixels'
    if isempty(screen)
        error('If unit is "pixels", a screen configuration must be supplied')
    end
    [target_locations(:,1),target_locations(:,2)] = screen.pix_to_deg(target_locations(:,1),target_locations(:,2));
elseif ~strcmp(unit, 'degrees')
    error('unit should be "pixels" or "degrees"')
end

% now, per target, compute data quality metrics
vars = {'eye','target_id','offset','offset_x','offset_y','rms_s2s','rms_s2s_x','rms_s2s_y','std','std_x','std_y','bcea','bcea_orientation','bcea_ax1','bcea_ax2','bcea_aspect_ratio'};
if include_data_loss
    vars = [vars {'data_loss','effective_frequency'}];
end
dq = table('Size',[length(targets),length(vars)], 'VariableTypes', [{'char'} repmat({'double'},1,length(vars)-1)], 'VariableNames', vars);
eyes = {'left','right'};
for e=1:length(eyes)
    if ~any(strcmp([eyes{e} '_x'],gaze.Properties.VariableNames))
        continue
    end
    for t=1:length(targets)
        oi = (e-1)*length(targets)+t;
        is_target = gaze.target_id==targets(t);
        dq_calc = DataQuality(gaze.([eyes{e} '_x'])(is_target), gaze.([eyes{e} '_y'])(is_target), gaze.timestamp(is_target)/1000, unit, screen);    % timestamps are in ms in the file

        dq.eye(oi) = eyes(e);
        dq.target_id(oi) = targets(t);
        [dq.offset(oi),dq.offset_x(oi),dq.offset_y(oi)]    = dq_calc.accuracy(target_locations(t,1), target_locations(t,2));
        [dq.rms_s2s(oi),dq.rms_s2s_x(oi),dq.rms_s2s_y(oi)] = dq_calc.precision_RMS_S2S();
        [dq.std(oi),dq.std_x(oi),dq.std_y(oi)] = dq_calc.precision_STD();
        [dq.bcea(oi),dq.bcea_orientation(oi),dq.bcea_ax1(oi),dq.bcea_ax2(oi),dq.bcea_aspect_ratio(oi)] = dq_calc.precision_BCEA();
        if include_data_loss
            dq.data_loss(oi) = dq_calc.data_loss_percentage();
            dq.effective_frequency(oi) = dq_calc.effective_frequency();
        end
    end
end

if ~advanced
    to_drop = ~ismember(dq.Properties.VariableNames,{'eye', 'target_id', 'offset', 'rms_s2s', 'std', 'bcea', 'data_loss', 'effective_frequency'});
    dq = removevars(dq,dq.Properties.VariableNames(to_drop));
end
