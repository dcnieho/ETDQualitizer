function dq = compute_data_quality_from_validation(gaze, unit, screen, include_data_loss)
arguments
    gaze                           {mustBeA(gaze,'table')}
    unit                (1,:) char {mustBeMember(unit,{'pixels','degrees'})}
    screen                         {mustBeScalarOrEmpty,mustBeA(screen,{'ScreenConfiguration','double'})} = []
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
vars = {'target_id','offset','offset_x','offset_y','rms_s2s','rms_s2s_x','rms_s2s_y','std','std_x','std_y','bcea','bcea_orientation','bcea_ax1','bcea_ax2','bcea_aspect_ratio'};
if include_data_loss
    vars = [vars {'data_loss','effective_frequency'}];
end
dq = table('Size',[length(targets),length(vars)], 'VariableTypes', repmat({'double'},1,length(vars)), 'VariableNames', vars);
for t=1:length(targets)
    is_target = gaze.target_id==targets(t);
    dq_calc = DataQuality(gaze.left_x(is_target), gaze.left_y(is_target), gaze.timestamp(is_target)/1000, unit, screen);    % timestamps are in ms in the file

    dq.target_id(t) = targets(t);
    [dq.offset(t),dq.offset_x(t),dq.offset_y(t)]    = dq_calc.accuracy(target_locations(t,1), target_locations(t,2));
    [dq.rms_s2s(t),dq.rms_s2s_x(t),dq.rms_s2s_y(t)] = dq_calc.precision_RMS_S2S();
    [dq.std(t),dq.std_x(t),dq.std_y(t)] = dq_calc.precision_STD();
    [dq.bcea(t),dq.bcea_orientation(t),dq.bcea_ax1(t),dq.bcea_ax2(t),dq.bcea_aspect_ratio(t)] = dq_calc.precision_BCEA();
    if include_data_loss
        dq.data_loss(t) = dq_calc.data_loss_percentage();
        dq.effective_frequency(t) = dq_calc.effective_frequency();
    end
end
