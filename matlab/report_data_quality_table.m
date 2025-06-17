function [txt,measures] = report_data_quality_table(dq_table)
arguments
    dq_table {mustBeA(dq_table,'table')}
end

% average over targets and eyes
measures.all  = removevars(groupsummary(removevars(dq_table, {'eye','target_id'}), 'file', 'mean'),'GroupCount');
measures.all.Properties.VariableNames(2:end) = dq_table.Properties.VariableNames(4:end);

% do summary statistics
measures.mean = mean(measures.all(:,2:end)   ,1);
measures.std  = std (measures.all(:,2:end),[],1);
measures.min  = min (measures.all(:,2:end),[],1);
measures.max  = max (measures.all(:,2:end),[],1);

% make text. A little overcomplete, user can trim what they don't want
% N.B.: do not include data loss/effective frequency, nor bcea. Bcea is
% niche, user who wants it can get that themselves. Data loss you'd really
% want to report for all the analyzed data, not just this validation
% procedure.
n_target = length(unique(dq_table.target_id));
n_subj   = size(measures.all,1);
version  = ETDQ_version();
txt = sprintf('For %d participants, the average inaccuracy in the data determined from a %d-point validation procedure using ETDQualitizer v%s (Niehorster et al., in prep) was %.2f° (STD %.2f°, range %.2f°--%.2f°). Average RMS-S2S precision was %.3f° (STD %.3f°, range %.3f°--%.3f°) and STD precision %.3f° (STD %.3f°, range %.3f°--%.3f°).',n_subj,n_target,version,measures.mean.offset,measures.std.offset,measures.min.offset,measures.max.offset,measures.mean.rms_s2s,measures.std.rms_s2s,measures.min.rms_s2s,measures.max.rms_s2s,measures.mean.std,measures.std.std,measures.min.std,measures.max.std);
