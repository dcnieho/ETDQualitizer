p = fileparts(mfilename("fullpath"));
if isfolder(fullfile(p,'..','matlab'))
    addpath(fullfile(p,'..','matlab'))
else
    addpath(fullfile(p,'..'))
end


% Run all tests
testClasses = {
    'TestScreenConfiguration'
    'TestFickConversions'
    'TestDataQualityMetrics'
    'TestDataQualityClass'
};
results = runtests(testClasses);

% Display results
disp(results);
