clear
close all
clc

projectRoot = '/Users/manon/Documents/Uni/Trinity4/FYP/CNSP-resources/CNSP';
datasetFolder = 'diliBach/';

% Load envelope TRF results
load(fullfile( ...
    projectRoot, ...
    'resultFigures', ...
    datasetFolder, ...
    '1-8Hz', ...
    'TRF_nSub20_cond1_stimIdx1.mat'));

fprintf('Number of subjects: %d\n', length(modelAll));

fprintf('\nmodelAll(1).w size:\n');
disp(size(modelAll(1).w));

fprintf('\nmodelAll(1) fields:\n');
disp(fieldnames(modelAll(1)));

nSubs = length(modelAll);
nLags = size(modelAll(1).w, 2);
nChans = size(modelAll(1).w, 3);

% TRF weights:
% subject × lag × channel
W = zeros(nSubs, nLags, nChans);

for s = 1:nSubs
    W(s,:,:) = squeeze(modelAll(s).w);
end

% e.g. W(3,8,20) -> TRF weights for subject 3 at lag 8 on electrode 20

fprintf('W size:\n');
disp(size(W));

% reproducible subject split
rng(42);   % fixed seed

subjectOrder = randperm(nSubs); % shuffling 

nA = floor(nSubs/2);

groupA = subjectOrder(1:nA);
groupB = subjectOrder(nA+1:end);

fprintf('Group A:\n');
disp(groupA);

fprintf('Group B:\n');
disp(groupB);

% average TRF for each group
meanA = squeeze(mean(W(groupA,:,:), 1));
meanB = squeeze(mean(W(groupB,:,:), 1));
size(meanA)


% sanity check
time = modelAll(1).t;

meanA_allChannels = mean(meanA, 2);
meanB_allChannels = mean(meanB, 2);

figure;

plot(time, meanA_allChannels, 'LineWidth', 2);
hold on;
plot(time, meanB_allChannels, 'LineWidth', 2);

xlabel('Latency (ms)');
ylabel('TRF weight');
legend('Group A', 'Group B');
title('Envelope TRF - Split-half comparison');
grid on;


% Component Analysis

W_A = W(groupB,:,:);

nDiscovery = length(groupB);

% Compute t-statistics
meanW = mean(W_A, 1);
stdW  = std(W_A, 0, 1);

tVals = squeeze(meanW ./ (stdW ./ sqrt(nDiscovery)));

% tVals:
% lag × channel


% Permutation test
% Sign-flipping creates data expected under H0: TRF = 0

% We are testing: How large could a t-value get purely by chance?

nPerm = 5000;
rng(42);
maxAbsT = zeros(nPerm,1);

for p = 1:nPerm

    % Randomly flip each subject's TRF
    signs = 2 * randi([0 1], nDiscovery, 1) - 1;

    permW = W_A .* reshape(signs, nDiscovery, 1, 1);

    permMean = mean(permW, 1);
    permStd  = std(permW, 0, 1);

    permT = squeeze(permMean ./ ...
        (permStd ./ sqrt(nDiscovery)));

    % Store the biggest statistic anywhere in the entire TRF
    maxAbsT(p) = max(abs(permT(:)));

end


% Corrected significance threshold

alpha = 0.05;
tThreshold = prctile(maxAbsT, 100*(1-alpha));
fprintf('\nCorrected t threshold: %.3f\n', tThreshold);


% Separate positive and negative significant points

positiveMask = tVals >= tThreshold;
negativeMask = tVals <= -tThreshold;


% Find latencies containing significant responses

positiveTimes = any(positiveMask, 2);
negativeTimes = any(negativeMask, 2);

time = modelAll(1).t;


% Automatically find consecutive temporal components

positiveComponents = findRuns(positiveTimes);
negativeComponents = findRuns(negativeTimes);

% Require components to persist across multiple time samples

minComponentSamples = 2;

% Length of each detected component in number of consecutive samples
positiveLengths = positiveComponents(:,2) - positiveComponents(:,1) + 1;
negativeLengths = negativeComponents(:,2) - negativeComponents(:,1) + 1;

% Remove components that occur at only one isolated time sample
positiveComponents = positiveComponents( ...
    positiveLengths >= minComponentSamples, :);

negativeComponents = negativeComponents( ...
    negativeLengths >= minComponentSamples, :);


% Characterise positive components

components = struct([]);

componentNum = 0;

for c = 1:size(positiveComponents,1)

    componentNum = componentNum + 1;

    startIdx = positiveComponents(c,1);
    endIdx   = positiveComponents(c,2);

    sectionT = tVals(startIdx:endIdx,:);

    % strongest positive point
    [peakT, linearIdx] = max(sectionT(:));

    [localLag, peakChannel] = ind2sub(size(sectionT), linearIdx);

    peakLag = startIdx + localLag - 1;

    % corrected p-value from permutation distribution
    correctedP = ...
        (1 + sum(maxAbsT >= abs(peakT))) / ...
        (nPerm + 1);

    components(componentNum).polarity = "Positive";
    components(componentNum).startTime = time(startIdx);
    components(componentNum).endTime = time(endIdx);
    components(componentNum).peakTime = time(peakLag);
    components(componentNum).peakChannel = peakChannel;
    components(componentNum).peakT = peakT;
    components(componentNum).pCorrected = correctedP;

    % total significant t-statistic magnitude
    componentMask = positiveMask(startIdx:endIdx,:);

    components(componentNum).mass = ...
        sum(sectionT(componentMask));

end


% Characterise negative components

for c = 1:size(negativeComponents,1)

    componentNum = componentNum + 1;

    startIdx = negativeComponents(c,1);
    endIdx   = negativeComponents(c,2);

    sectionT = tVals(startIdx:endIdx,:);

    % strongest negative point
    [peakT, linearIdx] = min(sectionT(:));

    [localLag, peakChannel] = ind2sub(size(sectionT), linearIdx);

    peakLag = startIdx + localLag - 1;

    correctedP = ...
        (1 + sum(maxAbsT >= abs(peakT))) / ...
        (nPerm + 1);

    components(componentNum).polarity = "Negative";
    components(componentNum).startTime = time(startIdx);
    components(componentNum).endTime = time(endIdx);
    components(componentNum).peakTime = time(peakLag);
    components(componentNum).peakChannel = peakChannel;
    components(componentNum).peakT = peakT;
    components(componentNum).pCorrected = correctedP;

    componentMask = negativeMask(startIdx:endIdx,:);

    components(componentNum).mass = ...
        sum(abs(sectionT(componentMask)));

end

% Print detected components

if isempty(components)

    fprintf('\nNo significant components found.\n');

else

    componentTable = struct2table(components);

    % strongest component first
    componentTable = sortrows(componentTable, ...
                              'mass', ...
                              'descend');

    fprintf('\nDetected components:\n');
    disp(componentTable);

end

function runs = findRuns(binaryVector)

    binaryVector = binaryVector(:);

    changes = diff([false; binaryVector; false]);

    starts = find(changes == 1);
    ends   = find(changes == -1) - 1;

    runs = [starts ends];

end