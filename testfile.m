p = ppg_anal;

%(SamlingFrequency, ResampleFrequency, ResampleBool, FilterFloor, FilterHigh, SegmentBool, SegmentLength(ms), InvertSignal)
%The data will be resampled to ResampleFrequency if ResampleBool is true.
%The data will be segmented into segments of SegmentLength if SegmentBool is true. If SegmentBool is false, the data will not be segmented and SegmentLength will be ignored.

%p.LoadPPG(1000, 300, true, 0.4, 8, true, 10000, false);
p.LoadDirectory(1000, 500, false, 0.4, 8, false, 10000, false);


%Calculates fiducial points and features for the current and every following segment.
%p.ProcessDataset();

%Finds the best cycle in the current segment.
%[Min1, Min2, cycle_corr_quality, cycle_skew_quality, segment_corr_quality, segment_skew_quality, cycle_quality, segment_quality]
[min1, min2, ~, ~, ~, ~] = p.FindBestCycle();

%Calculates fiducial points and features for the cycle between min1 and min2.
p.CalculateFiducial(min1, min2);

%The current cycle and its derivetives are stored in p.seg, p.VPG, p.APG and p.JPG
%The index of the current segment is stored in p.total_seg_idx.
%All fiducial points stored in p.fiducials and all features stored in p.features.

%Selects next segment in file/directory.
p.Next();

[min1, min2, ~, ~, ~, ~] = p.FindBestCycle();

p.CalculateFiducial(min1, min2);

%Generates output to a .mat file if set to true or three .xlsx files if set to false.
p.GenerateOutput(true);
