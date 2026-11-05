=== STIMULUS-ONLY CHECKS ===

Structure Integrity:
- stimIsStruct: Stimulus is a structure
- stimHasData: Has 'data' field
- stimHasFs: Has 'fs' (sampling frequency) field

Data Format:
- stimDataIsCell: data field is a cell array
- stimDataNotEmpty: data cell array is not empty
- stimDataNumeric: Contains numeric values
- stimFsPositive: fs is positive number (>0)
- stimNamesOptional: 'names' field optional but must be cell array if present

Data Consistency:
- stimDimsConsistent: Feature dimensions consistent across trials
- stimNoNaN: No NaN values in data
- stimNoInf: No Inf values in data

=== NEURAL-ONLY CHECKS ===

Structure Integrity:
- neuralIsStruct: Neural is a structure
- neuralHasData: Has 'data' field
- neuralHasFs: Has 'fs' field

Data Format:
- neuralDataIsCell: data field is a cell array
- neuralDataNotEmpty: data cell array is not empty
- neuralDataNumeric: Contains numeric values
- neuralFsPositive: fs is positive number (>0)
- neuralChanlocsOptional: 'chanlocs' field optional but must be structure if present

Data Consistency:
- neuralChansConsistent: Number of channels consistent across trials
- neuralNoNaN: No NaN values in data
- neuralNoInf: No Inf values in data

=== COMPATIBILITY CHECKS (BOTH FILES REQUIRED) ===

Sampling & Trials:
- fsMatch: Sampling frequencies match (stim.fs == neural.fs)
- nTrialsMatch: Number of trials match (size(stim.data,2) == length(neural.data))

Data Alignment:
- trialLengthsCompat: Trial lengths compatible (within 1% tolerance)

=== VALIDATION CATEGORIES ===
1. Structure Integrity
2. Required Fields  
3. Data Format
4. Data Consistency
5. Compatibility

=== CRITICAL FAILURES (Will fail validation) ===

Stimulus:
- Missing 'data' or 'fs' fields
- data not a cell array
- fs ≤ 0
- Empty data cell array

Neural:
- Missing 'data' or 'fs' fields
- data not a cell array  
- fs ≤ 0
- Empty data cell array

Compatibility:
- Mismatched sampling rates
- Different number of trials
- Severe trial length mismatches (>1%)