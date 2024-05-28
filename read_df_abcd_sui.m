function read_df_abcd_sui = importfile(filename, dataLines)
%IMPORTFILE Import data from a text file

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Setup the Import Options
opts = delimitedTextImportOptions("NumVariables", 31);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["VarName1", "X", "src_subject_id", "eventname", "reshist_addr1_adi_wsum", "ADI_pc1", "anytlfbSU", "mypiSU", "NeighSafety", "NeighCrime", "FamilyConflict", "site_id_l", "interview_age", "rel_family_id", "reshist_addr1_pb", "reshist_addr1_no2_2016_aavg", "Income", "HighestEd", "Male_bin", "White", "Black", "Hispanic", "Asian", "Other", "reshist_addr1_pm252016aa", "famhx_ss_momdad_dg_p", "famhx_ss_momdad_alc_p", "lead_1", "air_pollution_no2", "air_pollution_pm25", "ADI"];
opts.VariableTypes = ["double", "double", "string", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
opts = setvaropts(opts, 3, "WhitespaceRule", "preserve");
opts = setvaropts(opts, [7, 8, 12], "TrimNonNumeric", true);
opts = setvaropts(opts, [7, 8, 12], "ThousandsSeparator", ",");
opts = setvaropts(opts, [3, 4], "EmptyFieldRule", "auto");
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Import the data
read_df_abcd_sui = readtable(filename, opts);

end