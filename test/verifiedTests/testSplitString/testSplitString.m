% The COBRAToolbox: testSplitString.m
%
% Purpose:
%     - testSplitString tests the functionality of splitString()
%       and checks solution against a known solution.
%
% Authors:
%     - Lemmer El Assal February 2017
%

% define the path to The COBRAToolbox
pth = which('initCobraToolbox.m');
CBTDIR = pth(1:end-(length('initCobraToolbox.m') + 1));

initTest([CBTDIR, filesep, 'test', filesep, 'verifiedTests', filesep, 'testTools']);

ref_fields = {'Some';'Strings';'Delimited'};

testString1 = 'Some Strings Delimited';
testString2 = 'Some|Strings|Delimited';
fields1 = splitString(testString1);
fields2 = splitString(testString2,'\|');

assert(isequal(ref_fields,fields1))
assert(isequal(ref_fields,fields2))

% change the directory
cd(CBTDIR)
