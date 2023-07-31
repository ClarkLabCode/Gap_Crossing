% ALLCROSSSTATSNAMESMAKER Lets user create AllCrossStatsNames
%  
%  This script generates the cell array AllCrossStatsNames that holds the
%  names of all the experiments being held in AllCrossStats.
%  
%  Structure of the cell array:
%  Column 1: All generic controls
%  Column 2: [Gal4] > +; shts(5); shR(5)
%  Column 3: [Gal4] > +; +; +
%  Columns 4+: Other experiments related to [Gal4] in that row

% Initialize All_CrossStats_Names (41 drivers, up to 15 types of exp per driver
AllCrossStatsNames = cell(41,15);

% All generic controls
AllCrossStatsNames{1,1} = 'IsoD1';
AllCrossStatsNames{2,1} = 'IsoD1, dark';
AllCrossStatsNames{3,1} = 'old split empty > shts (2&3)';
AllCrossStatsNames{4,1} = 'split empty DBD > shts (2&3)';
AllCrossStatsNames{5,1} = 'old split empty > shts (2&3), dark';
AllCrossStatsNames{6,1} = 'split empty DBD > +;tshGal80;shR(5) (tsh screen)';
AllCrossStatsNames{7,1} = 'split empty DBD > +;tshGal80;shR(5) (1st test run)';
AllCrossStatsNames{8,1} = 'old single empty > shts (2&3)';
AllCrossStatsNames{9,1} = 'new single empty > shts (2&3)';
AllCrossStatsNames{10,1} = 'new split empty > shts (2&3)';
AllCrossStatsNames{11,1} = 'old split empty > +;tshGal80;shR(5)';
AllCrossStatsNames{12,1} = 'old split empty > TNT(2)';
AllCrossStatsNames{13,1} = 'shts (2&3)';
AllCrossStatsNames{14,1} = 'shts (2&3) > IsoD1';
AllCrossStatsNames{15,1} = 'IsoD1, bottom lit';
AllCrossStatsNames{16,1} = 'IsoD1, top lit';
AllCrossStatsNames{17,1} = 'TNT (2) > IsoD1';
AllCrossStatsNames{18,1} = 'IsoD1, one eye painted';
AllCrossStatsNames{19,1} = 'IsoD1, two eyes painted';
AllCrossStatsNames{20,1} = 'IsoD1 (second screen)';
AllCrossStatsNames{21,1} = 'IsoD1, dark (second screen)';
AllCrossStatsNames{22,1} = 'split empty DBD > +;tshGal80;shR(5) (tsh screen), dark';

% All first round screen silencing
AllCrossStatsNames{1,2} = 'split T2a > shts (2&3)';
AllCrossStatsNames{2,2} = 'split T3 > shts (2&3)';
AllCrossStatsNames{3,2} = 'split T4+T5 > shts (2&3)';
AllCrossStatsNames{4,2} = 'split C2+C3 > shts (2&3)';
AllCrossStatsNames{5,2} = 'split C3 > shts (2&3)';
AllCrossStatsNames{6,2} = 'LLPC1 > shts (2&3)';
AllCrossStatsNames{7,2} = 'split LPC1 > shts (2&3)';
AllCrossStatsNames{8,2} = 'split LPLC1 > shts (2&3)';
AllCrossStatsNames{9,2} = 'split LPLC2 > shts (2&3)';
AllCrossStatsNames{10,2} = 'LPLC3 > shts (2&3)';
AllCrossStatsNames{11,2} = 'LPLC4 > shts (2&3)';
AllCrossStatsNames{12,2} = 'CH > shts (2&3)';
AllCrossStatsNames{13,2} = 'FD1+3 > shts (2&3)';
AllCrossStatsNames{14,2} = 'H1 > shts (2&3)';
AllCrossStatsNames{15,2} = 'H2 > shts (2&3)';
AllCrossStatsNames{16,2} = 'Hx > shts (2&3)';
AllCrossStatsNames{17,2} = 'HS+VS > shts (2&3)';
AllCrossStatsNames{18,2} = 'V1 > shts (2&3)';
AllCrossStatsNames{19,2} = 'split LC4 > shts (2&3)';
AllCrossStatsNames{20,2} = 'split LC6 > shts (2&3)';
AllCrossStatsNames{21,2} = 'split LC9 > shts (2&3)';
AllCrossStatsNames{22,2} = 'split LC10 > shts (2&3)';
AllCrossStatsNames{23,2} = 'split LC11 > shts (2&3)';
AllCrossStatsNames{24,2} = 'split LC12 > shts (2&3)';
AllCrossStatsNames{25,2} = 'split LC13 > shts (2&3)';
AllCrossStatsNames{26,2} = 'split LC14 > shts (2&3)';
AllCrossStatsNames{27,2} = 'split LC14b > shts (2&3)';
AllCrossStatsNames{28,2} = 'split LC15 > shts (2&3)';
AllCrossStatsNames{29,2} = 'split LC16 > shts (2&3)';
AllCrossStatsNames{30,2} = 'split LC17 > shts (2&3)';
AllCrossStatsNames{31,2} = 'split LC18 > shts (2&3)';
AllCrossStatsNames{32,2} = 'split LC19 > shts (2&3)';
AllCrossStatsNames{33,2} = 'split LC20 > shts (2&3)';
AllCrossStatsNames{34,2} = 'split LC21 > shts (2&3)';
AllCrossStatsNames{35,2} = 'split LC22 > shts (2&3)';
AllCrossStatsNames{36,2} = 'split LC23 > shts (2&3)';
AllCrossStatsNames{37,2} = 'split LC24 > shts (2&3)';
AllCrossStatsNames{38,2} = 'split LC25 > shts (2&3)';
AllCrossStatsNames{39,2} = 'split LC26 > shts (2&3)';
AllCrossStatsNames{40,2} = 'Eye Movements 1 > shts (2&3)';
AllCrossStatsNames{41,2} = 'Eye Movements 2 > shts (2&3)';

% All first round screen IsoD1 controls
AllCrossStatsNames{1,3} = 'split T2a > IsoD1';
AllCrossStatsNames{2,3} = 'split T3 > IsoD1';
AllCrossStatsNames{3,3} = 'split T4+T5 > IsoD1';
AllCrossStatsNames{4,3} = 'split C2+C3 > IsoD1';
AllCrossStatsNames{5,3} = 'split C3 > IsoD1';
AllCrossStatsNames{6,3} = 'LLPC1 > IsoD1';
AllCrossStatsNames{7,3} = 'split LPC1 > IsoD1';
AllCrossStatsNames{8,3} = 'split LPLC1 > IsoD1';
AllCrossStatsNames{9,3} = 'split LPLC2 > IsoD1';
AllCrossStatsNames{10,3} = 'LPLC3 > IsoD1';
AllCrossStatsNames{11,3} = 'LPLC4 > IsoD1';
AllCrossStatsNames{12,3} = 'CH > IsoD1';
AllCrossStatsNames{13,3} = 'FD1+3 > IsoD1';
AllCrossStatsNames{14,3} = 'H1 > IsoD1';
AllCrossStatsNames{15,3} = 'H2 > IsoD1';
AllCrossStatsNames{16,3} = 'Hx > IsoD1';
AllCrossStatsNames{17,3} = 'HS+VS > IsoD1';
AllCrossStatsNames{18,3} = 'V1 > IsoD1';
AllCrossStatsNames{19,3} = 'split LC4 > IsoD1';
AllCrossStatsNames{20,3} = 'split LC6 > IsoD1';
AllCrossStatsNames{21,3} = 'split LC9 > IsoD1';
AllCrossStatsNames{22,3} = 'split LC10 > IsoD1';
AllCrossStatsNames{23,3} = 'split LC11 > IsoD1';
AllCrossStatsNames{24,3} = 'split LC12 > IsoD1';
AllCrossStatsNames{25,3} = 'split LC13 > IsoD1';
AllCrossStatsNames{26,3} = 'split LC14 > IsoD1';
AllCrossStatsNames{27,3} = 'split LC14b > IsoD1';
AllCrossStatsNames{28,3} = 'split LC15 > IsoD1';
AllCrossStatsNames{29,3} = 'split LC16 > IsoD1';
AllCrossStatsNames{30,3} = 'split LC17 > IsoD1';
AllCrossStatsNames{31,3} = 'split LC18 > IsoD1';
AllCrossStatsNames{32,3} = 'split LC19 > IsoD1';
AllCrossStatsNames{33,3} = 'split LC20 > IsoD1';
AllCrossStatsNames{34,3} = 'split LC21 > IsoD1';
AllCrossStatsNames{35,3} = 'split LC22 > IsoD1';
AllCrossStatsNames{36,3} = 'split LC23 > IsoD1';
AllCrossStatsNames{37,3} = 'split LC24 > IsoD1';
AllCrossStatsNames{38,3} = 'split LC25 > IsoD1';
AllCrossStatsNames{39,3} = 'split LC26 > IsoD1';
AllCrossStatsNames{40,3} = 'Eye Movements 1 > IsoD1';
AllCrossStatsNames{41,3} = 'Eye Movements 2 > IsoD1';

% Extra experiments from hits:
% T4+T5
AllCrossStatsNames{3,4} = 'split T4+T5 > 2nd screen UAS';
AllCrossStatsNames{3,5} = 'split T4+T5 > IsoD1 (2nd screen)';
AllCrossStatsNames{3,6} = 'split T4+T5 > 2nd screen UAS, dark';
AllCrossStatsNames{3,7} = 'split T4+T5 > +;tshGal80;shR(5)';
AllCrossStatsNames{3,8} = 'split T4+T5 > IsoD1 (tsh ctrl)';
AllCrossStatsNames{3,9} = 'split T4+T5 > shts (2&3) (2nd time)';
AllCrossStatsNames{3,10} = 'split T4+T5 > IsoD1 (2nd time)';
AllCrossStatsNames{3,11} = 'split T4+T5 > shts (2&3), bottom lit';
AllCrossStatsNames{3,12} = 'split T4+T5 > shts (2&3), top lit';
AllCrossStatsNames{3,13} = 'split T4+T5 > TNT (2)';

% C3
AllCrossStatsNames{5,4} = 'split C3 > 2nd screen UAS';
AllCrossStatsNames{5,5} = 'split C3 > IsoD1 (2nd screen)';
AllCrossStatsNames{5,6} = 'split C3 > 2nd screen UAS, dark';
AllCrossStatsNames{5,7} = 'split C3 > +;tshGal80;shR(5)';
AllCrossStatsNames{5,8} = 'split C3 > IsoD1 (tsh ctrl)';

% Hx
AllCrossStatsNames{16,4} = 'Hx > 2nd screen UAS';
AllCrossStatsNames{16,5} = 'Hx > IsoD1 (2nd screen)';
AllCrossStatsNames{16,6} = 'Hx > 2nd screen UAS, dark';
AllCrossStatsNames{16,7} = 'Hx > +;tshGal80;shR(5)';
AllCrossStatsNames{16,8} = 'Hx > IsoD1 (tsh ctrl)';

% LC10
AllCrossStatsNames{22,4} = 'split LC10 > 2nd screen UAS';
AllCrossStatsNames{22,5} = 'split LC10 > IsoD1 (2nd screen)';
AllCrossStatsNames{22,6} = 'split LC10 > 2nd screen UAS, dark';
AllCrossStatsNames{22,7} = 'split LC10 > +;tshGal80;shR(5)';
AllCrossStatsNames{22,8} = 'split LC10 > IsoD1 (tsh ctrl)';

% LC14b
AllCrossStatsNames{27,4} = 'split LC14b > 2nd screen UAS';
AllCrossStatsNames{27,5} = 'split LC14b > IsoD1 (2nd screen)';
AllCrossStatsNames{27,6} = 'split LC14b > 2nd screen UAS, dark';
AllCrossStatsNames{27,7} = 'split LC14b > +;tshGal80;shR(5)';
AllCrossStatsNames{27,8} = 'split LC14b > IsoD1 (tsh ctrl)';

% LC25
AllCrossStatsNames{38,4} = 'split LC25 > 2nd screen UAS';
AllCrossStatsNames{38,5} = 'split LC25 > IsoD1 (2nd screen)';
AllCrossStatsNames{38,6} = 'split LC25 > 2nd screen UAS, dark';
AllCrossStatsNames{38,7} = 'split LC25 > +;tshGal80;shR(5)';
AllCrossStatsNames{38,8} = 'split LC25 > IsoD1 (tsh ctrl)';

% LC26
AllCrossStatsNames{39,4} = 'split LC26 > 2nd screen UAS';
AllCrossStatsNames{39,5} = 'split LC26 > IsoD1 (2nd screen)';
AllCrossStatsNames{39,6} = 'split LC26 > 2nd screen UAS, dark';
AllCrossStatsNames{39,7} = 'split LC26 > +;tshGal80;shR(5)';
AllCrossStatsNames{39,8} = 'split LC26 > IsoD1 (tsh ctrl)';

% tsh controls for eye movements
% Eye Move 1
AllCrossStatsNames{40,7} = 'Eye Movement 1 > +;tshGal80;shR(5)';
AllCrossStatsNames{40,8} = 'Eye Movement 1 > IsoD1 (tsh ctrl)';

% Eye Move 2
AllCrossStatsNames{41,7} = 'Eye Movement 2 > +;tshGal80;shR(5)';
AllCrossStatsNames{41,8} = 'Eye Movement 2 > IsoD1 (tsh ctrl)';