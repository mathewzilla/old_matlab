README for code to generate figures in MH Evans et al "State of the Art in Whisker-Based Object Localisation: from beam equations to Gaussian Processes"

NOTE: Used as a github learning area 12/13 by MHE

Code is run from AR_master_code.m

Data (in own folder) consists of 4 separate sets from 3 different robots, to date only 3 are in a processable format:

1. data_XY.mat - XY positioning robot data, previously published in Evans et al Frontiers in Neurorobotics 2013 and Proc. ROBIO/SAB 2010 . Cell array. 4 sets x 26 speeds x 101 radii. Files are 750 x 2 (at 1kHz).

2. ScratchPeaksXY.mat - Scratchbot data. Previously published in Evans et al Frontiers in Neurorobotics 2013 and Proc. ROBIO 2010.  Cell array. 3 speeds, 3 radii, 8 contacts. Files are 2001 x 2. Probably at 1KHz (need to check)

3. roombaRadiusXY.mat - Crunchbot data.   4 whiskers, 6 contacts at different radii (90-130), 5 repeats. Files are 4000 x 2 (at 2kHz). Code had to be fiddled with during the loop as some data wasn't quite in the right place. roombaRadius{3,3,5} was truncated, so filled with noise from the start of the trial

4. Collated.mat - Scratchbot data. Old 'expt1', saw/sin whisker movement, 11 radii [90-190mm], previously unpublished.


HOW THE CODE WORKS ATM

User defines dataset to use (1-4) or lets the code loop through datasets

Data is loaded

User chooses whether to smooth or differentiate the data (will run all conditions and select best for paper)

PREPROCESSING seperates the data into training and test sets, providing (optionally) smooth or diff'd data in seperated and concatenated files (diff classifiers need different things)

ATM the code runs once, using an arbitrary file as the test data. Need to change code to run 4 times (leave one out) then combine/average scores

Training and testing occurs, set up/called using the old code. 

TO DO: This could be editted so a single classifier is called with a single script. User could then decide which classifier to try. This would return scores, times for training/testing.
    


TO DO: 

i) Get one whole classifier working on one dataset again TICK

ii) Get timer measurement working TICK

iii) Get other classifiers working again TICK

iv) Get one classifier working on another dataset TICK

v) Get all classifiers working on another dataset TICK

vi) Get all classifiers working on all datasets TICK

vii) Share code/data with collaborators to add different classifiers

viii) Decide on preprocessing/ data handling for paper, inc. 2 axes or just one?

ix) Decide on paper figures/ figure styles and code up

ix) Publish code and data with paper
