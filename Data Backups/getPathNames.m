
function paths = getPathNames
% paths = getPathNames

% folder location endings
paths.courtship = 'Courtship Videos/';
paths.single_trial = 'Trial Data/';
paths.raw_data = 'Raw Data/';
paths.grouped_trials = 'Data structures/';
paths.group_comparision = 'Grouped Data Structures/';
paths.NOAA = 'NOAA Data/';

% variable drive names
paths.storageDrive = 'Data Storage';
paths.removableDrive = 'OnTheGoData';

% general server paths
mac_server_path = '/Volumes/shared/Evyn/DATA/';
second_server_path = '\\syn-jeannelab-1.its.yale.edu\nas1-jeannelab-cc0928-Neuro-ps\Evyn\DATA\';

% root directory locations
paths.fixedDriveLocations  = {'ACADIA','EVYNPC','evyndickinson','rebeccaray'}; % computer names that have permanent local drives

% TODO: make one giant switchboard option where it only gives the location
% paths for the computer in question. [1/30/2026]


% ACADIA:
paths.acadiaServerPath = 'S:\Evyn\DATA\';
paths.acadiaLocalPath = 'D:\Evyn Lab Data\';
paths.acadiaServerTwoPath = second_server_path;

% TOGIAK:
paths.togiakServerPath = 'S:\Evyn\DATA\';
paths.togiakServerTwoPath =  second_server_path;

% EVYN PC:
paths.EvynPCServerPath = '\\svalbard.med.yale.internal\shared\Evyn\DATA\';
% paths.EvynPCServerPath = 'S:\Evyn\DATA\';
paths.EvynPCLocalPath = 'K:\DATA\';
paths.EvynPCServerTwoPath = second_server_path;

% EVYN M3 MBP:
paths.EvynMacServerPath = mac_server_path;
paths.EvynMacLocalPath = '/Users/evyndickinson/Documents/Jeanne Lab/DATA/';

% BECCA MAC AIR:
paths.BeccaLocalPath = '/Users/rebeccaray/Desktop/Jeanne Lab Data/';
paths.BeccaServerPath = mac_server_path;

% DENALI:
paths.denaliServerPath = 'Z:\Evyn\DATA\';

% CHILKAT: 
paths.chilkatServerPath = 'S:\Evyn\DATA\';
paths.chilkatServerTwoPath = second_server_path;

% SLEEPING GIANT:
paths.SGServerPath = 'S:\Evyn\DATA\';
paths.SGServerTwoPath = 'R:\Evyn\DATA\';

% Server Paths: 
paths.MacServerPath = mac_server_path;






