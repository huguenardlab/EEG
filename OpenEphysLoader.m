
function [data,header] = OpenEphysLoader(directory,varargin)
% [data,header] = OpenEphysLoader( directory, (chans,type,decimate) )
%
% Wrapper for "load_open_ephys_faster.m" provided by the OpenEphys group to
% make loading data easier, while also providing more thorough checks on
% dropped data packets and channels recorded
% -----------------------------------------------------------------------
%               >>> INPUTS >>>
% directory:
%   The name of the folder containing the data to load
%
% (chans):
%   A vector of channel numbers to load. Default = "Inf"... i.e. all channels 
%   available for the given "type" (see below)
%
% (type):
%   A string specifying the data format to load. Can be
%   "continuous", "spikes", "ADC", or "AUX". Note that, although "AUX" &
%   "ADC" channels are saved as .continous files, here we specify them as
%   separate "types" since AUX channels usually contain stimuli, etc.
%   Default is to load "continous" files
%
% (desiredfs)
%   an integer that gives the desired sampling rate (Hz).  A decimation factor
%   will then be calculated based on this
%
% (decimate)
%   An integer specifying the decimation value to reduce the memory load
%   when forming the data matrix. If "decimate" = 0 or 1, no decimation
%   will occur. Decimation will only be applied to "continuous" files.
%   Default = 0.
%
%               <<< OUTPUTS <<<
% data:
%   Column-oriented data matrix containing the channels specified by the
%   user as each column.
%
% header:
%   A structure containing header information from the recording, including 
%   channels pulled out of the folder, the sampling rate, etc.
% 
% -----------------------------------------------------------------------
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% -----------------------------------------------------------------------
% By JMS, 2/14/2017

global p

% Parse inputs
p = check_inputs( varargin );

% collect file names
[paths,filenames] = get_files( directory, p.type ); % extract all files of specified type
[data,header] = extract_header( paths{1} );
header.rootDir = directory;

% extract the possible channels available in the input directory
availableChans = cellfun( @cell2mat, cellfun( @(c)(regexp(c,'CH\d*', 'match')),filenames,'un',0 ), 'un',0 );
availableChans = cell2mat( cellfun( @(c)(str2double(c(3:end))),availableChans,'un', 0 ) );

% correct for "Inf" channels if necessary
nfiles = numel( paths );
nchans = numel( p.chans );
if isinf( p.chans )
    nchans = nfiles;
    p.chans = availableChans;
end 

% check if channels provided do not match those requested
if nchans > nfiles
    warning( 'More channels requested than files available...defaulting to available channels' );
    nchans = nfiles;
    p.chans = availableChans(ismember( availableChans,p.chans ));
elseif ~isequal( p.chans(ismember( p.chans,availableChans )),p.chans )
    warning( 'Not all channels requested are available...removing absent channels' );
    p.chans = p.chans(ismember( p.chans,availableChans ));
    nchans = numel( p.chans );
end
p.chans = sort( p.chans ); % sort for good measure 

% loop over files specified by requested channels 
for f = 1:nchans
    fileToOpen = paths{availableChans == p.chans(f)};
    header.chans(f)=p.chans(f);
    % load the file and store into "data"
    if max( strcmp( p.type, {'continuous','AUX','ADC'} ) ) > 0
        if p.decimate > 1
            [tmp,~,~] = load_open_ephys_data_faster( fileToOpen );
            if p.decimate > 10 
                fd = factor(p.decimate);
                for i = 1:numel( fd )
                    tmp = decimate( tmp,fd(i),'FIR' );
                end
                data(1:size(tmp,1),f) = tmp;
            else
                tmp=decimate( tmp,p.decimate,'FIR' );
                data(1:size(tmp,1),f) = tmp;
            end
        else
            [data(:,f),~,~] = load_open_ephys_data_faster( fileToOpen );
        end
    else
        [data{f},header.timeStamps{f},info] = load_open_ephys_data_faster( fileToOpen ); % no decimation for "spikes" or "events"
        if ~isempty( data{f} )
            header.gain(:,f) = info.gain(1,:);
            header.thresh(:,f) = info.thresh(1,:);
            if isfield(header,'sortedID')
                header.sortedID{f} = info.sortedId;
            end
        end
    end
end
end

%% Functions
function p = check_inputs( inputs )
    pnames = {'chans','type','decimate','desiredfs'};
    defaults = {Inf,'continuous',1,0};
                 
    p = inputParser;             
    for j = 1:numel(pnames)
        p.addParameter(pnames{j},defaults{j});
    end
     
    p.parse(inputs{:});
    p = p.Results;
    
    % check decimate value
    if isprime( p.decimate ) && p.decimate > 10
        p.decimate = p.decimate - 1; % makes it even
    end
end

function [paths,filenames] = get_files( directory,type )
    list = dir( directory );
    filenames = {list(3:end).name}; % store into cell array
    
    % filter for file with specified "type"
    filenames = filenames(~cellfun( @isempty,strfind( filenames,type ) ));    
    if strcmp( type,'continuous' )
        filenames = filenames(cellfun( @isempty,strfind( filenames,'AUX' ) )); % eliminate "AUX" channels
        filenames = filenames(cellfun( @isempty,strfind( filenames,'Continuous_Data' ) )); % eliminate the "continuous_data" file
        filenames = filenames(cellfun( @isempty,strfind( filenames,'ADC' ) )); % eliminate ADC channels
    end
    if isempty( filenames )
        error( ['No files found of specified type: .',type] );
    end
    
    % form full paths for the files
    paths = fullfile( directory,filenames ); 
end
 
function [data,header] = extract_header( file )
    global p
    
    [ts,info] = load_open_ephys_data_faster_headeronly( file );

    header.fs = info.header.sampleRate;
    header.allChans = p.chans;
    header.date = info.header.date_created;
    header.type = p.type;    
    
    if p.desiredfs > 0
        p.decimate=floor(header.fs / p.desiredfs);
    end
    
    switch p.type
        
        case {'continuous','AUX','ADC'}
            
            % pull out specific .continuous header info
            header.fs = header.fs / max( p.decimate,1 );
            header.bitVolts = info.header.bitVolts;
            header.decimateVal = p.decimate;
            header.timeStamps = ts(1:max( 1,p.decimate ):end); % to align with resampled data
            
            % check for dropped packets
            ts = ts/info.header.sampleRate;
            dt = diff( ts );
            timePerPacket = info.header.blockLength / info.header.sampleRate;
            dropped = dt > timePerPacket;
            nDropped = round( sum( dt(dropped) ) / timePerPacket );
            totalDropped = sum( dt(dropped) );
            if nDropped > 0
                fprintf( 'Warning: %i packets dropped (total: %04f sec) during recording session\n',...
                    nDropped,totalDropped );
            end        
            header.numPacketsDropped = nDropped;
            header.timeDropped = totalDropped;
            header.droppedTimeStamps = ts([false;dropped]);
            
            % preallocate data matrix
            N = numel( header.timeStamps );
            data = zeros( N,numel( p.chans ) );
            
        case 'spikes'
            header.chansPerElectrode = info.header.num_channels;
            header.gain = zeros( header.chansPerElectrode,numel( p.chans ) );
            header.thresh = header.gain;
            
            % preallocate data cell
            data = cell( 1,numel( p.chans ) );
            header.timeStamps = data;
            
            % check for sorted IDs
            if isfield(info,'sortedId')
                header.sortedID = data;
            end
    end
end
    