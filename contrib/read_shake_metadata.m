function [config,info,rupture,stations] = read_shake_metadata(fname)
% read_shake_metadata  Read metadata from a ShakeMap container.
%   [config,info,rupture,stations] = read_shake_metadata(fname);
%   Inputs:
%     - fname ShakeMap HDF file results.
%   Outputs:
%     - config Structure containing configuration data.
%     - info Structure ShakeMap processing parameters and map summary information.
%     - rupture Structure containing rupture coordinates and properties.
%     - stations Structure containing observed and instrument ground motion
%     data. (Currently empty pending further development.)
%    
    hdfinfo = h5info(fname);
    dicts = {'config'};
    json_strings = {'info.json','rupture','stations'};
    for i=1:length(hdfinfo.Groups)
        group = hdfinfo.Groups(i);
        gname = group.Name;
        if contains(gname,'config')
            config = get_dictionary(group,'config');
            x = 1;
        elseif contains(gname,'info.json')
            info = get_string(group);
        elseif contains(gname,'rupture')
            rupture = get_string(group);
        elseif contains(gname,'stations')
            %stations = get_string(group);
            stations = struct;
        end
    end    
end

function structure = get_string(group)
% Read and decode JSON string from HDF.
    jsonstr = group.Attributes(1).Value{1};
    structure = jsondecode(jsonstr);
end

function outstruct = get_dictionary(group,name)
% Recursively unpack structure from HDF groups.
    outstruct = struct;
    groups = group.Groups;
    if ~isempty(groups)
        for i=1:length(groups)
            ggroup = groups(i);
            gname = ggroup.Name;
            parts = strsplit(gname,'/');
            gname = parts{end};
            tstruct = get_dictionary(ggroup,gname);
            outstruct = setfield(outstruct,gname,tstruct);
        end
    else
        for i=1:length(group.Attributes)
            key = group.Attributes(i).Name;
            value = group.Attributes(i).Value;
            if iscell(value)
                if length(value) == 1
                    value = value{1};
                end
            else
                value = double(value);
            end
            outstruct = setfield(outstruct,key,value);
        end
    end
end