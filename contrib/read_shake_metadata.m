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
%     data.
%    
    hdfinfo = h5info(fname);
    dicts = {'config'};
    json_strings = {'info.json','rupture','stations'};
    for i=1:length(hdfinfo.Groups)
        group = hdfinfo.Groups(i);
        gname = group.Name;
        if contains(gname,'dictionaries')
            for j=1:length(group.Datasets)
                dataset = group.Datasets(j);
                if strcmpi(dataset.Name,'info.json')
                    dsetname = sprintf('%s/%s',gname,dataset.Name);
                    info_string = h5read(fname,dsetname);
                    info = jsondecode(info_string{1});
                elseif strcmpi(dataset.Name,'rupture')
                    dsetname = sprintf('%s/%s',gname,dataset.Name);
                    rupture_string = h5read(fname,dsetname);
                    rupture = jsondecode(rupture_string{1});
                elseif strcmpi(dataset.Name,'stations_dict')
                    dsetname = sprintf('%s/%s',gname,dataset.Name);
                    stations_string = h5read(fname,dsetname);
                    stations = jsondecode(stations_string{1});
                elseif strcmpi(dataset.Name,'config')
                    dsetname = sprintf('%s/%s',gname,dataset.Name);
                    config_string = h5read(fname,dsetname);
                    config = jsondecode(config_string{1});
                end
            end
        end
    end    
end

