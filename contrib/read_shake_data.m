function [dmean,mean_metadata,dstd,std_metadata] = read_shake_data(fname,imt,component)
% read_shake_data  Read IMT data from a ShakeMap container.
%   [dmean,mean_metadata,dstd,std_metadata] = read_shake_data(fname,imt,component);
%   Inputs:
%     - fname ShakeMap HDF file results.
%     - imt  Intensity Measure Type (MMI,PGA,SA(0.3),etc)
%     - component Ground motion component (Larger,Rotd50,etc.)
%   Outputs:
%     - dmean Mean IMT array.
%     - mean_medatata Structure containing metadata for mean array.
%       - units	Physical units of IMT or standard deviation.
%       - digits	Number of significant digits to use for the values.
%       - xmin	The eastern boundary of the grid (degrees longitude)
%       - xmax	The eastern boundary of the grid (degrees longitude)
%       - ymin	The southern boundary of the grid (degrees latitude)
%       - ymax	The northern boundary of the grid (degrees latitude)
%       - nx	The number of grid points in the x dimension
%       - ny	The number of grid points in the y dimension
%       - dx	The grid interval in the x dimension
%       - dy	The grid interval in the y dimension
%     - dstd Standard IMT array.
%     - std_medatata Structure containing metadata for mean array
%     (identical to mean_metadata).
%

    info = h5info(fname);
    for i=1:length(info.Groups)
        group = info.Groups(i);
        gname = group.Name;
        if contains(gname,imt,'IgnoreCase',1) && contains(gname,component,'IgnoreCase',1)
            for j=1:length(group.Datasets)
                dataset = group.Datasets(j);
                dname = dataset.Name;
                parts = strsplit(dname,'_');
                quantity = parts(2);
                %compstr = parts(4);
                if strcmpi(quantity,'mean')
                    [dmean,mean_metadata] = get_quantity(fname,dataset,gname);
                elseif strcmpi(quantity,'std')
                    [dstd,std_metadata] = get_quantity(fname,dataset,gname);
                end
            end
        end
    end    
end

function [quant,metadata] = get_quantity(fname,dataset,groupname)
    metadata = struct;
    for k=1:length(dataset.Attributes)
        key = dataset.Attributes(k).Name;
        value = dataset.Attributes(k).Value;
        if iscell(value)
            value = value{1};
        else
            value = double(value);
        end
        metadata = setfield(metadata,key,value);
    end
    dsetname = sprintf('%s/%s',groupname,dataset.Name);
    quant = h5read(fname,dsetname);
end

