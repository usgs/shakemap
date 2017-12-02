function [mean,mean_metadata,std,std_metadata] = read_shake_data(hdf_file,imt,component)
% read_shake_data  Read ShakeMap HDF file data.
%   [mean,mean_metadata,std,std_metadata] = read_shake_data(hdf_file,imt,component)
%   Inputs:
%     - hdf_file Path to a ShakeMap HDF5 results file.
%     - imt String intensity measure type (IMT)
%     ('MMI','PGA','SA(1.0)',etc.)
%     - component IMT component (usually 'Larger').
%   Outputs:
%     - mean Mean IMT value grid.
%     - mean_metadata Mean IMT value metadata.
%     - std Standard deviation of IMT value grid.
%     - std_metadata Standard deviation IMT value metadata.
    info = h5info(hdf_file);
    for i=1:length(info.Groups)
        group = info.Groups(i);
        gname = group.Name;
        if ~isempty(strfind(gname,'__imt'))
            parts = strsplit(gname,'_');
            imtstr = parts{3};
            compstr = strip(parts{4},'_');
            if strcmpi(imt,imtstr) && strcmpi(component,compstr)
                
                for j=1:length(group.Datasets)
                    dataset = info.Groups(i).Datasets(j);
                    dname = dataset.Name;
                    parts = strsplit(dname,'_');
                    quantity = strip(parts{2},'_');
                    dsetname = sprintf('%s/%s',gname,dname);
                    if strcmpi(quantity,'mean')
                        mean = h5read(hdf_file,dsetname);
                        mean_metadata = struct;
                        for k=1:length(dataset.Attributes)
                            attribute = dataset.Attributes(k);
                            key = attribute.Name;
                            if iscell(attribute.Value)
                                value = attribute.Value{1};
                            else
                                value = double(attribute.Value);
                            end
                            mean_metadata = setfield(mean_metadata,key,value);
                        end
                    end
                    if strcmpi(quantity,'std')
                        std = h5read(hdf_file,dsetname);
                        std_metadata = struct;
                        for k=1:length(dataset.Attributes)
                            attribute = dataset.Attributes(k);
                            key = attribute.Name;
                            if iscell(attribute.Value)
                                value = attribute.Value{1};
                            else
                                value = double(attribute.Value);
                            end
                            std_metadata = setfield(std_metadata,key,value);
                        end
                    end
                end
            else
                continue;
            end
            
        end
    end
end