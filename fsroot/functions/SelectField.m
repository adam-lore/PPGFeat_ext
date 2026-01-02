function [field, res] = SelectField(data)
    field = {};
    res = true;

    if ~isscalar(data) || ~isa(data, 'struct')
        return
    end

    while 1
        if isempty(field)
            mat_field_arr = fieldnames(data);
        else
            mat_field_arr = fieldnames(getfield(data, field{:}));
        end
        
        if (size(mat_field_arr) == 1)
            field_index = 1;
        else
            [field_index, tf] = listdlg('SelectionMode', 'single', 'ListString', mat_field_arr);
            if (tf == false)
                res = false;
                return
            end
            
        end

        field{end + 1} = mat_field_arr{field_index};

        if ~isscalar(getfield(data, field{:})) || ~isa(getfield(data, field{:}), 'struct')
            break
        end
    end
end