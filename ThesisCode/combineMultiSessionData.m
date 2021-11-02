function [combinedData] = combineMultiSessionData(dataField,analysisDir,type,previousData,sessionName,animal)

cd(analysisDir)
yes=0;
al=0;
if isempty(type)
    data = load(dataField);
else
    data = load([dataField,'_',type]);
end
field = fields(data);
if strcmp(dataField,'tSp')
    if isempty(previousData)
        combinedData.(sessionName).(animal){1} = [data.(field{:})];
        return
    else
        variables = fields(previousData);
        for v=1:length(variables)
            combinedData.(variables{v}) = [previousData.(variables{v})];
            if strcmp(variables(v),sessionName)
                yes=1;
            end
        end
        if yes==1
            Afield = fields(previousData.(sessionName));
            for a=1:length(Afield)
                if strcmp(Afield(a),animal)
                    al=1;
                end
            end
            if al==1
                prior=previousData.(sessionName).(animal);
                combinedData.(sessionName).(animal){length(prior)+1} = [data.(field{:})];
            else
                combinedData.(sessionName).(animal){1} = [data.(field{:})];
            end
        else
            combinedData.(sessionName).(animal){1} = [data.(field{:})];
        end
        return
    end
end
if strcmp(dataField,'tsHE')
    if isempty(previousData)
        for f=1:length(field)
            combinedData.(sessionName).(animal).(field{f}){1} = [data.(field{f})];
        end
        return
    else
        variables = fields(previousData);
        for v=1:length(variables)
            combinedData.(variables{v}) = [previousData.(variables{v})];
            if strcmp(variables(v),sessionName)
                yes=1;
            end
        end
        if yes==1
            Afield = fields(previousData.(sessionName));
            for a=1:length(Afield)
                if strcmp(Afield(a),animal)
                    al=1;
                end
            end
            if al==1
                for f=1:length(field)
                    prior=previousData.(sessionName).(animal).(field{f});
                    combinedData.(sessionName).(animal).(field{f}){length(prior)+1} = [data.(field{f})];
                end
            else
                for f=1:length(field)
                    combinedData.(sessionName).(animal).(field{f}){1} = [data.(field{f})];
                end
            end
        else
            for f=1:length(field)
                combinedData.(sessionName).(animal).(field{f}){1} = [data.(field{f})];
            end
        end
        return
    end
end
if strcmp(dataField,'Trials') 
    if isempty(previousData)
        for f=1:length(field)
            field2 = fields(data.(field{f}));
            for f2=1:length(field2)
                combinedData.(sessionName).(field{f}).(animal).(field2{f2}){1} = [data.(field{f}).(field2{f2})];
            end
        end
        return
    else
        variables = fields(previousData);
        for v=1:length(variables)
            combinedData.(variables{v}) = [previousData.(variables{v})];
            if strcmp(variables(v),sessionName)
                yes=1;
            end
        end
        if yes==1
            Afield = fields(previousData.(sessionName).(field{1}));
            for a=1:length(Afield)
                if strcmp(Afield(a),animal)
                    al=1;
                end
            end
            if al==1
                for f=1:length(field)
                    field2 = fields(data.(field{f}));
                    for f2=1:length(field2)
                        prior=previousData.(sessionName).(field{f}).(animal).(field2{f2});
                        combinedData.(sessionName).(field{f}).(animal).(field2{f2}){length(prior)+1} = [data.(field{f}).(field2{f2})];
                    end
                end
            else
                for f=1:length(field)
                    field2 = fields(data.(field{f}));
                    for f2=1:length(field2)
                        combinedData.(sessionName).(field{f}).(animal).(field2{f2}){1} = [data.(field{f}).(field2{f2})];
                    end
                end
            end
        else
            for f=1:length(field)
                field2 = fields(data.(field{f}));
                for f2=1:length(field2)
                    combinedData.(sessionName).(field{f}).(animal).(field2{f2}){1} = [data.(field{f}).(field2{f2})];
                end
            end
        end
        return
    end
end


if isempty(previousData)
    for f=1:length(field)
        combinedData.(sessionName).(animal).(field{f}){1} = [data.(field{f})];
    end
    return
else
    variables = fields(previousData);
    for v=1:length(variables)
        combinedData.(variables{v}) = [previousData.(variables{v})];
        if strcmp(variables(v),sessionName)
            yes=1;
        end
    end
    if yes==1
        Afield = fields(previousData.(sessionName));
        for a=1:length(Afield)
            if strcmp(Afield(a),animal)
                al=1;
            end
        end
        if al==1
            for f=1:length(field)
                prior=previousData.(sessionName).(animal).(field{f});
                combinedData.(sessionName).(animal).(field{f}){length(prior)+1} = [data.(field{f})];
            end
        else
            for f=1:length(field)
                combinedData.(sessionName).(animal).(field{f}){1} = [data.(field{f})];
            end
        end
    else
        for f=1:length(field)
            combinedData.(sessionName).(animal).(field{f}){1} = [data.(field{f})];
        end
    end
    return
end
