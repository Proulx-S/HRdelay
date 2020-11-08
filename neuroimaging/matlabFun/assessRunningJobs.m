function jobStr = assessRunningJobs(user)
% Output number of jobs running

[status,result] = unix(['showq -u ' user ' -r']);
ind=strfind(result,'active jobs');
i = 1;
foundIt1 = false;
foundIt2 = false;
jobStr = [];
while ~foundIt2
    tmp = result(ind(2)-i);
    if ~foundIt1
        if ~isnan(str2double(tmp))
            foundIt1 = true;
            jobStr = [tmp jobStr];
        end
    elseif foundIt1
        if ~isnan(str2double(tmp))
            jobStr = [tmp jobStr];
        else
            foundIt2 = true;
        end
    end
    i = i+1;
end
jobStr = str2double(jobStr);