function [tr,te] = newNormalization(tr,te,p)
if ~isfield(p,'rotTo')
    p.rotTo = 45;
end


%Center phase
if p.normToPlaid
    if p.filterData
        switch p.keepInfo
            case 'p'
                %Train
                phase = angle(tr);
                phase_shift = circ_mean(phase,[],1)-(p.rotTo/180*pi);
                phase = phase-repmat(phase_shift,size(phase,1),1);
                %             %%%%%%%%%%%%%%%%%%%
                %             phase = phase + repmat(rand(size(phase,1),1)*2*pi,1,size(phase,2));
                %             phase = phase + repmat(rand(1,size(phase,2))*2*pi,size(phase,1),1);
                %             %%%%%%%%%%%%%%%%%%%
                [trR,trI] = pol2cart(phase,1);
                %Test
                phase = angle(te);
                phase = phase-repmat(phase_shift,size(phase,1),1);
                [teR,teI] = pol2cart(phase,1);
            case 'm'
                trR = tr;
                teR = te;
        end
    else
        error('double-check that')
    end
else
    if p.filterData
        switch p.keepInfo
            case 'p'
                %Train
                phase = angle(tr);
                phase_shift = circ_mean(phase,[],1)-(p.rotTo/180*pi);
                phase = phase-repmat(phase_shift,size(phase,1),1);
                %             %%%%%%%%%%%%%%%%%%%
                %             phase = phase + repmat(rand(size(phase,1),1)*2*pi,1,size(phase,2));
                %             phase = phase + repmat(rand(1,size(phase,2))*2*pi,size(phase,1),1);
                %             %%%%%%%%%%%%%%%%%%%
                [trR,trI] = pol2cart(phase,1);
                %Test
                phase = angle(te);
                phase = phase-repmat(phase_shift,size(phase,1),1);
                [teR,teI] = pol2cart(phase,1);
            case 'm'
                %Train
                trR = abs(tr);
                %Test
                teR = abs(te);
        end
    else
        error('that is weird, phase_shift based on training applied to testing')
        %Train
        phase = angle(tr); mag = abs(tr);
        phase_shift = circ_mean(phase,[],1)-(p.rotTo/180*pi);
        phase = phase-repmat(phase_shift,size(phase,1),1);
        [trR,trI] = pol2cart(phase,mag);
        %Test
        phase = angle(te); mag = abs(te);
        phase = phase-repmat(phase_shift,size(phase,1),1);
        [teR,teI] = pol2cart(phase,mag);
    end    
end

if ~p.doCrossVal
    %Center
    %real
    trR_shift = mean(trR,1);
    trR = trR-repmat(trR_shift,size(trR,1),1);
    teR = teR-repmat(trR_shift,size(teR,1),1);
    %imag
    if exist('trI','var')
        trI_shift = mean(trI,1);
        trI = trI-repmat(trI_shift,size(trI,1),1);
        teI = teI-repmat(trI_shift,size(teI,1),1);
    end
    
    %Scale
    if exist('trI','var')
        tr_scale = mean([std(trR,[],1); std(trI,[],1)]);
        trR = trR./repmat(tr_scale,size(trR,1),1); trI = trI./repmat(tr_scale,size(trI,1),1);
        teR = teR./repmat(tr_scale,size(teR,1),1); teI = teI./repmat(tr_scale,size(teI,1),1);
        
        tr = complex(trR,trI);
        te = complex(teR,teI);
    else
        tr_scale = std(trR,[],1);
        trR = trR./repmat(tr_scale,size(trR,1),1);
        teR = teR./repmat(tr_scale,size(teR,1),1);
        
        tr = trR;
        te = teR;
    end
end

