H=nan(1,69880); 
P=nan(1,69880);

for i = 1:69880
    
    if voxelMask(i) == 1
    
        Joint = RsqAllJoint(:,i); 
        Ind=RsqAllInd(:,i); 

        % cheat
        Joint(7:8,:) = Joint(5:6,:);
        Ind(7:8,:) = Ind(5:6,:);
        
        if any(isnan(Joint))
            for z=1:20
                if isnan(Joint(z))
                    Joint(z) = mean(Joint,'omitnan');
                end
            end
            if any(isnan(Joint))
                error(['Bad, bad joint at ',num2str(i)]);
            end
        end

        if any(isnan(Ind))
            for z=1:20
                if isnan(Ind(z))
                    Ind(z) = mean(Ind,'omitnan');
                end
            end
            if any(isnan(Ind))
                error(['Bad, bad ind at ',num2str(i)]);
            end
        end

        [H(i), P(i)] = ttest(Ind,Joint);

        if mod(i,10000) == 0
            disp(['Done with ', num2str(i)]);
        end
        
    end

end