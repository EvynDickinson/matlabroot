


function result = runBatch_step1_updates(filePath, plate,pix2mm,con_type,radii)
        tic 
        try load([filePath(1:end-7) '.mat'])
        catch 
            result = false;
            return
        end
        disp(['Loaded ' folder ' ' expName])

        % V2_UPDATE (5.21.25)
        for arena = 1:4 
            arenaData(arena).plate = plate;
            arenaData(arena).con_type = con_type;
            arenaData(arena).pix2mm = pix2mm;
        end

        for arena = 1:4
            % Pull fly coordinate position data for this arena 
            X = T.(['x' arenaIdx{arena}]);
            Y = T.(['y' arenaIdx{arena}]);
            
            % Get center position for each well
            centers = arenaData(arena).wellcenters;
        
            % Calculate distance to food and occupancy counts
            for well = 1:4
                % center of the well coordinates
                dX = X - centers(1,well);
                dY = Y - centers(2,well);
        
                % Find distance to well center
                dist2well = sqrt(dX.^2 + dY.^2);
                
                % Find points within arena:
                loc = dist2well<=radii; % tracked points within circle
                N = sum(loc,2); % how many points tracked within that circle
                dist2well_mm = dist2well ./ pix2mm;
                
                % store data
                arenaData(arena).dist2well(:,well) = mean(dist2well_mm,2,'omitnan');
                arenaData(arena).dist2well_err(:,well) = std(dist2well_mm,0,2,'omitnan');
                arenaData(arena).occ_N(:,well) = N;
                arenaData(arena).occ_P(:,well) = N/nflies(arena);
            end
        end
        
        initial_vars{end+1} = 'filePath';
        clearvars('-except',initial_vars{:})

        % save the new data trial information
        disp('Saving...')
        save(filePath,'-v7.3')
        if exist(filePath, 'file')
            result = true;
        end
        toc
end























        