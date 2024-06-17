

function getFlyCount(img,savePath,xlrows)

baseFolder = getCloudPath;

% load the outlines of the arena and their coordinates
load([savePath ' arena coordinates.mat'])

[excelfile, Excel, xlFile] = load_QuadBowlExperiments;

% Number of flies:
nframes = length(img);
nflies = nan(1,4); 

for arena = 1:4
    arenaData(arena).nflies = excelfile{xlrows(arena),Excel.numflies};
    nflies(arena) = arenaData(arena).nflies;
end

if any(isnan(nflies)) 
    nflies = [];
    % manual count of flies
    fprintf('\nCount the number of flies in the picture by clicking them\n then hit ENTER\n')
    T = true;
    while T        
        fprintf('Number of flies counted \n     A      B      C      D\n    ---    ---    ---    ---\n')
        for jj = 1:nframes
            demoImg = img(jj).data; %display the arena frame image
            PT(jj).frame = readPoints(demoImg);
            for arena = 1:4
                centre = arenaData(arena).centre;
                X = PT(jj).frame(1,:);
                Y = PT(jj).frame(2,:);
                % find points within each arena sphere
                loc = (((X-centre(1)).^2 + (Y-centre(2)).^2).^0.5)<=arenaData(arena).r;
                nflies(arena,jj) = sum(loc);
            end
            disp(nflies(:,jj)')
        end
        disp('Number of flies:')
        A = nflies(1,:)'; B = nflies(2,:)'; C = nflies(3,:)'; D = nflies(4,:)';
        numFlies = [A,B,C,D];
        disp(table(A,B,C,D))
        
        % Skip Frame Speed Hack
        % if there is a row of zeros (skipped tracking frame) remove the empty frame
        skipFrameloc = sum(numFlies==0,2)==4; %skipped counting frame locations
        mismatchLoc = diff(numFlies)==0; % across frame mismatched numbers
        % if we skipped labeling a frame but the first two counts match up,
        % then auto add the number from the first two counts to the list
        if any(skipFrameloc) && all(mismatchLoc(1,:)) 
                numFlies(skipFrameloc,:) = numFlies(1,:);
        end
        % Resume normal miscounted fly error catching protocol
        count_match = (diff(numFlies,(nframes-1),1)==0);
        if all(count_match)
            nflies = [];
            for arena = 1:4
                nflies(arena) = median(numFlies(:,arena)); 
                arenaData(arena).nflies = nflies(arena);
            end
            T = false;
        else
%             disp('    A     B     C     D')
%             disp(nflies')
            switch questdlg('Nonmatching fly counts, manually select number?:')
                case 'Yes'
                    for arena = 1:4
                        arenaNums{arena} = num2str(median(numFlies(:,arena)));
                        if ~count_match(arena)
                            arenaNums{arena} = 'NaN';
                        end
                    end
                    prompt = {'Arena A','Arena B', 'Arena C', 'Arena D'}; dlgtitle = 'Input';
                    answer = inputdlg(prompt,dlgtitle,[1 25],arenaNums);
                    nflies = [];
                    for arena = 1:4
                        arenaData(arena).nflies = str2double(answer{arena});
                        nflies(arena) = arenaData(arena).nflies;
                    end
                    T = false;
                case 'No'
                    T = true;
                case 'Cancel'
                    return
            end
        end
    end
    % write number of flies into the excel sheet
    try
        for arena = 1:4
            writematrix(nflies(arena),xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.numflies) num2str(xlrows(arena))]);
        end
    catch
        h = warndlg('Close Experiment Summary excel file and then close this warning box');
        uiwait(h)
        for arena = 1:4
            writematrix(nflies(arena),xlFile,'Sheet','Exp List','Range',[Alphabet(Excel.numflies) num2str(xlrows(arena))]);
        end
    end
end
disp('Number of flies written to excel:')
disp('     A     B     C     D')
disp(nflies)

end