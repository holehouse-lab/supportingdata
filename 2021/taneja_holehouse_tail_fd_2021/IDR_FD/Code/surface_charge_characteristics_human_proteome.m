
%in order to run this script:
%go to 'Import Data' in Matlab, then select ../Processed_Data/matlab_metadata.csv
%you then need to convert the column type to a Text (not numeric)

%to generate the mean electrostatic potential and patchiness with respect
%to the FD:IDR junction, run surface_charge_characteristics_human_proteome(matlabmetadata, 1)

%to generate the mean electrostatic potential and patchiness with respect
%to random residues, run surface_charge_characteristics_human_proteome(matlabmetadata, 0)

function output = surface_charge_characteristics_human_proteome(matlabmetadata, idr_nearest_only)

apbs_input_dir = '../PDB_Structure_FD_candidates_within10/APBS_Output/Ion_Concentration_10';

parfor (i = 1:height(matlabmetadata), 8)  %run in parallel
%for i=1:height(matlabmetadata) -- run this sequentially 
    disp(i);
    type = string(matlabmetadata{i,'type'});
    if (idr_nearest_only && strcmp('IDR_nearest',type)) || (~idr_nearest_only && ~strcmp('IDR_nearest',type))
        pdb_id = string(matlabmetadata{i,'PDB_ID'});
        chain = string(matlabmetadata{i,'chain'});
        pdb_id_chain = strcat(pdb_id, '_', chain); 
        fp = matlabmetadata{i,'fp'};
        fn = matlabmetadata{i,'fn'};

        disp(pdb_id_chain);

        idr_x = matlabmetadata{i,'x_loc'};
        idr_y = matlabmetadata{i,'y_loc'};
        idr_z = matlabmetadata{i,'z_loc'};

        disp(type);
        gen_surface_charge_characteristics(apbs_input_dir, pdb_id_chain, idr_x, idr_y, idr_z, type);        
    end 
    
end 

end
    


function output = gen_surface_charge_characteristics(apbs_input_dir, pdb_id_chain, idr_x, idr_y, idr_z, type)

    dx_file = strcat(apbs_input_dir, '/', pdb_id_chain, '.dx');
    ply_file = strcat(apbs_input_dir,'/', pdb_id_chain, '.ply');
        
    if isfile(dx_file) && isfile(ply_file)    
        output1=gen_surf_pot_and_coords(dx_file, ply_file);
        X = output1;
        X(any(isnan(X), 2),:)=[];
        pcutoff = 2;
        ncutoff = -2;
        correl_l = 2; 
        radius_start = correl_l; 
        patchiness_all = patchcorr_optimized(X,pcutoff,ncutoff,correl_l); 
        q_mean_all = mean(X(:,4)); 
        patchiness_q_f_radius = patchcorr_optimized_wrt_idr(X,pcutoff,ncutoff,radius_start, correl_l, idr_x, idr_y, idr_z); 
        
        if (strcmp(type,'IDR_nearest'))
            csvwrite(strcat('../Processed_Data/FD_Surface_Characteristics/', pdb_id_chain, '_patchiness_all', '.csv'),patchiness_all); 
            csvwrite(strcat('../Processed_Data/FD_Surface_Characteristics/', pdb_id_chain, '_q_all', '.csv'),q_mean_all); 
            csvwrite(strcat('../Processed_Data/FD_Surface_Characteristics/', pdb_id_chain, '_patchiness_q_f_radius', '.csv'),patchiness_q_f_radius);
        else
            csvwrite(strcat('../Processed_Data/FD_Surface_Characteristics_Random_Residue/', pdb_id_chain, '_type_', type, '_patchiness_all', '.csv'),patchiness_all); 
            csvwrite(strcat('../Processed_Data/FD_Surface_Characteristics_Random_Residue/', pdb_id_chain, '_type_', type, '_q_all', '.csv'),q_mean_all); 
            csvwrite(strcat('../Processed_Data/FD_Surface_Characteristics_Random_Residue/', pdb_id_chain, '_type_', type, '_patchiness_q_f_radius', '.csv'),patchiness_q_f_radius);
        end
        
    end
end
    

function output1 = gen_surf_pot_and_coords(dx_file, ply_file)

    PlyData = pcread(ply_file);
    PlyCoor = PlyData.Location;

    dxFileId = fopen(dx_file,'r');
    for j = 1:4
        fgetl(dxFileId); % Irrelevant lines
    end
    gridCountsLine = fgetl(dxFileId);
    nxyz = sscanf(gridCountsLine, 'object 1 class gridpositions counts %d %d %d');
    nx = nxyz(1); 
    ny = nxyz(2); 
    nz = nxyz(3);
    originLine = fgetl(dxFileId);
    origin = sscanf(originLine, 'origin %f %f %f');
    dxLine = fgetl(dxFileId);
    dx = sscanf(dxLine, 'delta %f %f %f'); dx = dx(1);
    dyLine = fgetl(dxFileId);
    dy = sscanf(dyLine, 'delta %f %f %f'); dy = dy(2);
    dzLine = fgetl(dxFileId);
    dz = sscanf(dzLine, 'delta %f %f %f'); dz = dz(3);
    fclose(dxFileId);
    % Get potential information
    dxFileObj = importdata(dx_file, ' ', 11);
    dxFileData = dxFileObj.data';
    dxFileData = reshape(dxFileData, size(dxFileData,2)*3, 1); %put each potential value on a single row
    dxFileData(isnan(dxFileData)) = []; % Remove NaN entries
    potDX = reshape(dxFileData, nz, ny, nx); % nz * ny * nx
    potDX = permute(potDX, [2,3,1]); % ny * nx * nz -> more natural plotting
    %potDX = reshape(dxFileData, nx, ny, nz); 
    % The grid on x,y,z
    xGrid = origin(1) + (0:nx-1)*dx;
    yGrid = origin(2) + (0:ny-1)*dy;
    zGrid = origin(3) + (0:nz-1)*dz;

    surfPot = interp3(xGrid, yGrid, zGrid, potDX, PlyCoor(:,1),PlyCoor(:,2), PlyCoor(:,3), 'linear');

    output1=[PlyCoor,surfPot];
 
end
    
  

    




function [patchiness]=patchcorr_optimized(X,pcutoff,ncutoff,correl_l)

%q is the vector of potentials
%X is the matrix of coordinates on the surface to which q is mapped
%cutoff is the cutoff potential below which points are considered neutral
%correl_l is the search radius
Idx = rangesearch(X(:,1:3),X(:,1:3),correl_l,'SortIndices',false);
[total, count] = patchiness_helper_all(Idx,X, pcutoff, ncutoff); 
patchiness = total/count; 

end 
    





function [patchiness_q_f_radius]=patchcorr_optimized_wrt_idr(X,pcutoff,ncutoff,radius_start, correl_l, idr_x, idr_y, idr_z)
%numpatch is the integer number of patches
%area is a vector of the areas of each of the patches
%q is the vector of potentials
%X is the matrix of coordinates on the surface to which q is mapped
%cutoff is the cutoff potential below which points are considered neutral
%correl_l is the search radius

patchiness_q_f_radius = []; 

total_protein_length_approximation = norm(max(X(:,1:3)) - min(X(:,1:3)));

Idx = rangesearch(X(:,1:3),X(:,1:3),correl_l,'SortIndices',false);
[total, count] = patchiness_helper_all(Idx,X, pcutoff, ncutoff); 
patchiness_all = total/count; 


A = [idr_x, idr_y, idr_z];
A = A';
A = reshape(A,1,[]);

patchiness_curr_radius = 100; 
i = 0; 
while patchiness_curr_radius ~= patchiness_all
    
    nearby_points_idr = rangesearch(A,X(:,1:3),radius_start+i,'SortIndices',false);
    nearby_points_idr_idx = find(~cellfun(@isempty,nearby_points_idr)); %gets row numbers in X that are within radius_start+i within A 
    
    if ~isempty(nearby_points_idr_idx)
    
        %find points within patchiness search radius 
        %Idx{i}  --> points in X that are near nearby_points_idr_idx
        X_nearby = X(nearby_points_idr_idx,:);
        Idx = rangesearch(X_nearby(:,1:3),X(:,1:3),correl_l,'SortIndices',false);
        
        [total, count] = patchiness_helper_f_radius(Idx,nearby_points_idr_idx,X, pcutoff, ncutoff); 
        patchiness_curr_radius = total/count; 

        net_q_curr_radius = mean(X(nearby_points_idr_idx,4)); 

        matrix_val = [radius_start+i, patchiness_curr_radius, net_q_curr_radius, size(nearby_points_idr_idx,1), count, total];
        patchiness_q_f_radius = [patchiness_q_f_radius; matrix_val]; 

    end 
    i = i+1; 
end

end


