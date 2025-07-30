function visualise_nii(file)

data = load_nii(file);
origin = abs(data.hdr.hist.originator(1:3));
img = data.img;
[Nx,Ny,Nz] = size(img);
if false
    for i = 1:Nx
        sliced = abs(img(i,:,:));
        if sum(sliced,[1,2,3])> 0
            min_x = i -1;
            break
        end
    end
    for j = 0:(Nx-1)
        i = Nx - j;
        sliced = abs(img(i,:,:));
        if sum(sliced,[1,2,3])> 0
            max_x = i +1;
            break
        end
    end
    for i = 1:Ny
        sliced = abs(img(:,i,:));
        if sum(sliced,[1,2,3])> 0
            min_y = i -1;
            break
        end
    end
    for j = 0:(Ny-1)
        i = Ny - j;
        sliced = abs(img(:,i,:));
        if sum(sliced,[1,2,3])> 0
            max_y = i+1;
            break
        end
    end
    for i = 1:Nz
        sliced = abs(img(:,:,i));
        if sum(sliced,[1,2,3])> 0
            min_z = i -1;
            break
        end
    end
    for j = 0:(Nz-1)
        i = Nz - j;
        sliced = abs(img(:,:,i));
        if sum(sliced,[1,2,3])> 0
            max_z = i+1;
            break
        end
    end
    
    clear options
    options.cut_from_L = min_x;
    options.cut_from_R = max_x;
    options.cut_from_P =min_y*(min_y < origin(2)) + 0*(min_y >= origin(2));
    options.cut_from_A = max_y;
    options.cut_from_I = min_z;
    options.cut_from_S = max_z;
    
    data = clip_nii(data,options);
end


view_nii(data);