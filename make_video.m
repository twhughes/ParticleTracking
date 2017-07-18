function [] = make_video(field_array)

    figure(); clf; colormap(redblue);
    scale = max(abs(field_array(:)));
    scale_factor = 0.2;    
    
    for i = (1:400)
        E = real(field_array*exp(-1i*i/10));
        imagesc(flipud(transpose([E;E;E;E;E;E;E;E])),[-scale*scale_factor,scale*scale_factor]);   
        pause(0.01); clf;
    end
    
end