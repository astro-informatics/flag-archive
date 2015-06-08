function [] = flagmovie( f, file, L, P, layers, zoomfactor, caxisvals, showslices )

fig = figure('Renderer','zbuffer', 'Position',[0 0 650 550],'Color',[1 1 1]);
set(gca,'NextPlot','replaceChildren');
writerObj = VideoWriter( file );
writerObj.FileFormat
writerObj.FrameRate = 4;
writerObj.Quality = 100;
open(writerObj);

nblayers = length(layers)
for ilayer = 1:nblayers
    layer = layers(ilayer)
    fb = f;
    flaglet_plot_f(fb, 'L', L, 'P', P, 'layer', layer, 'ShowSlices', showslices )
    if ilayer == 1
        zoom(zoomfactor)
    end
    colormap(jet)
    %colormap(hot(256))
    caxis(caxisvals)
    v = axis;
    axis([-1 1 -1 1 -1 1]);
    F(ilayer) = getframe;
    frame = getframe;
    writeVideo(writerObj,frame);    
end
close(writerObj);

end