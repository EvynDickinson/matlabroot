
function separateAxes

    ax = gca;
    % set the X axis vertex start to its the orignial point
    origin = get(ax.XAxis.Axle,'VertexData');
    origin(1,1) = ax.XTick(1);
    set(ax.XAxis.Axle,'VertexData',origin);

    % set the Y axis vertex start to its the orignial point
    origin = get(ax.YAxis.Axle,'VertexData');
    origin(2,1) = ax.YTick(1);
    set(ax.YAxis.Axle,'VertexData',origin);

end
