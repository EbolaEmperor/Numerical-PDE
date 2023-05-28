B = load("result.txt");
t=-2:0.05:25;
m = size(B,2);
f=figure;
for j=1:size(B,1)
    %% 将y=2*cos(x+n(idx))的图像绘制出来
    plot(t, B(j, 2:m));
    xlim([-2,25]);
    ylim([-0.4,1.1]);
    drawnow
    %% 将所有的图像显示在同一个图窗中
    frame=getframe(f);
    i{j}=frame2im(frame);
    [A,map]=rgb2ind(i{j},256);
    %% 将上面合成后的图像导出为gif文件
    filename='advection7.gif';
    if j==1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.04);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
    end
end