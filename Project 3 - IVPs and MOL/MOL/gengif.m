B = load("result-dense.txt");
t=0.01:0.01:0.99;
m = size(B,2);
f=figure;
for j=1:size(B,1)
    plot(t, B(j, 2:m));
    xlim([0,1]);
    ylim([-0.6,1.6]);
    drawnow
    %% 将所有的图像显示在同一个图窗中
    frame=getframe(f);
    i{j}=frame2im(frame);
    [A,map]=rgb2ind(i{j},256);
    %% 将上面合成后的图像导出为gif文件
    filename='heat.gif';
    if j==1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.04);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.04);
    end
end