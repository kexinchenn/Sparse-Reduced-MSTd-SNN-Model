clear variables
%%
% 8 direction and 5 speed
numSpeeds = 5;
vecSpeeds = 2.^(linspace(log2(0.5), log2(16), numSpeeds));
numDirs = 8;
vecDirs = (0:numDirs-1)/numDirs*2*pi;

[flow, MT, MSTd] = setupModel(vecSpeeds, vecDirs);
%% direction tuning
x = [0:0.01:2*pi];

for d=vecDirs
    resp = MT.fncDir(x,d);
    plot(x,resp,'LineWidth',2)
    hold on
end
hold off
xlim([0 2*pi])
% yticks([0:0.2:1])
xticks(linspace(0,2*pi,8))
xticklabels({'0','pi/4','pi/2','3pi/4','pi','5pi/4','3pi/2','2pi'})
xlabel('Direction')
ylabel('MT Response')
%% speed tuning
x = [0:0.01:100];
for s=vecSpeeds
    resp = MT.fncSpeed(x,s);
    plot(log2(x),resp,'LineWidth',2)
    hold on
% pause
end
hold off
xlim([log2(0.25), log2(32)])
xticks(linspace(log2(0.5), log2(16), numSpeeds))
xticklabels(linspace(0.5, 16, numSpeeds))
yticks([0:0.2:1])
xlabel('Speed (deg/sec)')
ylabel('MT Response')