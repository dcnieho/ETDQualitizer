close all

% Often, eye trackers report gaze position as a location on the screen.
% When converting this position to an eye orientation, this position is
% first transformed to a vector with respect to a reference position and
% this on-screen vector is then transformed to an eye orientation using
% trigonometry. This calculation has the assumption that the axis from the
% observer's eye to the reference position on the screen is perpendicular
% to the screen. These assumptions in practice often do not hold, for
% instance because the reference position is not right ahead of the eye but
% shifted vertically, or because the screen is tilted around a horizontal
% axis. With this simulation we calculate the error that violation of
% this assumption introduces in the estimated gaze angle.
%
% We perform this simulation in 1D (imagine vertical) as that is
% sufficient to illustrate the magnitude of errors introduced.
%
% What we compute is the following. We take a fixed on-screen gaze position
% on the screen as indicated by an eye tracker. The position corresponds to
% a 10 deg gaze angle with respect to the reference position (screen
% center) in the condition where the assumption holds, i.e., the eye is
% positioned right in front of the reference position and the screen is not
% tilted. We then either shift or tilt the screen and compute what the
% required _actual_ gaze angle is that would produce the same on-screen
% gaze position. The difference between the gaze angle that would be
% computed given the screen position (always 10 deg) and the actual gaze
% position that would produce this gaze position is the error due to
% violation of the given assumption.


screenDistance = 650;   % the screen is positioned 650 mm aways from the eye
gazeDirection = 10;     % the gaze direction is 10 deg below the reference position on the screen

shiftRange = 200;       % the screen is shifted vertically by up to 200 mm
tiltRange = 20;         % the screen is titled along a horizontal axis by up to 20 degrees

nstep = 151;
ylims = [-1.3 .2];

% compute the gazePosition on the screen w.r.t. the reference position if
% the assumptions hold
gazePosition = screenDistance*tand(gazeDirection);


% When the screen is shifted, compute the gaze rotation from the reference
% position that yields the same on-screen gaze position
shifts = linspace(-shiftRange,shiftRange,nstep);
shiftAngles = atand((gazePosition+shifts)/screenDistance)-atand(shifts/screenDistance);


% when the screen is tilted, compute the gaze rotation from the reference
% position that yields the same on-screen gaze position
tilts = linspace(-tiltRange,tiltRange,nstep);
viewingDistance = sqrt(screenDistance^2+gazePosition^2-2*screenDistance*gazePosition*cosd(90+tilts));   % length of gazeVector to screen, using law of cosines
tiltAngles = asind(gazePosition*sind(90+tilts)./viewingDistance);   % using law of sines


% plot actual gaze directions when the screen is shifted
figure
plot(shifts([1 end]),[0 0],'Color',[.7 .7 .7],'LineStyle','--','LineWidth',1.2)
hold on
plot(shifts,shiftAngles-gazeDirection,'k','LineWidth',2)
axis tight
ax = gca;
box off
ax.YAxis.TickLabelFormat = '%.1f';
ax.XLabel.FontSize = 16;
ax.XAxis.FontSize = 13;
ax.YLabel.FontSize = 16;
ax.YAxis.FontSize = 13;
xlabel('Screen shift (mm)')
ylabel('Error in estimated gaze angle (deg)')
ylim(ylims)
yyaxis right
ylabel('Error in estimated gaze angle (%)')
ax.YLabel.FontSize = 16;
ax.YAxis(2).TickLabelFormat = '% .0f';
ax.YAxis(2).FontSize = 13;
ax.YAxis(2).Color = [0 0 0];
ylim(ylims/gazeDirection*100)
print('error_shift.png','-dpng','-r300')

% illustrate the geometry of the two extreme cases and the case when the
% assumption holds
figure
screen = 150;
subplot(3,1,1), hold on, axis ij
s = 1;
title(sprintf('screen shift = %d mm',shifts(s)))
plot(0,0,'ko')
plot(screenDistance*[1 1],[-screen screen]+shifts(s),'k-')
plot([0 screenDistance],[0 shifts(s)],'b-')
plot(screenDistance*[1 1],[0 gazePosition]+shifts(s),'c-')
plot(screenDistance,gazePosition+shifts(s),'ro')
plot([0 screenDistance],[0 gazePosition+shifts(s)],'r-')
xlim([0 screenDistance])
ylim([-screen+min(shifts) screen+max(shifts)])

subplot(3,1,2), hold on, axis ij
s = ceil(length(shifts)/2);
title(sprintf('screen shift = %d mm',shifts(s)))
plot(0,0,'ko')
plot(screenDistance*[1 1],[-screen screen]+shifts(s),'k-')
plot([0 screenDistance],[0 shifts(s)],'b-')
plot(screenDistance*[1 1],[0 gazePosition]+shifts(s),'c-')
plot(screenDistance,gazePosition+shifts(s),'ro')
plot([0 screenDistance],[0 gazePosition+shifts(s)],'r-')
xlim([0 screenDistance])
ylim([-screen+min(shifts) screen+max(shifts)])

subplot(3,1,3), hold on, axis ij
s = length(shifts);
title(sprintf('screen shift = %d mm',shifts(s)))
plot(0,0,'ko')
plot(screenDistance*[1 1],[-screen screen]+shifts(s),'k-')
plot([0 screenDistance],[0 shifts(s)],'b-')
plot(screenDistance*[1 1],[0 gazePosition]+shifts(s),'c-')
plot(screenDistance,gazePosition+shifts(s),'ro')
plot([0 screenDistance],[0 gazePosition+shifts(s)],'r-')
xlim([0 screenDistance])
ylim([-screen+min(shifts) screen+max(shifts)])


% plot actual gaze directions when the screen is tilted
figure
plot(tilts([1 end]),[0 0],'Color',[.7 .7 .7],'LineStyle','--','LineWidth',1.2)
hold on
plot(tilts,tiltAngles-gazeDirection,'k','LineWidth',2)
axis tight
ax = gca;
box off
ax.YAxis.TickLabelFormat = '%.1f';
ax.XLabel.FontSize = 16;
ax.XAxis.FontSize = 13;
ax.YLabel.FontSize = 16;
ax.YAxis.FontSize = 13;
xlabel('Screen tilt (deg)')
ylabel('Error in estimated gaze angle (deg)')
ylim(ylims)
yyaxis right
ylabel('Error in estimated gaze angle (%)')
ax.YLabel.FontSize = 16;
ax.YAxis(2).TickLabelFormat = '% .0f';
ax.YAxis(2).FontSize = 13;
ax.YAxis(2).Color = [0 0 0];
ylim(ylims/gazeDirection*100)
print('error_tilt.png','-dpng','-r300')

% illustrate the geometry of the two extreme cases and the case when the
% assumption holds
figure
subplot(3,1,1), hold on, axis ij
s = 1;
title(sprintf('screen tilt = %d deg',tilts(s)))
plot(0,0,'ko')
plot(screenDistance+[-1 1]*screen*sind(tilts(s)), [-1 1]*screen*cosd(tilts(s)), 'k-')
plot(screenDistance+gazePosition*sind(tilts(s)),gazePosition*cosd(tilts(s)),'ro')
plot([0 screenDistance+gazePosition*sind(tilts(s))],[0 gazePosition*cosd(tilts(s))],'r-')
plot([0 screenDistance],[0 0],'b-')
plot(screenDistance+[0 gazePosition*sind(tilts(s))],[0 gazePosition*cosd(tilts(s))],'c-')
ylim([-1 1]*screen);
xlim([0 screenDistance+screen*sind(tiltRange)])

subplot(3,1,2), hold on, axis ij
s = ceil(length(shifts)/2);
title(sprintf('screen tilt = %d deg',tilts(s)))
plot(0,0,'ko')
plot(screenDistance+[-1 1]*screen*sind(tilts(s)), [-1 1]*screen*cosd(tilts(s)), 'k-')
plot(screenDistance+gazePosition*sind(tilts(s)),gazePosition*cosd(tilts(s)),'ro')
plot([0 screenDistance+gazePosition*sind(tilts(s))],[0 gazePosition*cosd(tilts(s))],'r-')
plot([0 screenDistance],[0 0],'b-')
plot(screenDistance+[0 gazePosition*sind(tilts(s))],[0 gazePosition*cosd(tilts(s))],'c-')
ylim([-1 1]*screen);
xlim([0 screenDistance+screen*sind(tiltRange)])

subplot(3,1,3), hold on, axis ij
s = length(shifts);
title(sprintf('screen tilt = %d deg',tilts(s)))
plot(0,0,'ko')
plot(screenDistance+[-1 1]*screen*sind(tilts(s)), [-1 1]*screen*cosd(tilts(s)), 'k-')
plot(screenDistance+gazePosition*sind(tilts(s)),gazePosition*cosd(tilts(s)),'ro')
plot([0 screenDistance+gazePosition*sind(tilts(s))],[0 gazePosition*cosd(tilts(s))],'r-')
plot([0 screenDistance],[0 0],'b-')
plot(screenDistance+[0 gazePosition*sind(tilts(s))],[0 gazePosition*cosd(tilts(s))],'c-')
ylim([-1 1]*screen);
xlim([0 screenDistance+screen*sind(tiltRange)])
