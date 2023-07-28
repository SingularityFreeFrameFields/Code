function [] = VisualizeInitialization(Initialization, narrowBand)
%% Docstring
%{
 Takes in a field initialization and plots it!
%}
complexVecs = Initialization(1:(length(Initialization) / 2)) + 1i .* Initialization(((length(Initialization) / 2) + 1) : end);

realPart = Initialization(1:(length(Initialization) / 2));
imPart = Initialization(((length(Initialization) / 2) + 1) : end);

FrameField = zeros(length(complexVecs), 4);
for i = 1:length(complexVecs)
    FrameField(i, :) = transpose(roots([1 0 0 0 -complexVecs(i)]));
end
%FrameField = FrameField ./ (100 .* abs(FrameField));
[i,j] = find(narrowBand);

%%visualization
% hold on; axis equal;
% q1 = quiver(j, i, real(FrameField(:,1)), -imag(FrameField(:,1)), 0.2, 'blue');
% quiver(j, i, real(FrameField(:,2)), -imag(FrameField(:,2)), 0.2, 'blue');
% quiver(j, i, real(FrameField(:,3)), -imag(FrameField(:,3)), 0.2, 'blue');
% quiver(j, i, real(FrameField(:,4)), -imag(FrameField(:,4)), 0.2, 'blue');
% quiver(j, i, 1.5*realPart, -1.5*imPart, 0.4, 'magenta');
% %q1.ShowArrowHead = 'off';

% plot singularities
%plot(initial_singsxy(:,1),initial_singsxy(:,2),'mp')

