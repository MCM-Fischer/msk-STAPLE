function [BCS, JCS, PelvisBL] = Fischer2019_pelvis(Pelvis, side_raw, ...
    result_plots, debug_plots, in_mm)
%FISCHER2019_pelvis Calculation of the pelvic CS based on [Fischer 2019].
%
%   This is a wrapper function to integrate the MATLAB function: 
%       https://github.com/RWTHmediTEC/PelvicLandmarkIdentification
%   
%   The method is described in detail in:
%       Fischer, M. C. M. et al. A robust method for automatic 
%       identification of landmarks on surface models of the pelvis. 
%       Sci. Rep. 9, 13322 (2019)
%       https://doi.org/10.1038/s41598-019-49573-4
%
%   ATTENTION: Consider the requirements for the mesh of the pelvis 
%   described in the document header of pelvicLandmarkID.
%
% AUTHOR: Maximilian C. M. Fischer
% VERSION: 1.0.0
% DATE: 2024-08-19
% COPYRIGHT (C) 2024 Maximilian C. M. Fischer
% LICENSE: CC BY-NC 4.0
%

if nargin<2;    side_raw = 'r';    end
if nargin<3;    result_plots = 1;  end
if nargin<4;    debug_plots = 0;   end
if nargin<5;    in_mm = 1;         end
if in_mm == 1;  dim_fact = 0.001;  else;  dim_fact = 1; end

% get side id correspondent to body side (used for hip joint parent)
% no need for sign, left and right rf are identical
[~, side_low] = bodySide2Sign(side_raw);

disp('---------------------')
disp('  FISCHER - PELVIS   '); 
disp('---------------------')
disp(['* Hip Joint   : ', upper(side_low)]);
disp(['* Method      : ', 'https://doi.org/10.1038/s41598-019-49573-4']);
disp(['* Result Plots: ', convertBoolean2OnOff(result_plots)]);
disp(['* Debug  Plots: ', convertBoolean2OnOff(debug_plots)]);
disp(['* Triang Units: ', 'mm']);
disp('---------------------')
disp('Initializing method...')

% inertial axes
[~, CenterVol, InertiaMatrix] =  TriInertiaPpties(Pelvis);

pelvisMesh.vertices = Pelvis.Points;
pelvisMesh.faces = Pelvis.ConnectivityList;

disp('Landmarking...')

[~, Landmarks] = pelvicLandmarkID(pelvisMesh, 'visu',0, 'CS','APP' ,'debug',debug_plots);

LASIS = Landmarks.ASIS(1,:);
RASIS = Landmarks.ASIS(2,:);
LPSIS = Landmarks.PSIS(1,:);
RPSIS = Landmarks.PSIS(2,:);
SYMP = Landmarks.PS;

% Sanity check if bone landmarks were correctly identified
if norm(RASIS-LASIS)<norm(RPSIS-LPSIS)
    % inform user
    disp('GIBOK_pelvis.')
    warndlg('Inter-ASIS distance is shorter than inter-PSIS distance. Better check manually.')
end

% ISB reference system
PelvisOr = (RASIS+LASIS)'/2.0;

% segment reference system
BCS.CenterVol = CenterVol;
BCS.Origin = PelvisOr;
BCS.InertiaMatrix = InertiaMatrix;
BCS.V = CS_pelvis_ISB(RASIS, LASIS, RPSIS, LPSIS);
% CS.V = RotPseudoISB2Glob;

% storing joint details
JCS.ground_pelvis.V = BCS.V;
JCS.ground_pelvis.Origin = PelvisOr;
JCS.ground_pelvis.child_location    = PelvisOr'*dim_fact; % [1x3] as in OpenSim
JCS.ground_pelvis.child_orientation = computeXYZAngleSeq(BCS.V); % [1x3] as in OpenSim

% define hip parent
hip_name = ['hip_', side_low];
JCS.(hip_name).parent_orientation   = computeXYZAngleSeq(BCS.V);

% Export bone landmarks: [3x1] vectors
PelvisBL.RASI     = RASIS'; 
PelvisBL.LASI     = LASIS'; 
PelvisBL.RPSI     = RPSIS'; 
PelvisBL.LPSI     = LPSIS'; 
PelvisBL.SYMP     = SYMP';

% debug plot
label_switch = 1;
if result_plots == 1
    figure('Name', ['STAPLE | bone: pelvis | side: ', side_low])
    plotTriangLight(Pelvis, BCS, 0); hold on
    quickPlotRefSystem(BCS);
    quickPlotRefSystem(JCS.ground_pelvis);
    trisurf(triangulation([1 2 3],[RASIS; SYMP; LASIS]),...
        'facealpha',0.4,'facecolor','y','edgecolor','k');
    
    % plot markers and labels
    plotBoneLandmarks(PelvisBL, label_switch);
end

% final printout
disp('Done.');

end