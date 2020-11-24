% GETJOINTPARAMS Assemble a structure with all the information required to
% create a CustomJoint of a specified lower limb joint. Normally this
% function is used after the geometrical analyses on the bones and before
% generating the joint for the automatic OpenSim model using 
% createCustomJointFromStruct. It is assumed that the inputs will contain
% enough information (location and orientation) to define the joint
% reference system. NOTE: A body is connected to ground with a 
% free_to_ground joint if no other specifics are provided, please see
% examples on partial models for a practical example.
% IMPORTANT: modifying the values of the fields of JointParamsStruct output
% structure allows to modify the joint model according to the preferences
% of the researcher.
%
% JointParamsStruct = getJointParams(joint_name, parentBodyStruct, childBodyStruct)
%
% Inputs:
%   joint_name - name of the lower limb joint for which we want to create
%       the structure containing all parameters (string).
%
%   parentBodyStruct - a MATLAB structure collecting the information
%       derived from analysing the parent bone. This structure normally is
%       organised as
%       parentBodyStruct.aJointName.parent/child_location/orientation and
%       it is produced by using the algorithms to analyse the bone
%       geometries.
%
%   childBodyStruct - exactly the same as parentBodyStruct, but for the
%       structure obtained analysing the bone geometry that serves as child
%       of the joint.
%
% Outputs:
%   JointParamsStruct - a structure collecting all the information required
%       to define an OpenSim CustomJoint. The typical fields of this
%       structure are the following: name, parent, child, coordsNames,
%       coordsTypes, ROM and rotationAxes. An example of JointParamsStruct is
%       the following:
%         JointParamsStruct.name                = 'hip_r';
%         JointParamsStruct.parent              = 'pelvis';
%         JointParamsStruct.child               = 'femur_r';
%         JointParamsStruct.coordsNames         = {'hip_flexion_r','hip_adduction_r','hip_rotation_r'};
%         JointParamsStruct.coordsTypes         = {'rotational', 'rotational', 'rotational'};
%         JointParamsStruct.coordRanges         = {[-120 120], [-120 120], [-120 120]};% in degrees 
%         JointParamsStruct.rotationAxes        = 'zxy';
%
% See also CREATECUSTOMJOINTFROMSTRUCT, CREATELOWERLIMBJOINTS.
%
%-------------------------------------------------------------------------%
%  Author:   Luca Modenese
%  Copyright 2020 Luca Modenese
%-------------------------------------------------------------------------%

function JointParamsStruct = getJointParams(joint_name, parentBodyStruct, childBodyStruct, side)

%% Ensure that location and orientation are set on both parent and child
% empty input as parentBodyCS means all set to zero
if isempty(parentBodyStruct)
    parentBodyStruct.(joint_name).parent_location     = [0.0000	0.0000	0.0000];
    parentBodyStruct.(joint_name).parent_orientation  = [0.0000	0.0000	0.0000];
end
% if joint_name is unavailable set it to empty (child)
if ~isfield(childBodyStruct, joint_name)
    childBodyStruct.(joint_name) = [];
end
% if joint_name is unavailable set it to empty (parent)
if ~isfield(parentBodyStruct, joint_name)
    parentBodyStruct.(joint_name) = [];
end

% if child_location is unavailable use parent
if isfield(childBodyStruct.(joint_name), 'child_location')
    JointParamsStruct.child_location   = childBodyStruct.(joint_name).child_location;
else
    JointParamsStruct.child_location   = parentBodyStruct.(joint_name).parent_location;
end

% if parent_location is unavailable use child
if isfield(parentBodyStruct.(joint_name), 'parent_location')
    JointParamsStruct.parent_location   = parentBodyStruct.(joint_name).parent_location;
else
    JointParamsStruct.parent_location   = childBodyStruct.(joint_name).child_location;
end

% if child_orientation is unavailable use parent
if isfield(childBodyStruct.(joint_name), 'child_orientation')
    JointParamsStruct.child_orientation   = childBodyStruct.(joint_name).child_orientation;
else
    JointParamsStruct.child_orientation   = parentBodyStruct.(joint_name).parent_orientation;
end

% if parent_orientation is unavailable use child
if isfield(parentBodyStruct.(joint_name), 'parent_orientation')
    JointParamsStruct.parent_orientation  = parentBodyStruct.(joint_name).parent_orientation;
else
    JointParamsStruct.parent_orientation  = childBodyStruct.(joint_name).child_orientation;
end

%% assign the parameters required to create a CustomJoint
switch joint_name
    case 'ground_pelvis'
        JointParamsStruct.name                = 'ground_pelvis';
        JointParamsStruct.parent              = 'ground';
        JointParamsStruct.child               = 'pelvis';
        JointParamsStruct.coordsNames         = {'pelvis_tilt','pelvis_list','pelvis_rotation', 'pelvis_tx','pelvis_ty', 'pelvis_tz'};
        JointParamsStruct.coordsTypes         = {'rotational', 'rotational', 'rotational', 'translational', 'translational','translational'};
        JointParamsStruct.coordRanges         = {[-90 90], [-90 90] , [-90 90], [-10, 10] , [-10, 10] , [-10, 10]};
        JointParamsStruct.rotationAxes        = 'zxy';
    case 'free_to_ground'
        cb                                    = childBodyStruct.free_to_ground.child; % cb = current bone (for brevity)
        JointParamsStruct.name                = ['ground_',cb];
        JointParamsStruct.parent              = 'ground';
        JointParamsStruct.child               = cb;
        JointParamsStruct.coordsNames         = {['ground_', cb,'_rz'],['ground_', cb,'_rx'],['ground_', cb,'_ry'],...
                                                 ['ground_', cb,'_tx'],['ground_', cb,'_ty'],['ground_', cb,'_tz']};
        JointParamsStruct.coordsTypes         = {'rotational', 'rotational', 'rotational', 'translational', 'translational','translational'};
        JointParamsStruct.coordRanges         = {[-120 120], [-120 120] , [-120 120], [-10, 10] , [-10, 10] , [-10, 10]};
        JointParamsStruct.rotationAxes        = 'zxy';  
    case ['hip_', side]
        JointParamsStruct.name                = ['hip_', side];
        JointParamsStruct.parent              = 'pelvis';
        JointParamsStruct.child               = ['femur_', side];
        JointParamsStruct.coordsNames         = {['hip_flexion_', side],['hip_adduction_', side],['hip_rotation_', side]};
        JointParamsStruct.coordsTypes         = {'rotational', 'rotational', 'rotational'};
        JointParamsStruct.coordRanges         = {[-120 120], [-120 120], [-120 120]};
        JointParamsStruct.rotationAxes        = 'zxy';
    case ['knee_', side]
        JointParamsStruct.name               = ['knee_', side];
        JointParamsStruct.parent             = ['femur_', side];
        JointParamsStruct.child              = ['tibia_', side];
        JointParamsStruct.coordsNames        = {['knee_angle_', side]};
        JointParamsStruct.coordsTypes        = {'rotational'};
        JointParamsStruct.coordRanges        = {[-120 10]};
        JointParamsStruct.rotationAxes       = 'zxy';   
    case ['ankle_', side]
        JointParamsStruct.name               = ['ankle_', side];
        JointParamsStruct.parent             = ['tibia_', side];
        JointParamsStruct.child              = ['talus_', side];
        JointParamsStruct.coordsNames        = {['ankle_angle_', side]};
        JointParamsStruct.coordsTypes        = {'rotational'};
        JointParamsStruct.coordRanges        = {[-90 90]};
        JointParamsStruct.rotationAxes       = 'zxy';
    case ['subtalar_', side]
        JointParamsStruct.name               = ['subtalar_', side];
        JointParamsStruct.parent             = ['talus_', side];
        JointParamsStruct.child              = ['calcn_', side];
        JointParamsStruct.coordsNames        = {['subtalar_angle_', side]};
        JointParamsStruct.coordsTypes        = {'rotational'};
        JointParamsStruct.coordRanges        = {[-90 90]};
        JointParamsStruct.rotationAxes       = 'zxy';
    case 'patellofemoral_r'
        JointParamsStruct.name               = ['patellofemoral_', side];
        JointParamsStruct.parent             = ['femur_', side];
        JointParamsStruct.child              = ['patella_', side];
        JointParamsStruct.coordsNames        = {['knee_angle_',side,'_beta']};
        JointParamsStruct.coordsTypes        = {'rotational'};
%         JointParamsStruct.coordRanges            = {[-90 90]};
        JointParamsStruct.rotationAxes       = 'zxy';
    otherwise
        error(['getJointParams.m Unsupported joint ',joint_name ,'.']);
end