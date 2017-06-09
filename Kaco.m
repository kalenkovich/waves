% PARAMS
PARAMS.subject_name = 'Kaco';
PARAMS.protocol_name = 'waves';
PARAMS.freesurfer_folder = 'E:\freesurfer\subjects\kaco';
PARAMS.recordings_folder = 'E:\subjects_original_data\Kaco\kaco_spont01.ds';
PARAMS.max_distance = 0.04;

% Add Kaco's Freesurfer folder to the Brainstorm database
% Then add the recordings folder

from_bst_initialize(PARAMS);