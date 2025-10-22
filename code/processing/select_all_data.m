function [user_path, base_path, func_nii, t1_nii, t1_nii2, ho_gz, json_path] = select_all_data()
% Inputs:
%   (none; all paths are chosen interactively via GUI)
%
% Outputs:
%   user_path (char)   — chosen user root
%   base_path (char)   — dataset root
%   func_nii (char)    — path to selected functional EPI NIfTI
%   t1_nii (char)      — path to T1-weighted anatomical, run-01
%   t1_nii2 (char)     — path to T1-weighted anatomical, run-02
%   ho_gz (char)       — path to Harvard–Oxford atlas NIfTI
%   json_path (char)   — path to the BIDS JSON sidecar for the selected EPI task (TR, slice timing, etc.)
%
% Purpose:
%   Interactive selector that gathers all paths needed for the pipeline:
%   (1) user/session root, (2) dataset base, (3) EPI (functional), (4) two T1 runs,
%   (5) Harvard–Oxford atlas, and (6) the EPI’s BIDS JSON sidecar. 

    % --- USER (start at OS home) ---                                           
    if ispc, home_dir = getenv('USERPROFILE'); else, home_dir = getenv('HOME'); end 
    if isempty(home_dir) || ~isfolder(home_dir), home_dir = pwd; end                
    user_path = uigetdir(home_dir, 'Select USER path');                             % interactive selection of paths.
    if isequal(user_path,0), error('No USER path selected.'); end                   % fail downstream paths that depend on this root.

    % --- BASE (inside USER) ---                                                    % define dataset root that contain the EPI task
    base_path = uigetdir(user_path, 'Select BASE dataset folder');                  
    if isequal(base_path,0), error('No BASE path selected.'); end                   

    % --- EPI (inside BASE) ---                                                     % select the functional run that the SW will run analysis for
    [fE,pE] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...           
                         'Select functional (EPI) file', base_path);
    if isequal(fE,0), error('No EPI selected.'); end                               
    func_nii = fullfile(pE,fE);                                                     

    % --- Infer SUBJECT directory from EPI path (walk up to a 'sub-*' folder) --- 
    cur = fileparts(func_nii);                                                      % folder containing EPI. start from EPI’s directory to search upwards.
    for i = 1:6                                                                    
        [parent, name] = fileparts(cur);                                           
        if startsWith(name,'sub-','IgnoreCase',true)                                
            subject_dir = cur; break;                                               % once found, this anchors subsequent searches.
        end
        if isempty(parent) || strcmp(parent,cur)                                    % if reached filesystem root or stuck without finding the folder; stop to avoid infinite loop.
            subject_dir = fileparts(func_nii); break;                               % fall back to EPI directory
        end
        cur = parent;                                                              
    end

    % get the subject root if it exists                                            
    anat_root = fullfile(subject_dir,'anat');                                    
    if ~isfolder(anat_root), anat_root = subject_dir; end                           

    % --- T1 picks ---                                                              %selecting both runs supports later decisions (choose best or compute average).
    [f1,p1] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, 'Select T1 run-01', anat_root);  % Notes: start browsing at anat_root to reduce navigation time.
    if isequal(f1,0), error('No T1 run-01 selected.'); end                          
    t1_nii = fullfile(p1,f1);                                                       

    [f2,p2] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, 'Select T1 run-02', anat_root);  
    if isequal(f2,0), error('No T1 run-02 selected.'); end                         
    t1_nii2 = fullfile(p2,f2);                                                    

    % --- Harvard–Oxford atlas ---    
    [ho_file, ho_path] = uigetfile({'*.nii;*.nii.gz','NIfTI (*.nii, *.nii.gz)'}, ...
                                   'Select Harvard-Oxford atlas', user_path);
    if isequal(ho_file,0), error('No atlas selected.'); end                         
    ho_gz = fullfile(ho_path, ho_file);                                            

    % --- Derive BIDS JSON sidecar from the functional file ---                   % Notes: JSON carries TR, slice timing, phase encoding, etc., essential for STC & preprocessing.
    [func_dir, func_base] = fileparts(func_nii);                
    if endsWith(func_base, '.nii')                               
        [~,func_base] = fileparts(func_nii);                     
        [func_dir, ~, ~] = fileparts(func_nii);                                                                         
    end
    json_path = fullfile(func_dir, [func_base '.json']);                               

    % If not found, let user pick manually (needed e.g. for slice timing)         
    if ~isfile(json_path)
        [jfile, jpath] = uigetfile('*.json', 'Select the BIDS JSON for this run', func_dir);  
        if isequal(jfile,0)
            error('JSON file is required. Aborting.');                              
        end
        json_path = fullfile(jpath, jfile);                                         
    end
end
