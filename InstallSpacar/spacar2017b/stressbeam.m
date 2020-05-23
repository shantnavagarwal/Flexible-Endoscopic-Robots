function [stress, mesh, ind, stressExtrema] = stressbeam(sbd, elem_nrs, time_steps, varargin)

%ELEMSTRESS Computes the stresses in 2D or 3D beam type element(s) with 
%   cross-sectional data read from the DAT file.
%
%WHEN CALLING STRESSBEAM DIRECTLY, USE ONE OF THE FOLLOWING INPUT SYNTAX's:
%
%   [STRESS, MESH] = STRESSBEAM(FILENAME.SBD, ELEM_NUMBERS, TIME_STEPS) 
%   Computes the stresses in beam type element(s) defined by ELEM_NUMBERS, 
%   on time steps defined by TIME_STEPS. If TIME_STEPS is empty, the last 
%   time step or load step will be used to compute the stresses. 
%   The outputs MESH and STRESS are structure arrays, where the indices 
%   correspond to the n'th component of ELEM_NUMBERS. The structure for the
%   n'th supplied element number is organized as follows:
%   MESH(n).Interior.x        (Matrix containing the x-positions)
%   MESH(n).Interior.y        (Matrix containing the y-positions)
%   MESH(n).Interior.z        (Matrix containing the z-positions)
%   STRESS(n).Interior.       (Structure containing the stress components)
%                     .Mx     (Torsional shear stresses about x-axis)
%                     .My     (Bending normal stresses about y-axis)
%                     .Mz     (Bending normal stresses about z-axis)
%                     .Fx     (Normal stresses in x-dir induced by Fx)
%                     .Fy     (Shear stresses in y-dir induced by Fy)
%                     .Fz     (Shear stresses in z-dir induced by Fz) 
%                     .userOutput (user defined stress output, defaults to
%                             to the mises stresses)
%   Here MESH(n).Interior corresponds to the locations of the cross-
%   sectional stresses in STRESS(n).Interior, along the length of the 
%   element. The mesh and stresses on the external surfaces, using the same
%   format as for the interior cross-sectional stresses, are output in: 
%   MESH(n).Surf.Front,     STRESS(n).Surf.Front,
%   MESH(n).Surf.Back,      STRESS(n).Surf.Back,
%   MESH(n).Surf.Outer,     STRESS(n).Surf.Outer,
%   MESH(n).Surf.Inner,     STRESS(n).Surf.Inner.
%   Where MESH(n).Surf.Inner is empty for rectangular or circular cross-
%   sections because these do not have an inner external surface, unlike 
%   there hollow variants. 
%   Incase elements with no cross-section (cross-section type: 'line') are 
%   considered, the output is as follows:
%   MESH(n).x, MESH(n).y, MESH(n).z, STRESS(n).Mx STRESS(n).My ...etc
%   Where STRESS contains the stress resultants and not the actual
%   stresses since a line does not have a cross-section!
%
%   [STRESS, MESH] = STRESSBEAM(FILENAME.SBD, ELEM_NUMBERS, TIME_STEPS, OPTIONS) 
%   Computes the stresses based on the options supplied in the structure 
%   OPTIONS. 
%   Valid fields and inputs for OPTIONS are:
%       OPTIONS.type - Changes the data written to the STRESS.(...).userOutput 
%           field, to: 
%           'mises' (outputs the von Mises stresses, default)
%           'sigX'  (outputs the sum normal stresses in x-dir)
%           'tauXY' (outputs the sum shear stresses in y-dir)
%           'tauXZ' (outputs the sum stear stresses in z-dir)
%           'Mx'    (outputs the torsional shear stresses about x-dir only)
%           'My'    (outputs the normal bending stresses about y-dir only)
%           'Mz'    (outputs the normal bending stresses about z-dir only)
%           'Fx'    (outputs the normal stress in x-dir induced by Fx only)
%           'Fy'    (outputs the shear stress in y-dir induced by Fy only) 
%           'Fz'    (outputs the shear stress in z-dir induced by Fz only)
%       OPTIONS.exterior   - Computes stresses at the external surfaces 
%           only: true (1) or false (0, default)
%       OPTIONS.density.tot - The total number of stress evaluations per 
%           element. It can also be one of the following strings:
%           'coarse' (~1e2 evaluations)
%           'normal' (~1e3 evaluations, default)
%           'fine'   (~1e4 evaluations)
%       OPTIONS.density.x   - The number of stress evaluations in the x-dir
%       OPTIONS.density.y   - The number of stress evaluations in the y-dir
%       OPTIONS.density.z   - The number of stress evaluations in the z-dir
%       OPTIONS.density.r   - The number of stress evaluations in the 
%                            radial direction (for circular c. section)
%       OPTIONS.density.phi - The number of stress evaluations in the
%                            circular direction (for circular c. section)
%   The fields density.tot, density.x, density.y and density.z can be used 
%   together as can density.tot, density.r and density.phi.
%
%   [STRESS, MESH, STRESSIND] = STRESSBEAM(FILENAME.SBD, ELEM_NUMBERS, ...)
%   Outputs the STRESSIND array which links the element numbers to the
%   positions in the STRESS and MESH structure arrays. 
%
%   [STRESS, MESH, STRESSIND, STRESSEXTREMA] = STRESSBEAM(FILENAME.SBD, ...)
%   Outputs the stress extrema's in the structure STRESSEXTREMA, with
%   fields:
%       STRESSEXTREMA.max = maximum occuring stress
%       STRESSEXTREMA.min = minimum occuring stress

%SPAVISUAL USAGE OF STRESSBEAM:
%   
%   [STRESS, MESH, ...] = STRESSBEAM(FILENAME.SBD, ..., CTYPE, CDIM)
%   Computes the stresses with cross-sectional type CTYPE and dimensions
%   CDIM. CTYPE and CDIM are cell arrays containing the cross-sectional
%   type data and dimensions determined from GETINFO.m.

% Created by S.E. Boer

%% Load necessary data
if isa(sbd,'struct') % calling syntax from within SPAVISUAL/SPAGUI
    it = sbd.it;
    time_steps = sbd.time;
%     if isempty(time_steps) || time_steps==0
%         time_steps = 1;
%     end
    elemLength = sbd.lzero(elem_nrs);
    opts = varargin{1};
    cType = varargin{2};
    cDim = varargin{3};
    if numel(varargin)>3
        hWaitDialog = varargin{4};
        progressAmount = varargin{5};
    else
        hWaitDialog = [];
        progressAmount = [];
    end
else
    it = getfrsbf(sbd,'it');
    % Determine input
    if nargin < 4
        opts = [];
        [cType, cDim] = getCsectData(sbd(1:end-4),it);
        cType = cType(elem_nrs);
        cDim = cDim(elem_nrs);
    elseif nargin == 4
        opts = varargin{1};
        [cType, cDim] = getCsectData(sbd(1:end-4),it);
        cType = cType(elem_nrs);
        cDim = cDim(elem_nrs);
    elseif nargin == 5
        opts = varargin{1};
        cType = {varargin{2}.CrossSection};
        cDim = {varargin{2}.Dimensions};
        if ~(length(cType)==length(elem_nrs) && length(cType)==length(elem_nrs))
           error('ERROR: wrong cross-sectional data supplied') 
        end
    elseif nargin == 6
        % These options should only be used by Spavisual
        opts = varargin{1};
        cType = varargin{2};
        cDim = varargin{3};
    elseif nargin > 6
        error('incorrect input')
    end
    hWaitDialog = [];
    progressAmount = [];
end

% Check cross-sectional data
if any(cellfun(@isempty,cType))
    error('Cross sectional data is not supplied for all elements where stress output is requested')
end

if all(strcmpi(cType,'rect') | strcmpi(cType,'recthol') | strcmpi(cType,'circ') | strcmpi(cType,'circhol'))
    elem_dimension = 3;
elseif all(strcmpi(cType,'line'))
    elem_dimension = 2;
else
    error('Possible mixing of ''line'' type cross-sections with 3D type cross-sections')
end

if ~isfield(opts, 'density'), opts.density = []; end
if ~isfield(opts, 'exterior'), opts.exterior = 0; end
if ~isfield(opts, 'type')
    if elem_dimension == 2
        [opts.type] = deal('fres');
    elseif elem_dimension == 3
        [opts.type] = deal('mises');
    end
else
    if ( elem_dimension==2 ) && any(strcmpi(opts(1).type, ...
            {'mises','sigx','tauxy','tauxz'}))
        warning('SPAVISUAL:InputWarning','%s\n%s',...
            ['Stress component ''' opts(1).type ''' not valid for cross-sectional type ''line'''],...
            'Stress component set to ''Fres''.');
    [opts.type] = deal('fres');
    end
end

if numel(opts)==1
    opts = repmat(opts,length(elem_nrs),1);
end

% Get necessary data from Spacar
if ~isa(sbd,'struct')
    if isempty(time_steps)
        time_steps = getfrsbf(sbd,'tdef');
    end
    elemLength = getfrsbf(sbd,'rl0');
    elemLength = elemLength(elem_nrs);
end

%% Determining the mesh
if elem_dimension == 2
    mesh(length(elem_nrs)) = struct('Line',[], 'Type',[],'Dimensions',[],...
    'Color',[],'DrawMode',[],'Edge',[]);
elseif elem_dimension == 3
    mesh(length(elem_nrs)) = struct('DrawMode',[], 'Line',[], ...
        'Type',[], 'Dimensions',[],'Edge',[],'Surf',[],'Interior',[]);
end
lineMesh(length(elem_nrs),1) = {[]};

for i = 1 : length(elem_nrs)
    density = meshDensity(it(elem_nrs(i)), cType{i}, cDim{i}, elemLength(i), opts(i).density);
    mesh(i) = meshCreate(cType{i}, cDim{i}, density, [0 0 0], opts(i).exterior);
    lineMesh{i} = mesh(i).Line;
end

%% Computing the stress resultants and stresses
% Compute the stress resultants
if ~isempty(hWaitDialog)
    [N, M, ind] = stressResultants(sbd, elem_nrs, time_steps, lineMesh, hWaitDialog, progressAmount);
else
    [N, M, ind] = stressResultants(sbd, elem_nrs, time_steps, lineMesh);
end
Nnames = fieldnames(N);
Mnames = fieldnames(M);
stress(length(elem_nrs)) = struct;

for i = 1 : length(elem_nrs)
    for j = 1 : length(Mnames)
        stress(i).Line.(['M' Mnames{j}]) = M(i).(Mnames{j});
    end
    for j = 1 : length(Nnames)
        stress(i).Line.(['F' Nnames{j}]) = N(i).(Nnames{j});
    end
    
    if elem_dimension == 3
        % Compute the stresses from the stress resultants based on specific
        % cross-sectional data
        dumStruct = stressCrossSect(mesh(i), N(i), M(i));
        Snames = fieldnames(dumStruct);
        for j = 1 : length(Snames)
            stress(i).(Snames{j}) = dumStruct.(Snames{j});
        end
    end
end

[stress, maxStress, minStress] = stressOutput(stress, opts(1).type);

stressExtrema.max = maxStress;
stressExtrema.min = minStress;

%% Determine input parameters
function varargout = getCsectData(filenm,it)

fid = fopen([filenm '.dat'], 'rt');
if fid<0, error(['Error: unable to open ' filenm '.dat']); end
try
    fstrm = char(fread(fid, '*uint8')');
catch errormsg % the try-catch combination makes sure the file is always closed
    fclose(fid);
    error(errormsg);
end
fclose(fid);

comment_tokens = '%#;';
keywords = {'beamvis','trussvis'};
for i=1:numel(keywords)
    search_string = sprintf(...
        '^\\s*%s[ \\t]+([^%s\\n\\r]+)|[\\n\\r]\\s*%s[ \\t]+([^%s\\n\\r]+)', ...
        keywords{i},comment_tokens,keywords{i},comment_tokens);
    tokens = regexpi(fstrm,search_string,'tokens');
    
    for j=1:numel(tokens)
        cType = strread(tokens{j}{1},'%s');
        if isnan(str2double(cType{1})) && any(strcmpi(cType{1}, {'rect' 'recthol' 'circ' 'circhol' 'line'}))
            cType = cType{1};
            par = str2num(tokens{j}{1}(length(cType)+1:end)); %#ok<ST2NM>
        elseif ~isnan(str2double(cType{1}))
            cType = 'rect';
            par = str2num(tokens{j}{1}); %#ok<ST2NM>
        else
            continue;
        end
        
        np = numel(par);
        if strcmpi(cType,'rect') && np >= 2
            dims = [abs(par(end-1)), abs(par(end))];    % correct for negative inputs
            if strcmpi(keywords{i},'trussvis')
                dims = sqrt((dims(1)*dims(2))/pi);            % convert to circular cross-section
                cType = 'circ';
            end
            par = par(1:end-2);
        elseif strcmpi(cType,'recthol') && np >= 4
            dims = [abs(par(end-3)), abs(par(end-2)), abs(par(end-1)), abs(par(end))];
            if strcmpi(keywords{i},'trussvis')
                dims = sqrt( (4*dims(3)*dims(4)*(dims(1)*dims(2)-1)) /pi); % convert to circular cross-section
                cType = 'circ';
            end
            par = par(1:end-4);
        elseif strcmpi(cType,'circ') && np >= 1
            dims = abs(par(end));
            par = par(1:end-1);
        elseif strcmpi(cType,'circhol') && np >= 2
            dims = [abs(par(end-1)), abs(par(end))];
            par = par(1:end-2);
        elseif strcmpi(cType,'line')
            dims = [];
        else
            continue;
        end
        
        if isempty(par)
            if strcmpi(keywords{i},'beamvis')
                sel = (it==1 | it==5 | it==15 | it==16);
            elseif strcmpi(keywords{i},'trussvis')
                sel = (it==2 | it==7);
            end
            CrossSection(sel,1) = {cType}; %#ok<AGROW>
            Dimensions(sel,1) = {dims}; %#ok<AGROW>
        else
            CrossSection(par,1) = {cType}; %#ok<AGROW>
            Dimensions(par,1) = {dims}; %#ok<AGROW>
        end
    end
end

try
    varargout{1} = CrossSection;
    varargout{2} = Dimensions;
catch %#ok<CTCH>
    CrossSection(1:numel(it)) = {'line'};
    Dimensions(1:numel(it)) = {[]};
    varargout{1} = CrossSection;
    varargout{2} = Dimensions;
end
