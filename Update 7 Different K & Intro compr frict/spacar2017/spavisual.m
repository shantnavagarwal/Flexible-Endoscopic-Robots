function varargout = spavisual(filename,varargin)

% SPAVISUAL(FNAME) visualizes SPACAR results.
% FNAME:    SPACAR file name (without extension)
%
% SPAVISUAL(FNAME,SELMODE) visualizes the SELMODE-th vibration mode after a 
% mode 4, 7, 9 or 10 run or the SELMODE-th buckling mode after a mode 8 run. 
% SELMODE must be a scalar.
%
% SPAVISUAL(FNAME,ANIM_OPTS,VISU_OPTS) allows overriding the settings in the
% VISUALIZATION section in the SPACAR DAT-file. 
% ANIM_OPTS and VISU_OPTS are structures with one or more of fields, as 
% listed below. Each field is used to override the SPACAR default settings
% and/or the settings in the SPACAR DAT-file. Pass in an empty matrix for 
% ANIM_OPTS if you wish to override VISU_OPTS only.
%
% H = SPAVISUAL(...) returns the handle of the figure in which the results 
% are visualized.
%
% See also the SPAVISUAL documentation for more information.

%   rewrite of original SPAVISUAL code (by Jan Bennink/Ronald Aarts)
%   current version:    >0.4 (alpha)
%   date last modified: 2017-09-21
%   author(s):          Tjeerd van der Poel, Steven Boer

narginchk(1,4);
nargoutchk(0,1);

%% check command line input arguments (for correct data type)______________
% possible calling syntax:
% 1: spavisual(fname)
% 2: spavisual(fname,selmode)
% 3: spavisual(fname,modeopts,visopts,graphopts) (any of the last three
%    structures may be omitted or specified as empty)

if nargin>1
    counter = {'first','second','third','fourth'};
    % parse second through to last variable argument (in reversed order)
    for i=nargin-1:-1:2
        if ~isempty(varargin{i}) && ~isa(varargin{i},'struct')
            error('SPAVISUAL:InputError',...
                'The %s argument to SPAVISUAL must be a structure or empty.',...
                counter{i});
        end
    end
    % parse first variable argument
    if ~isempty(varargin{1}) && ~isa(varargin{1},'struct')
        if nargin>2
            error('SPAVISUAL:InputError',...
                'The %s argument to SPAVISUAL must be a structure or empty.',...
                counter{1});
        elseif isempty(varargin{1}) || ~(isnumeric(varargin{1}) && numel(varargin{1})==1)
            error('SPAVISUAL:InputError',...
                'When calling SPAVISUAL(FNAME,SELMODE), SELMODE must be a scalar.');
        end
    end
end

% check existence and date of DAT and SBD files
if ~(exist([filename '.dat'],'file')==2)
    error('SPAVISUAL:FileNotFound',['Unable to find the SPACAR DAT-file ' filename '.dat.']);
end
sbd   = sprintf('%s.sbd',filename);    % set SBD file
if ~(exist(sbd,'file')==2)
    error('SPAVISUAL:FileNotFound','Unable to find the SBD-file %s.', sbd);
end
if isOlderThanDatFile(sbd)
    error('SPAVISUAL:FileDateError', ...
        'The SBD-file (%s) is older than the SPACAR DAT-file (%s).\n%s',...
        sbd, [filename '.dat'], 'Please rerun the SPACAR analysis.');
end


%% Read SPAVISUAL settings_________________________________________________
% set defaults, read options from DAT file and parse command line arguments
[elemSettings,animSettings,axisSettings,figSettings] = getInfo(filename,varargin);

%% Load data from SBD file_________________________________________________
lnp   = getfrsbf(sbd,'lnp');    % location matrix for the nodes
ln    = getfrsbf(sbd,'ln');     % connection matrix for the nodes (to elements)
le    = getfrsbf(sbd,'le');
Lzero = getfrsbf(sbd,'rl0');    % initial length of elements
mode  = getfrsbf(sbd,'mode');   % mode of last SPACAR run
it    = getfrsbf(sbd,'it');     % list of element types for each element in mechanism
% number of fixed, calculable, input, dynamic and kinematic coordinates
nxp   = getfrsbf(sbd,'nxp');
range = animSettings.Range;

% the following information is required for stress computations
try
    stress_elems = find((it==1 | it==2 | it==5 | it==7 | it==15 | it==16));
    nldeform = getfrsbf(sbd,'kdef')';
    Ct = getfrsbf(sbd,'dr0');
    Ct = Ct(:,1);
catch %#ok<CTCH>
    [nldeform, Ct] = getMiscData(filename, stress_elems);
end
clear stress_elems

% NOTE: each column is a timestep, nodal coordinates along rows
x      = getfrsbf(sbd,'x',range)';
time   = getfrsbf(sbd,'t',range);
sig    = getfrsbf(sbd,'sig',range)';
if length(x(:,1)) == 1
    x = x';
    sig = sig';
end

%% precompute initial rotation matrices of elements _______________________
% with respect to global coordinate system

%   this is rather cumbersome, as the rotation information is
%   stored differently in rxyz for the various elements,
%   see GETROTMTX for details
Rip = cell(numel(it),1);
Riq = cell(numel(it),1);
rxyzp = getfrsbf(sbd,'rxyz');   % initial orientation of p-nodes of elements
rxyzq = getfrsbf(sbd,'rq');     % initial orientation of q-nodes of elements
for i=1:numel(it)
    Rip{i} = getRotMtx(rxyzp(i,:),it(i));
    if any([15 16]==it(i)) % PLTUBE or TUBE
        Riq{i} = getRotMtx(rxyzq(i,:),it(i));
    else
        Riq{i} = Rip{i};
    end
end
clear rxyzp rxyzq i

%% Determining the shear coefficients _____________________________________
try                             
    shearCoef = getfrsbf(sbd,'stif');
    if animSettings.ModelDimension == 2
        shearCoef = shearCoef(:,3);
        shearCoef(:,2) = 0;
        shearCoef = fliplr(shearCoef);
    elseif animSettings.ModelDimension == 3
        shearCoef = shearCoef(:,5:6);
    end
catch %#ok<CTCH>
    % Deze warning nog even commenten totdat de nieuwe SPACAR uit komt
%     warning('SPAVISUAL:InputWarning','%s\n%s',...
%         'Old SPACAR version. Beam elements with large shear effects will be',...
%         'drawn incorrectly.');
    shearCoef = zeros(length(it),2);
end

%% create SPACARDATA structure ____________________________________________
spacardata = struct('filename', filename, 'x', x, 'le', le, 'ln', ln, ...
    'lnp', lnp, 'mode', mode, 'it', it, 'nxp', nxp, 'rip', [], 'riq', [], ...
    'lzero', Lzero, 'shearcoef', shearCoef, 'time', time, 'modes',[], ...
    'sig', sig, 'nldeform', nldeform, 'ct', Ct);
spacardata.rip = Rip;
spacardata.riq = Riq;
clear Rip Riq shearCoef Lzero Ct sig x le ln lnp it nxp time nldeform

%% preprocessing for the various modes ____________________________________
switch mode
    case 0          % DOF Test
        % For the time being don't expect an SVD analysis to be stored in
        % the SBD file, so repeat it.
        % nep: [fixed, calculable, rheonomic, dynamic, kinematic] deformations
        nep = getfrsbf(sbd,'nep');
        % nxp: [fixed, calculable, rheonomic, dynamic and kinematic] coordinates
        nxp = getfrsbf(sbd,'nxp');
        % BigD: [D0D0 DcD0 DmD0; D0Dm DcDm DmDm; D0Dc DcDc DmDc]
        BigD = getfrsbf(sbd,'bigd',1);
        % partioning of rows of BigD: [fixed, rheonomic, dynamic, kinematic, 
        % calculable] deformations
        % partioning of columns of BigD: [fixed, kinematic, calculable, 
        % rheonomic, dynamic] coordinates
        
        % select only [DcD0;DcDm]
        Dcc = BigD( 1:(nep(1)+nep(3)+nep(4)) , nxp(1)+(1:nxp(2)) );
        % perform SVD on Dcc, note that singular vector U is in terms of 
        % fixed, rheonomic and dynamic deformations and singular vector V 
        % is in terms of calculable coordinates
        [ U, s, V ] = svd(Dcc);
        s = diag(s);
        % Number of near-zero singular values
        nsing = length(find(s<sqrt(eps)*s(1)));
        if length(U)>length(V)
            nover  = length(U)-length(V)+nsing;
            nunder = nsing;
        elseif length(U)<length(V)
            nover  = nsing;
            nunder = length(V)-length(U)+nsing;
        else
            nover  = nsing;
            nunder = nsing;
        end
        animSettings.NumOverconstraints = nover;
        animSettings.NumUnderconstraints = nunder;

        % extend underconstraint to all coordinates
        Vext = [zeros(nxp(1),nunder); ...
            V(:,end-nunder+1:end); zeros(sum(nxp(3:5)),nunder) ];
        
        %%% Marijn edit: scale eigenvectors for visualization
        sc = calculateEigvScale(spacardata.x(:,1),Vext,spacardata.lnp,spacardata.ln,spacardata.it);
        
        % store relevant SVD results for later use
        spacardata.underconstraint = Vext*diag(sc);
        spacardata.overconstraint = U(:,end-nover+1:end);
    case {1,2}    % Deformation, Inverse Dynamics (2x)
        % nothing to be done?
    case {3,4,7,9,10}    % State Space, Eigenfrequencies, State Space in equilibrium
        % compute mode shapes and frequencies
        sbm = sprintf('%s.sbm',filename);
        if ~exist(sbm,'file')==2
            error('SPAVISUAL:FileNotFound','Unable to find the SBM-file %s.', sbm);
        end
        if isOlderThanDatFile(sbm)
            error('SPAVISUAL:FileDateError', ...
                'The SBM-file (%s) is older than the SPACAR DAT-file (%s).\n%s',...
                sbm, [filename '.dat'], 'Please rerun the SPACAR analysis.');
        end
        
        if mode==3 % in mode 3, NDDOF is zero in SBD file, use NDOF instead
            nddof = getfrsbf(sbd,'ndof');
        else
            nddof   = getfrsbf(sbd,'nddof');
        end
        animSettings.Modes.Maximum = min(animSettings.Modes.Maximum,nddof);
        nmodes  = animSettings.Modes.Maximum;
        
        if nmodes>0 % (no point in computing modes otherwise)
            % nep: [fixed, calculable, rheonomic, dynamic, kinematic] deformations
            nep = getfrsbf(sbd,'nep');
            % nxp: [fixed, calculable, rheonomic, dynamic, kinematic] coordinates
            nxp = getfrsbf(sbd,'nxp');
            
            % check if first order geometric transfer functions are present
            try
                computeModes = true;
                dum = getfrsbf(sbd,'dx',1); %#ok<NASGU>
                clear dum
            catch ME
                if strcmpi(ME.message,'ID not found')
                    % first order geometric transfer functions are not present
                    errordlg(sprintf('%s %s (%s).\n%s.\n%s.', ...
                        'The first order geometric transfer functions', ...
                        'were not found in the SBD-file', sbd, ...
                        'Please set the SPACAR output level to an appropriate value', ...
                        'See also the OUTLEVEL keyword in the SPACAR manual'), ...
                        'SPAVISUAL: Missing geometric transfer functions', ...
                        'modal');
                    computeModes = false;
                else
                    rethrow(ME);
                end
            end
            
            if computeModes
                % compute modes for all time/load steps
                for i=1:length(range)
                    if mode~=3
                        % read mass and stiffness matrices
                        M = getfrsbf(sbm,'M0',range(i));
                        K = getfrsbf(sbm,'K0',range(i)) + getfrsbf(sbm,'N0',range(i)) + ...
                            getfrsbf(sbm,'G0',range(i));
                        C = getfrsbf(sbm,'C0',range(i)) + getfrsbf(sbm,'D0',range(i));
                    else % use LTV file
                        % NOTE: DOF ordering may be different!!!
                        % CHECK WITH RONALD AARTS/JAAP MEIJAARD
                        ltv = sprintf('%s.ltv',filename);
                        if ~exist(ltv,'file')==2
                            error('SPAVISUAL:FileNotFound','Unable to find the LTV-file %s.', ltv);
                        end
                        if isOlderThanDatFile(ltv)
                            error('SPAVISUAL:FileDateError', ...
                                'The LTV-file (%s) is older than the SPACAR DAT-file (%s).\n%s',...
                                ltv, [filename '.dat'], 'Please rerun the SPACAR analysis.');
                        end
                        M = getfrsbf(ltv,'M0',range(i));
                        K = getfrsbf(ltv,'K0B',range(i));
                        C = getfrsbf(ltv,'C0B',range(i));
                    end
                    
                    % Calculate eigenfrequencies and modes_____________________________
                    % use only dynamic DOF part of stiffness/mass/damp matrix
                    
                    if animSettings.Modes.Damping
                        % Include damping matrix in computation of
                        % eigenmodes
                        A = -[C(1:nddof,1:nddof) K(1:nddof,1:nddof); -eye(nddof) zeros(size(K(1:nddof,1:nddof)))];
                        B = blkdiag(M(1:nddof,1:nddof), eye(nddof) );
                        
                        if nddof > nmodes
                            eopts.disp = 0;  % do not display anything
                            [V,D] = eigs(A, B, 2*nmodes, 'sm', eopts);
                            %[V,D] = eigs(K(1:nddof,1:nddof),M(1:nddof,1:nddof),nmodes,'sm',eopts);
                        else
                            % use only dynamic DOF part of stiffness matrix
                            [V,D] = eig(A, B);
                            %[V,D] = eig(K(1:nddof,1:nddof),M(1:nddof,1:nddof));
                        end
                        D  = diag(D);  % eigenvalues
                        % sort eigenvalues and eigenvectors by eigenvalue magnitude
                        [d,o] = sort(D(:)); %#ok<ASGLU>
                        V = V(length(M)+1:end, o(1:2:end));
                        d = abs(D(o(1:2:end))*1i);
                        % compute the eigenfrequencies
                        spacardata.modes.frequencies{i} = d*[1 1/(2*pi)]; % [rad/s Hz]
                        
                        % Make modeshapes orthonormal with respect to M
                        V = V*diag(1./sqrt(diag(V.'*M*V)));
                        
                        if animSettings.Modes.ScaleEq   % Equalize norm of real and complex part
                            V = real(V)+imag(V)*diag(1i*sqrt(diag(real(V)'*real(V)) ./ diag(imag(V)'*imag(V))));
                        end
                    else
                        if nddof > nmodes
                            eopts.disp = 0;  % do not display anything
                            [V,D] = eigs(K(1:nddof,1:nddof),M(1:nddof,1:nddof),nmodes,'sm',eopts);
                        else
                            % use only dynamic DOF part of stiffness matrix
                            [V,D] = eig(K(1:nddof,1:nddof),M(1:nddof,1:nddof));
                        end
                        D  = diag(D);  % eigenvalues
                        % sort eigenvalues and eigenvectors by eigenvalue magnitude
                        [d,o] = sort(abs(D(:))); %#ok<ASGLU>
                        V = V(:,o);
                        d = D(o);
                        % compute the eigenfrequencies
                        spacardata.modes.frequencies{i} = sqrt(d)*[1 1/(2*pi)]; % [rad/s Hz]
                        
                        % normalize mode shapes to have norm 1 (and store them)
                        V = V*diag(1./sqrt(diag(V'*V)));
                    end
                    % first order geometric transfer function of coordinates with
                    % respect to DOFs, with partitioning of DOF vector:
                    %   q = [xp' xm' ep' em' xk' ek']':
                    %       xp: rheonomic coordinates
                    %       xm: dynamic coordinates
                    %       ep: rheonomic deformations
                    %       em: dynamic coordinates
                    %       xk: kinematic coordinates
                    %       ek: kinematic deformations
                    dx = getfrsbf(sbd,'dx',range(i));
                    
                    % select only dynamic DOFs (mode shapes are computed with rheonomic
                    % deformations/coordinates as static inputs, kinematic DOFs are
                    % also ignored)
                    if mode~=3
                        dx = dx(:,[ nxp(3)+(1:nxp(4)) nxp(3)+nxp(4)+nep(3)+(1:nep(4)) ]);
                    else % after a SPACAR mode 3 run, dynamic dofs are
                        % "renamed" as rheonomic
                       dx = dx(:,1:nxp(3)+nxp(4)+nep(3)+nep(4));
                    end
                    
                    %%% Marijn edit: scale eigenvectors for visualization
                    sc = calculateEigvScale(spacardata.x(:,i),dx*V,spacardata.lnp,spacardata.ln,spacardata.it);
                    
                    % mode shapes in terms of nodal coordinates
                    spacardata.modes.vibration{i}  = dx*V*diag(sc);
                end
                
                clear D d o V nmodes nddof K M sbm ltv dx
            end
            clear computeModes nep nxp i tdef
        end
    case 8          % Buckling
        sbm = sprintf('%s.sbm',filename);
        if ~exist(sbm,'file')==2
            error('SPAVISUAL:FileNotFound','Unable to find the SBM-file %s.', sbm);
        end
        if isOlderThanDatFile(sbm)
            error('SPAVISUAL:FileDateError', ...
                'The SBM-file (%s) is older than the SPACAR DAT-file (%s).\n%s',...
                sbm, [filename '.dat'], 'Please rerun the SPACAR analysis.');
        end

        % nep: [fixed, calculable, rheonomic, dynamic, kinematic] deformations
        nep = getfrsbf(sbd,'nep');  
        % nxp: [fixed, calculable, rheonomic, dynamic, kinematic] coordinates
        nxp = getfrsbf(sbd,'nxp');

        % compute modes for all load steps
        
        for i=1:length(range)
            k0  = getfrsbf(sbm,'k0',range(i));
            g0  = getfrsbf(sbm,'g0',range(i));
            dx  = getfrsbf(sbd,'dx',range(i));
    
            [V,D] = eig(-k0,g0);    % solve buckling eigenvalue problem
            D = diag(D);            % eigenvalues: load multipliers
            [f,o] = sortrows([abs(D(:)) V']);   % sort by load multiplier absolute value
            % sort mode shapes and load multipliers
            V = f(:,2:end)';      % sorted mode shapes
            f = D(o);             % sorted load multipliers
        
            animSettings.Modes.Maximum = min(length(f),...
                animSettings.Modes.Maximum); % store maximum number of modes
            nmodes = animSettings.Modes.Maximum;
            
            if nmodes>0 % (no point in storing non-existing information)
                % store load multipliers
                spacardata.modes.loadmultipliers{i} = f(1:nmodes);
                
                % select only dynamic degrees of freedom
                dx = dx(:,[ nxp(3)+(1:nxp(4)) nxp(3)+nxp(4)+nep(3)+(1:nep(4)) ]);
                
                %%% Marijn edit: scale eigenvectors for visualization
                sc = calculateEigvScale(spacardata.x(:,i),dx*V(:,1:nmodes),spacardata.lnp,spacardata.ln,spacardata.it);
                    
                % store deformation for each buckling mode
                spacardata.modes.buckling{i} = dx*V(:,1:nmodes)*diag(sc);
            end
            clear dx f V D sel k0 g0 nmodes
        end
        clear i tdef nxp nep sbm
    otherwise       % mode not implemented => exit
        error('SPAVISUAL:UnknownSpacarMode','Unsupported SPACAR mode %d', mode);
end

%% start the GUI
h_fig = spagui(spacardata,figSettings,axisSettings,animSettings,elemSettings);

% Return the figure handle if it has been requested
if nargout>0, varargout{1} = h_fig; end

end

function isOlder = isOlderThanDatFile(inputfile)

    [p,name,ext] = fileparts(inputfile);

    dirinfo = struct2cell(dir([p pathsep name '.*']));
    extensions = cellfun(@(x) regexpi(x,'[.](.*)','tokens'),dirinfo(1,:),'UniformOutput',false);
    extensions = cellfun(@(x) x{1}{1}, extensions, 'UniformOutput', false);
    datenums = [dirinfo{5,:}];
    idxdat = strcmpi('dat',extensions);
    idxfile = strcmpi(ext,extensions);

    if datenums(idxdat)>datenums(idxfile)
        isOlder = true;
    else
        isOlder = false;
    end
    
end