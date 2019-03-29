% -*- matlab -*- (enables emacs matlab mode)

DATA_TYPE = 'IPNL';
DATA_FILE ='../analysis_output/17-01-06_18-00_10MGamma_keV1099_MLEM_TIMING_200MBq.txt';
SM_DIR = '/Users/mattia.fontana/PhD/Work/Backup_CC/sps/SimulationOutput/Geant4/CC_NuclearMedicine/NM_ComparisonAnger/reconstruction_MLEM/Matrix/';        % directory to write the system matrix
SENSITIVITY_FILE = '';
%% directory with the results
RESULTS_DIR = '../results/AngerComparison_withTIMING_forPaper/2614keV/';


%%   Data management
SAMPLES =1;                            % for statistical studies
COUNTS_PER_SAMPLE=500000;                  % nb of events to treat
PRESELECT=0;                           % COUNTS_PER_SAMPLE is:
                                       % if 0, nb of detected events
                                       % if 1, nb of useful events
FIRST =1;                              % first event to treat (no preselect)
                                       % start at 1; only for _CALC_


%%   Iterations
FIRST_ITERATION=1;                     % continue iterations;
                                       % 1 means from the beginning
ITERATIONS = 15;
LAST_ITERATION=FIRST_ITERATION+ITERATIONS;


%%   System matrix management and calculation
TYPE_SM='EV_CALC_RAM';             % 'SM_PRECALC' (re-used matrix), 'SM_CALC' (write on disk), 'EV_PRECALC_DISK'
                                   % 'EV_PRECALC_RAM' , 'EV_CALC_DISK', 'EV_CALC_RAM' (write in the RAM memory)
CLEAR_SM=1;                        % Delete at the end =1 or not =0 (SM or events)
COUNTS_PER_PARTITION = 50;         % Only for TYPE=='SM_CALC'

ALGORITHM = 'CV';          % CV, RTS, RTV
MODEL = 'cos0rho0';       % cos0rho0 que le K et la gaussienne,
                          % cos1rho1 cos/rho, cos1rho2 cos/rho^2, cos0rho2 1/rho^2
WIDTH_FACTOR=1;         % thicker cones



%%   Energy a priori
ENERGY_FLAG = 2;               % ANY=0, RANGE=1, KNOWN=2
if (ENERGY_FLAG==1)            % total energy in some range
  ENERGY_MIN = 1147.5;         % lower bound, keV
  ENERGY_MAX = 1402.5;         % upper bound, keV
end
if (ENERGY_FLAG==2)            % known total energy
  ENERGY_TOTAL=1099;           % the total energy, keV
end


%% Hodoscope a priori
HODOSCOPE_FLAG=0;                           % ON=1, OFF =0;
if (HODOSCOPE_FLAG > 0)
  BEAM_SIGMA = 2 ;                  % Gaussian hodoscope, cm
  BEAM_WIDTH_FACTOR =3 ;            % number of sigma
  BEAM_ENTRY_POINT =[-10, 0, 0];         % a point on the beam line
  BEAM_DIRECTION =[1, 0, 0] ;            % direction of the beam
  shift=0;
  BEAM_FIRST_ITERATION=FIRST_ITERATION+shift;
  BEAM_ITERATIONS=ITERATIONS-shift;
  BEAM_INCLUSION='ONLYHODO';         % CONSTANT, LINEAR, FORCE, FORCEINV, ALTERNATE
end


%%   Reconstructed volume
VOLUME_DIMENSIONS = [ 5, 5, 5] ;   % in cm
%VOXELS = [20, 20, 2] ;
VOXELS = [51, 51, 51] ;              % voxels in volume
VOLUME_CENTRE = [0, 0, 0];             % volume to reconstruct, cm


%%   Camera properties
CAMERA_CENTRE = [0, 0, 0];             % camera position, cm
CAMERA_NORMAL = [0, 0, -1];            % vector orthogonal to the
                                       % camera, pointing towards
                                       % the source
%% Dresden camera
NB_LAYERS = 7;
LAY_CENTRE_1= [0,0,-10];   % supposed to be on a line || Oz, if >1
LAY_CENTRE_2= [0,0,-11];   % supposed to be on a line || Oz, if >1
LAY_CENTRE_3= [0,0,-12];   % supposed to be on a line || Oz, if >1
LAY_CENTRE_4= [0,0,-13];   % supposed to be on a line || Oz, if >1
LAY_CENTRE_5= [0,0,-14];   % supposed to be on a line || Oz, if >1
LAY_CENTRE_6= [0,0,-15];   % supposed to be on a line || Oz, if >1
LAY_CENTRE_7= [0,0,-16];   % supposed to be on a line || Oz, if >1
ABS_CENTRE = [0,0,-31];
LAY_SIZE = [9,9,0.2];
ABS_SIZE = [28,19,3.5];
LAY_VOXELS= [90,90,2];
ABS_VOXELS= [280,190,35];
% considering spatial uncertainty
SPATIAL_UNCERTAINTY=0;   % 1=yes, 0=no
if (SPATIAL_UNCERTAINTY==0)
  LAY_VOXEL_SAMPLING = [1,1,1];     % nb points to sample a voxel
  ABS_VOXEL_SAMPLING = [2,2,5];    % nb points to sample a voxel
end
% Camera frame (devrait-il être direct ?)
Ox = [1, 0, 0];     % parallel to scatterer edge
Oy = [0, 1, 0];     % parallel to scatterer edge
Oz = [0, 0, -1];    % orthogonal to the camera, tw the source
