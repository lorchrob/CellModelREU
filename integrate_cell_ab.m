function integrate_cell_ab(varargin)

  %%  Various notes
  %  case of files/variables:  We adopt the convention that all file and
  %    variable names use lowercase.  Windows does not differentiate
  %    between upper and lower case (for instance Windows considers
  %    file.txt, FILE.txt, and FIle.txt to all be the same file) while
  %    Matlab and Linux does (file.txt, FILE.txt, and FIle.txt are
  %    considered to be three different files in linux and VAR VAr and var
  %    are considered to be three different variables in Matlab).  For the
  %    sake of the latter two situations, we try to always use lower case.
  %  underscores:  When using multiple words to name files/variables, I try
  %    to use underscores to help maintain readability of the code.
  %  file structure:
  %    cancer_cell_XXXX_XX_XX:  Where current version (with version
  %        corresponding to the date XXXX_XX_XX) of code is stored
  %      matlab_files:  All the m/matlab files for the code
  %      system_specific_files:  various simulations are possible for
  %          different systems
  %        cancer_cell_SYSTEM_NAME:  various systems include flow through 
  %            a tapered microchannel (cancer_cell_microchannel), a tube
  %            (cancer_cell_tube), or even just a resting cell in a square
  %            fluid region (cancer_cell_square)
  %          flexpde_files:  ALL files needed for solving the coupled
  %              system at a single time step including flexpde
  %              templates/pde files and any include files those templates
  %              use.  Also, a couple outputs from flexpde are stored here
  %          data:  Data from each time step are stored in mat-files in
  %              this directory
  %  varargin:  By feeding in "parameter pairs" we are allowed to change
  %    settings in the program.  For instance:
  %    integrate_cell_ab('system_name','cancer_cell_in_square') changes
  %      default system to the cancer cell in a square
  %    integrate_cell_ab('system_name','cancer_cell_in_square','n_nodes',
  %      40) changes the default system and the default number of external
  %      nodes to 40
  %    Using varargin, one should be able to use the program below without
  %    changing it in order to consider multiple situations (but we can't
  %    use it to consider all desired situations yet...still some work on
  %    the program to do).
  
  %%  Major changes (vs 2018 version)
  %  fade_in_globs eliminated...simulations suggested one could tease out a
  %    solution in particularly difficult situations if one used this
  %    variable carefully in flexpde "stages" but recent additions to
  %    flexpde (initial equations section, for instance) have rendered such
  %    potential not significant enough to keep around
  %  tangential vectors used...instead of using cosines and sines, the code
  %    has been rewritten to use just tangential components of the
  %    elements.  These components are labeled ndxis, ndyis, ndxos, ndyos
  %    for normalized dx and normalized dy for the inner and outer nodes.
  %  xm/ym eliminated (except for use in initialization)...to enable
  %    simpler indexing, we store all nodes in a single location/list.
  %    There are n_nodes external nodes which are always listed first in
  %    clockwise order and then n_inodes internal nodes which are then
  %    listed (in random order).

  %%  Various simulation specific definitions
  %  p-structure that holds "program-specific" variable values
  %  Type of system...typical examples of different systems listed below
%   p.system_name = 'cancer_cell_in_square';
%   p.system_name = 'cancer_cell_in_square_no_flow';
%   p.system_name = 'cancer_cell_in_circle';
%   p.system_name = 'cancer_cell_in_microchannel';
%   p.system_name = 'cancer_cell_in_tube';
%   p.system_name = 'cancer_cell_into_wall';
%   p.system_name = 'cancer_cell_into_corner';
%   p.system_name = 'cancer_cell_in_box';
%   p.system_name = 'cancer_cell_in_tube';
%   p.system_name = 'red_blood_cell_in_shear_flow';
%   p.system_name = 'red_blood_cell_in_tube_fluid_interior';
  p.system_name = 'cancer_cell_in_microchannel_rounded';
  
  %  Which flexpde to use (6 or 7)
  p.flexpde_version = 7;

  %  There are various degrees of "loudness" that one can use when running
  %  flexpde.
  %  'n -s ' use the Nongui version of flexpde to run flexpde silently in
  %  the background
  %  ' -s ' use the gui version in "silent mode"...runs in background
  %  without showing any window (but still makes pictures and stores those
  %  pictures afterwards)
  %  ' -r -x ' use the gui version with window, run it, and exit afterwards
  %  ' -r ' gui version run only, user must exit the program by hand
  %  See flexpde help for more options including -q and others.  Note the
  %  importance of spaces in the above options.
  p.flexpde_args = 'n -s ';
  
  %  With Jose's work, we started using different time integrators.  Choose
  %  your favorite from the time_integrator programs included in this file
  %  including ab_advance_cell and adapt_fe
  p.time_integrator = @ab_advance_cell;
  
  %  When true, don't try and continue a previous simulation, start from
  %  scratch instead and write over any old simulation information present
  p.start_from_scratch = true;
  
  %  Factor to refine entire fluid mesh using "ngrid"
  p.ngrid_factor = 1;
  
  %  Unique data suffix identifier/name for a particular simulation/set of
  %  runs.  For example, if one runs a simulation for the
  %  cancer_cell_in_microchannel system using 20 external nodes but then
  %  decides they want to try using 40 external nodes, one might consider
  %  using the program call integrate_cell_ab('n_nodes',40,'sim_suffix',
  %  '_40_nodes') to use 40 external nodes and to store the resulting
  %  matlab data in the folder 'data_40_nodes' instead of overwriting data
  %  from the previous 20 external node run already stored in the 'data'
  %  folder.
  p.sim_suffix = '';
  
  %  Figure numbers to use
  p.cell_fig_num = 1;
  p.area_mon_fig_num = 2;
  
  %  FlexPDE options
  p.directlimit = 15000;
  p.regrid = 'on';
  p.find_cond_num = false;

  %  Flexpde complains when the distance from a boundary (or cell) node to
  %  another boundary segment is < 0.0001*the length of the adjacent
  %  boundary segment
  %  (http://www.pdesolutions.com/forum/phpBB3/viewtopic.php?f=2&t=242)
  %  This stores that value of 0.0001 for later use to try to make sure
  %  such situations do not occur, seen later
  p.dist_to_bound_frac = 0.0001;
  
  %  For debugging (set to true)
  p.randomize = false;
  
  %  Typical loop used to replace default values with user fed in values
  %  via the varargin capabilities of Matlab
  p_var_names = fieldnames(p);
  for vac = 1:2:numel(varargin)
    if any(strcmp(varargin{vac},p_var_names))
      p.(varargin{vac}) = varargin{vac+1};
    end
  end
  %  Time integrator seems to need a "locally defined" function handle
  if ischar(p.time_integrator)
    eval(['p.time_integrator = @',p.time_integrator]);
  end
    
  %  FlexPDE7 stores stuff differently vs FlexPDE6 so we need extra
  %  directory specific machinations to read and write flexpde files
  %  smoothly
  if p.flexpde_version == 7
    p.add_dir = ['temp_output',my_slash];
  else
    p.add_dir = '';
  end
  %  Include location of velocity info
%   switch p.system_name
%     case {'cancer_cell_in_microchannel','cancer_cell_in_microchannel_slower'}
%       p.vel_filename = 'uinfo_2.dat';
%     otherwise
      p.vel_filename = 'uinfo.dat';
%   end
  
  %%  Some initializations  
  %  System parameters that we are currently thinking about changing.
  %  There are other system parameters defined in the pde-file but those
  %  are not changed here.  Any needed parameter not defined here must be
  %  defined in the pde file and vice versa.  See program for specifics.
  s = get_system_parameters(p,varargin{:});
  if s.dt_min == s.dt_max, p.time_integrator = @adapt_rk2; end
  if s.enforce_symmetry, warning('Symmetry being enforced!'); end
  
  %  Initial cell shape/information
  [cnm2,cnm1,s,ts,areas,volests] = get_initial_cell_info(s,p,varargin{:});
  
  %  For debugging, stretch out the cell and hope it comes back to its
  %  original (almost) circular shape
%   c0.xn = 2*c0.xn;
%   %  Bump only the centermost node
%   cnm2.xn(41) = cnm2.xn(41)+0.5;
%   cnm1.xn(41) = cnm1.xn(41)+0.5;
  
  %%  Write various include files that we only need once per simulation
  %  Flexpde allows us to write external text files that we can include in
  %  any flexpde script.  We write some of these include files below.
  %  These particular include files need be written only once per
  %  simulation.
  %  The various system parameters mentioned above
  write_system_parameters(s);
  %  Write the global variables include file, a list of the names assigned
  %  to the velocity components for each node
  write_global_variables(cnm2,s);
  %  Write the boundaries for the cell
  write_cell_boundaries(cnm2,s,p);  
  %  Write pressure constraint
  write_pressure_constraint(cnm2,s);
  %  Write "friction" constraint...see cancer_cell_in_square_no_flow
  write_friction_constraint(cnm2,s);
  %  Write the report for the cell velocities
  write_summary_report(cnm2,s,p);
  
  %%  For plots
  %  Close figure that we want to use to plot the cell trajectory
  close(figure(p.cell_fig_num));
  %  Open a fresh figure window for it
  figure(p.cell_fig_num);
  %  Reshape the figure window
  set(p.cell_fig_num,'Units','inches','position',[3,6,10,4]);
  %  Make an axis (for later plotting)
  p.cell_ax_hand = axes;
  %  Do the same for the area monitor figure window
  close(figure(p.area_mon_fig_num)); figure(p.area_mon_fig_num);
  set(p.area_mon_fig_num,'Units','inches','position',[3,1,10,4]);
  p.area_mon_ax_hand = axes;
  
  %%  Initial forward euler integration to get AB2 going
  %  We first comment on the AB2 algorithm we are using.  It is a 2nd order
  %  Adams-Bashforth method with adaptive time stepping.  It uses the two
  %  previous time steps to estimate the location of the cell at the
  %  current time step.  It also chooses a time step that makes sure all of
  %  several criterion hold.  See ab_advance_cell for details.
  %
  %  We note, however, that because of the way that choosing the time step
  %  works, the time step size associated with a current cell shape
  %  estimate corresponds to the difference in time between the previous
  %  cell shape estimate and the current cell shape estimate, not the
  %  difference in time between the current cell shape and next cell shape
  %  estimate
  %
  %  Below, c_{n-2} = cnm2 and c_{n-1} = cnm1 are used in the eventual AB2
  %  procedure.
  
  %  Add a semi-random perturbation using cos and sin functions (for
  %  debugging)
  if p.randomize
    coef1 = 0.21; coef2 = 0.3; ths = 2*pi*(1:cnm2.n_enodes)'./...
      cnm2.n_enodes;
    cnm2.xn(1:cnm2.n_enodes) = cnm2.xn(1:cnm2.n_enodes)+...
      coef1*cos(2*ths)+coef2*cos(3*ths);
    cnm2.yn(1:cnm2.n_enodes) = cnm2.yn(1:cnm2.n_enodes)+...
      coef1*sin(2*ths)+coef2*sin(3*ths);
    cnm2.xn(end) = cnm2.xn(end)+0.13;
    cnm2.yn(end) = cnm2.yn(end)-0.14;
    cnm2 = node_info(cnm2,s);
    cnm1 = cnm2;
  end
  
  %  Get initial velocities for the cell (if needed)
  if isequal(cnm2,cnm1)
    cnm2 = get_vels(cnm2,s,p,true);
    %  Initialization of monitors to check if time stepping seems stable
    areas(end+1) = cnm2.area;
    volests(end+1) = cnm2.volest;
    cnm2.t = 0; ts = cnm2.t;
    %  Plot this initial cell
    cnm2 = plot_cell(cnm2,s,'r',p,ts,areas,volests);
    %  Store initial cell structure
    store_solution(cnm2,s);

    %  Take an adaptive rk2 time step forward to initialize multistep
    %  method (or just take a regular step if already using adapt_rk2).
    cnm1 = adapt_rk2(cnm2,cnm2,s,p);
    cnm1.t = cnm2.t+cnm1.dt;
    ts(2) = cnm1.t;
    areas(2) = cnm1.area;
    volests(2) = cnm1.volest;
  end
  
  %  Plot the cell again
  cnm1 = plot_cell(cnm1,s,'r',p,ts,areas,volests);
  %  Store cell positions and velocities at next time step
  store_solution(cnm1,s);
  c = cnm1;
  
  %  For the microfluidic channel, goal is to get cell to this point
%   while max(c.xn(:)) < (max(s.xb)-0.1)
  while (min(c.xn(:)) < 170) | ...
      isempty(strfind(p.system_name,'microchannel'))
    
    %  Find new position and velocities of cell...the new velocities in
    %  such a way that conditions are met to make sure nothing non-physical
    %  happens and that the cell velocities aren't changing too rapidly
    %  (see ab_advance_cell for details)
    fprintf('time step %d: \n',cnm1.time_step+1);
    [c,s] = p.time_integrator(cnm2,cnm1,s,p);
    %  Store monitors and time step stuff
    areas(end+1) = c.area;
    ts(end+1) = ts(end)+c.dt;
    c.t = ts(end);
    volests(end+1) = c.volest;
    %  Plot the cell again
    c = plot_cell(c,s,'r',p,ts,areas,volests);
    %  Store cell positions and velocities at next time step
    store_solution(c,s);
    
    %  Prepare for next AB2 step
    cnm2 = cnm1;
    cnm1 = c;
    
    if ts(end) >= s.t_max
      return
    end
    
  end

end

function s = get_system_parameters(p,varargin)

  %  Various system-specific (not cell-specific) information
  
  %  Boundary of system
  %  Used to rotate domain of system and, later on, the cell initial
  %  position
  s.domain_rot_ang = 0;
  %  Used for tube widths and lengths...and maybe other things
  s.char_width = 20;
  s.char_length = 20;
  %  For gravitational force (default is 0)
  s.grav = 0;
  %  Used to decide whether or not to explicitly force constant
  %  cross-sectional area of the cell.  Multiple options,
  %  'Aref','Vref','constA','constV','constmix','none'
  s.area_constraint = 'Aref';
  %  Only used when calculating volume estimates
  %    xaxi-assume approximate axisymmetry about x axis
  %    yaxi-assume approximate asisymmetry about y axis
  %    trap-use trapezoidal cross-sectional assumption
  %    none-just returns zero...avoids costly calculations
  s.vol_est_type = 'none';
  %  Only used for 'constmix' above...area_frac*dAdt+(1-area_frac)*dVdt = 0
  s.area_frac = 0.5;
  %  Whether or not to enforce symmetry wrt y-axis
  s.enforce_symmetry = false;
  %  For square simulations, how hard to try to make sure that the cell
  %  center stays near (0,0)
  if strcmp(p.system_name,'cancer_cell_in_square')
    s.rest_force_mag = 100;
  else
    s.rest_force_mag = 0;
  end
  %  Makes it so the cell stays inside an "imaginary" microchannel boundary
  %  for "cancer_cell_in_square" simulations.  Hence the default value is 0
  %  cause we are not usually interested in this case.  When we are,
  %  however, it mimics what a cell would look like (according to our
  %  model) when it gets stuck in the microchannel array
  s.imag_bnd_force_mag = 0;
  %  Even if we select the velocities just right so that darea/dt = dAdt =
  %  exactly 0, the time stepping is imperfect and while the area at the
  %  next time step after time stepping will be close to the area before
  %  the time step, it will differ by O(dt) (or maybe O(dt^2)).  This
  %  restoration factor keeps the area from wandering too far away from the
  %  initial area value.  In addition, if we really want, we could start
  %  out with an initial area that is too large and this restoration factor
  %  would cause the area to drift back towards a realistic (Aref) value.
  s.area_constraint_restoration_factor = 1e-3;
  %  Allow user to change any above values
  fn = fieldnames(s);
  for vac = 1:2:numel(varargin)
    if any(strcmp(varargin{vac},fn))
      s.(varargin{vac})= varargin{vac+1};
    end
  end
  %  Assign pright if need be
  if strncmpi(s.area_constraint,'const',5)
    warning('Constant area not yet set up for multiple cells');
%     s.pright = 0; 
  end
  switch p.system_name
    case 'cancer_cell_wgravwall_fict'
      s.xb = [-30,30,30,-30]; s.yb = [-30,-30,30,30];
      %  Note the "wall" actually forms a polygon, useful in some
      %  instances...for "gravwall", it is theoretically needed but the
      %  algorithm assumes polygon, so we make an oversized polygon
      s.xw = [-40,40,40,-40]; s.yw = [-10,-10,40,40];
    case 'cancer_cell_wgravwall'
      s.xb = [-30,30,30,-30]; s.yb = [-7,-7,30,30];
      s.xw = [-40,40,40,-40]; s.yw = [-40,-40,40,40];
    case {'cancer_cell_in_microchannel','cancer_cell_in_microchannel_slower'}
      %  Tapered microchannel from TruongVo et al
      s.xb = [-190,-150, 150,150,190,190,150,150,-150,-190];
      s.yb = [-15,  -15,-2.5,-15,-15, 15, 15,2.5,  15,  15];
%       %  NOT USED YET.  We found another temporary (but not as potentially
%       %  powerful) work-around in flexpde.
%       %  Strings used for boundary condition declaration in flexpde (ith
%       %  condition corresponds to line between (xb(i),yb(i)) to
%       %  (xb(i+1),yb(i+1))       
%       %  Slip conditions (on bottom of domain before the channel starts)
%       s.bc_strs{1} = ' nobc(u) value(v) = 0 nobc(p) ';
%       %  No-slip conditions (on bottom and top walls)
%       [s.bc_strs{2:3}] = deal(' value(u) = 0 value(v) = 0 nobc(p) ');
%       s.bc_strs{4} = s.bc_strs{1};
%       %  Outlet conditions
%       s.bc_strs{5} = ' value(p) = pright ';
%       s.bc_strs{6} = s.bc_strs{1};
%       [s.bc_strs{7:8}] = deal(s.bc_strs{2});
%       s.bc_strs{9} = s.bc_strs{1};
%       %  Inlet conditions
%       s.bc_strs{10} = ' value(p) = pright+deltap ';
    case 'cancer_cell_in_microchannel_rounded'
      tvec2 = [300,12.5]; that2 = tvec2./norm(tvec2);
      nhat2 = [that2(2),-that2(1)];
      tmpr = 3;
      nhat4 = [-1,0];
      th34 = atan2(abs(nhat2(2)*nhat4(1)-nhat2(1)*nhat4(2)),...
        sum(nhat2.*nhat4));
      th910 = th34;
      s.xb = [-190,-150,150-tmpr*tan(th34/2)*that2(1),150,150,190,...
        190,150,150,150-tmpr*tan(th910/2)*that2(1),-150,-190];
      s.yb = [-15,-15,-2.5-tmpr*tan(th34/2)*that2(2),-2.5-tmpr*tan(th34/2),-15,-15,...
        15,15,2.5+tmpr*tan(th34/2),2.5+tmpr*tan(th34/2)*that2(2),15,15];
      %  For use with arcs, when present.  For now we make two arcs at the
      %  top and bottom of the outlet of size 1 micron.  There is one entry
      %  for each boundary element (not for each node)
      s.ra = [0,0,tmpr,0,0,0,0,0,tmpr,0,0,0];
      s.th = [0,0,th34,0,0,0,0,0,th34,0,0,0];
      %  We also put the arc centers.  Note that for straight segments,
      %  these values (which we set to zero) are never used.
      s.arcx = [0,0,150-tmpr,0,0,0,0,0,150-tmpr,0,0,0];
      s.arcy = [0,0,-2.5-tmpr*tan(th34/2),0,0,0,0,0,2.5+tmpr*tan(th34/2),0,0,0];
%       %  For debugging
%       th = linspace(0,2*pi);
%       plot(s.xb,s.yb,'b',s.arcx(3)+tmpr*cos(th),s.arcy(3)+tmpr*sin(th),...
%         s.arcx(9)+tmpr*cos(th),s.arcy(9)+tmpr*sin(th));
%       keyboard
    case 'cancer_cell_in_circle'
      %  Circle for pressure tests...just radius for now
      s.xb = [20];
      s.yb = [20];
    case {'cancer_cell_in_square','cancer_cell_in_square_no_flow'}
      %  Square for pressure tests
      s.xb = [-15,15,15,-15];
      s.yb = [-15,-15,15,15];
      %  This is the x-location at which we believe the cell gets
      %  stuck...in the simulations the cell's center of mass should
      %  start and eventually end up at/very close to this location.  Note
      %  xloc = 0 corresponds to the center of the microchannel (see xw and
      %  yw above)
      s.xloc = 140;
      %  When using "wall forces" to make the shape conform to microchannel
      %  walls, these "wall boundaries" are used...note they never appear
      %  explicitly in flexpde boundary declarations but only implicitly
      %  appear in "wall forces" (see fxnwf below)
      s.xw = [-190,-150, 150,150,190,190,150,150,-150,-190]-s.xloc;
      s.yw = [-15,  -15,-2.5,-15,-15, 15, 15,2.5,  15,  15];
    case {'cancer_cell_chemotaxing_no_flow','cancer_cell_chemotaxing'}
      %  "Huge" square to allow for migration
      s.xb = [-30,30,30,-30];
      s.yb = [-30,-30,30,30];
    case {'cancer_cell_in_tube','cancer_cell_in_tube_fluid_interior',...
        'cancer_cell_in_box','red_blood_cell_in_shear_flow',...
        'red_blood_cells_in_shear_flow','red_blood_cell_in_tube',...
        'red_blood_cell_in_tube_fluid_interior'}
      %  Cancer cell flowing in rectilinear tube or just floating in a box
      %  or red blood cell in shear flow
      s.xb = 0.5*s.char_length*[-1,1,1,-1];
      s.yb = 0.5*s.char_width*[-1,-1,1,1];%0.5*(30+5)/2*[-1,-1,1,1];
      s.ra = [0,0,0,0];
    case 'cancer_cell_into_wall'
      %  Cancer cell impinged on a wall
      %  Upside down T
      s.xb = [-30,30,30,10,10,-10,-10,-30];
      s.yb = [0,0,10,10,30,30,10,10];
    case 'cancer_cell_into_corner'
      %  Cancer cell impinged on a corner
      s.xb = [-12,0,12,20,6,6,-6,-6,-20];
      s.yb = [-16,0,-16,-10,-16+8+50/3,30,30,-16+8+50/3,-10];
  end
  %  Rotate domain...boundary conditions may or may not need to be adjusted
  [s.xb,s.yb] = deal(...
    cos(s.domain_rot_ang)*s.xb-sin(s.domain_rot_ang)*s.yb,...
    sin(s.domain_rot_ang)*s.xb+cos(s.domain_rot_ang)*s.yb);
  
  %  Lubrication theory information
  %  When nodes are within this distance of other objects (nodes or walls)
  %  we will begin to add in lubrication terms to keep them from hitting
  %  the objects while still allowing them to move with a reasonable 
  %  trajectory
  s.lub_threshold_frac = 0.01;
  %  When we are not within the lubrication threshold above, we choose a
  %  small enough time step so that dist does not change by more than this
  %  amount
  s.dist_max_rel_change = 0.1;
  %  When we are within the lubrication threshold above (and hence some
  %  lubrication terms are included in our force estimates), we choose a
  %  small enough time step so that the dist does not change by more than
  %  this amount
  s.lub_dist_max_rel_change = 0.01;
  %  Limits how quickly area of the cell can change over a time step
  s.area_max_rel_change = 0.01;
  %  This is the maximum size of the lubrication element that is
  %  superimposed over the top of the node in question.  The value actually
  %  used in lubrication calculations can actually be smaller (see
  %  lubrication special cases pg 35-36...in my dissertation).
  s.lub_max_little_l_frac = 10;
  %  Fraction associated with the distance at which we start using a
  %  variable incompressibility penalty value (k)
  s.var_pen_dist_frac = 50;
  %  Strength of elastic-like interaction between solids in our system
  %  For lubrication, an additional deltass*f(dist)*normal gets added
  %  (later)
  s.cr = 1e6;
  %  Strength of lubrication terms...default should be 1 as this
  %  corresponds to the exact estimates from lubrication theory
  s.lub_coef = 1;
  %  If "linear" then we use hooke's law.  If 1/r, we use 1/r "elasticity".
  %  We also introduce "quadratic" which has slope zero (i.e. "gently
  %  introduced") at lub_threshold_frac
  s.ss_type = '1/r_zero_slope';
  %  Strength of lubrication interaction (should be between 0-off and 1-on
  %  from the way we've done things)
  s.lub_elements = 1;
  
  %  Imposed force coefficients
  
  
  %  Maximum and minimum allowable time step and initial timestep
  s.dt_max = 10;
  s.dt_min = s.dt_max*1e-6;
  if isequal(func2str(p.time_integrator),'adapt_rk2')
    s.dt_0 = 0.1;
  else
    s.dt_0 = s.dt_min;
  end
  %  We allow the time step to grow by at most this amoount each step
  s.dt_grow = 2;
  %  Maximum time...1e6 effectively unlimited
  s.t_max = 1e6;
  
  %  We try to choose a time step small enough so that our nodal velocity
  %  estimates change by no more than this amount.  Because flexpde is not
  %  perfect nor is its mesh and corresponding solution process, we may not
  %  be able to choose a time step small enough to guarantee this...but
  %  hopefully (see time stepping later on).
  s.v_change_max = 0.01;
  s.dt_safety = 0.8;
  %  Characteristic velocity.  All we need is a decent guess to get things
  %  started.
  %  For red blood cells...this becomes the flow speed on the top/bottom of
  %  the shear flow regon giving a shear rate of 2*v_char/char_width
  s.dv_tol = 0.001;
  s.v_char = (7.4+1.4)/2;
  s.v_char = 30*3.5/((30+5)/2);%(11+3.5)/2;
  %  For flow in a tube...average velocity of 1 micron/ms across width
  if ~isempty(strfind(p.system_name,'in_tube'))
    s.v_char = 1;
  else ~isempty(strfind(p.system_name,'chemotaxis'))
    s.v_char = 0.1;
  end
  
  %  Use of staged lubrication effects.  In flexpde one can slowly turn on
  %  terms in equations by "staging" them.  If we use "s.lub_on" as defined
  %  below, eventually flexpde will calculate results using 0^2, use those
  %  results as initial guesses in order to calculate results using
  %  (1/9)^2, use those results as initial guesses in order to calculate
  %  results using (2/9)^2 all the way to s.lub_on = 1.  lub_on multiplies
  %  all additional lubrication terms that are added when the distance
  %  between the cell and other objects is < lub_threshold.  Hence if
  %  lub_on = 0 then all such lubrication terms are zeroed out and the
  %  solution is found as if there were no such lubrication terms.  Staging
  %  the "lub_on" variable then effectively slowly introduces these
  %  additional lubrication terms.
  s.lub_on = linspace(0,1,10).^2;
  %  Turns off staging.  Again, if we let lub_on = 0, no additional
  %  lubrication terms are ever added.  In the past, doing this led to
  %  physically unrealistic situations where cell nodes eventually
  %  penetrate nearby walls.
  s.lub_on = 1;
  
  %  How long we let flexpde run for before we think about intervening (in
  %  case flexpde did something wrong)...30 minutes
  s.max_run_time_coupled = 30*60;
  %  How much longer we let flexpde run the first time (we don't have any
  %  good guesses yet for the nodal velocities)
  s.max_run_time_initial_add = 60*60;
  
  %  Grid limit to use for fpde runs (see gridlimit in select section for
  %  fpde)
  s.grid_limit = 15;

  %  Allow user to change any above values
  fn = fieldnames(s);
  for vac = 1:2:numel(varargin)
    if any(strcmp(varargin{vac},fn))
      s.(varargin{vac})= varargin{vac+1};
    end
  end

end

function [cnm2,cnm1,s,ts,areas,volests] = get_initial_cell_info(s,p,varargin)

  %  Various cell-specific information  
  %  To be consistent, we try to use dissertation names of objects
  %  (includes violation of "only use lowercase" convention)
  
  %%  Constants that are the same throughout any given simulation:
  %  Initial cell shape "type" (like a circle or a star or something
  %  different, see second section below)
  c.initial_cell_type = 'circle';
  
  %%  Allow multiple internal nodes; 2/1/2019
  %  In this version of the code, xm and ym are no longer active parts of
  %  the corresponding dynamical system, they are only used for
  %  initialization of the cell shape.
  c.mult_int_nodes = false;
  
  %  Different initial shapes:  Node positions
  %  Couple of options for central node positions (note: size(xm) should 
  %    == size(ym))
  switch p.system_name
    case {'cancer_cell_in_microchannel',...
        'cancer_cell_in_microchannel_slower',...
        'cancer_cell_in_microchannel_rounded'}
      c.xm = -170;
    otherwise
      c.xm = 0;
  end
  c.ym = 0;%-2.489;
  %  While we had multiple cells working at one point, that is not
  %  currently the case
%   c.xm = -140+[-4,4];
%   c.ym = [0,0];
  
  %  Number of outer nodes and # of cells (> 1 cell not yet possible)
  c.n_enodes = 20;
  c.n_cells = numel(c.xm);
  
  %  Length scales
  %  ls-a length scale.  We assume we start off with a circular cell with
  %    radius of size ls.  We currently adjust this so that ls produces
  %    liref below.  ls is only used to get liref and not used elsewhere.
  %  Normal cell
  c.ls = 11.2/2;
  %  Cancer cell
  c.ls = 12.4/2;
  
  %  Initialize times and corresponding areas
  ts = [];
  areas = [];
  volests = [];
  
  %  Redefine any of above using varargin
  fn = fieldnames(c);
  for vac = 1:2:numel(varargin)
    if any(strcmp(varargin{vac},fn))
      c.(varargin{vac}) = varargin{vac+1};
    end
  end
  
  %  More length scales
  %  We make it so that the reference length of the internal elements
  %  produce a polygon with the same area as a circle with radius ls.  Note
  %  that as n_nodes -> infinity, ls -> lref.
  c.liref = sqrt(2*pi*c.ls^2/c.n_enodes/sin(2*pi/c.n_enodes));
  %  Reference length of exterior segment.  We need some estimate of
  %  sphericity in order to estimate this accurately.  "Sphericity" =
  %  preferred perimeter of shape/perimeter of circle with preferred area
  %  of shape where here perimeter of actual shape = loref*n_nodes and
  %  preferred area = Aref (see below)
  %  "sphericity" taken from red blood cell model
  c.loref = 0.97*2*sin(2*pi/c.n_enodes/2)*c.liref/0.832231354014029;
  %  "sphericity" of 1
  c.loref = 2*sin(2*pi/c.n_enodes/2)*c.liref;
  %  If the elasticity is such that it keeps outer elements from extending
  %  beyond a certain point, the below max length becomes important
  %  ELASTICITY IN OUTER ELEMENTS MUST MATCH CONDITION BELOW!
  c.lomax = 2*c.loref;
  %  Effectively unbounded
  c.lomax = 1e6*c.loref;
  
  %  Reference cross-sectional area
  c.Aref = pi*c.ls^2;
  %  Pressure associated with cross-sectional area variations
  c.kp = 50;%*60;
  
  %  Viscosities of membrane elements with outer and inner elasticities
  c.mume = 200;
  %  Note the inner viscoelasticity is adjusted by the number of nodes to
  %  keep the average force per unit of "cell cytoplasm" constant
  c.mumi = 100*20/c.n_enodes;
  %  For if we make the interior a fluid.  Based off of "cell interior has
  %  a viscosity approximately 5X that of plasma"...may not be completely
  %  accurate.
  c.mufi = 5;
  
  %  Bending elasticity
  c.kbe = 0.9;
  c.kbi = 0;
  %  Bending moment (m) approximation "type"...newm-use curvature;
  %  oldmlin-use angle; oldmnonlin-use 2*tan(angle/2)
  c.mtype = 'oldmlin';
  %  Linear elasticity of external elements
  c.kte = 12;
  c.ktlin = true;
  %  Linear elasticity of internal elements (added 2018/07/09 for cancer
  %  cells
  if strncmpi(p.system_name,'cancer_cell',11)
    c.kti = 0;
  else
    c.kti = 0;
  end
  %  Area elasticity...based on "Three-dimensional numerical simulation of
  %  red blood cell motion in Poiseuille flows" Shi, Pan, Glowinski.  Their
  %  model is 3d, ours is 2d, so we can't really use their model to
  %  estimate this.  This is something we will need to find a way to
  %  estimate on our own (if we plan to go ahead with the 2d model).
  c.ka = 0;
  
  %  Imposed force on nodes.  If true, we interpolate the force using the
  %  imposed forces defined elsewhere (see FORCES ON NODES below and
  %  following)
  c.interp_imposed_force = false;
  
  %%  Calculate/find important directory names
  s.data_loc = ['..',my_slash,'system_specific_files',my_slash,...
    p.system_name,my_slash,'data',p.sim_suffix,my_slash];
  temp_flexpde_files_loc = ['..',my_slash,'system_specific_files',my_slash,...
    p.system_name,my_slash,'flexpde_files',my_slash];
  s.flexpde_files_loc = strrep(temp_flexpde_files_loc,'flexpde_files',...
    ['flexpde_files',p.sim_suffix]);
  if ~isequal(temp_flexpde_files_loc,s.flexpde_files_loc)
    if exist(s.flexpde_files_loc,'file')
      warning(sprintf(['Replacing \n',...
        strrep(s.flexpde_files_loc,my_slash,[my_slash,my_slash]),...
        '\nwith files from (if any exist)\n',...
        strrep(temp_flexpde_files_loc,my_slash,[my_slash,my_slash])]));
    end
    [~,~,~] = copyfile(temp_flexpde_files_loc,s.flexpde_files_loc);
  end
    
  %%  Initial cell shape setup
  %  Read in past efforts to see if we can restart an already existing
  %  simulation
  [~,~,~] = mkdir(s.data_loc);
  %  Get list of all files in the directory.  We assume everything is nice
  %  and ordered starting with time_step_0000, then time_step_0001, and
  %  consecutively until the most recent flexpde run, time_step_XXXX.
  data_file_list = dir([s.data_loc,'time_step_*.mat']);
  %  If we are told to start_from_scratch, we erase all previous results
  if p.start_from_scratch, delete([s.data_loc,'time_step_*.mat']); end
  %  Note if start_from_scratch is true, we advance the time step all the
  %  way to XXXX.  If not, we just set time_step = -1 which the program
  %  later uses to appropriately restart the given simulation
  c.time_step = -1+(~p.start_from_scratch)*numel(data_file_list);
  
  %  Pressure at the right hand side of the system...designed for
  %  microchannel flow but may come in handy for other systems as well
  c.pright_prev = 0;
  
  %  If we start from scratch or if we haven't taken any time steps yet,
  %  use one of our predefined initial shapes
  if c.time_step < 0
    %  Reset the time step value for use later
    c.time_step = 0;
    c.dt = s.dt_0;
    
    %  Angles often used for initial shapes
    ths = linspace(0,2*pi,c.n_enodes+1)';
    ths = ths(1:end-1);
    
    %  Outer node positions (based on central node location(s)).  Note xn 
    %  is an n_nodes x n_cells matrix
    switch c.initial_cell_type
      case 'circle'
        %  Circular.  Note xn is a n_nodes x n_cells matrix
        c.xn = bsxfun(@plus,c.xm,c.liref*cos(ths));
        c.yn = bsxfun(@plus,c.ym,c.liref*sin(ths));
      case 'crown'
        if c.n_enodes ~= 5, error('Use 5 nodes for crown!'); end
        a = fzero(@(a) polyarea(a*[-1,1,1,0,-1],[0,0,a,2*c.ym,a])-...
          c.Aref,[0,100]);
        c.xn = a*[-1;1;1;c.xm;-1];
        c.yn = [0;0;a;2*c.ym;a];
        %  Change loref to be a little more realistic (so that the initial
        %  lengths above are between 0 and 2 lorefs)
        c.loref = a*sqrt(2);
      case 'inflated_circle'
        %  Like circle above with just a slightly bigger area to start out
        %  with
        c.xn = bsxfun(@plus,c.xm,1.1*c.liref*cos(ths));
        c.yn = bsxfun(@plus,c.ym,1.1*c.liref*sin(ths));
      case 'star'
        %  Circular.  Note xn is a n_nodes x n_cells matrix
        c.xn = bsxfun(@plus,c.xm,c.liref*cos(ths));
        c.yn = bsxfun(@plus,c.ym,c.liref*sin(ths));
        %  Star-shape.  Especially useful for stability tests
        rx = c.xn./sqrt(c.xn.^2+c.yn.^2);
        ry = c.yn./sqrt(c.xn.^2+c.yn.^2);
        eveninds = mod(1:c.n_enodes,2) == 0;
        c.xn(eveninds) = c.xn(eveninds)+rx(eveninds)/2;
        c.yn(eveninds) = c.yn(eveninds)+ry(eveninds)/2;
        c.xn(~eveninds) = c.xn(~eveninds)-rx(~eveninds)/2;
        c.yn(~eveninds) = c.yn(~eveninds)-ry(~eveninds)/2;
      case 'crescent'
        %  Crescent shape.  Noncircular fat crescent-like shape.  May not be
        %  good when we enforce nonlinear elasticities (as lengths are never
        %  allowed to be > lomax but below may not enforce that)
        c.xn = bsxfun(@plus,c.xm,c.liref*cos(ths)+0.2*c.liref*cos(2*ths));
        c.yn = bsxfun(@plus,c.ym,c.liref*sin(ths));
        %  If one hands in xn and yn, the initial starting position is
        %  automatically written over.
      otherwise
        for vac = 1:2:numel(varargin)
          if strcmpi('xn',varargin{vac})
            c.xn = varargin{vac+1};
          elseif strcmpi('yn',varargin{vac})
            c.yn = varargin{vac+1};
          end
        end
        if ~isfield(c,'xn') | ~isfield(c,'yn')
          error(['You didn''t choose an appropriate "initial_cell_type"',...
            'or neglected to feed in xn and yn!']);
        end
    end
    %  If we want multiple internal nodes, add them now
    c = setup_int_nodes(c);
      
    %  Calculate important initial/reference configuration values
    c = node_info(c,s);
    
    %  Store reference lengths, angles, and areas for later use:
    c.lrefs = c.ls;
    c.alphrefs = c.alphs;
    c.earearefs = c.eareas;
    %  Currently our "reference angle" for external elements is actually
    %  assumed to be "pi"...that is a "flat" membrane is the natural
    %  configuration and any deviation from flat curvature is resisted.
    %  This may change at a later date and I don't believe this reference
    %  alpha is actually ever used.
    for nc = 1:c.n_enodes
      c.alphrefs{nc}(end) = pi;
    end
    
    %  Node velocities
    switch c.initial_cell_type
      case {'circle','crescent'}
        %  Rigid body motion with characteristic velocity to the right
        c.un = s.v_char*ones(size(c.xn));
        c.vn = zeros(size(c.xn));
%       case {'inflated_circle','star'}
      otherwise
        %  Assume nodes start at rest in this configuration
        c.un = zeros(size(c.xn));
        c.vn = zeros(size(c.xn));
    end
    %  Reference "volume", based on axisymmetric assumptions...could be
    %  based on something else...and initial cell configuration and other
    %  things (see est_3d_volume)
    c.Vref = est_3d_volume(c.xn,c.yn,s.vol_est_type);
  else
    %  Load previous times and areas/volume estimates
    for tc = 0:(c.time_step-1)
      run_info = load([s.data_loc,'time_step_',num2str(tc,'%04g')]);
      ts(tc+1) = run_info.c.t;
      areas(tc+1) = run_info.c.area;
      volests(tc+1) = run_info.c.volest;
    end
    %  Take shape, velocities, and time step from the last loaded data set
    run_info = load([s.data_loc,'time_step_',num2str(c.time_step,'%04g')]);
    c = run_info.c;
    s = run_info.s;
    c.pright_prev = run_info.c.pright_prev;
    if c.time_step > 0
      run_info = load([s.data_loc,'time_step_',...
        num2str(c.time_step-1,'%04g')]);
      cnm2 = run_info.c;
    end
  end
  
  %  If the domain has been rotated about the origin, rotate the cell's
  %  position about the origin as well
  [c.xm,c.ym] = deal(...
    cos(s.domain_rot_ang)*c.xm-sin(s.domain_rot_ang)*c.ym,...
    sin(s.domain_rot_ang)*c.xm+cos(s.domain_rot_ang)*c.ym);
  [c.xn,c.yn] = deal(...
    cos(s.domain_rot_ang)*c.xn-sin(s.domain_rot_ang)*c.yn,...
    sin(s.domain_rot_ang)*c.xn+cos(s.domain_rot_ang)*c.yn);
  
  %  Rotate cell about internal node by this angle
  c.initial_cell_rot_ang = 0;
  
  %  A variable that can be used to "phase in" the global equations for the
  %  global variables using flexpdes staging tools
%   c.fade_in_globs = 0;
  
  %  Redefine any of above using varargin (especially initial_cell_rot_ang
  %  here).  We use setdiff to avoid redefining parameters already
  %  redefined in a similar previous varargin statement.
  fn2 = fieldnames(c);
  fn = setdiff(fieldnames(c),union(fn,{'xn','yn'}));
  for vac = 1:2:numel(varargin)
    if any(strcmp(varargin{vac},fn))
      c.(varargin{vac}) = varargin{vac+1};
    end
  end
  
  %  Rotate the cell about xm and ym (this may mess up initial starting
  %  guesses for nodal velocities though usually that is ok)
  dx = c.xn-mean(c.xm);
  dy = c.yn-mean(c.ym);
  [dx,dy] = deal(...
    cos(c.initial_cell_rot_ang)*dx-sin(c.initial_cell_rot_ang)*dy,...
    sin(c.initial_cell_rot_ang)*dx+cos(c.initial_cell_rot_ang)*dy);
  c.xn = bsxfun(@plus,dx,mean(c.xm));
  c.yn = bsxfun(@plus,dy,mean(c.ym));
  
  %  Probably never used but just in case someone wants to redefine c.xn
  %  and c.yn we include below.
  fn = {'xn','yn'};
  for vac = 1:2:numel(varargin)
    if any(strcmp(varargin{vac},fn))
      c.(varargin{vac}) = varargin{vac+1};
    end
  end
  
  %  Adjust some of the lubrication values
  s.lub_threshold = s.lub_threshold_frac*1;20/c.n_enodes;
  s.lub_max_little_l = s.lub_max_little_l_frac*s.lub_threshold_frac;
  
  %  Initialize "previous time steps" data
  cnm1 = c;
  %  Technically redundant, but to be safe
  if c.time_step == 0
    cnm2 = c;
  end
  
end

function write_system_parameters(s)
  %  We first write any lub_on information.  This is done separately
  %  because lub_on can sometimes be stages (see earlier)
  [~,~,~] = mkdir(s.flexpde_files_loc);
  fid = fopen([s.flexpde_files_loc,'m_lubon.txt'],'w');
  if numel(s.lub_on) == 1
    fprintf(fid,'  lub_on = %g\r\n',adj(s.lub_on));
  else
    fprintf(fid,'  lub_on = staged(');
    for lc = 1:(numel(s.lub_on)-1)
      fprintf(fid,'%g,',adj(s.lub_on(lc)));
    end
    fprintf(fid,'%g)\r\n',adj(s.lub_on(end)));
  end
  s = rmfield(s,'lub_on');
  fclose(fid);
  
  %  Various fields which we ignore
  if isfield(s,'bc_strs'), s = rmfield(s,'bc_strs'); end
  
  %  We then write out other system parameters.  Note that some of these
  %  are not actually used by flexpde, only matlab.
  fid = fopen([s.flexpde_files_loc,'m_system_parameters.txt'],'w');
  fn = fieldnames(s);
  for fnc = 1:numel(fn)
    name = fn{fnc};
    val = s.(name);
    if isscalar(val)
      if isequal(val,round(val))
        fprintf(fid,'  %s = %d\r\n',name,val);
      else
        fprintf(fid,'  %s = %20.20g\r\n',name,adj(val));
      end
    elseif ~ischar(val)
      fprintf(fid,'  %s = array(',name);
      if isequal(val,round(val))
        for vc = 1:(numel(val)-1)
          fprintf(fid,'%d,',val(vc));
        end
        fprintf(fid,'%d',val(end));
      else
        for vc = 1:(numel(val)-1)
          fprintf(fid,'%20.20g,',adj(val(vc)));
        end
        fprintf(fid,'%20.20g',adj(val(end)));
      end
      fprintf(fid,')\r\n');
    end
  end
  fclose(fid);
end

function write_global_variables(c,s)
  my_cell = {'un','vn'};
 
  fid = fopen([s.flexpde_files_loc,'m_global_variables.txt'],'w');
  for mycc = 1:numel(my_cell)
    name = my_cell{mycc};
    for cc = 1:c.n_cells
      for nc = 1:c.n_nodes
        fprintf(fid,'  %s_%s_%s(0.15)\r\n',name,num2str(nc,'%02g'),...
          num2str(cc,'%02g'));
      end
    end
  end
  fclose(fid);
  
  fid = fopen([s.flexpde_files_loc,'m_p_global_variables.txt'],'w');
  if strncmpi(s.area_constraint,'const',5) && (c.n_cells > 1)
    error('Constant area constraint not set up for multiple cells');
%   if s.area_constraint
%     for cc = 1:c.n_cells
%       fprintf(fid,'  pin_%s\r\n',num2str(cc,'%02g'));
%     end
%   else
  elseif ~isequal(s.area_constraint,'none')
    fprintf(fid,'  pright');
  end  
  fclose(fid);
  
end

function write_cell_boundaries(c,s,p)

  fid = fopen([s.flexpde_files_loc,'m_cell_boundaries.txt'],'w');
  for cc = 1:c.n_cells
    for nc = 1:c.n_enodes
      nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
      nnextsuf = sprintf('_%s_%s',num2str(mod(nc,c.n_enodes)+1,'%02g'),...
        num2str(cc,'%02g'));
      fprintf(fid,'    start ''cell%s'' (xn%s,yn%s)\r\n',nsuf,nsuf,nsuf);
      %  This will, in general, give better results but since I do not
      %  recall if it is necessary, I have commented it out for the time
      %  being to maximize code efficiency.
%       fprintf(fid,'      mesh_spacing = lo%s/2\r\n',nsuf);
%       fprintf(fid,'      mesh_spacing = min(lo%s,loref)/2\r\n',nsuf);
      if isempty(strfind(p.system_name,'no_flow'))
        fprintf(fid,'        value(u) = un%s+s_norm%s*(un%s-un%s)\r\n',...
          nsuf,nsuf,nnextsuf,nsuf);
        fprintf(fid,'        value(v) = vn%s+s_norm%s*(vn%s-vn%s)\r\n',...
          nsuf,nsuf,nnextsuf,nsuf);
        if ~isempty(strfind(p.system_name,'fluid_interior'))
          fprintf(fid,'        nobc(pin)\r\n');
        else
          fprintf(fid,'        nobc(p)\r\n');
        end
      end
      fprintf(fid,'      line to (xn%s,yn%s)\r\n',nnextsuf,nnextsuf);
    end
  end
  fclose(fid);

end

%  Not currently used, we found an alternative flexpde approach (not as
%  good perhaps, but good enough for now)
function write_system_boundary(s,max_length)

  fid = fopen([s.flexpde_files_loc,'m_system_boundaries.txt'],'w');
  for bc = 1:2

  end
  fclose(fid);

end

function write_pressure_constraint(c,s)
  
  fid = fopen([s.flexpde_files_loc,'m_pressure_constraint.txt'],'w');
  if strncmpi(s.area_constraint,'const',5)
    %  Chooses the internal pressures in order to conserve cross-sectional
    %  area.
    if c.n_cells > 1
      error('Constant area constraint not set up for multiple cells yet.');
%       for cc = 1:c.n_cells
%         fprintf(fid,'pin_%02g:  dAdt_%02g = (%g)*(1-A_%02g/Aref)',cc,cc,...
%           s.area_constraint_restoration_factor,cc);
%       end
    else
      if isequal(s.area_constraint,'constA')
        fprintf(fid,'pright:  dAdt_01 = (%g)*(1-A_01/Aref)',...
          adj(s.area_constraint_restoration_factor));
      elseif isequal(s.area_constraint,'constV')
        fprintf(fid,'pright:  dVdt_01 = (%g)*(1-Vest_01/Vref))',...
          adj(s.area_constraint_restoration_factor));
      else
        fprintf(fid,['pright:  ((%g)*dAdt_01/Aref+(1-%g)*dVdt_01/Vref-',...
          '(%g)*(1-(%g)*A_01/Aref-(1-%g)*Vest_01/Vref)) = 0'],...
          adj(s.area_frac),adj(s.area_frac),...
          adj(s.area_constraint_restoration_factor),...
          adj(s.area_frac),adj(s.area_frac));
      end
    end
  elseif ~isequal(s.area_constraint,'none')
    %  Chooses the pressure at the left side to be such that the average
    %  pressure(s) on the cell boundary(ies) is zero
    fprintf(fid,'pright:  0 =');
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
        if mod(nc-1,5) == 0, fprintf(fid,'  '); end
        if (mod(nc-1,5) == 0)&&((cc ~= 1)||(nc ~= 1)), fprintf(fid,'  '); end
        fprintf(fid,'line_integral(p,''cell%s'')',nsuf);
        if (cc < c.n_cells) || (nc < c.n_enodes), fprintf(fid,'+'); end
        if (mod(nc-1,5) == 4), fprintf(fid,'\r\n'); end
      end
    end    
  end
  fclose(fid);
  
end

function write_friction_constraint(c,s)
  
  fid = fopen([s.flexpde_files_loc,'m_friction_constraint.txt'],'w');
  fprintf(fid,'friction:  ');
  for cc = 1:c.n_cells
    for nc = 1:c.n_enodes
      nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
      if mod(nc-1,5) == 0, fprintf(fid,'  '); end
      if (mod(nc-1,5) == 0)&&((cc ~= 1)||(nc ~= 1)), fprintf(fid,'  '); end
      fprintf(fid,'line_integral(normal(sigma11,sigma12),''cell%s'')',nsuf);
      if (cc < c.n_cells) || (nc < c.n_enodes), fprintf(fid,'+'); end
      if (mod(nc-1,5) == 4), fprintf(fid,'\r\n'); end
    end
    fprintf(fid,'    =0\r\n');
  end
  fclose(fid);
  
end

function write_summary_report(c,s,p)

  fid = fopen([s.flexpde_files_loc,'m_summary_report.txt'],'w');
  tmp_filename = p.vel_filename;
  tmp_filename = strrep(tmp_filename,'_','');
  tmp_filename = strrep(tmp_filename,'0','');
  tmp_filename = strrep(tmp_filename,'1','');
  tmp_filename = strrep(tmp_filename,'2','');
  tmp_filename = strrep(tmp_filename,'3','');
  tmp_filename = strrep(tmp_filename,'4','');
  tmp_filename = strrep(tmp_filename,'5','');
  tmp_filename = strrep(tmp_filename,'6','');
  tmp_filename = strrep(tmp_filename,'7','');
  tmp_filename = strrep(tmp_filename,'8','');
  tmp_filename = strrep(tmp_filename,'9','');
  fprintf(fid,['  summary export file = ''',...
    tmp_filename,'''\r\n']);
  fprintf(fid,'    report(n_nodes)\r\n    report(n_cells)\r\n');
  for cc = 1:c.n_cells
    csuf = sprintf('_%s',num2str(cc,'%02g'));
    for nc = 1:c.n_nodes
      nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
      fprintf(fid,'    report(un%s) as ''un(%d,%d)''\r\n',nsuf,nc,cc);
      fprintf(fid,'    report(vn%s) as ''vn(%d,%d)''\r\n',nsuf,nc,cc);
    end
  end
  if isfield(c,'pright_prev')
    fprintf(fid,'    report(pright) as ''pright_prev''\r\n');
  end
  fclose(fid);

end

function c = get_vels(c,s,p,initial_run)

    %  Update cell information.  (Note, these defns depend only on
    %  position, not velocity...if they depended on velocity, things might
    %  get tricky.)
    c = node_info(c,s);
    
    %  Find cell/wall nodes that are close to objects (cell/walls), label
    %  them, and calculate useful info for them for later on
    [l,c.d] = get_lubrication_info(c,s);

    %  Write down variable (or constant if no two objects are close to
    %  each other) penalty function for compressibility condition:  
    write_incompressibility_penalty_func(c,s);

    %  Rewrite the boundary of the domain if so desired for a particular
    %  system.  The way we know it is so desired is if s.bc_strs is
    %  actually defined.  Not yet used...see flexpde workaround.
    if isfield(s,'bc_strs')
      write_system_boundary(s,min(c.d(:))/p.dist_to_bound_frac);
    end

    %  Write the global equations
    write_global_equations(c,s,p);

    %  Write out the cell information
    write_cell_info(c,s,l,p);

    %  Write some ics to help speed up calcs
    write_initial_conds(c,s);
    
    %  Run fpde a couple of times to get nodal velocities
    run_flexpde(c,s,p,initial_run);
    
    %  Read the flexpde nodal velocities
    [c.un,c.vn,c.pright_prev] = read_nodal_velocities(s,p);

    %  If requested, estimate the condition number
    if p.find_cond_num
      if p.flexpde_version == 7
        [c.cond_num_matrix,c.cond_num_eigs] = find_cond_nums(s,p);
      else
        warning('Can''t find condition number using version 6 yet!');
      end
    end

end

function write_incompressibility_penalty_func(c,s)

  of1 = fopen([s.flexpde_files_loc,'m_incompressibility.txt'],'w');
 
  %  This turns on stiffer "incompressibility" if the distance between the
  %  two cells is less than lub_threshold.  This is done by increasing the
  %  value of the penalty K in the div(grad(P)) = K*div([u,v]) equation in
  %  regions where the fluid is close to an object (cell or wall). The
  %  higher K is in this equation, the more the solution tries to maintain
  %  incompressibility.
  %
  %  See the text file above for the resulting formula but know that it
  %  finds the distance from a point in the fluid to the closest object
  %  (wall or cell)
  %
  %  Don't turn on the variable penalty function unless we are within this
  %  distance of the wall
  if (c.d < s.var_pen_dist_frac*s.lub_threshold)
    %  Write some useful lengths, magnitudes, dot, and cross products
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        ncp1 = mod(nc,c.n_enodes)+1;
	      fprintf(of1,['lvmag_%02d_%02d = ',...
          'sqrt((xn_%02d_%02d-xn_%02d_%02d)^2+',...
          '(yn_%02d_%02d-yn_%02d_%02d)^2)\r\n'],...
          nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc);
	      fprintf(of1,['pvmag_%02d_%02d = sqrt((x-xn_%02d_%02d)^2+',...
          '(y-yn_%02d_%02d)^2)\r\n'],nc,cc,nc,cc,nc,cc);
        fprintf(of1,['pvdotlv_%02d_%02d = ',...
          '((x-xn_%02d_%02d)*(xn_%02d_%02d-xn_%02d_%02d)+',...
          '(y-yn_%02d_%02d)*(yn_%02d_%02d-yn_%02d_%02d))\r\n'],nc,cc,...
          nc,cc,ncp1,cc,nc,cc,nc,cc,ncp1,cc,nc,cc);
        fprintf(of1,['pvcrosslv_%02d_%02d = ',...
          '((x-xn_%02d_%02d)*(yn_%02d_%02d-yn_%02d_%02d)-',...
          '(y-yn_%02d_%02d)*(xn_%02d_%02d-xn_%02d_%02d))\r\n'],nc,cc,...
          nc,cc,ncp1,cc,nc,cc,nc,cc,ncp1,cc,nc,cc);
      end
    end
    fprintf(of1,'\r\nk = ');
    %  Use those lengths, magnitudes, dot, and cross products to produce a
    %  long if statement giving end formula for the variable penatly
    %  function
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        ncp1 = mod(nc,c.n_enodes)+1;
	      fprintf(of1,'min(if (0 < pvdotlv_%02d_%02d)',nc,cc);
	      fprintf(of1,'and (pvdotlv_%02d_%02d < lvmag_%02d_%02d^2)\r\n',...
          nc,cc,nc,cc);
	      fprintf(of1,...
          'then abs(pvcrosslv_%02d_%02d)/lvmag_%02d_%02d\r\n',...
          nc,cc,nc,cc);
	      fprintf(of1,['else min(pvmag_%02d_%02d,',...
          'pvmag_%02d_%02d),\r\n'],nc,cc,ncp1,cc);
      end
    end
    %  Right sided parentheses to match all the previous left-sided
    %  parentheses
    fprintf(of1,'1\r\n');
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        fprintf(of1,')');
      end
    end
    fprintf(of1,'\r\n');
  else
    %  Just make K = 1 if no objects are close to any other objects
    fprintf(of1,'k=1\r\n');
  end
  
  %  This is an "incompressibility measure" for the cell itself. It is
  %  easy to obtain a formula for the area of the cell (a trapezoid-like
  %  integration rule)...the rate at which the area changes can then be
  %  found using the chain rule on the summation formula.  The resulting
  %  formula involves only nodal positions and nodal velocities (as seen
  %  below/in the txt file)
  %  It can (and has in the past) be used to guarantee the cross-sectional
  %  area of the cell does not change over time steps
  %  It can also be used, if desired, to estimate the rate at which area is
  %  changing with respect to time
  %  Additional quantities added for volume constraints
  cinds = 1:c.n_enodes;
  cindsp1 = [2:c.n_enodes,1];
  for cc = 1:c.n_cells
    A = polyarea(c.xn(:,cc),c.yn(:,cc));
    yc = sum((c.yn(cinds,cc)+c.yn(cindsp1,cc)).*...
      (c.xn(cinds,cc).*c.yn(cindsp1,cc)-...
      c.xn(cindsp1,cc).*c.yn(cinds,cc)))/6/A;
    fprintf(of1,'yc_%02d = %g\r\n',cc,yc);
    fprintf(of1,'dAdt_%02d = ',cc);
    for nc = 1:c.n_enodes
      ncp1 = mod(nc,c.n_enodes)+1;
      fprintf(of1,['(-(vn_%02d_%02d+vn_%02d_%02d)/2)*',...
        '(xn_%02d_%02d-xn_%02d_%02d)+',...
        '(-(yn_%02d_%02d+yn_%02d_%02d)/2)*',...
        '(un_%02d_%02d-un_%02d_%02d)'],...
        ncp1,cc,nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc);
      if nc < c.n_enodes,
        fprintf(of1,'+\r\n');
      else
        fprintf(of1,'\r\n');
      end
    end
    fprintf(of1,['dycdt_%02d = (-1/A_%02d)*dAdt_%02d*yc_%02d',...
      '+1/(6*A_%02d)*('],cc,cc,cc,cc,cc);
    for nc = 1:c.n_enodes
      ncp1 = mod(nc,c.n_enodes)+1;
      fprintf(of1,['(vn_%02d_%02d+vn_%02d_%02d)*',...
        '(xn_%02d_%02d*yn_%02d_%02d-xn_%02d_%02d*yn_%02d_%02d)+\r\n',...
        '(yn_%02d_%02d+yn_%02d_%02d)*',...
        '(un_%02d_%02d*yn_%02d_%02d-un_%02d_%02d*yn_%02d_%02d+',...
        'xn_%02d_%02d*vn_%02d_%02d-xn_%02d_%02d*vn_%02d_%02d)'],...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc,ncp1,cc,nc,cc,...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc,ncp1,cc,nc,cc,...
        nc,cc,ncp1,cc,ncp1,cc,nc,cc);
      if nc < c.n_enodes
        fprintf(of1,'+\r\n');
      else
        fprintf(of1,')\r\n');
      end
    end
    %  If we are enforcing symmetry, we can use this which estimates the
    %  derivative of the volume of the cell if we assume a cell that is
    %  axisymmetric about the y-axis.
    if s.enforce_symmetry
      if c.n_enodes/2 ~= round(c.n_enodes/2)
        error('When using symmetry must have even number of outer nodes!');
      end
      fprintf(of1,'dVdt_%02d = pi/6*((',cc);
      for nc = 1:c.n_enodes
        ncp1 = mod(nc,c.n_enodes)+1;
        fprintf(of1,['(un_%02d_%02d-un_%02d_%02d)*',...
        '(yn_%02d_%02d^2+yn_%02d_%02d^2+',...
        'yn_%02d_%02d*yn_%02d_%02d)+\r\n',...
        '(xn_%02d_%02d-xn_%02d_%02d)*',...
        '(2*yn_%02d_%02d*vn_%02d_%02d+2*yn_%02d_%02d*vn_%02d_%02d+',...
        'vn_%02d_%02d*yn_%02d_%02d+yn_%02d_%02d*vn_%02d_%02d)'],...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc,ncp1,cc,...
        nc,cc,ncp1,cc,nc,cc,nc,cc,ncp1,cc,ncp1,cc,...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc);
        if nc == c.n_enodes/2
          fprintf(of1,')-\r\n(');
        elseif nc == c.n_enodes
          fprintf(of1,'))\r\n');
        else
          fprintf(of1,'+\r\n');
        end
      end
    end
  end
  
  fclose(of1);

end

function write_incompressibility_penalty_func_old(c,s)

  of1 = fopen([s.flexpde_files_loc,'m_incompressibility.txt'],'w');
 
  %  This turns on stiffer "incompressibility" if the distance between the
  %  two cells is less than lub_threshold.  This is done by increasing the
  %  value of the penalty K in the div(grad(P)) = K*div([u,v]) equation in
  %  regions where the fluid is close to an object (cell or wall). The
  %  higher K is in this equation, the more the solution tries to maintain
  %  incompressibility.
  %
  %  See the text file above for the resulting formula but know that it
  %  finds the distance from a point in the fluid to the closest object
  %  (wall or cell)
  %
  %  Don't turn on the variable penalty function unless we are within this
  %  distance of the wall
  if (c.d < s.var_pen_dist_frac*s.lub_threshold)
    %  Write some useful lengths, magnitudes, dot, and cross products
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        ncp1 = mod(nc,c.n_enodes)+1;
	      fprintf(of1,['lvmag_%02d_%02d = ',...
          'sqrt((xn_%02d_%02d-xn_%02d_%02d)^2+',...
          '(yn_%02d_%02d-yn_%02d_%02d)^2)\r\n'],...
          nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc);
	      fprintf(of1,['pvmag_%02d_%02d = sqrt((x-xn_%02d_%02d)^2+',...
          '(y-yn_%02d_%02d)^2)\r\n'],nc,cc,nc,cc,nc,cc);
        fprintf(of1,['pvdotlv_%02d_%02d = ',...
          '((x-xn_%02d_%02d)*(xn_%02d_%02d-xn_%02d_%02d)+',...
          '(y-yn_%02d_%02d)*(yn_%02d_%02d-yn_%02d_%02d))\r\n'],nc,cc,...
          nc,cc,ncp1,cc,nc,cc,nc,cc,ncp1,cc,nc,cc);
        fprintf(of1,['pvcrosslv_%02d_%02d = ',...
          '((x-xn_%02d_%02d)*(yn_%02d_%02d-yn_%02d_%02d)-',...
          '(y-yn_%02d_%02d)*(xn_%02d_%02d-xn_%02d_%02d))\r\n'],nc,cc,...
          nc,cc,ncp1,cc,nc,cc,nc,cc,ncp1,cc,nc,cc);
      end
    end
    fprintf(of1,'\r\nk = ');
    %  Use those lengths, magnitudes, dot, and cross products to produce a
    %  long if statement giving end formula for the variable penatly
    %  function
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        ncp1 = mod(nc,c.n_enodes)+1;
	      fprintf(of1,'min(if (0 < pvdotlv_%02d_%02d)',nc,cc);
	      fprintf(of1,'and (pvdotlv_%02d_%02d < lvmag_%02d_%02d^2)\r\n',...
          nc,cc,nc,cc);
	      fprintf(of1,...
          'then abs(pvcrosslv_%02d_%02d)/lvmag_%02d_%02d\r\n',...
          nc,cc,nc,cc);
	      fprintf(of1,['else min(pvmag_%02d_%02d,',...
          'pvmag_%02d_%02d),\r\n'],nc,cc,ncp1,cc);
      end
    end
    %  Right sided parentheses to match all the previous left-sided
    %  parentheses
    fprintf(of1,'1\r\n');
    for cc = 1:c.n_cells
      for nc = 1:c.n_enodes
        fprintf(of1,')');
      end
    end
    fprintf(of1,'\r\n');
  else
    %  Just make K = 1 if no objects are close to any other objects
    fprintf(of1,'k=1\r\n');
  end
  
  %  This is an "incompressibility measure" for the cell itself. It is
  %  easy to obtain a formula for the area of the cell (a trapezoid-like
  %  integration rule)...the rate at which the area changes can then be
  %  found using the chain rule on the summation formula.  The resulting
  %  formula involves only nodal positions and nodal velocities (as seen
  %  below/in the txt file)
  %  It can (and has in the past) be used to guarantee the cross-sectional
  %  area of the cell does not change over time steps
  %  It can also be used, if desired, to estimate the rate at which area is
  %  changing with respect to time
  %  Additional quantities added for volume constraints
  cinds = 1:c.n_enodes;
  cindsp1 = [2:c.n_enodes,1];
  for cc = 1:c.n_cells
    A = polyarea(c.xn(:,cc),c.yn(:,cc));
    yc = sum((c.yn(cinds,cc)+c.yn(cindsp1,cc)).*...
      (c.xn(cinds,cc).*c.yn(cindsp1,cc)-...
      c.xn(cindsp1,cc).*c.yn(cinds,cc)))/6/A;
    fprintf(of1,'yc_%02d = %g\r\n',cc,yc);
    fprintf(of1,'dAdt_%02d = ',cc);
    for nc = 1:c.n_enodes
      ncp1 = mod(nc,c.n_enodes)+1;
      fprintf(of1,['(-(vn_%02d_%02d+vn_%02d_%02d)/2)*',...
        '(xn_%02d_%02d-xn_%02d_%02d)+',...
        '(-(yn_%02d_%02d+yn_%02d_%02d)/2)*',...
        '(un_%02d_%02d-un_%02d_%02d)'],...
        ncp1,cc,nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc);
      if nc < c.n_enodes,
        fprintf(of1,'+\r\n');
      else
        fprintf(of1,'\r\n');
      end
    end
    fprintf(of1,['dycdt_%02d = (-1/A_%02d)*dAdt_%02d*yc_%02d',...
      '+1/(6*A_%02d)*('],cc,cc,cc,cc,cc);
    for nc = 1:c.n_enodes
      ncp1 = mod(nc,c.n_enodes)+1;
      fprintf(of1,['(vn_%02d_%02d+vn_%02d_%02d)*',...
        '(xn_%02d_%02d*yn_%02d_%02d-xn_%02d_%02d*yn_%02d_%02d)+\r\n',...
        '(yn_%02d_%02d+yn_%02d_%02d)*',...
        '(un_%02d_%02d*yn_%02d_%02d-un_%02d_%02d*yn_%02d_%02d+',...
        'xn_%02d_%02d*vn_%02d_%02d-xn_%02d_%02d*vn_%02d_%02d)'],...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc,ncp1,cc,nc,cc,...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc,ncp1,cc,nc,cc,...
        nc,cc,ncp1,cc,ncp1,cc,nc,cc);
      if nc < c.n_enodes
        fprintf(of1,'+\r\n');
      else
        fprintf(of1,')\r\n');
      end
    end
    %  If we are enforcing symmetry, we can use this which estimates the
    %  derivative of the volume of the cell if we assume a cell that is
    %  axisymmetric about the y-axis.
    if s.enforce_symmetry
      if c.n_enodes/2 ~= round(c.n_enodes/2)
        error('When using symmetry must have even number of outer nodes!');
      end
      fprintf(of1,'dVdt_%02d = pi/6*((',cc);
      for nc = 1:c.n_enodes
        ncp1 = mod(nc,c.n_enodes)+1;
        fprintf(of1,['(un_%02d_%02d-un_%02d_%02d)*',...
        '(yn_%02d_%02d^2+yn_%02d_%02d^2+',...
        'yn_%02d_%02d*yn_%02d_%02d)+\r\n',...
        '(xn_%02d_%02d-xn_%02d_%02d)*',...
        '(2*yn_%02d_%02d*vn_%02d_%02d+2*yn_%02d_%02d*vn_%02d_%02d+',...
        'vn_%02d_%02d*yn_%02d_%02d+yn_%02d_%02d*vn_%02d_%02d)'],...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc,nc,cc,ncp1,cc,...
        nc,cc,ncp1,cc,nc,cc,nc,cc,ncp1,cc,ncp1,cc,...
        nc,cc,ncp1,cc,nc,cc,ncp1,cc);
        if nc == c.n_enodes/2
          fprintf(of1,')-\r\n(');
        elseif nc == c.n_enodes
          fprintf(of1,'))\r\n');
        else
          fprintf(of1,'+\r\n');
        end
      end
    end
  end
  
  fclose(of1);

end

function write_global_equations(c,s,p)
  %  To save some space, we loop through the following variables rather
  %  than explicitly writing code for each separately.
  my_cell = {'un','vn'};
  
  %  We write the nodal equations for the nodal velocities.  Each equation
  %  effectively includes fluid forces, viscoelastic forces from the
  %  external elements, viscous forces from internal elements, external
  %  bending elasticity, and, occassionally, additional supplemental
  %  lubrication force terms.  The nasty details, however, are hidden in
  %  the expressions for the integrals for the fluid forces are found in
  %  write_cell_info.m/m_cell_info.txt.
  %
  %  The corresponding equations we are writing are given on pg 29 of my
  %  thesis.  There they are written in vector form while here they are
  %  written in component form and sines and cosines are used to extract
  %  the proper components.  Look at m_global_equations.txt to see them.
  fid = fopen([s.flexpde_files_loc,'m_global_equations.txt'],'w');
  %  Loop through un, vn
  for mycc = 1:numel(my_cell)
    name = my_cell{mycc};
    xy = char(name(1)+3);
    %  Loop through all cells
    for cc = 1:c.n_cells
      csuf = sprintf('_%s',num2str(cc,'%02g'));
      for nc = 1:c.n_nodes
        nsuf = sprintf('_%s%s',num2str(nc,'%02g'),csuf);
        nprev = mod(nc-2,c.n_enodes)+1;
        nnext = mod(nc,c.n_enodes)+1;
        fprintf(fid,'  %s%s:  fixed_node%s*%s%s+(1-fixed_node%s)*',...
          name,nsuf,nsuf,name,nsuf,nsuf);
        fprintf(fid,'((1-interp_imposed_force)*f%sn%s',xy,nsuf);
        for connc = 1:numel(c.nc{nc})
          if (c.mumi == 0) & (nc > c.n_enodes)
            if connc == 1
              fprintf(fid,'*0+%s%s',name,nsuf);
            end
          else
            cn = c.nc{nc}(connc);
            cnsuf = sprintf('_%s%s',num2str(cn,'%02g'),csuf);
            fprintf(fid,'+t%ss%s%s*tb%s%s',xy,nsuf,cnsuf,...
              nsuf,cnsuf);
            if nc <= c.n_enodes
              if (cn == nnext) || (cn == nprev)
                fprintf(fid,'-n%ss%s%s*qb%s%s',xy,...
                  nsuf,cnsuf,nsuf,cnsuf);
%               elseif cn == nprev
%                 fprintf(fid,'-n%ss%s%s*qb%s%s',xy,...
%                   nsuf,cnsuf,nsuf,cnsuf);
              end
            end
          end
        end
        if nc <= c.n_enodes
          nsuf = ['_',char(name(1)+3),nsuf];
          nprevsuf = ['_',char(name(1)+3),sprintf('_%s%s',...
            num2str(mod(nc-2,c.n_enodes)+1,'%02g'),csuf)];
          fprintf(fid,'+int%s-int_s_norm%s+int_s_norm%s',nsuf,nsuf,...
            nprevsuf);
        end
        fprintf(fid,') = 0\r\n');
      end
    end
  end
%         %  un and vn stuff
%         for nc = [c.n_enodes,1:c.n_enodes-1]
%           nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
%           nnextsuf = sprintf('_%s_%s',num2str(mod(nc,c.n_enodes)+1,...
%             '%02g'),num2str(cc,'%02g'));
%           if name(1) == 'u'
%             fprintf(fid,['  un%s: (1-interp_imposed_force)*fxn%s+',...
%               'tob%s*cos(th%s)+tob%s*(-cos(th%s))+',...
%               'qb%s*(-sin(th%s))+qb%s*sin(th%s)+int_x%s-',...
%               'int_s_norm_x%s+int_s_norm_x%s+tib%s*cos(ph%s))+',...
%               'fade_in_globs*(un%s-(%20.20g)) = 0\r\n'],...
%               nnextsuf,nnextsuf,...
%               nnextsuf,nnextsuf,nsuf,nsuf,...
%               nnextsuf,nnextsuf,nsuf,nsuf,...
%               nnextsuf,nnextsuf,nsuf,nnextsuf,nnextsuf,nnextsuf,...
%               adj(c.(name)(mod(nc,c.n_enodes)+1,cc)));
%           else
%             if s.enforce_symmetry && ...
%                 ((nc == c.n_enodes) || (nc == c.n_enodes/2))
% %               warning(['Assuming symmetric flow for now ',...
% %                 '(with even # of nodes).']);
%               fprintf(fid,'  vn%s: vn%s = 0\r\n',nnextsuf,nnextsuf);
%             else
%               fprintf(fid,['  vn%s: (1-fade_in_globs)*(fyn%s*(1-interp_imposed_force)+',...
%                 'tob%s*sin(th%s)+tob%s*(-sin(th%s))+qb%s*(cos(th%s))+',...
%                 'qb%s*(-cos(th%s))+int_y%s-int_s_norm_y%s+',...
%                 'int_s_norm_y%s+tib%s*sin(ph%s))',...
%                 '+fade_in_globs*(vn%s-(%20.20g)) = 0\r\n'],...
%                 nnextsuf,nnextsuf,...
%                 nnextsuf,nnextsuf,nsuf,nsuf,...
%                 nnextsuf,nnextsuf,nsuf,nsuf,...
%                 nnextsuf,nnextsuf,nsuf,nnextsuf,nnextsuf,nnextsuf,...
%                 adj(c.(name)(mod(nc,c.n_enodes)+1,cc)));
%             end
%           end
%         end
%       end
%     end
%   end
  fclose(fid);
  
  %  For the square_no_flow situation, we readjust the
  %  m_global_equations.txt file.  Note that we use the Linux specific
  %  "sed" command
  if isequal(p.system_name,'cancer_cell_in_square_no_flow') &&...
      s.enforce_symmetry
    %  Check that n_nodes is even.  Rounding and use of a 1e-6 tolerance is
    %  probably not necessary, but we are just being safe.
    if abs(c.n_enodes/2-round(c.n_enodes/2)) < 1e-6
      opp_nod_ind = round(c.n_enodes/2)+1;
      onsuf = ['_',num2str(opp_nod_ind,'%02g'),'_01'];
      %  Replace the un_01_01 equation with un_01_01 = -un_??_01 where ??
      %  should be the node exactly opposite node 1.  This requires there
      %  is an even number of external nodes
      tmpstr = '';
      for nc = 1:c.n_enodes
        nsuf = ['_',num2str(nc,'%02g'),'_01'];
        tmpstr = [tmpstr,sprintf('un%s+',nsuf)];
      end
      system(['sed -i ''/um_01/c\  um_01: ',tmpstr(1:end-1),'+um_01 = 0',...
        ''' ',s.flexpde_files_loc,'m_global_equations.txt']);
      system(['sed -i ''/vn_01_01/c\  vn_01_01: vn_01_01 = 0'' ',...
        s.flexpde_files_loc,'m_global_equations.txt']);
      system(['sed -i ''/vn',onsuf,'/c\  vn',onsuf,': vn',onsuf,...
        ' = 0'' ',s.flexpde_files_loc,'m_global_equations.txt']);
    end
  end
end

function write_cell_info(c,s,l,p)

  %  Write information
  fid = fopen([s.flexpde_files_loc,'m_cell_info.txt'],'w');
  fn = fieldnames(c);
  
  %  Most node-node info--note that while we "grant" flexpde knowledge of
  %  all these values, flexpde does not necessarily use them all as Matlab
  %  does the necessary calculations (e.g. eareas, earearefs, see
  %  "triareainfo.m" stuff).
  for fnc = 1:numel(fn)
    name = fn{fnc};
    val = c.(name);
    if isequal(class(val),'cell')
      for cc = 1:c.n_cells
        csuf = sprintf('_%s',num2str(cc,'%02g'));
        for nc = 1:c.n_nodes
          nsuf = sprintf('_%s%s',num2str(nc,'%02g'),csuf);
          switch name
            %  These guys don't interact with other cells.
            case {'anfs','eareas','earearefs'}
              for elec = 1:numel(c.eareas{nc})
                elsuf = sprintf('_%02g%s',elec,csuf);
                fprintf(fid,'%s%s_%02g = %30.20g\r\n',name,nsuf,elec,...
                  val{nc}(elec));
              end
            otherwise
              for connc = 1:numel(c.nc{nc})
                cn = c.nc{nc}(connc);
                cnsuf = sprintf('_%s%s',num2str(cn,'%02g'),csuf);
                fprintf(fid,'%s%s%s = %30.20g\r\n',name,nsuf,cnsuf,...
                  val{nc}(connc));
              end
          end
        end
      end
    else
      switch name
        case {'xn','yn','un','vn','dt','nlist','enlist','inlist'}
        otherwise
          val = c.(name);
          if 1 == 1%~ischar(val)
            if ~ischar(val)
              if isequal(val,round(val))
                fprintf(fid,'  %s = %d\r\n',name,val);
              else
                fprintf(fid,'  %s = %20.20g\r\n',name,adj(val));
              end
            end
          end
      end
    end
  end
    
  %  Do some geom calcs
  %  i+1 and i-1 indices for calculating various geometric quantities
  ie = 1:c.n_enodes;
  ip1 = mod(ie,c.n_enodes)+1;
  im1 = mod(ie-2,c.n_enodes)+1;
  %  Below we assumed there may be more than 1 cell
  %  For more than 1 cell, the structure v will have vectors as fields
  %  For more than 1 cell, the structure m will have matrices as fields
  %  Internal node location(s)
%   v.xm = c.xm;
%   v.ym = c.ym;
  %  Current area of cell
  v.A = polyarea(c.xn(1:c.n_enodes),c.yn(1:c.n_enodes));
  %  Current "volume estimate" for cell (see est_3d_volume)
  v.Vest = est_3d_volume(c.xn(1:c.n_enodes),c.yn(1:c.n_enodes),...
    s.vol_est_type);
  %  Background pressure in cell...always set to 0 for now
  v.back_pres = 0*v.A;
  %  Interior pressure in cell (background pressure can be used to adjust
  %  this but we never do as there are better ways)
  if ~(strncmpi(s.area_constraint,'const',5) && (c.n_cells > 1))
    if isequal(s.area_constraint,'Aref')
      v.pin = c.kp*(1-v.A/c.Aref)+v.back_pres;
    else
      v.pin = c.kp*(1-v.Vest/c.Vref)+v.back_pres;
    end
  end
  
  %  External node computations (we may or may not generalize this later)
  m.xn = c.xn;
  m.yn = c.yn;
  
  %  Allows us to fix the nodes in place
  m.fixed_node = false(size(m.xn));
  %  I believe this fixes the central node in place for our initial single
  %  cell setup
%   m.fixed_node(:) = true;
%   m.fixed_node(41) = false;
  
  %  External line segment vector components in x and y direction
  m.dxe = c.xn(ip1,:)-c.xn(ie,:);
  m.dye = c.yn(ip1,:)-c.yn(ie,:);
  %  Length of external and internal line segments
  m.le = sqrt(m.dxe.^2+m.dye.^2);
  %  Corresponding angles of orientation for external and internal segments
  m.th = atan2(m.dye,m.dxe);
  %  Use cross and dot prod to get angle between successive ext segments
  m.mycross = m.dxe(im1,:).*m.dye(ie,:)-m.dye(im1,:).*m.dxe(ie,:);
  m.mydot = m.dxe(ie,:).*m.dxe(im1,:)+m.dye(ie,:).*m.dye(im1,:);
  m.al = atan2(m.mycross,m.mydot);
  for cc = 1:c.n_cells
    for nc = 1:c.n_enodes
      nnext = mod(nc,c.n_enodes)+1;
      for connc = 1:numel(c.nc{nc})
        cn = c.nc{nc}(connc);
        if (cn == nnext)
          m.lrefe(nc,cc) = c.lrefs{nc}(connc);
          m.le(nc,cc) = c.ls{nc}(connc);
        end
      end
    end
  end
  m.dxxe = (m.dxe(ie,:)-m.dxe(im1,:))./m.lrefe.^2;
  m.dyye = (m.dye(ie,:)-m.dye(im1,:))./m.lrefe.^2;
  m.magcurv = sqrt(m.dxxe.^2+m.dyye.^2);
  m.newm = sign(m.al).*c.kbe.*m.magcurv;
  m.oldmlin = c.kbe./c.loref.*m.al;
  m.oldmnonlin = c.kbe./c.loref.*2*tan(m.al/2);
  if isfield(c,'kblin'); error('kblin option no longer used!'); end
  m.m = m.(c.mtype);
  m.qb = (m.m-m.m(ip1,:))./m.le;
  
  %  For debugging only
  m.dxxxxe = (m.dxxe(im1,:)-2*m.dxxe+m.dxxe(ip1,:))./c.loref^2;
  m.dyyyye = (m.dyye(im1,:)-2*m.dyye+m.dyye(ip1,:))./c.loref^2;
  m.fperlen = c.kbe*sqrt(m.dxxxxe.^2+m.dyyyye.^2);
  m.dal = (m.al-m.al(ip1,:));
  m.qbold = c.kbe./c.loref./m.le.*m.dal;
    
%   if c.mtype
%     %  Used when calculating shear stress along external elements due to
%     %  bending elasticity...linear so allows sharp corners
%     m.dal = (m.al-m.al(ip1,:));
%   else
%     %  Nonlinear bending elasticity-no sharp corners because bending moment
%     %  becomes infinite in such cases
%     m.dal = 2*tan(m.al/2)-2*tan(m.al(ip1,:)/2);
%   end
  %  Corresponding shear stress along external elements due to bend elast
%   m.qb = c.kbe./c.loref./m.lo.*m.dal;
  
  %  FORCES ON NODES (2018/06/21)
  if isfield(s,'xw')
    %  We start with a restorative force that brings our cell back closer
    %  to (0,0).  This force should be neglible vs everyone else.
    fxnrest = -s.rest_force_mag*mean(c.xn)*ones(c.n_nodes,c.n_cells)./...
      c.n_enodes;
    fynrest = -s.rest_force_mag*mean(c.yn)*ones(c.n_nodes,c.n_cells)./...
      c.n_enodes;
    %  We now try a force that forces the cell to conform to the
    %  microchannel.  In the square simulations the microchannel does not
    %  explicitly appear in the boundary declarations in flexpde.  It is,
    %  instead, implemented implicitly by these "wall forces" (aka wf's)
    %  below
    %
    %  To implement first we find distances from nodes to the imaginary
    %  microchannel boundary (given by xw and yw) at a certain x-location
    %  along the microchannel (given by xloc)
    %
    %  This finds distances from every cell node to every line segment on
    %  the (imaginary) microchannel boundary (numel(xw)-1) by n_nodes
    [dists,dot_dist,overlap,norm_x,norm_y,norm_arc,type] = ...
      dist_from_pt_to_line_segs(c.xn(1:c.n_enodes,:),...
      c.yn(1:c.n_enodes,:),s.xw,s.yw);
    %  This selects out only the segments on the boundary that the cell
    %  nodes are closest to: mindists-vector of n_nodes with distances and
    %  inds-vector of indices corresponding to the line segment on the
    %  boundary that is closest to each node
    [mindists,inds] = min(dists,[],1);
    %  Convert these indices into "linear indices" for easier use with
    %  matrices
    lininds = sub2ind(size(dists),inds,1:c.n_enodes);
    %  Decide if the cell nodes are inside (1) or outside (0) of the
    %  boundary
    my_eps = 0.1;
    in_or_out = inpolygon(c.xn(1:c.n_enodes,:),c.yn(1:c.n_enodes,:),...
      s.xw,s.yw);
    %  This defines a force that pushes nodes back towards the wall if they
    %  lie outside the imaginary boundary (xw and yw).  It is proportional
    %  to the distance the node is from the channel and 0 if the node is
    %  already in the channel.  Note that the center node will ideally
    %  always lie inside the boundary so we have effectively made fxmwf =
    %  fymwf = 0.
    %   fxnwf = (s.imag_bnd_force_mag*(~in_or_out').*norm_x(lininds).*...
    %     (max(min(my_eps,mindists+my_eps),0)./my_eps).^(1))';
    %   fynwf = (s.imag_bnd_force_mag*(~in_or_out').*norm_y(lininds).*...
    %     (max(min(my_eps,mindists+my_eps),0)./my_eps).^(1))';
    %   fxnwf = (s.imag_bnd_force_mag*(~in_or_out').*norm_x(lininds).*mindists)';
    %   fynwf = (s.imag_bnd_force_mag*(~in_or_out').*norm_y(lininds).*mindists)';
    %   ramp_func = @(x,e) (x > 0).*(x < e).*x.^2./e.^3.*(3*e-2.*x)+(x >= e);
    ramp_func = @(x,e) (x > -e).*(x < 0).*(x+e).^2.*(e-2*x)./e.^3+(x >= 0);
    fxnwf = zeros(c.n_nodes,c.n_cells);
    fynwf = zeros(c.n_nodes,c.n_cells);
    fxnwf(1:c.n_enodes,:) = (s.imag_bnd_force_mag*(norm_x(lininds).*...
      ramp_func((1-2*in_or_out').*mindists,my_eps)))';
    fynwf(1:c.n_enodes,:) = (s.imag_bnd_force_mag*(norm_y(lininds).*...
      ramp_func((1-2*in_or_out').*mindists,my_eps)))';
    %  For debugging only
    %   plot(c.xn,c.yn);
    %   xlims = xlim; ylims = ylim;
    %   hold on;
    %   plot(s.xw-s.xloc,s.yw);
    %   quiver(c.xn,c.yn,norm_x(lininds)',norm_y(lininds)')
    %   axis equal
    %   xlim(xlims); ylim(ylims);
    fxngf = zeros(c.n_nodes,c.n_cells);
    fyngf = zeros(c.n_nodes,c.n_cells);
    if ~isempty(strfind(p.system_name,'grav'))
      for nc = 1:c.n_nodes
        area_ave = 0;
        for cc = 1:numel(c.eareas{nc})
          area_ave = area_ave+c.eareas{nc}(cc);
        end
        area_ave = area_ave/numel(c.eareas{nc});
        area_ave = 1;
        fyngf(nc) = -s.grav*area_ave;
      end
    end
    
    m.fxn = fxnrest+fxnwf+fxngf;
    m.fyn = fynrest+fynwf+fyngf;
    
  else
    m.fxn = 0*c.xn; m.fyn = 0*c.yn;
  end
  
  %  We added on areal forces (2019/05/19)
  if c.n_cells > 1
    error('stuff needs to be fixed for multiple cells');
  end
  for cc = 1:c.n_cells
    for nc = 1:c.n_nodes
      fxnafs(nc,cc) = 0;
      fynafs(nc,cc) = 0;
      for ec = 1:size(c.anfs{nc},1)
        fxnafs(nc,cc) = fxnafs(nc,cc)+c.ka*c.anfs{nc}(ec,1);
        fynafs(nc,cc) = fynafs(nc,cc)+c.ka*c.anfs{nc}(ec,2);
      end
    end
  end
  m.fxn = m.fxn-fxnafs;
  m.fyn = m.fyn-fynafs;
    
  %  Write the "vector" info
  fn = fieldnames(v);
  for fnc = 1:numel(fn)
    name = fn{fnc};
    val = v.(name);
    for vc = 1:numel(val)
      fprintf(fid,'  %s_%s = %20.20g\r\n',name,num2str(vc,'%02g'),...
        adj(val(vc)));
    end
  end
  
  %  Write the "matrix" info
  fn = fieldnames(m);
  for fnc = 1:numel(fn)
    name = fn{fnc};
    val = m.(name);
    for nc = 1:size(val,1)
      for cc = 1:c.n_cells
        fprintf(fid,'  %s_%s_%s = %20.20g\r\n',....
          name,num2str(nc,'%02g'),num2str(cc,'%02g'),...
          adj(val(nc,cc)));
      end
    end
  end
  
  %  Tensions, external arclengths and shears
  for cc = 1:c.n_cells
    csuf = sprintf('_%s',num2str(cc,'%02g'));
    for nc = 1:c.n_nodes
      nsuf = sprintf('_%s%s',num2str(nc,'%02g'),csuf);
      nprev = mod(nc-2,c.n_enodes)+1;
      nnext = mod(nc,c.n_enodes)+1;
      for connc = 1:numel(c.nc{nc})
        cn = c.nc{nc}(connc);
        cnsuf = sprintf('_%s%s',num2str(cn,'%02g'),csuf);
        %  Do tension calculations and print them to file.  Shear
        %  calculations involve a little more work and are found below
        if (nc <= c.n_enodes) && ((cn == nnext) || (cn == nprev))
          tmp = 'e';  else, tmp = 'i'; end
        fprintf(fid,'tb%s%s = kt%s*(ls%s%s/lrefs%s%s-1)+',...
          nsuf,cnsuf,tmp,nsuf,cnsuf,nsuf,cnsuf);
        fprintf(fid,['mum%s*((un%s-un%s)*dxs%s%s+',...
          '(vn%s-vn%s)*dys%s%s)/ls%s%s^2\r\n'],tmp,...
          cnsuf,nsuf,nsuf,cnsuf,cnsuf,nsuf,nsuf,cnsuf,nsuf,cnsuf);
        if tmp == 'e'
          if cn == nnext
            fprintf(fid,'qb%s%s = qb%s\r\n',nsuf,cnsuf,nsuf);
            fprintf(fid,['s_norm%s = sqrt((x-xn%s)^2+(y-yn%s)^2)/',...
              'ls%s%s\r\n'],nsuf,nsuf,nsuf,nsuf,cnsuf);
          elseif cn == nprev
            fprintf(fid,'qb%s%s = qb%s\r\n',nsuf,cnsuf,cnsuf);
          end
        end
      end
    end
  end
  
  %  The above definitions are in terms of just real numbers while the
  %  below definitions are made in terms of other variables.  Again, look
  %  at the corresponding m_cell_info.txt to understand better.
  
  %  Write "vector" definitions that depend on other variables
  %  Here it's just the fluid stress along an external element with f
  %  corresponding to normal stress and g corresponding to shear stress
  %  along an element
  for cc = 1:c.n_cells
    csuf = sprintf('_%s',num2str(cc,'%02g'));
    fprintf(fid,'  f%s = normal(sn)+pin%s\r\n',csuf,csuf);
    fprintf(fid,'  g%s = tangential(sn)\r\n',csuf);
  end
  
  %  Finally, even more complicated definitions involving integrals of 
  %  variables.  These integrals are along line segments.
  %  Loop through the cells
  for cc = 1:c.n_cells
    %  Loop through the line segments on a cell
    for nc = 1:c.n_enodes
      nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
      nnextsuf = sprintf('_%s_%s',num2str(mod(nc,c.n_enodes)+1,'%02g'),...
        num2str(cc,'%02g'));
      csuf = sprintf('_%s',num2str(cc,'%02g'));
      %  Find segments and nodes that we will add lubrication forces to.
      %  Most of the time, these will be empty lists.
      %  On very rare occassion these may actually be relatively
      %  complicated lists corresponding to, for instance,
      %  three nodes close together.  Because these are so rare, we do not
      %  yet have any code implemented to deal with such special cases.
      s_inds = find(l.seglist(cc).inds == nc);
      n_inds = find(l.nodelist(cc).inds == nc);
      %  For now, just take one node for each line segment, the one closest
      %  to that line segment, and just take one segment for each node, the
      %  one closest to that node.  I can't think of ever seeing an example
      %  where we didn't do this in the end.  Note, we often pick up two
      %  segments for a given node when that given node is closest to and
      %  endpoint of those two segments.
      if numel(s_inds) > 1
        keyboard
      end
      if numel(n_inds) > 1
        [~,ind] = min(l.nodelist(cc).dists(n_inds));
        n_inds = n_inds(ind);
      end
      %  Find the "coefficient matrix" in formula 2.19 in dissertation
      %  Loop through cell line segments that are close to other nodes
      if ~isempty(s_inds)
        %  quantity 1 = l_l/h_l in 2.19.  Usually l_l = lub_max_little_l
        %  but when the node opposite the segment is near a node of the
        %  segment, l_l may be smaller (down to 0.5 lub_max_little_l, see
        %  "special cases" in dissertation, pgs 35-36).
        q1s = s.lub_max_little_l*(1/2+...
          l.seglist(cc).overlap(s_inds)/s.lub_max_little_l)./...
          l.seglist(cc).dists(s_inds);
        %  Realizing some possible complications with regards to continuity
        %  with respect to segment-segment angle, we adjust q1s to be
        %  simpler for now (possibly changed later):
        q1s = s.lub_max_little_l*(1/2+1/2)./l.seglist(cc).dists(s_inds);
        %  Ramping function that slowly introduces lubrication forces with
        %  the function being 0 when the distance from the line segment to
        %  the node in question is s.lub_threshold and 1 when the distance
        %  is 0
        r = 1-l.seglist(cc).dists(s_inds)/s.lub_threshold;
        %  For solid-solid interactions
        switch s.ss_type
          case 'linear'
            q2s = s.lub_max_little_l*r;
          case '1/r_zero_slope'
            q2s = s.lub_max_little_l*s.lub_threshold./...
              l.seglist(cc).dists(s_inds).*r.^2;
          case '1/r_nonzero_slope'
            q2s = s.lub_max_little_l*r./l.seglist(cc).dists(s_inds);
          case 'quadratic'
            q2s = s.lub_max_little_l*r.^2;
          otherwise
            error('ss_type not appropriately defined!');
        end
        
        %  The coefficient matrix in 2.19 (there, h_l =
        %  -mu*matrix*vector(u)_r).  s denotes segment
        scxxs = s.lub_coef*r.*(q1s.^3.*l.seglist(cc).norm_x(s_inds).^2+...
          q1s.*l.seglist(cc).norm_y(s_inds).^2);
        scxys = s.lub_coef*r.*(q1s.^3.*l.seglist(cc).norm_x(s_inds).*...
          l.seglist(cc).norm_y(s_inds)-...
          q1s.*l.seglist(cc).norm_x(s_inds).*...
          l.seglist(cc).norm_y(s_inds));
        scyxs = scxys;
        scyys = s.lub_coef*r.*(q1s.^3.*l.seglist(cc).norm_y(s_inds).^2+...
          q1s.*l.seglist(cc).norm_x(s_inds).^2);
%         sssx = (-1).^(l.seglist(cc).cinds(s_inds) ~= cc).*...
%           q2s.*l.seglist(cc).norm_x(s_inds);
%         sssy = (-1).^(l.seglist(cc).cinds(s_inds) ~= cc).*q2s.*...
%           l.seglist(cc).norm_y(s_inds);
        ss_sx = -q2s.*l.seglist(cc).norm_x(s_inds);
        ss_sy = -q2s.*l.seglist(cc).norm_y(s_inds);
      end
      %  Loop through cell nodes that are close to other line segments
      %  Quantities below correspond to same quantities as above.  One
      %  note, for two cells that are close to each other, the below
      %  calculations will be redoing the above calculations but that is ok
      %  for now.  This is because the matrix corresponding to the
      %  lubrication forces for a segment in a node-segment pairing (where
      %  the node and segment are very close to each other) will be the
      %  same as the matrix corresponding to the lubrication forces for the
      %  node in the same node-segment pairing.
      if ~isempty(n_inds)
        q1s = s.lub_max_little_l*(1/2+...
          min(l.nodelist(cc).overlap(n_inds)/s.lub_max_little_l,1/2))./...
          l.nodelist(cc).dists(n_inds);
        %  Realizing some possible complications with regards to continuity
        %  with respect to segment-segment angle, we adjust q1s to be
        %  simpler for now (possibly changed later):
        q1s = s.lub_max_little_l*(1/2+1/2)./l.nodelist(cc).dists(n_inds);
        r = 1-l.nodelist(cc).dists(n_inds)/s.lub_threshold;
        %  For solid-solid interactions
        q2s = s.lub_max_little_l*s.lub_threshold./...
          l.nodelist(cc).dists(n_inds).*r.^2;
        %  n denotes nodes
        ncxxs = r.*(q1s.^3.*l.nodelist(cc).norm_x(n_inds).^2+...
          q1s.*l.nodelist(cc).norm_y(n_inds).^2);
        ncxys = r.*(q1s.^3.*l.nodelist(cc).norm_x(n_inds).*...
          l.nodelist(cc).norm_y(n_inds)-...
          q1s.*l.nodelist(cc).norm_x(n_inds).*...
          l.nodelist(cc).norm_y(n_inds));
        ncyxs = ncxys;
        ncyys = r.*(q1s.^3.*l.nodelist(cc).norm_y(n_inds).^2+...
          q1s.*l.nodelist(cc).norm_x(n_inds).^2);
%         nssx = (-1).^(l.nodelist(cc).cinds(n_inds) ~= cc).*q2s.*...
%           l.nodelist(cc).norm_x(n_inds);
%         nssy = (-1).^(l.nodelist(cc).cinds(n_inds) ~= cc).*q2s.*...
%           l.nodelist(cc).norm_y(n_inds);
        nssx = -q2s.*l.nodelist(cc).norm_x(n_inds);
        nssy = -q2s.*l.nodelist(cc).norm_y(n_inds);
      end
      %  Below we have the line integrals that tell us about the fluid
      %  forces acting near a node.  Integration by parts can be used to
      %  tell us that the important integrals are those found in 2.11 in
      %  the dissertation.  Ultimately, to get those integrals for a given
      %  segment, we 1) integrate the x and y components of fluid stress
      %  ,((sigma+pin*identity matrix)*normal vector to cell line segment),
      %  on the cell to get int_x and int_y values (see txt file) and 2)
      %  integrate the normalized arclength times the x and y components of
      %  fluid stress to get int_s_norm_x and int_s_norm_y.  After figuring
      %  out angles, these integrals end up going into the global/nodal
      %  equations as seen in the "write_global_equations" program and
      %  those equations are basically 2.12 and 2.13 in the disseration.
      fprintf(fid,['  int_x%s = line_integral(',...
        '(sigma11+pin%s)*nxs%s%s+sigma21*nys%s%s+',...
        'interp_imposed_force*((1-s_norm%s)*fxn%s+s_norm%s*fxn%s),',...
        '''cell%s'',1)'],nsuf,csuf,nsuf,nnextsuf,nsuf,nnextsuf,...
        nsuf,nsuf,nsuf,nnextsuf,nsuf);
      if isequal(p.system_name,'cancer_cell_in_tube_fluid_interior')
        %  Use opposite sign (-) because normal points in opposite
        %  direction on the inside of the cell
        fprintf(fid,['-line_integral(',...
          '(sigma11+pin%s)*nxs%s%s+sigma21*nys%s%s,',...
          '''cell%s'',2)'],csuf,nsuf,nnextsuf,nsuf,nnextsuf,nsuf);
      end
      %  x-lubrication forces on the current segment due to a cell/wall
      %  node.  Note that we already calculated the matrix from 2.19 above
      %  but we must be careful...in general to get the lubrication force
      %  acting on a point P due to cell/wall structure near point Q we
      %  must find the relative velocity between the two points: u_r =
      %  u_P-u_Q.  In the case of nodes below, P is a point on a segment
      %  closest to a particularly close node while Q is that node.
      for ac = 1:numel(s_inds)
        %  coefficients out front of 2.19 expression
        fprintf(fid,'-lub_on*mu*(%g*(',scxxs(ac));
        %  Force arising from x-component of velocity = u
        %  Find horizontal velocity u of the point on the segment that is
        %  closest to the node using linear interpolation (point P)
        fprintf(fid,'((1-(%g))*un_%02d_%02d+(%g)*un_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        %  Horizontal velocity u on the nearby node (point Q)
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          %  Nearby node corresponds to a wall, 0 velocity
          fprintf(fid,'-0)');
        else
          %  Nearby node corresponds to cell, use the cell's node's
          %  velocity
          fprintf(fid,'-un_%02d_%02d)',l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)));
        end
        %  Force arising from y-component of vertical velocity = v
        fprintf(fid,'*lub_elements+(%g)*(',adj(scxys(ac)));
        %  Interpolation to find vertical velocity = v of point P
        fprintf(fid,'((1-(%g))*vn_%02d_%02d+(%g)*vn_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        %  Estimate velocity of the corresponding nearby node (see above)
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)*lub_elements+cr*(%g)/mu)',adj(ss_sx));
        else
          fprintf(fid,...
            '-vn_%02d_%02d)*lub_elements+cr*(%g)/mu)',...
            l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)),adj(ss_sx));
        end
      end
      %  x-lubrication forces on a node in the current line segment
      %  resulting from a different segment acting on that node
      %  Here we use the usual +lub_on
      for ac = 1:numel(n_inds)
        %  Force arising from u
        %  Plug in u at the node
        fprintf(fid,'-lub_on*mu*(%g*(un_%02d_%02d-',adj(ncxxs(ac)),nc,cc);    
        %  Find u for the point on the different segment that is closest to
        %  the node
        if l.nodelist(cc).cinds(n_inds(ac)) == 0
          fprintf(fid,'0)');
        else
          fprintf(fid,'((1-(%g))*un_%02d_%02d+(%g)*un_%02d_%02d))',...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            l.nodelist(cc).sinds(n_inds(ac)),...
            l.nodelist(cc).cinds(n_inds(ac)),...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            mod(l.nodelist(cc).sinds(n_inds(ac))+1,c.n_enodes),...
            l.nodelist(cc).cinds(n_inds(ac)));
        end
        %  Force arising from v
        %  Velocity of the node
        fprintf(fid,'*lub_elements+(%g)*(vn_%02d_%02d-',adj(ncxys(ac)),...
          nc,cc);  
        %  Find v for the point on the different segment that is closest to
        %  the node
        if l.nodelist(cc).cinds(n_inds(ac)) == 0
          fprintf(fid,'0)*lub_elements+cr*(%g)/mu)',adj(nssx));
        else
          fprintf(fid,['((1-(%g))*vn_%02d_%02d+(%g)*vn_%02d_%02d))*',...
            'lub_elements+cr*(%g)/mu)'],...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            l.nodelist(cc).sinds(n_inds(ac)),...
            l.nodelist(cc).cinds(n_inds(ac)),...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            mod(l.nodelist(cc).sinds(n_inds(ac))+1,c.n_enodes),...
            l.nodelist(cc).cinds(n_inds(ac)),adj(nssx));
        end
      end
      fprintf(fid,'\r\n');
      %  Integral along segment of (sigma\dot yhat)
      fprintf(fid,['  int_y%s = line_integral(',...
        'sigma21*nxs%s%s+(sigma22+pin%s)*nys%s%s+',...
        'interp_imposed_force*((1-s_norm%s)*fyn%s+s_norm%s*fyn%s),',...
        '''cell%s'',1)'],nsuf,nsuf,nnextsuf,csuf,nsuf,nnextsuf,...
        nsuf,nsuf,nsuf,nnextsuf,nsuf);
      if isequal(p.system_name,'cancer_cell_in_tube_fluid_interior')
        %  Use opposite sign (-) again
        fprintf(fid,['-line_integral(',...
          'sigma21*nxs%s%s+(sigma22+pin%s)*nys%s%s,',...
          '''cell%s'',2)'],nsuf,nnextsuf,csuf,nsuf,nnextsuf,nsuf);
      end
      %  y-lubrication forces acting on current segment due to another
      %  node...see above for more details
      for ac = 1:numel(s_inds)
        fprintf(fid,'-lub_on*mu*(%g*(',adj(scyxs(ac)));
        fprintf(fid,'((1-(%g))*un_%02d_%02d+(%g)*un_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)');
        else
          fprintf(fid,'-un_%02d_%02d)',l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)));
        end
        fprintf(fid,'*lub_elements+(%g)*(',adj(scyys(ac)));
        fprintf(fid,'((1-(%g))*vn_%02d_%02d+(%g)*vn_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)*lub_elements+cr*(%g)/mu)',adj(ss_sy));
        else
          fprintf(fid,...
            '-vn_%02d_%02d)*lub_elements+cr*(%g)/mu)',...
            l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)),adj(ss_sy));
        end
      end
      %  y-lubrication forces acting on a node for the current segment
      %  due to another segment
      for ac = 1:numel(n_inds)
        fprintf(fid,'-lub_on*mu*(%g*(un_%02d_%02d-',adj(ncyxs(ac)),nc,cc);
        if l.nodelist(cc).cinds(n_inds(ac)) == 0
          fprintf(fid,'0)');
        else
          fprintf(fid,'((1-(%g))*un_%02d_%02d+(%g)*un_%02d_%02d))',...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            l.nodelist(cc).sinds(n_inds(ac)),...
            l.nodelist(cc).cinds(n_inds(ac)),...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            mod(l.nodelist(cc).sinds(n_inds(ac))+1,c.n_enodes),...
            l.nodelist(cc).cinds(n_inds(ac)));
        end
        fprintf(fid,'*lub_elements+(%g)*(vn_%02d_%02d-',adj(ncyys(ac)),...
          nc,cc);
        if l.nodelist(cc).cinds(n_inds(ac)) == 0
          fprintf(fid,'0)*lub_elements+cr*(%g)/mu)',adj(nssy));
        else
          fprintf(fid,...
            ['((1-(%g))*vn_%02d_%02d+(%g)*vn_%02d_%02d))*',...
            'lub_elements+cr*(%g)/mu)'],...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            l.nodelist(cc).sinds(n_inds(ac)),...
            l.nodelist(cc).cinds(n_inds(ac)),...
            adj(l.nodelist(cc).normarc(n_inds(ac))),...
            mod(l.nodelist(cc).sinds(n_inds(ac))+1,c.n_enodes),...
            l.nodelist(cc).cinds(n_inds(ac)),adj(nssy));
        end
      end
      fprintf(fid,'\r\n');
      %  integral along segment of (s_norm*sigma\dot xhat)
      fprintf(fid,['  int_s_norm_x%s = line_integral(',...
        's_norm%s*((sigma11+pin%s)*nxs%s%s+sigma21*nys%s%s+',...
        'interp_imposed_force*((1-s_norm%s)*fxn%s+s_norm%s*fxn%s)),',...
        '''cell%s'',1)'],nsuf,nsuf,csuf,nsuf,nnextsuf,nsuf,nnextsuf,...
        nsuf,nsuf,nsuf,nnextsuf,nsuf);
      if isequal(p.system_name,'cancer_cell_in_tube_fluid_interior')
        %  Use opposite sign (-) again
        fprintf(fid,['-line_integral(s_norm%s*(',...
          '(sigma11+pin%s)*nxs%s%s+sigma21*nys%s%s),''cell%s'',2)'],...
          nsuf,csuf,nsuf,nnextsuf,nsuf,nnextsuf,nsuf);
      end
      %  x-lubrication forces on current segment from a different node (see
      %  above).  Note there are no contributions to this integral from
      %  lubrication forces on nodes for the current segment arising from
      %  different segments because at those nodes, s_norm = 0.
      for ac = 1:numel(s_inds)
        fprintf(fid,'-lub_on*mu*(%g)*(%g*(',...
          adj(l.seglist(cc).normarc(s_inds(ac))),adj(scxxs(ac)));
        fprintf(fid,'((1-(%g))*un_%02d_%02d+(%g)*un_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)');
        else
          fprintf(fid,'-un_%02d_%02d)',l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)));
        end
        fprintf(fid,'*lub_elements+(%g)*(',adj(scxys(ac)));
        fprintf(fid,'((1-(%g))*vn_%02d_%02d+(%g)*vn_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)*lub_elements+cr*(%g)/mu)',adj(ss_sx));
        else
          fprintf(fid,['-vn_%02d_%02d)*',...
            'lub_elements+cr*(%g)/mu)'],...
            l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)),adj(ss_sx));
        end
      end
      fprintf(fid,'\r\n');
      %  integral along segment of (s_norm*sigma\dot yhat)
      fprintf(fid,['  int_s_norm_y%s = line_integral(',...
        's_norm%s*(sigma21*nxs%s%s+(sigma22+pin%s)*nys%s%s+',...
        'interp_imposed_force*((1-s_norm%s)*fyn%s+s_norm%s*fyn%s)),',...
        '''cell%s'',1)'],nsuf,nsuf,nsuf,nnextsuf,csuf,nsuf,nnextsuf,...
        nsuf,nsuf,nsuf,nnextsuf,nsuf);
      if isequal(p.system_name,'cancer_cell_in_tube_fluid_interior')
        %  Use opposite sign (-) again
        fprintf(fid,['-line_integral(s_norm%s*(',...
          'sigma21*nxs%s%s+(sigma22+pin%s)*nys%s%s),''cell%s'',2)'],...
          nsuf,nsuf,nnextsuf,csuf,nsuf,nnextsuf,nsuf);
      end
      %  y-lubrication forces acting on current segment from a different
      %  node (see above)
      for ac = 1:numel(s_inds)
        fprintf(fid,'-lub_on*mu*%g*(%g*(',...
          adj(l.seglist(cc).normarc(s_inds(ac))),adj(scyxs(ac)));
        fprintf(fid,'((1-(%g))*un_%02d_%02d+(%g)*un_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)');
        else
          fprintf(fid,'-un_%02d_%02d)',l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)));
        end
        fprintf(fid,'*lub_elements+(%g)*(',adj(scyys(ac)));
        fprintf(fid,'((1-(%g))*vn_%02d_%02d+(%g)*vn_%02d_%02d)',...
          adj(l.seglist(cc).normarc(s_inds(ac))),nc,cc,...
          adj(l.seglist(cc).normarc(s_inds(ac))),mod(nc,c.n_enodes)+1,cc);
        if l.seglist(cc).cinds(s_inds(ac)) == 0
          fprintf(fid,'-0)*lub_elements+cr*(%g)/mu)',adj(ss_sy));
        else
          fprintf(fid,...
            '-vn_%02d_%02d)*lub_elements+cr*(%g)/mu)',...
            l.seglist(cc).ninds(s_inds(ac)),...
            l.seglist(cc).cinds(s_inds(ac)),adj(ss_sy));
        end
      end
      fprintf(fid,'\r\n');
    end
  end
  
  fclose(fid);

end

function write_initial_conds(c,s)

  fid = fopen([s.flexpde_files_loc,'m_glob_initial_vals.txt'],'w');
  if ~isfield(c,'um')
    for cc = 1:c.n_cells
      c.um(cc) = 0; c.vm(cc) = 0;
      for nc = 1:c.n_enodes
        c.un(nc,cc) = 0;
        c.vn(nc,cc) = 0;
      end
    end
  end
  
  for cc = 1:c.n_cells
    csuf = sprintf('_%s',num2str(cc,'%02g'));
    for nc = 1:c.n_nodes
      nsuf = sprintf('_%s_%s',num2str(nc,'%02g'),num2str(cc,'%02g'));
      fprintf(fid,'  un%s = %20.20g\r\n',nsuf,adj(c.un(nc,cc)));
      fprintf(fid,'  vn%s = %20.20g\r\n',nsuf,adj(c.vn(nc,cc)));
    end
  end
  fclose(fid);

end

function run_flexpde(c,s,p,initial_run)

  %  FlexPDE7 introduces a few new tools which we assume we always have
  %  access to.  This includes the INITIAL EQUATIONS capability and the
  %  fact that it seems that FlexPDE7n (the nongui version) does not hang
  %  when there are mesh issues.  Rather it quits without writing any data
  %  to file (which we utilize below).  At the same time, we still
  %  implement a "temp.pde" file which allows us, in cases where flexpde
  %  does not work ideally, to try multiple changes in temp.pde to get
  %  things working without messing up the original pde template.
  filenames = {[s.flexpde_files_loc,p.system_name,'.pde'],...
    [s.flexpde_files_loc,'temp.pde']};

  %  Maximum time that we allow fpde to run for the coupled solve.  If this
  %  time is exceeded we assume that bad things happened and try to
  %  "intervene"
  max_run_time = s.max_run_time_coupled+...
    initial_run*s.max_run_time_initial_add;
  
  %  While there used to be other options, with our staging taking place
  %  within flexpde, there is now really only 2 feasible options:
  %    1)  Run code normal-like
  %    2)  Resort to vandenberg iteration
  success = false;
  
  %  For direct solving...to "turn off" and use iterative solver only just
  %  set directlimit to 0
  if (p.flexpde_version == 7) && (p.directlimit > 0)
    direct_select_args = {'direct','on','directlimit',p.directlimit};
  else
    direct_select_args = {};
  end
  
  %  Select parameters needed to find condition number
  if p.find_cond_num
    cond_select_args = {'debug','mxcond'};
  else
    cond_select_args = {};
  end
  
  %  Try the normal route first:
  [~,~,~] = copyfile(filenames{1},filenames{2},'f');
  select_params = {cond_select_args{:},direct_select_args{:},...
    'gridlimit',s.grid_limit,'regrid',p.regrid};
  %  Actually run it.  Flexpde occassionally fails on startup (crash on
  %  launch) for reasons other than mesh construction.  This loop fixes
  %  that.  It also uses 10*max_time_normal in case flexpde is having a lot
  %  of trouble converging.
  success = run_fpdefile(s.flexpde_files_loc,filenames{2},...
    max_run_time,numel(s.lub_on),s,p,select_params);
  
  %  If for some reason the above didn't work, we try vandenberg
  %  iteration...a last resort, takes significantly longer.  Also, when
  %  using vandenberg, the normal tolerances should be adjusted (to be
  %  smaller I believe) to get the same accuracy levels.  I haven't worked
  %  out exactly how to adjust those tolerances quite yet.
  if ~success
    select_params = {cond_select_args{:},direct_select_args{:},...
      'gridlimit',s.grid_limit,'vandenberg','on','regrid',p.regrid};
    success = run_fpdefile(s.flexpde_files_loc,filenames{2},...
      max_run_time,numel(s.lub_on),s,p,select_params);
  end

  %  None of the three options worked!  Despair!
  if ~success
    warning('FPDE doesn''t seem to be converging!');
    keyboard
  end
  
end

function run_flexpde_old(c,s,p,initial_run)

  %  FlexPDE must successfully complete both mesh construction and the
  %  solution finding process.  In complex situations (e.g. cells are close
  %  to walls or to themselves) its good to proceed in two different steps:
  %  1)  Find ngrid and aspect values (FlexPDE "select" parameters) that
  %  make the mesh construction work and 2)  Find the solution.  Finding
  %  the solution can be further split into two stages (literally "stages",
  %  the term used by FlexPDE).  In the first stage a solution is found by
  %  prescribing the nodal velocities on the cell(s).  In the second stage,
  %  the nodal velocities are then solved for with the previous solution's
  %  details serving to assist with an initial guess for the solution for
  %  this second stage.  Both stages are handled internally in FlexPDE.
  %  We put a bunch of file names in a structure for easy access later.
  %  filenames{1} (with _grid_only)-pde file used to figure out the proper
  %  ngrid and aspect values to use to make the mesh construction process
  %  work.  This file must be separately defined in the flexpde_files
  %  location inside
  %  system_specific_files/system_name/flexpde_filesmaybeasuffix
  %  If it is not defined, the code will proceed.  This may be suitable
  %  when no mesh construction issues are expected but may result in
  %  failure if even just one mesh construction issue is encountered.
  %  filenames{2}-the normal pde file.  Ideally this is constructed with
  %  the two stages mentioned above though it need not be.
  %  filenames{3}-We copy the normal pde file here so that way if something
  %  goes wrong, we can adjust temp.pde as much as we like trying to make
  %  it work before making permanent changes to the normal pde file.
  %  filenames{4}-location of the mesh made by the _grid_only pde file
  %  filenames{5}-ocassionally we copy 
  filenames = {[s.flexpde_files_loc,p.system_name,'_grid_only.pde'],...
    [s.flexpde_files_loc,p.system_name,'.pde'],...
    [s.flexpde_files_loc,'temp.pde'],...
    [s.flexpde_files_loc,p.system_name,p.add_dir,'transfer_grid_only.dat'],...
    [s.flexpde_files_loc,p.system_name,p.add_dir,'transfer.dat']};
  
  %  While usually the default grid settings work just fine, occassionally
  %  they don't.  In that case, we try other grid parameters/settings
  %  (ngrid-controls mesh density with higher numbers more dense;
  %  aspect-max side length on finite element triangle/min side length on
  %  finite element triangle) until hopefully it works.
  %  There is no best way to choose this.
  ngrid = p.ngrid_factor*[20,30,40,60,80,100,150,200];
  aspect = [2,2,2,2,2,4,8,16,32,64];
  %  Maximum time that we allow fpde to run for the mesh construction step.
  %  If this time is exceeded, we assume that fpde had troubles making the
  %  mesh using default grid settings and we try other grid settings.
  max_time_grid_only = s.max_run_time_grid_only+...
    initial_run*s.max_run_time_initial_add;
  
  %  Mesh construction step (only do if a corresponding grid_only exists)
  if exist(filenames{1},'file')
    %  We copy the "template" for the prescribed velocities pde file over
    %  to a new file before running it.  This allows us to play with
    %  temp.pde without having to worry about corrupting the template.
    [~,~,~] = copyfile(filenames{1},filenames{3},'f');
    
    %  Assume no success
    success = false;
    
    for ngc = 1:numel(ngrid)
      write_addl_select(s,'ngrid',ngrid(ngc),'aspect',aspect(ngc));
      %  Run fpde with prescribed velocities
      ptemp = p;
      ptemp.vel_filename = strrep(strrep(strrep(strrep(strrep(...
        strrep(strrep(strrep(strrep(strrep(strrep(p.vel_filename,...
        '_',''),'0',''),'1',''),'2',''),'3',''),'4',''),'5',''),...
        '6',''),'7',''),'8',''),'9','');
      success = run_fpdefile(s.flexpde_files_loc,filenames{3},...
        max_time_grid_only,numel(s.lub_on),s,ptemp);
      %  If successful
      if success
        %  It appears these values for the parameters ngrid and aspect
        %  allow flexpde to produce a successful mesh.  Store them for
        %  later use.
        working_ngrid = ngrid(ngc);
        working_aspect = aspect(ngc);
        break;
      end
    end
    
    %  No success in making the mesh.
    if ~success
      warning('FPDE doesn''t seem to be converging!');
      keyboard
    end
  else
    working_ngrid = ngrid(1);
    working_aspect = aspect(1);
  end
      
  %  Maximum time that we allow fpde to run for the coupled solve.  If this
  %  time is exceeded we assume that bad things happened and try to
  %  "intervene"
  max_time_normal = s.max_run_time_coupled*(1+...
    initial_run*s.max_run_time_initial_add/s.max_run_time_grid_only);
  
  %  While there used to be other options, with our staging taking place
  %  within flexpde, there is now really only 2 feasible options:
  %    1)  Run code normal-like
  %    2)  Resort to vandenberg iteration
  success = false;
  
  %  For direct solving...to "turn off" and use iterative solver only just
  %  set directlimit to 0
  if (p.flexpde_version == 7) && (p.directlimit > 0)
    direct_select_args = {'direct','on','directlimit',p.directlimit};
  else
    direct_select_args = {};
  end
  
  %  Select parameters needed to find condition number
  if p.find_cond_num
    cond_select_args = {'debug','mxcond'};
  else
    cond_select_args = {};
  end
  
  %  Try the normal route first:
  [~,~,~] = copyfile(filenames{2},filenames{3},'f');
  if isequal(p.system_name,'cancer_cell_in_tube_fluid_interior') && ...
      initial_run
    write_addl_select(s,'ngrid',working_ngrid,'aspect',working_aspect,...
      cond_select_args{:},'direct','on','directlimit',15000,...
      'gridlimit',s.grid_limit,'regrid',p.regrid);
  else
    write_addl_select(s,'ngrid',working_ngrid,'aspect',working_aspect,...
      cond_select_args{:},direct_select_args{:},'gridlimit',s.grid_limit,...
      'regrid',p.regrid);
  end
  %  Actually run it.  For direct solving, FlexPDE7 may take a
  %  long time if the matrices involved are big enough...hence
  %  the extra time buffer.
  %  For flexpde7 and windows 8.1 flexpde occassionally fails
  %  for reasons other than mesh construction.  This loop fixes that.
  maxits = 10;
  waittimes = [10*60,max_time_normal];
  waittimes(3:maxits) = 10*max_time_normal;
  itc = 0;
  while (~success) && (itc < maxits)
    itc = itc+1;
    success = run_fpdefile(s.flexpde_files_loc,filenames{3},...
      waittimes(itc),numel(s.lub_on),s,p);
  end
  
  %  If for some reason the above didn't work, we try vandenberg
  %  iteration...a last resort, takes significantly longer
  if ~success
    write_addl_select(s,'ngrid',working_ngrid,'aspect',working_aspect,...
      cond_select_args{:},direct_select_args{:},...
      'gridlimit',s.grid_limit,'vandenberg','on','regrid',p.regrid);
  end
  while (~success) && (itc < maxits)
    itc = itc+1;
    success = run_fpdefile(s.flexpde_files_loc,filenames{3},...
      2*waittimes(itc),numel(s.lub_on),s,p);
  end

  %  None of the three options worked!  Despair!
  if ~success
    warning('FPDE doesn''t seem to be converging!');
    keyboard
  end
  
end

function deltavmag = find_deltavmag(cprev,ccurr,vchar)
  vminfrac = 0.01;
  deltavmagn = ((ccurr.un-cprev.un).^2+...
    (ccurr.vn-cprev.vn).^2)./...
    max((cprev.un.^2+cprev.vn.^2),vminfrac*vchar^2);
  deltavmag = sqrt(max(deltavmagn(:)));
end

function reldvmag = find_deltavmag2(cprev,ccurr,vchartol)
  vprevmagn = sqrt(cprev.un.^2+cprev.vn.^2);
  %  eps command defends against 0/dividing by zero.
  vprevmagn = max(vprevmagn,eps(vprevmagn));
  dvmagn = sqrt((ccurr.un-cprev.un).^2+(ccurr.vn-cprev.vn).^2);
  %  To defend against the cell's velocities becoming very small (it is
  %  very possible that one of the cell's nodes may briefly come to rest
  %  depending on the situation), we have an "absolute tolerance" in the
  %  denominator as well
  reldvmagn = dvmagn./max(vprevmagn,vchartol);
    
  reldvmag = max(reldvmagn(:));
end

function [cn,s] = ab_advance_cell(cnm2,cnm1,s,p)
 
  %  If flexpde gives us some bogus velocities, we do the forward euler
  %  step again with a smaller time step in the hopes that the next time
  %  forward euler doesn't give us the bogus velocities
  dtnm3h = cnm1.dt;
  dtbig = min(s.dt_max,cnm1.dt*s.dt_grow);
  my_rel_tol = (1-s.dt_safety)*s.dt_min;
  dtsma = 0;
  cn = cnm1;
  %  Distance and area tests to eliminate bumping into walls and
  %  oscillating areas
  dminbef = find_distance(cnm1,s);
  if numel(dminbef) > 1
    error('Not ready for multiple cells here!');
  end
  ie = 1:cn.n_enodes;
  ip1 = mod(ie,cn.n_enodes)+1;
  Abef = polyarea(cnm1.xn(ie,:),cnm1.yn(ie,:));
  
  while (dtbig-dtsma) > my_rel_tol
    
    dtnm1h = (dtbig+dtsma)/2;
    
    cn = ab_step_approx(cn,cnm1,cnm2,dtnm1h,dtnm3h,true);
    dminaft = find_distance(cn,s);
    Aaft = polyarea(cn.xn(ie,:),cn.yn(ie,:));
    dxo = cn.xn(ip1,:)-cn.xn(ie,:);
    dyo = cn.yn(ip1,:)-cn.yn(ie,:);
    loaft = sqrt(dxo.^2+dyo.^2);
    deltavmag = find_deltavmag2(cnm1,cn,s.v_char*s.dv_tol);
    
    if any((abs(dminaft-dminbef)./dminbef > s.dist_max_rel_change) ||...
      ((abs(dminaft-dminbef)./dminbef > s.lub_dist_max_rel_change) && ...
      (dminaft < s.lub_threshold)) || ...
      (abs(Abef-Aaft)./Abef > s.area_max_rel_change) || ...
      (((dminaft-s.lub_threshold).*(dminbef-s.lub_threshold) < 0) ...
      && ((dminaft < (1-s.lub_threshold)*s.lub_threshold) || ...
      (dminaft > (1+s.lub_threshold)*s.lub_threshold)))) || ...
      (max(loaft(:)) >= cn.lomax) || ...
      (deltavmag >= s.v_change_max*s.dt_safety)
      dtbig = dtnm1h;
    else
      dtsma = dtnm1h;
    end
  end
  dtnm1h = max(dtsma,s.dt_min);
  cn = ab_step_approx(cn,cnm1,cnm2,dtnm1h,dtnm3h,false);
  
  deltavmag = 2*s.v_change_max;
  while true
    cn = get_vels(cn,s,p,false);
    deltavmag = find_deltavmag2(cnm1,cn,s.v_char*s.dv_tol);
    fprintf('dv = %g; dt = %g; distaft = %g\n.',deltavmag,dtnm1h,dminaft);
    if (deltavmag >= s.v_change_max) && (dtnm1h > s.dt_min/s.dt_safety)
      dtnm1h = fminbnd(@(dtnm1hnew) (find_deltavmag2(cnm1,...
        get_vels_approx(cn,cnm1,cnm2,dtnm1h,dtnm3h,dtnm1hnew),...
        s.v_char*s.dv_tol)-...
        s.v_change_max*s.dt_safety)^2,...
        s.dt_min,dtnm1h);
      cn = ab_step_approx(cn,cnm1,cnm2,dtnm1h,dtnm3h,false);
    else
      break;
    end
  end
  
  if (deltavmag >= s.v_change_max)
    warning('Current solution doesn''t seem to be converging');
%     keyboard
  end
  
  cn.dt = dtnm1h;
  cn.time_step = cnm1.time_step+1;
  tmp = find_distance(cn,s);
  cn.d = min(tmp(:));
  if cn.mumi == 0
    cn.xm = mean(cn.xn);
    cn.ym = mean(cn.yn);
  end
  %  Readjust tube/cell for tube flow
  if ~isempty(strfind(p.system_name,'in_tube'))
    [cn,s] = recenter(cn,s,p);
%     tmp_xmean = mean([cn.xn(:);cn.xm(:)]);
%     s.xb = s.xb-mean(s.xb(:))+tmp_xmean;
%     cn.xn(:) = cn.xn(:)-tmp_xmean;
%     cn.xm(:) = cn.xm(:)-tmp_xmean;
  end
  if s.enforce_symmetry
    %  Also, preserve symmetry for now...not good for off-centered cells!
%     warning(['Assuming symmetric flow for now (with even # of nodes).']);
    cn.yn(1) = cn.yn(1)-0.05*cn.yn(1);
    cn.yn(1+cn.n_nodes/2) = cn.yn(1+cn.n_nodes/2)-...
      0.05*cn.yn(1+cn.n_nodes/2);
    cn.ym = cn.ym-0.05*cn.ym;
    tmpy = mean([cn.yn(2:cn.n_nodes/2),...
      abs(cn.yn(cn.n_nodes:-1:(2+cn.n_nodes/2)))],2);
    cn.yn(2:cn.n_nodes/2) = tmpy;
    cn.yn(cn.n_nodes:-1:(2+cn.n_nodes/2)) = -tmpy;
    tmpx = mean([cn.xn(2:cn.n_nodes/2),...
      cn.xn(cn.n_nodes:-1:(2+cn.n_nodes/2))],2);
    cn.xn(2:cn.n_nodes/2) = tmpx;
    cn.xn(cn.n_nodes:-1:(2+cn.n_nodes/2)) = tmpx;
  end
  
end

function [cn,s] = fe_advance_cell_nonadapt(cnm1,s,p)
 
  %  If flexpde gives us some bogus velocities, we do the forward euler
  %  step again with a smaller time step in the hopes that the next time
  %  forward euler doesn't give us the bogus velocities
  cn = cnm1;

  cn.xm = cnm1.xm+cnm1.dt*cnm1.um;
  cn.xn = cnm1.xn+cnm1.dt*cnm1.un;
  cn.ym = cnm1.ym+cnm1.dt*cnm1.vm;
  cn.yn = cnm1.yn+cnm1.dt*cnm1.vn;
  
  %  Enforces symmetry about x-axis
  if s.enforce_symmetry
    %  Also, preserve symmetry for now...not good for off-centered cells!
    %     warning(['Assuming symmetric flow for now (with even # of nodes).']);
    cn.yn(1) = cn.yn(1)-0.05*cn.yn(1);
    cn.yn(1+cn.n_nodes/2) = cn.yn(1+cn.n_nodes/2)-...
      0.05*cn.yn(1+cn.n_nodes/2);
    cn.ym = cn.ym-0.05*cn.ym;
    tmpy = mean([cn.yn(2:cn.n_nodes/2),...
      abs(cn.yn(cn.n_nodes:-1:(2+cn.n_nodes/2)))],2);
    cn.yn(2:cn.n_nodes/2) = tmpy;
    cn.yn(cn.n_nodes:-1:(2+cn.n_nodes/2)) = -tmpy;
    tmpx = mean([cn.xn(2:cn.n_nodes/2),...
      cn.xn(cn.n_nodes:-1:(2+cn.n_nodes/2))],2);
    cn.xn(2:cn.n_nodes/2) = tmpx;
    cn.xn(cn.n_nodes:-1:(2+cn.n_nodes/2)) = tmpx;
  end
  
  %  In tube, readjust for y-axis symmetry (we shift the cell rather than
  %  shifting the tube each time step)
  [cn,s] = recenter(cn,s,p);
  
  cn = get_vels(cn,s,p,false);
  
  cn.time_step = cnm1.time_step+1;
  tmp = find_distance(cn,s);
  cn.d = min(tmp(:));
  
end

function cnnew = get_vels_approx(cn,cnm1,cnm2,dtnm1h,dtnm3h,dtnm1hnew)
  %  Uses quadratic interpolation to estimate values of velocity
  cnnew = cn;
  coef1 = (dtnm3h^2*dtnm1hnew+dtnm3h*dtnm1hnew^2)/...
    ((dtnm1h+dtnm3h)*dtnm1h*dtnm3h);
  coef2 = ((-dtnm3h-dtnm1h)*dtnm1hnew^2+(dtnm1h^2-dtnm3h^2)*dtnm1hnew+...
    (dtnm1h+dtnm3h)*dtnm1h*dtnm3h)/((dtnm1h+dtnm3h)*dtnm1h*dtnm3h);
  coef3 = (-dtnm1h^2*dtnm1hnew+dtnm1h*dtnm1hnew^2)/...
    ((dtnm1h+dtnm3h)*dtnm1h*dtnm3h);
  cnnew.un = coef1*cn.un+coef2*cnm1.un+coef3*cnm2.un;
  cnnew.vn = coef1*cn.vn+coef2*cnm1.vn+coef3*cnm2.vn;
end

function cn = ab_step_approx(cn,cnm1,cnm2,dtnm1h,dtnm3h,do_vels)
  %  These are exact for ab2 adaptive
  cn.xn = cnm1.xn+dtnm1h*(cnm1.un+dtnm1h/dtnm3h/2*(cnm1.un-cnm2.un));
  cn.yn = cnm1.yn+dtnm1h*(cnm1.vn+dtnm1h/dtnm3h/2*(cnm1.vn-cnm2.vn));
  %  These are estimates for velocities using ab2 adaptive underlying
  %  approximating function
  if do_vels
    cn.un = cnm1.un+dtnm1h/dtnm3h*(cnm1.un-cnm2.un);
    cn.vn = cnm1.vn+dtnm1h/dtnm3h*(cnm1.vn-cnm2.vn);
  end
end

function [caft,s] = adapt_rk2(~,cbef,s,p)

%  Effectively moved to end of function...this makes it so that the data
%  stored is more accurate and, if appropriately implemented, does not
%  change the algorithm (adaptive runge-kutta 2)
%   %  Get velocities
%   cbef = get_vels(cbef,s,p,false);

  %  Initialize
  caft = cbef;
  %  Take time step according to current dt...assumes no issues with
  %  crossing walls/boundaries (see below to fix this if so desired)
  caft.xn = cbef.xn+cbef.dt*cbef.un;
  caft.yn = cbef.yn+cbef.dt*cbef.vn;
  
  %  Get new velocity estimates
  caft = get_vels(caft,s,p,false);
  
  %  Find estimated change in velocities
  reldvmag = find_deltavmag2(cbef,caft,s.v_char*s.dv_tol);
  
  %  Readjust dt using linear interpolation between velocities
  %  Note that it appears that the dt assigned to a given c structure seems
  %  to correspond to the time step taken to get to that current cell
  %  shape, not the time step that occurs after that shape is encountered
  caft.dt = min(s.v_change_max*caft.dt./reldvmag,s.dt_grow*caft.dt);
  caft.dt = min(max(caft.dt,s.dt_min),s.dt_max);
  
  %  Readjust position using an RK2 step
  caft.xn = cbef.xn+caft.dt.*((2*cbef.dt-caft.dt)./2./cbef.dt.*cbef.un+...
    caft.dt./2./cbef.dt*caft.un);
  caft.yn = cbef.yn+caft.dt.*((2*cbef.dt-caft.dt)./2./cbef.dt.*cbef.vn+...
    caft.dt./2./cbef.dt*caft.vn);
  
  %  Interpolate velocities for this new dt...doesn't actually affect time
  %  stepping, but data stored is better/accurate
  caft.un = cbef.un+(caft.un-cbef.un)*caft.dt./cbef.dt;
  caft.vn = cbef.vn+(caft.vn-cbef.vn)*caft.dt./cbef.dt;
  
%   %  Readjust position using this new velocity...this is a FE step though
%   %  we will probably adjust it to be a RK2 step in the near future
%   caft.xm = cbef.xm+caft.dt*cbef.um;
%   caft.xn = cbef.xn+caft.dt*cbef.un;
%   caft.ym = cbef.ym+caft.dt*cbef.vm;
%   caft.yn = cbef.yn+caft.dt*cbef.vn;
  
  %  Find estimated change in velocities (should now equal s.v_change_max)
  reldvmag = find_deltavmag2(cbef,caft,s.v_char*s.dv_tol);

  %  Correction in case in middle of a tube, to match Jose's code
  if ~isempty(strfind(p.system_name,'in_tube'))
    [caft,s] = recenter(caft,s,p);
  end
  
  %  Augment the time step
  caft.time_step = cbef.time_step+1;
  tmp = find_distance(caft,s);
  caft.d = min(tmp(:));
  
  %  Find the velocity (to be used in the next time step procedure...at the
  %  same time, it is the most accurate velocity estimate at the current
  %  time)
  caft = get_vels(caft,s,p,false);
  
  %  Some output
  fprintf('dv = %g; dt = %g; distaft = %g\n.',reldvmag,caft.dt,caft.d);

end

function [caft,s] = fe_advance_cell(cbef,s,dt,filenames)
 
  %  If flexpde gives us some bogus velocities, we do the forward euler
  %  step again with a smaller time step in the hopes that the next time
  %  forward euler doesn't give us the bogus velocities
  while true
    %  Distance and area tests to eliminate bumping into walls and
    %  oscillating areas
    dminbef = find_distance(cbef,s);
    Abef = polyarea(cbef.xn,cbef.yn);
    
    lengthscale = max(cbef.xn(:,1))-min(cbef.xn(:,1));
    %   dt = 0.02*lengthscale./max(max(abs(c.un(:))),max(abs(c.vn(:))));
    caft = cbef;
    caft.xm = cbef.xm+dt*cbef.um;
    caft.xn = cbef.xn+dt*cbef.un;
    caft.ym = cbef.ym+dt*cbef.vm;
    caft.yn = cbef.yn+dt*cbef.vn;
    
    dminaft = find_distance(caft,s);
    Aaft = polyarea(caft.xn,caft.yn);
    if numel(dminaft) > 1
      error('Not ready for multiple cells here!');
    end
    while any(abs(dminaft-dminbef)./dminbef > 0.1) || ...
        (abs(Abef-Aaft)/Abef > 0.01) || ...
        (((dminaft-s.lub_threshold)*(dminbef-s.lub_threshold) < 0) ...
        && ((dminaft < (1-s.lub_threshold)*s.lub_threshold) || ...
        (dminaft > (1+s.lub_threshold)*s.lub_threshold)))
      dt = dt/s.dt_grow;
      caft.xm = cbef.xm+dt*cbef.um;
      caft.xn = cbef.xn+dt*cbef.un;
      caft.ym = cbef.ym+dt*cbef.vm;
      caft.yn = cbef.yn+dt*cbef.vn;
      dminaft = find_distance(caft,s);
      Aaft = polyarea(caft.xn,caft.yn);
    end
    caft = get_vels(caft,s,filenames,false);
    deltavmag = find_deltavmag2(cbef,caft,s.v_char*s.dv_tol);
    fprintf('dv = %g; dt = %g; distaft = %g.\n',deltavmag,dt,dminaft);
    if deltavmag == 0
      dt = min(dt,s.dt_min);
    elseif deltavmag < s.v_change_max
      break;
    else
      dt = 0.5*dt;
      if dt <= s.dt_min
        warning('Current solution doesn''t seem to be converging');
%         keyboard
        break;
      end
    end
  end
  caft.dt = dt;
  caft.d = min(dminaft(:));
  
end

function ccurr = correct_advance_cell(cprev,ccurr,s)
 
  %  A Forward Euler step with some recentering to make sure we don't bump
  %  into walls
  dminbef = find_distance(cprev,s);
  Abef = polyarea(cprev.xn,cprev.yn);
  
  %  ccurr has the estimates for the new velocities but the time step from
  %  the previous predictor step
  dtp = ccurr.dt;
  dtc = dtp;
  caft = ccurr;
  caft.xm = cprev.xm+dtc*((2*dtp-dtc)/2/dtp*cprev.um+dtc/2/dtp*ccurr.um);
  caft.xn = cprev.xn+dtc*((2*dtp-dtc)/2/dtp*cprev.un+dtc/2/dtp*ccurr.un);
  caft.ym = cprev.ym+dtc*((2*dtp-dtc)/2/dtp*cprev.vm+dtc/2/dtp*ccurr.vm);
  caft.yn = cprev.yn+dtc*((2*dtp-dtc)/2/dtp*cprev.vn+dtc/2/dtp*ccurr.vn);
  
  dminaft = find_distance(caft,s);
  Aaft = polyarea(caft.xn,caft.yn);
  while any(abs(dminaft-dminbef)./dminbef > 0.1) || ...
      (abs(Abef-Aaft)/Abef > 0.01)
    dtc = dtc/s.dt_grow;
    caft.xm = cprev.xm+dtc*((2*dtp-dtc)/2/dtp*cprev.um+dtc/2/dtp*ccurr.um);
    caft.xn = cprev.xn+dtc*((2*dtp-dtc)/2/dtp*cprev.un+dtc/2/dtp*ccurr.un);
    caft.ym = cprev.ym+dtc*((2*dtp-dtc)/2/dtp*cprev.vm+dtc/2/dtp*ccurr.vm);
    caft.yn = cprev.yn+dtc*((2*dtp-dtc)/2/dtp*cprev.vn+dtc/2/dtp*ccurr.vn);
    dminaft = find_distance(caft,s);
    Aaft = polyarea(caft.xn,caft.yn);
  end
  ccurr = caft;
  ccurr.dt = dtc;
  
  ccurr.un = (dtc*ccurr.un+(dtp-dtc)*cprev.un)./dtp;
  ccurr.vn = (dtc*ccurr.vn+(dtp-dtc)*cprev.vn)./dtp;
  
end

function store_solution(c,s)
  save([s.data_loc,'time_step_',num2str(c.time_step,'%04g')],'c','s');
end

function [success] = run_fpdefile(parent_dir,filename,max_time,...
  numstages,s,p,select_params)

  %  Linux often, somehow, allows programs to exit and control to return to
  %  matlab/the calling program before all of the files that the called
  %  program is trying to write upon exit are completely written.  I
  %  believe linux stores the "to be written" data in a temporary buffer
  %  upon program exit and later (after program exit) gets around to
  %  writing the data to its final location on the hard drive.  Because of
  %  this, we have to be careful when trying to subsequently read those
  %  files and check to make sure they have been fully written to the hard
  %  drive before trying to proceed.  This may be unnecessary in
  %  windows/mac.
  
  %  There are two files we will eventually read from.  grid_filename has
  %  info regarding the most recent (successful) run of flexpde including,
  %  most importantly, level of grid refinement info.
  if p.flexpde_version == 7
    grid_fileloc = [parent_dir,p.add_dir,'temp_log.txt'];
    vel_fileloc = [parent_dir,p.add_dir,p.vel_filename];
  else
    grid_fileloc = [parent_dir,'temp.log'];
    vel_fileloc = [parent_dir,p.add_dir,p.vel_filename];
  end
  %  Wipe the slate
  if exist(grid_fileloc,'file'), delete(grid_fileloc); end
  if exist(vel_fileloc,'file'), delete(vel_fileloc); end
  
  %  Custom call for Mac (should work for now but this bypasses a lot of
  %  the error checking below)***need to work with Jose for this
  if ismac
    flexpde_location = '/Applications/FlexPDE6/flexpde6.app/Contents/MacOS/flexpde6';
    [status,result] = system([flexpde_location,' -q ',filename]);
    %  Assume success
    success = true;
    return;
  end
  
  %  How often we check to see if flexpde has been completed:
  my_pause = 0.1;
  
  %  Assume failure until proven wrong
  success = false;
  grid_ref = 0;
  
  %  While usually the default grid settings work just fine, occassionally
  %  they don't.  In that case, we try other grid parameters/settings
  %  (ngrid-controls mesh density with higher numbers more dense;
  %  aspect-max side length on finite element triangle/min side length on
  %  finite element triangle) until hopefully it works.
  %  There is no best way to choose this.
  ngridprev = [20,25,30,35,40,60,80,100,150,200];
  aspectprev = [2,2,2,2,2,4,8,16,32,64];
  ngrid = 20:200;
  aspect = interp1(ngridprev,aspectprev,ngrid,'previous');
  
  for ngc = 1:numel(ngrid)
    %  Write flexpde "select" parameters
    write_addl_select(s,'ngrid',ngrid(ngc),'aspect',aspect(ngc),...
      select_params{:});

    %  Call SimpleSystem in an effort to run "filename"
    pid = SimpleSystem(filename,p,s);
    
    %  Keep on checking that the files have been 1) initiated 2) completed
    for myc = 0:my_pause:(max_time-my_pause)
      %  Pause briefly before proceeding (avoids continual checking which
      %  can occupy cpu time/power)
      pause(my_pause);
      %  Check if files exist yet...may not have been fully completed yet if
      %  this check is passed
      if exist(grid_fileloc,'file')
        %  Read off the number of grid refinements
        if ~ispc
          [status,result] = system(['tail -n 100 ',grid_fileloc]);
        else
          [status,result] = system(['PowerShell.exe -windowstyle ',...
            'hidden "Get-Content ',grid_fileloc,' -Tail 100']);
        end
        %  First check if there was a grid failure...hopefully corresponds
        %  always to "Mesh Generation Failed"
        ks = strfind(result,'Mesh Generation Failed');
        if ~isempty(ks)
          fprintf('Grid generation failure.\n');
          break;
        end
        %  Next check if solver finished.  Hopefully corresponds to
        %  "--DONE--" A similar thing is done for the velocity file in
        %  read_nodal_velocities later.
        if exist(vel_fileloc,'file')
          ks = strfind(result,'--DONE--');
          if ~isempty(ks)
            ks = strfind(result,' GRID ');
            %  Assume that if, for some reason, flexpde doesn't print out '
            %  GRID ', then we still completed successfully.
            if ~isempty(ks)
              grid_ref = sscanf(result(ks(end):end),' GRID %d');
              success = true;
              break
            end
          end
        end
        %  No success yet, if unix, try and force system to write more of
        %  the files.  Idea:  file info may be stored in a temporary buffer
        %  and not yet written to hard drive.
        if isunix
          [~,~] = system('sync');
        end
      end
    end
    
    %  Run wasn't successful!  Kill flexpde.
    if ~success
      if isunix
        system(['kill -9 ',num2str(pid)]);
      else
        system(['taskkill /PID ',num2str(pid),' /F']);
      end
    end
    
    %  Run involved too many grid refinements for our taste, try again
    if grid_ref == s.grid_limit
      success = false;
    end
    
    if success
      break;
    end
    
  end
  
  pause(my_pause);

end

function write_addl_select(s,varargin)
  
  fid = fopen([s.flexpde_files_loc,'m_addl_select.txt'],'w');
  for vac = 1:2:numel(varargin)
    name = varargin{vac};
    val = varargin{vac+1};
    switch class(val)
      case 'double'
        fprintf(fid,'  %s = %g\r\n',name,adj(val));
      otherwise
        switch name
          case {'debug','print'}
            fprintf(fid,'  %s(%s)\r\n',name,val);
          otherwise
            fprintf(fid,'  %s %s\r\n',name,val);
        end
    end
  end
  fclose(fid);

end

function pid = SimpleSystem(pdefilename,p,s)
  
  %  FlexPDE call via "system" function-allows us to grab the process ID 
  %  and kill flexpde if it doesn't seem to be converging to a solution for
  %  our system (as measured by amount of time flexpde has spent trying to
  %  solve it, see calling program)
  if ~ispc
    flexpde_location = sprintf(['/home/jarobarb/FlexPDE%d/flexpde%d',...
      p.flexpde_args],p.flexpde_version,p.flexpde_version);
    [status,result] = system([flexpde_location,p.flexpde_args,...
      pdefilename,' & echo $!']);
    pid = str2num(result(1,1:end));
  else
    if (p.flexpde_version == 7) && (s.directlimit > 0)
      flexpde_location = sprintf(['C:\\FlexPDE%d\\FlexPDE%d',...
        strrep(p.flexpde_args(1),' ','')],p.flexpde_version,...
        p.flexpde_version);      
      [status,result] = system(['PowerShell.exe -windowstyle hidden ',...
        '"Start-process ''',flexpde_location,''' -ArgumentList ''-s'','...
        '''',pdefilename,'''" -PassThru']);
    else
      flexpde_location = sprintf('C:\\FlexPDE%d\\FlexPDE%d',...
        p.flexpde_version,p.flexpde_version);      
      [status,result] = system(['PowerShell.exe -windowstyle hidden ',...
        '"Start-process ''',flexpde_location,''' -ArgumentList ''-q'','...
        '''',pdefilename,'''" -PassThru']);      
    end
    
    %  Reads powershell output to obtain process id for killing if needed
    k = strfind(lower(result),lower('FlexPDE'));
    %  Windows 8.1 pid
    %  pid = str2num(result(k-8:k-1));
    %  Windows 10 pid
    pid = str2num(result(k-11:k-4));
  end
%   fprintf('%s\n',result);
  fprintf('%s-pid %g. ',datestr(rem(now,1)),pid);
  
end

function [cond_num_matrix,cond_num_eigs] = find_cond_nums(s,p)
  fid = fopen([s.flexpde_files_loc,p.add_dir,'temp_debug.txt'],'r');
  
  %  This uses past tools which I think are now obsolete
%   matrix_names = {'II','IB','IG','IK','BI','BB','BG','BK',...
%     'GI','GB','GG','GK','KI','KB','KG','KK'};
%   
%   mystr = fgetl(fid);
%   while ~strncmpi(' Matrix MATRIX ** length=',mystr,25)
%     mystr = fgetl(fid);
%   end
%   n = sscanf(mystr,' Matrix MATRIX ** length=%d');
%   
%   for mnc = 1:numel(matrix_names)
%     curr_matrix_name = matrix_names{mnc};
%     while ~strncmpi([curr_matrix_name,' matrix'],mystr,9)
%       mystr = fgetl(fid);
%     end
%     mystr = fgetl(fid);
%     matrix_struc.(curr_matrix_name) = sparse(n,n);
%     while strncmpi(mystr,'ROW',3)
%       [row_vals] = sscanf(mystr,'ROW %d OF %d');
%       i = row_vals(1);
%       m = row_vals(2); n = row_vals(2);
%       mystr = fgetl(fid);
%       while strncmpi(mystr,'(',1)
%         tmp2 = sscanf(mystr,'(%d) %g ');
%         inds = tmp2(1:2:end-1);
%         vals = tmp2(2:2:end);
%         matrix_struc.(curr_matrix_name)(i,inds) = vals;
%         mystr = fgetl(fid);
%       end
%     end
%     if ~isfield(matrix_struc,curr_matrix_name)
%       matrix_struc.(curr_matrix_name) = sparse(m,n);
%     end
%       
%   end
%   
%   %%  Now read in the eigenvalues
%   for ric = 1:2
%     lc = 0;
%     while ~contains(mystr,'Lambda')
%       mystr = fgetl(fid);
%     end
%     mystr = fgetl(fid);
%     k = strfind(mystr,'|');
%     while ~isempty(k)
%       mystr = mystr(k+1:end);
%       tmp = sscanf(mystr,'%g');
%       nls = numel(tmp);
%       lambda(lc+[1:nls],ric) = tmp;
%       lc = lc+nls;
%       mystr = fgetl(fid);
%       k = strfind(mystr,'|');
%     end
%   end

  %  New stuff that uses "mxcond".  This tool (dgecon which uses dgetrf)
  %  does not take advantage of the fact that the underlying matrices are
  %  sparse so it takes a while.  This might be fixed at a later date.
  mystr = fgetl(fid);
  while ~isequal(mystr,-1)
    mystr = fgetl(fid);
    if strncmpi('Condition number = ',mystr,19)
      cond_num_matrix = sscanf(mystr,'Condition number = %g');
    end
  end
  
  fclose(fid);
  
%   BigMatrix = sparse(m,n);
%   for mnc = 1:numel(matrix_names)
%     BigMatrix = BigMatrix+matrix_struc.(matrix_names{mnc});
%   end
%  
%   %  For debugging only
%   close(figure(1)); figure(1);
%   for mnc = 1:numel(matrix_names)
%     subplot(4,4,mnc);
%     spy(matrix_struc.(matrix_names{mnc}));
%   end  
%   close(figure(2)); figure(2);
%   spy(BigMatrix);
%  
%   cond_num_matrix = condest(BigMatrix);
%   if isinf(cond_num_matrix)
%     cond_num_matrix = 9e99;
%   end
%   
%   eigmags = sum(lambda.^2,2);
%   cond_num_eigs = max(eigmags)/min(eigmags);
  cond_num_eigs = cond_num_matrix;
    
end

function [un,vn,pright_prev] = read_nodal_velocities(s,p)

  uinfo_filename = [p.add_dir,p.vel_filename];
  if numel(s.lub_on) > 1
    uinfo_filename = strrep(uinfo_filename,'.',...
      ['_',num2str(numel(s.lub_on)),'.']);
  end
  
  %  We're hoping that if flexpde is still writing stuff to this file then
  %  fopen won't work and fid == -1
  maxits = 10; j = 1;
  while j < maxits
    fid = fopen([s.flexpde_files_loc,uinfo_filename],'r');
    if fid == -1      
      [~,~] = system('sync');
      pause(0.1);
      j = j+1;
    else
      break
    end
  end
  if j == maxits
    keyboard;
  end
  
  my_str = fgetl(fid);
  while ~any(ismember(my_str,'}'))
    my_str = fgetl(fid);
  end
  fgetl(fid);
  
  %  n_cells and n
  eval([fgetl(fid),';'])
  eval([fgetl(fid),';'])
  %  um, vm, un, and vn
  for i = 1:n_cells
    for j = 1:n_nodes
      eval([fgetl(fid),';'])
      eval([fgetl(fid),';'])
    end
  end
  
  eval([fgetl(fid),';'])
  
  fclose(fid);

end

function x = adj(x)
  %  Function adjusts x so that fprintf doesn't print out 3 digits of
  %  accuracy for flexpde's sake...not an issue in linux
  maxval = 9.9999999e99;
  if abs(x) < 1e-99
    x = 0;
  elseif abs(x) > maxval
    x = sign(x)*maxval;
  end
end

function [c,s] = recenter(c,s,p)
  if c.mumi == 0
    c.xm = mean(c.xn(:)); c.ym = mean(c.yn(:));
  else
%    c.xm = cnm2.xm+c.dt*cnm2.um; c.ym = cnm2.ym+c.dt*cnm2.vm;
  end
  %  If we are in a tube, recenter the tube
  if ~isempty(strfind(p.system_name,'in_tube'))
    s.xb = s.xb-mean(s.xb(:))+mean([c.xn(:);c.xm(:)]);
    write_system_parameters(s);
%     cnm1.xn = cnm1.xn-mean(cnm1.xn);
  end
end

function c = node_info(c,s)
  %  Because this info is stored for each node, some info is redundant. For
  %  instance, each length will get stored twice because each edge has two
  %  nodes.  This is for simplicity and should have negligible impact on
  %  efficiency since the solving process should dominate cpu time.
  c.ls = c.nc;
  c.dxs = c.nc;
  c.dys = c.nc;
  c.nxs = c.nc;
  c.nys = c.nc;
  if isfield(c,'eareas'), firsttime = false; else, firsttime = true; end
  %  These should eventually be estimated in a better way (in particular,
  %  we want a circular shape to correspond to a steady configuration) but
  %  for now we just make them 1s.
  for nc = 1:numel(c.nc)
    cnis = c.nc{nc};
    c.dxs{nc} = c.xn(cnis)-c.xn(nc);
    c.dys{nc} = c.yn(cnis)-c.yn(nc);
    c.ls{nc} = sqrt(c.dxs{nc}.^2+c.dys{nc}.^2);
    c.txs{nc} = c.dxs{nc}./c.ls{nc};
    c.tys{nc} = c.dys{nc}./c.ls{nc};
    %  These make unit vectors perpendicular to the tangent vectors so that
    %  "normal" crossed with "tangent" yields 1.  The current
    %  counterclockwise orientation of external nodes makes the normal from
    %  external node i to i+1 point outwards (tangent vector points from
    %  node i to i+1)
    c.nxs{nc} = c.tys{nc};
    c.nys{nc} = -c.txs{nc};
    
%     %  Debugging
%     close(figure(5)); figure(5)
%     quiver(c.xn(nc)*ones(size(cnis)),c.yn(nc)*ones(size(cnis)),...
%       c.txs{nc}'.*c.ls{nc}',c.tys{nc}'.*c.ls{nc}',0);
%     hold on
%     midxs = (c.xn(nc)*ones(size(cnis))+c.xn(cnis)')./2;
%     midys = (c.yn(nc)*ones(size(cnis))+c.yn(cnis)')./2;
%     plot(c.xn(cnis),c.yn(cnis),'x');
%     quiver(midxs,midys,c.nxs{nc}',c.nys{nc}');
%     axis equal;
      
    % Between each pair of connections, there is an angle, which we
    % calculate:
    cnisp1 = circshift(1:numel(c.txs{nc}),-1);
    c.crosses{nc} = c.txs{nc}.*c.tys{nc}(cnisp1)-...
      c.tys{nc}.*c.txs{nc}(cnisp1);
    c.dots{nc} = c.txs{nc}.*c.txs{nc}(cnisp1)+...
      c.tys{nc}.*c.tys{nc}(cnisp1);
    c.alphs{nc} = atan2(c.crosses{nc},c.dots{nc});
    c.alphs{nc} = c.alphs{nc}+2*pi*(c.alphs{nc}<0);
    sum(c.alphs{nc});
    %  Between each pair of connections, there is also an element with
    %  associated info (like area of that element)
    if nc <= c.n_enodes
      axs = c.xn(cnis(1:end-1),:); ays = c.yn(cnis(1:end-1),:);
      bxs = c.xn(cnis(2:end),:); bys = c.yn(cnis(2:end),:);
    else
      axs = c.xn(cnis,:); ays = c.yn(cnis,:);
      bxs = c.xn(circshift(cnis,-1),:); bys = c.yn(circshift(cnis,-1),:);
    end
    xxi = c.xn(nc)*ones(size(axs)); xyi = c.yn(nc)*ones(size(bxs));
    a = [axs,ays]; b = [bxs,bys]; xi = [xxi,xyi];
    if firsttime
      c.eareas{nc} = triareainfo(xi,a,b);
    else
      %  c.anfs-"normalized forces".  This is the "geometric portion" of
      %  the area-related forces and multiplying by ka, the area elastic
      %  force modulus, will recover the actual forces
      [c.eareas{nc},c.anfs{nc}] = ...
        triareainfo(xi,a,b,c.earearefs{nc},1);
    end
  end
  %  Store area and "volume estimate"
  c.area = polyarea(c.xn(1:c.n_enodes),c.yn(1:c.n_enodes));
  c.volest = est_3d_volume(c.xn(1:c.n_enodes),...
    c.yn(1:c.n_enodes),s.vol_est_type);
end

function c = setup_int_nodes(c)

  if c.mult_int_nodes
    boungeominfo = [2;c.n_enodes;c.xn;c.yn];
    nameofbound = ['cell']';
    %  One can use set formulas to intersect and union bounded regions,
    %  but we just have one bounded region here
    setformulas = 'cell';
    [dl,bt] = decsg(boungeominfo,setformulas,nameofbound);
    pdegplot(dl,'EdgeLabels','on','FaceLabels','on');
    model = createpde;
    geometryFromEdges(model,dl);
    generateMesh(model,'Hmax',c.loref','Hmin',c.loref,'Hgrad',1,...
      'GeometricOrder','linear');
    pdemesh(model);
    
    tol = 1e-12;
    [a,b] = ismembertol([c.xn,c.yn],model.Mesh.Nodes',tol,'ByRows',true);
    if ~isequal(b,[1:c.n_enodes]')
      keyboard;
    end
    
    %  Redefine c.xn to correspond to both exterior and interior nodes
    c.xn = model.Mesh.Nodes(1,:)';
    c.yn = model.Mesh.Nodes(2,:)';
    
    %  Element matrix, used to extract reasonable interior elements...shorter
    %  name, easier to use.
    elemat = model.Mesh.Elements;
  else
    c.xn = [c.xn;c.xm];
    c.yn = [c.yn;c.ym];
    
    for enc = 1:(c.n_enodes-1)
      elemat(:,enc) = [enc;c.n_enodes+1;enc+1];
    end
    elemat(:,c.n_enodes) = [c.n_enodes;c.n_enodes+1;1];
  end
  c.n_nodes = numel(c.xn);
  c.n_inodes = c.n_nodes-c.n_enodes;

  %  For the "connection cells", for each connected point the first column 
  %  is its indices in terms of xs and ys, the second column is whether
  %  it is an external (false) or internal (true) node, the third column 
  %  is the index in terms of xm-ym or xn-yn depending on whichever is
  %  appropriate according to the second column.
  for nc = 1:numel(c.xn)
    tmp = elemat(:,find(any(nc == elemat)));
    tmp2 = setdiff(tmp(:),nc);
    if nc <= c.n_enodes
      nnext = mod(nc,c.n_enodes)+1;
      start_ind = nnext;
    else
      start_ind = tmp2(1);
    end
    inds = start_ind;
    for connc = 1:numel(tmp2)-1
      tri_ind = find(any(inds(connc) == tmp),1);
      inds(connc+1) = setdiff(tmp(:,tri_ind),[nc,inds(connc)]);
      [m,n] = size(tmp);
      tmp = tmp(:,setdiff(1:n,tri_ind));
    end
    if ispolycw(c.xn(inds),c.yn(inds))
      inds = circshift(fliplr(inds),1);
    end
    c.nc{nc} = inds;
  end
  
  %  For plotting later, a list of all the indices corresponding to the
  %  line segments
  inds = [1:c.n_enodes]';
  c.enlist = [inds,mod(inds,c.n_enodes)+1];
  tmp = sort(c.enlist,2);
  c.nlist = [];
  for nc = 1:c.n_nodes
    for cnc = 1:numel(c.nc{nc})
      c.nlist = [c.nlist;nc,c.nc{nc}(cnc)];
    end
  end
  %  There are twice as many indices as we actually need for our plots,
  %  sort and eliminate duplicates:
  c.nlist = unique(sort(c.nlist,2),'rows');
  c.inlist = setdiff(c.nlist,tmp,'rows');

  %  Check the connectivities
  do_connectivity_plot = false;
  if do_connectivity_plot
    myc = 'rgbmyck'; myc = [myc,myc]; myc = [myc,myc]; myc = [myc,myc];
    close(figure(1)); figure(1); plot(c.xn,c.yn,'k'); hold on;
    for nc = 1:numel(c.nc)
      for connc = 1:numel(c.nc{nc})
        plot([c.xn(nc),c.xn(c.nc{nc}(connc))],...
          [c.yn(nc),c.yn(c.nc{nc}(connc))],myc(nc),'LineWidth',2);
        pause
      end
    end
  end
            
  warning(['Mesh not yet perfected (e.g. equidistant nodes, ',...
    'equal connectivity, others?']);

end