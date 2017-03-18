function mpc2mat(mpc,mat_file)
%MPC2MAT Converts an "mpc" MATPOWER structure to a CSV file.
%   MPC2MAT(MPT, MAT_FILE)
%   Converts an "mpc" MATPOWER structure to a CSV file.
%
%   Input:
%      MPC : MATPOWER case structure (version 2).
%      MAT_FILE: name of output CSV file.
% 
%   Note: 
%      It only supports quadratic costs. If costs are not
%      quadratic, then random quadratic costs are assigned.

disp('Executing mpc2mat');

%% version
if mpc.version ~= '2'
  error('mpc2mat: Only MPC version 2 is supported')
end

%% gen length
sgen = size(mpc.gen);
scost = size(mpc.gencost);  
if scost(1) > 0 && scost(1) ~= sgen(1)
  error('mpc2mat: Bad generator cost data')
end

%% file parts
[path name ext] = fileparts(mat_file);
if ~strcmp(ext,'.mat')
   error('mpc2mat: File extension should be .mat');
end

%% open input file
[fid, msg] = fopen(mat_file,'w');
if fid < 0
   disp(msg);
   error('mpc2mat: Cannot open file:  %s',mat_file);
end

%% case name
fprintf(fid,'%s\n',name);

%% base power
fprintf(fid,'%.0f\n',mpc.baseMVA);

%% bus
fprintf(fid,'BUS\n');
fprintf(fid,'number,type,Pd (MW),Qd (MVAr),Gs (MW),Bs (MVAr),area,Vm (p.u.),Va (degress),');
fprintf(fid,'basekV (kV),zone,maxVm (p.u.),minVm (p.u.)\n');
sbus = size(mpc.bus);
for i=1:sbus(1)
    fprintf(fid,'%d,',mpc.bus(i,1));
    fprintf(fid,'%d,',mpc.bus(i,2));
    fprintf(fid,'%.10e,',mpc.bus(i,3));
    fprintf(fid,'%.10e,',mpc.bus(i,4));
    fprintf(fid,'%.10e,',mpc.bus(i,5));
    fprintf(fid,'%.10e,',mpc.bus(i,6));
    fprintf(fid,'%d,',mpc.bus(i,7));
    fprintf(fid,'%.10e,',mpc.bus(i,8));
    fprintf(fid,'%.10e,',mpc.bus(i,9));
    fprintf(fid,'%.10e,',mpc.bus(i,10));
    fprintf(fid,'%d,',mpc.bus(i,11));
    fprintf(fid,'%.10e,',mpc.bus(i,12));
    fprintf(fid,'%.10e\n',mpc.bus(i,13));
end
fprintf(fid,'END\n');

%% gen
fprintf(fid,'GEN\n');
fprintf(fid,'bus,Pg (MW),Qg (MVAr),Qmax (MVAr),Qmin (MVAr),Vg (p.u.),mBase (MVA),status,');
fprintf(fid,'Pmax (MW),Pmin (MW)\n')
for i=1:sgen(1)
    fprintf(fid,'%d,',mpc.gen(i,1));
    fprintf(fid,'%.10e,',mpc.gen(i,2));
    fprintf(fid,'%.10e,',mpc.gen(i,3));
    fprintf(fid,'%.10e,',mpc.gen(i,4));
    fprintf(fid,'%.10e,',mpc.gen(i,5));
    fprintf(fid,'%.10e,',mpc.gen(i,6));
    fprintf(fid,'%.10e,',mpc.gen(i,7));
    fprintf(fid,'%d,',mpc.gen(i,8));
    fprintf(fid,'%.10e,',mpc.gen(i,9));
    fprintf(fid,'%.10e\n',mpc.gen(i,10));
end
fprintf(fid,'END\n');

%% branch
fprintf(fid,'BRANCH\n');
fprintf(fid,'bus from,bus to,r (p.u.),x (p.u.),b (p.u.),rateA (MVA),rateB (MVA),rateC (MVA),');
fprintf(fid,'ratio,angle (degrees),status,min angle diff (degrees),max angle diff (degrees)\n')
sbranch = size(mpc.branch);
for i=1:sbranch(1)
    fprintf(fid,'%d,',mpc.branch(i,1));
    fprintf(fid,'%d,',mpc.branch(i,2));
    fprintf(fid,'%.10e,',mpc.branch(i,3));
    fprintf(fid,'%.10e,',mpc.branch(i,4));
    fprintf(fid,'%.10e,',mpc.branch(i,5));
    fprintf(fid,'%.10e,',mpc.branch(i,6));
    fprintf(fid,'%.10e,',mpc.branch(i,7));
    fprintf(fid,'%.10e,',mpc.branch(i,8));
    fprintf(fid,'%.10e,',mpc.branch(i,9));
    fprintf(fid,'%.10e,',mpc.branch(i,10));
    fprintf(fid,'%d,',mpc.branch(i,11));
    fprintf(fid,'%.10e,',mpc.branch(i,12));
    fprintf(fid,'%.10e\n',mpc.branch(i,13));
end
fprintf(fid,'END\n');

%% cost
fprintf(fid,'COST\n');
fprintf(fid,'Q2 ($/hr MW2),Q1 ($/hr MW),Q0 ($/hr)\n');
for i=1:scost(1)
    if mpc.gencost(i,1) == 2 && mpc.gencost(i,4) == 3
       fprintf(fid,'%.10e,',mpc.gencost(i,5));
       fprintf(fid,'%.10e,',mpc.gencost(i,6));
       fprintf(fid,'%.10e\n',mpc.gencost(i,7));
    else
       warning('unsupported cost function; using random quadratic costs')
       fprintf(fid,'%.10e,',0.04*rand()+0.01);
       fprintf(fid,'%.10e,',40.*rand()+10.);
       fprintf(fid,'%.10e\n',0.);
    end
end
fprintf(fid,'END\n');

%% close
disp('Done');
fclose(fid);

