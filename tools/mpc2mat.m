function mpc2mat(mpc,mat_file)
%MPC2MAT Converts an "mpc" MATPOWER structure to a CSV file.
%   MPC2MAT(MPT, MAT_FILE)
%   Converts an "mpc" MATPOWER structure to a CSV file.
%
%   Input:
%      MPC : MATPOWER case structure.
%      MAT_FILE: name of output CSV file.

disp('Executing mpc2mat');

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
for i=1:size(mpc.bus)(1)
    fprintf(fid,'%d,',mpc.bus(i,1));
    fprintf(fid,'%d,',mpc.bus(i,2));
    fprintf(fid,'%.5f,',mpc.bus(i,3));
    fprintf(fid,'%.5f,',mpc.bus(i,4));
    fprintf(fid,'%.5f,',mpc.bus(i,5));
    fprintf(fid,'%.5f,',mpc.bus(i,6));
    fprintf(fid,'%d,',mpc.bus(i,7));
    fprintf(fid,'%.5f,',mpc.bus(i,8));
    fprintf(fid,'%.5f,',mpc.bus(i,9));
    fprintf(fid,'%.5f,',mpc.bus(i,10));
    fprintf(fid,'%d,',mpc.bus(i,11));
    fprintf(fid,'%.5f,',mpc.bus(i,12));
    fprintf(fid,'%.5f\n',mpc.bus(i,13));
end
fprintf(fid,'END\n');

%% gen
fprintf(fid,'GEN\n');
fprintf(fid,'bus,Pg (MW),Qg (MVAr),Qmax (MVAr),Qmin (MVAr),Vg (p.u.),mBase (MVA),status,');
fprintf(fid,'Pmax (MW),Pmin (MW)\n')
for i=1:size(mpc.gen)(1)
    fprintf(fid,'%d,',mpc.gen(i,1));
    fprintf(fid,'%.5f,',mpc.gen(i,2));
    fprintf(fid,'%.5f,',mpc.gen(i,3));
    fprintf(fid,'%.5f,',mpc.gen(i,4));
    fprintf(fid,'%.5f,',mpc.gen(i,5));
    fprintf(fid,'%.5f,',mpc.gen(i,6));
    fprintf(fid,'%.5f,',mpc.gen(i,7));
    fprintf(fid,'%d,',mpc.gen(i,8));
    fprintf(fid,'%.5f,',mpc.gen(i,9));
    fprintf(fid,'%.5f\n',mpc.gen(i,10));
end
fprintf(fid,'END\n');

%% branch
fprintf(fid,'BRANCH\n');
fprintf(fid,'bus from,bus to,r (p.u.),x (p.u.),b (p.u.),rateA (MVA),rateB (MVA),rateC (MVA),');
fprintf(fid,'ratio,angle (degrees),status,min angle diff (degrees),max angle diff (degrees)\n')
for i=1:size(mpc.branch)(1)
    fprintf(fid,'%d,',mpc.branch(i,1));
    fprintf(fid,'%d,',mpc.branch(i,2));
    fprintf(fid,'%.5f,',mpc.branch(i,3));
    fprintf(fid,'%.5f,',mpc.branch(i,4));
    fprintf(fid,'%.5f,',mpc.branch(i,5));
    fprintf(fid,'%.5f,',mpc.branch(i,6));
    fprintf(fid,'%.5f,',mpc.branch(i,7));
    fprintf(fid,'%.5f,',mpc.branch(i,8));
    fprintf(fid,'%.5f,',mpc.branch(i,9));
    fprintf(fid,'%.5f,',mpc.branch(i,10));
    fprintf(fid,'%d,',mpc.branch(i,11));
    fprintf(fid,'%.5f,',mpc.branch(i,12));
    fprintf(fid,'%.5f\n',mpc.branch(i,13));
end
fprintf(fid,'END\n');

%% cost
fprintf(fid,'COST\n');
fprintf(fid,'gen index,Q2 ($/hr MW2),Q1 ($/hr MW),Q0 ($/hr)\n');
for i=1:size(mpc.gencost)(1)
    if mpc.gencost(i,1) == 2 && mpc.gencost(i,4) == 3
       fprintf(fid,'%d,',i-1);
       fprintf(fid,'%.8f,',mpc.gencost(i,5));
       fprintf(fid,'%.8f,',mpc.gencost(i,6));
       fprintf(fid,'%.8f\n',mpc.gencost(i,7));
    end
end
fprintf(fid,'END\n');

%% close
disp('Done');
fclose(fid);

