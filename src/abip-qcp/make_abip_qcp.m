mex_type = 'qcp';
    debug = 0;

platform = convertCharsToStrings(computer('arch'));
link = ' ';

if platform == "win64"
    mkl_lib_path = fullfile(mkl_path, "lib", "intel64");
    fprintf('Linking MKL in Windows \n');
    link = [link, ' -lmkl_intel_lp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "maci64"
    mkl_lib_path = fullfile(mkl_path, "lib");
    fprintf('Linking MKL in MacOS \n');
    link = [link, join(mkl_lib_path + ["/libmkl_intel_lp64.a","/libmkl_core.a", "/libmkl_sequential.a "])];
elseif platform == "glnxa64"
    mkl_lib_path = fullfile(mkl_path, "lib", "intel64");
    fprintf('Linking MKL in Linux \n');
    link = [link, ' -lmkl_intel_ilp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
else
    error('Unsupported platform.\n');
end

lib_path = join('-L' + mkl_lib_path);
mexfname = "abip_" + mex_type;
psrc = fullfile('.', 'source');  

pmex = fullfile('.', 'mex');

mex_file = fullfile(pmex, ['abip_'  mex_type  '_mex.c']);
src_files = dir( [psrc '/*.c'] );
srclist = [];

for i = 1:length(src_files)
    srclist = [srclist,convertCharsToStrings(src_files(i).name)];
end


src = fullfile(psrc, srclist);

pcs = fullfile('.', 'csparse', 'Source');  

cs_files = dir( [pcs '/*.c'] );
cslist = [];
for i = 1:length(cs_files)
    cslist = [cslist,convertCharsToStrings(cs_files(i).name)];
end
cs = fullfile(pcs, cslist);

pldl = fullfile('.', 'qdldl', 'src');  

ldl_files = dir( [pldl '/*.c'] );
ldllist = [];
for i = 1:length(ldl_files)
    ldllist = [ldllist,convertCharsToStrings(ldl_files(i).name)];
end
ldl = fullfile(pldl, ldllist);

pamd = fullfile('.', 'amd');  

amd_files = dir( [pamd '/*.c'] );
amdlist = [];
for i = 1:length(amd_files)
    amdlist = [amdlist,convertCharsToStrings(amd_files(i).name)];
end
amd = fullfile(pamd, amdlist);

src = [src,cs,ldl,mex_file,amd];

pinc = "include";

cs_include = fullfile("csparse", "Include");

ldl_include = fullfile("qdldl", "include");

mkl_include = fullfile(mkl_path, "include");

amd_include = "amd";

inc = [pinc,mkl_include,cs_include,ldl_include, amd_include];
inc = join("-I" + inc);

if(debug == 0)
    debugcommand = "-O";
else
    debugcommand = "-g ";
end
mexcommand = "mex " + debugcommand + " -output " + join([mexfname, lib_path, src, inc, link]);
fprintf("%s\n",mexcommand);
eval(replace(mexcommand, "Program Files (x86)", "'Program Files (x86)'"));

if(mex_type == "ml")
    fprintf("Successfull compiled mex function abip_ml\n");
    fprintf("Usage:[sol,info] = abip_ml(data,settings)\n ");
    fprintf("data struct contais X,y,lambda\n");
end
if(mex_type == "qp")
    fprintf("Successfull compiled mex function abip_qp\n");
    fprintf("Usage:[sol,info] = abip_qp(data,settings) \n ");
    fprintf("data struct contais A,Q,c,rl,ru,lb,ub and L(optional)\n");
end
if(mex_type == "socp")
    fprintf("Successfull compiled mex function abip_socp\n");
    fprintf("Usage:[sol,info] = abip_socp(data,cone,settings)\n ");
    fprintf("data struct contais A,b,c\n");
end
if(mex_type == "qcp")
    fprintf("Successfull compiled mex function abip_qcp\n");
    fprintf("Usage:[sol,info] = abip_qcp(data,cone,settings)\n ");
    fprintf("data struct contais A,Q,b,c\n");
end