using JLD2
using DelimitedFiles

function errorCalc(calc,ref,M) # calculate relative root mean square error

    eps = (sum(abs.(calc-ref).^2)/M).^(1/2)
    refNorm = (sum(abs.(ref).^2)/M).^(1/2);
    epsRel = eps/refNorm;

    return epsRel
end

#==========================================================================================
                            Export and plot results for specific test case
==========================================================================================#
data_file = "runtimes1x1_$(2642)DOFs_$(Int(1000))Hz.JLD2"
f = jldopen(data_file)

hist_pa = f["hist_pa"]
pa = f["pa"]
pasAN = f["pasAN"]
v_rAN_V = f["v_rAN_V"]
v_thetaAN_V = f["v_thetaAN_V"]
v_n0 = f["v_n0"]
vt_sum = f["vt_sum"]
ang_axis = f["ang_axis"]

perm = sortperm(ang_axis)
plt1 = scatter(ang_axis,real.(pa),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Pa)"); plot!(ang_axis[perm],real.(pasAN[perm]),label="Analytical",linewidth=2,color=:blue);
title!("Frequency = $(1000) Hz");
plt2 = scatter(ang_axis,real.(v_n0),label="BEM",marker=:cross,markersize=2,color=:black);
ylabel!("Re(Vn)"); plot!(ang_axis[perm],real.(v_rAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
plt3 = scatter(ang_axis,real.(vt_sum),label="BEM",marker=:cross,markersize=2,color=:black);
plot!(ang_axis[perm],real.(v_thetaAN_V[perm]),label="Analytical",linewidth=2,color=:blue);
xlabel!("Angle [deg]"); ylabel!("Re(Vt)");
plt4 = plot(plt1,plt2,plt3,layout=(3,1))


writedlm( "ang_axis.csv",  ang_axis, ',')
writedlm( "ang_axisPerm.csv",  ang_axis[perm], ',')
writedlm( "paReal.csv",  real.(pa), ',')
writedlm( "paImag.csv", imag.(pa), ',')
writedlm( "paAbs.csv",  abs.(pa), ',')
writedlm( "vnReal.csv",  real.(v_n0), ',')
writedlm( "vnImag.csv",  imag.(v_n0), ',')
writedlm( "vnAbs.csv",  abs.(v_n0), ',')
writedlm( "vtReal.csv",  real.(vt_sum), ',')
writedlm( "vtImag.csv",  imag.(vt_sum), ',')
writedlm( "vtAbs.csv",  abs.(vt_sum), ',')
writedlm( "paANReal.csv",  real.(pasAN[perm]), ',')
writedlm( "paANImag.csv",  imag.(pasAN[perm]), ',')
writedlm( "paANAbs.csv",  abs.(pasAN[perm]), ',')
writedlm( "vnANeal.csv",  real.(v_rAN_V[perm]), ',')
writedlm( "vnANImag.csv",  imag.(v_rAN_V[perm]), ',')
writedlm( "vnANAbs.csv",  abs.(v_rAN_V[perm]), ',')
writedlm( "vtANReal.csv",  real.(v_thetaAN_V[perm]), ',')
writedlm( "vtANImag.csv",  imag.(v_thetaAN_V[perm]), ',')
writedlm( "vtANAbs.csv",  abs.(v_thetaAN_V[perm]), ',')
close(f)

#==========================================================================================
                            Read condition numbers for specific test case
==========================================================================================#
data_file = "cond_$(2642)DOFs_freq$(Int(1000)).jld2"
f = jldopen(data_file)

condF1 = f["condF1"]
condF4 = f["condF4"]
condF10 = f["condF10"]
condGa = f["condGa"]
condGh = f["condGh"]
condGv = f["condGv"]
condHa = f["condHa"]
condHh = f["condHh"]
condHv = f["condHv"]
close(f)

#==========================================================================================
                            Calculate discretization and FMM approximation error for specific test case
==========================================================================================#
data_file = "resultsFMM1x1_$(494)DOFs_$(Int(100))Hz.JLD2"
f = jldopen(data_file)

paFMM = f["paFMM"]
pa = f["pa"]
pasAN = f["pasAN"]
perm = f["perm"]
ang_axis = f["ang_axis"]
M  = size(paFMM,1)

errorCalc(pa,pasAN,M)
errorCalc(paFMM,pa,M)

#==========================================================================================
                            Export memory for LGM, assembly and solution
==========================================================================================#
i = 0;
tLGM = zeros(Float64, 8,7); # allocations for setting up the building blocks + setting up system matrix (except 1n-iter since there is no system matrix)
tSol = zeros(Float64, 8,7); # allocations for solving system
varInfo = zeros(Float64, 8,7); # memory required for storing LGM + system matrix (except 1n-iter)
varInfoAlloc = zeros(Float64, 8,7); # allocations for setting up system matrix (except 1n-iter)

for freq in [100,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;
        data_file = "memory_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 1n iter
        tLGM[i,j] =  f["memory"].columns[6][1]; # tLGM
        varInfo[i,j] = f["varLGM"]; # size LGM
        #varInfoAlloc[i,j] = 0; # size LGM
    
        # 1n dense
        tLGM[i+2,j] =  f["memory"].columns[6][3] + f["memory"].columns[6][1];  # tLGM + tLGM_dense1x1
        varInfo[i+2,j] = f["varLGM"] + f["varLGMdense1x1"]; # size LGM + size LGM_dense1x1
        tSol[i+2,j] =  f["memory"].columns[6][4] # tDense1x1
        varInfoAlloc[i+2,j] =  f["memory"].columns[6][3]; # tLGM_dense1x1

        # 4n dense 
        tLGM[i+4,j] =  f["memory"].columns[6][5] + f["memory"].columns[6][1]; # tLGM + tLGM_dense4x4
        varInfo[i+4,j] = f["varLGM"] + f["varLGMdense4x4"]; # size LGM + size LGM_dense4x4
        tSol[i+4,j] =  f["memory"].columns[6][6] # tDense4x4
        varInfoAlloc[i+4,j] =  f["memory"].columns[6][5]; #  tLGM_dense4x4

        # 10n dense 
        tLGM[i+6,j] =  f["memory"].columns[6][7] + f["memory"].columns[6][1]; # tLGM + tLGM_dense10x10
        varInfo[i+6,j] = f["varLGM"] + f["varLGMdense10x10"]; # size LGM + size LGM_dense10x10
        tSol[i+6,j] =  f["memory"].columns[6][8] # tDense10x10
        varInfoAlloc[i+6,j] = f["memory"].columns[6][7]; #  tLGM_dense10x10

        close(f);

        if freq ==1000
            data_file = "memoryMaxIter_$(M)DOFs_$(Int(freq))Hz.JLD2";
            f = jldopen(data_file);
            tSol[i,j] =  f["memory"].columns[6][1]; # tIter for maxiter=20 at 1000Hz (from other results since maxiter was wrongly set to 1)
            close(f);
        else
            data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
            f = jldopen(data_file);
            tSol[i,j] = f["runtimes"].columns[6][4]; # tIter for maxiter=inf at 100Hz (inf because less mvps than restart required, from other results since maxiter was wrongly set to 1)
            close(f);
        end
    end
end

# export
writedlm( "tLGM.csv",  tLGM, ',')
writedlm( "tSol.csv",  tSol, ',')
writedlm( "varInfo.csv",  varInfo, ',')
writedlm( "varInfoAlloc.csv",  varInfoAlloc, ',')

#==========================================================================================
                            Export runtimes and accumulated allocations for LGM, assembly and solution for 1n-iter/dense (run on local machine)
==========================================================================================#
i = 0;
runtimesMedian1n = zeros(Float64, 4,7); # median runtime for setting up system matrix (except 1n-iter) + solving system
runtimesMin1n = zeros(Float64, 4,7); # minimum runtime for setting up system matrix (except 1n-iter) + solving system
runtimesMinLGM1n = zeros(Float64, 4,7); # runtime for setting up LGM
runtimesMemory1n = zeros(Float64, 4,7); # accumulated allocations for setting up system matrix (except 1n-iter) + solving system
runtimesMemoryLGM1n = zeros(Float64, 4,7); # allocations for setting up LGM

for freq in [100,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;
        data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 1n iter
        runtimesMedian1n[i,j] = f["runtimes"].columns[3][4];
        runtimesMin1n[i,j] = f["runtimes"].columns[2][4];
        runtimesMinLGM1n[i,j] = f["runtimes"].columns[2][1];
        runtimesMemoryLGM1n[i,j] = f["runtimes"].columns[6][1];
        runtimesMemory1n[i,j] = f["runtimes"].columns[6][4];
    
        # 1n dense
        runtimesMedian1n[i+2,j] = f["runtimes"].columns[3][3];
        runtimesMin1n[i+2,j] = f["runtimes"].columns[2][3];
        runtimesMinLGM1n[i+2,j] = f["runtimes"].columns[2][1];
        runtimesMemoryLGM1n[i+2,j] = f["runtimes"].columns[6][1];
        runtimesMemory1n[i+2,j] = f["runtimes"].columns[6][3];

        close(f);

    end
end

runtimesMedian1n /= 1e9; # from nano seconds to seconds
runtimesMin1n /= 1e9; # from nano seconds to seconds
runtimesMinLGM1n /= 1e9; # from nano seconds to seconds
runtimesMemory1n /= 1e9; # from bytes to gigabytes
runtimesMemoryLGM1n /= 1e9; # from bytes to gigabytes

# export
writedlm( "runtimesMedian1n.csv",  runtimesMedian1n, ',')
writedlm( "runtimesMin1n.csv",  runtimesMin1n, ',')
writedlm( "runtimesMinLGM1n.csv",  runtimesMinLGM1n, ',')
writedlm( "runtimesMemory1n.csv",  runtimesMemory1n, ',')
writedlm( "runtimesMemoryLGM1n.csv",  runtimesMemoryLGM1n, ',')

#==========================================================================================
                            Export runtimes and accumulated allocations for LGM, assembly and solution
==========================================================================================#
i = 0;
runtimesMedian = zeros(Float64, 12,7); # median runtime for setting up system matrix (except 1n-iter) + solving system
runtimesMin = zeros(Float64, 12,7); # minimum runtime for setting up system matrix (except 1n-iter) + solving system
runtimesMinLGM = zeros(Float64, 12,7); # runtime for setting up LGM
runtimesMemory = zeros(Float64, 12,7); # accumulated allocations for setting up system matrix (except 1n-iter) + solving system
runtimesMemoryLGM = zeros(Float64, 12,7); # allocations for setting up LGM

for freq in [100,500,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;
        data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 1n iter
        runtimesMedian[i,j] = f["runtimes"].columns[3][4];
        runtimesMin[i,j] = f["runtimes"].columns[2][4];
        runtimesMinLGM[i,j] = f["runtimes"].columns[2][1];
        runtimesMemoryLGM[i,j] = f["runtimes"].columns[6][1];
        runtimesMemory[i,j] = f["runtimes"].columns[6][4];
    
        # 1n dense
        runtimesMedian[i+3,j] = f["runtimes"].columns[3][3];
        runtimesMin[i+3,j] = f["runtimes"].columns[2][3];
        runtimesMinLGM[i+3,j] = f["runtimes"].columns[2][1];
        runtimesMemoryLGM[i+3,j] = f["runtimes"].columns[6][1];
        runtimesMemory[i+3,j] = f["runtimes"].columns[6][3];

        close(f);

        data_file = "runtimes4x4_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 4n dense 
        runtimesMedian[i+6,j] = f["runtimes"].columns[3][3];
        runtimesMin[i+6,j] = f["runtimes"].columns[2][3];
        runtimesMinLGM[i+6,j] = f["runtimes"].columns[2][1];
        runtimesMemoryLGM[i+6,j] = f["runtimes"].columns[6][1];
        runtimesMemory[i+6,j] = f["runtimes"].columns[6][3];

        close(f);

        data_file = "runtimes10x10_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 4n dense 
        runtimesMedian[i+9,j] = f["runtimes"].columns[3][3];
        runtimesMin[i+9,j] = f["runtimes"].columns[2][3];
        runtimesMinLGM[i+9,j] = f["runtimes"].columns[2][1];
        runtimesMemoryLGM[i+9,j] = f["runtimes"].columns[6][1];
        runtimesMemory[i+9,j] = f["runtimes"].columns[6][3];

        close(f);

    end
end

runtimesMedian /= 1e9; # from nano seconds to seconds
runtimesMin /= 1e9; # from nano seconds to seconds
runtimesMinLGM /= 1e9; # from nano seconds to seconds
runtimesMemory/= 1e9; # from bytes to gigabytes
runtimesMemoryLGM /= 1e9; # from bytes to gigabytes

# export
writedlm( "runtimesMedian.csv",  runtimesMedian, ',')
writedlm( "runtimesMin.csv",  runtimesMin, ',')
writedlm( "runtimesMinLGM.csv",  runtimesMinLGM, ',')
writedlm( "runtimesMemory.csv",  runtimesMemory, ',')
writedlm( "runtimesMemoryLGM.csv",  runtimesMemoryLGM, ',')

#==========================================================================================
                            Export runtimes and accumulated allocations for rhs
==========================================================================================#
i = 0;
runtimesMedianRHS = zeros(Float64, 12,7); # median runtime for setting up rhs
runtimesMinRHS = zeros(Float64, 12,7); # minimum runtime for setting up rhs
runtimesMemoryRHS = zeros(Float64, 12,7); # accumulated allocations for setting up rhs

for freq in [100,500,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;
        data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 1n iter
        runtimesMedianRHS[i,j] = f["runtimes"].columns[3][2];
        runtimesMinRHS[i,j] = f["runtimes"].columns[2][2];
        runtimesMemoryRHS[i,j] = f["runtimes"].columns[6][2];
    
        # 1n dense
        runtimesMedianRHS[i+3,j] = f["runtimes"].columns[3][2];
        runtimesMinRHS[i+3,j] = f["runtimes"].columns[2][2];
        runtimesMemoryRHS[i+3,j] = f["runtimes"].columns[6][2];
        close(f);

        data_file = "runtimes4x4_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 4n dense 
        runtimesMedianRHS[i+6,j] = f["runtimes"].columns[3][2];
        runtimesMinRHS[i+6,j] = f["runtimes"].columns[2][2];
        runtimesMemoryRHS[i+6,j] = f["runtimes"].columns[6][2];
        close(f);

        data_file = "runtimes10x10_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 10n dense 
        runtimesMedianRHS[i+9,j] = f["runtimes"].columns[3][2];
        runtimesMinRHS[i+9,j] = f["runtimes"].columns[2][2];
        runtimesMemoryRHS[i+9,j] = f["runtimes"].columns[6][2];
        close(f);

    end
end

runtimesMedianRHS /= 1e9; # from nano seconds to seconds
runtimesMinRHS /= 1e9; # from nano seconds to seconds
runtimesMemoryRHS/= 1e9; # from bytes to gigabytes

# export 
writedlm( "runtimesMedianRHS.csv",  runtimesMedianRHS, ',')
writedlm( "runtimesMinRHS.csv",  runtimesMinRHS, ',')
writedlm( "runtimesMemoryRHS.csv",  runtimesMemoryRHS, ',')

#==========================================================================================
                            Export runtimes and accumulated allocations for reconstruction of vv
==========================================================================================#
i = 0;
runtimesMedianREC = zeros(Float64, 12,7); # median runtime for setting up reconstruction
runtimesMinREC = zeros(Float64, 12,7); # minimum runtime for setting up reconstruction
runtimesMemoryREC = zeros(Float64, 12,7); # accumulated allocations for reconstruction

for freq in [100,500,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;  
        data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 1n iter
        runtimesMedianREC[i,j] = f["runtimes"].columns[3][5];
        runtimesMinREC[i,j] = f["runtimes"].columns[2][5];
        runtimesMemoryREC[i,j] = f["runtimes"].columns[6][5];
    
        # 1n dense
        runtimesMedianREC[i+3,j] = f["runtimes"].columns[3][5];
        runtimesMinREC[i+3,j] = f["runtimes"].columns[2][5];
        runtimesMemoryREC[i+3,j] = f["runtimes"].columns[6][5];
        close(f);

        data_file = "runtimes4x4_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 4n dense 
        runtimesMedianREC[i+6,j] = f["runtimes"].columns[3][4];
        runtimesMinREC[i+6,j] = f["runtimes"].columns[2][4];
        runtimesMemoryREC[i+6,j] = f["runtimes"].columns[6][4];
        close(f);

        data_file = "runtimes10x10_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        # 10n dense 
        runtimesMedianREC[i+9,j] = f["runtimes"].columns[3][4];
        runtimesMinREC[i+9,j] = f["runtimes"].columns[2][4];
        runtimesMemoryREC[i+9,j] = f["runtimes"].columns[6][4];
        close(f);

    end
end

runtimesMedianREC /= 1e9; # from nano seconds to seconds
runtimesMinREC /= 1e9; # from nano seconds to seconds
runtimesMemoryREC/= 1e9; # from bytes to gigabytes

# export
writedlm( "runtimesMedianREC.csv",  runtimesMedianREC, ',')
writedlm( "runtimesMinREC.csv",  runtimesMinREC, ',')
writedlm( "runtimesMemoryREC.csv",  runtimesMemoryREC, ',')

#==========================================================================================
                            Export GMRES logs for solution and reconstruction
==========================================================================================#
i = 0;
numMVPS = zeros(Float64, 3,7); # number matrix vector products performed in outer GMRES loop to solve 1n-iter system
relErrorIter = zeros(Float64, 3,7); # relative error in outer GMRES loop to solve 1n-iter system
numMVPSRec = zeros(Float64, 3,7); # number matrix vector products performed in reconstruction step to compute dpa (based on 1n-direct)
relErrorIterRec = zeros(Float64, 3,7); # relative error in outer GMRES loop in reconstruction step to compute dpa (based on 1n-direct)

for freq in [100,500,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;
        data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);
        hist_pa = f["hist_pa"]
        numMVPS[i,j] = hist_pa.mvps;
        relErrorIter[i,j] = hist_pa.data[:reltol];
        close(f);

        data_file = "results1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);
        hist_dpa = f["hist_dpaDense"]
        numMVPSRec[i,j] = hist_dpa.mvps;
        relErrorIterRec[i,j] = hist_dpa.data[:reltol];
        close(f);

    end
end

# export
writedlm( "numMVPS.csv",  numMVPS, ',');
writedlm( "relErrorIter.csv",  relErrorIter, ',');
writedlm( "numMVPSRec.csv",  numMVPS, ',');
writedlm( "relErrorIterRec.csv",  relErrorIter, ',');

#==========================================================================================
                            Export discretization and FMM approximation error wrt pa for quadratic physics discretization
==========================================================================================#
i = 0;
epsDiscr = zeros(Float64, 2,9);
epsFMM = zeros(Float64, 2,9);

for freq in [100,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062,69038,153978] # quadratic meshes
        j +=1;

        data_file = "resultsFMM1x1_$(M)DOFs_$(Int(freq))Hz.JLD2"
        f = jldopen(data_file);

        paFMM = f["paFMM"]
        pa = f["pa"]
        pasAN = f["pasAN"]

        # 1n-iter discretization error in surface nodes compared to analytical solution
        epsDiscr[i,j] =  errorCalc(pa,pasAN,M)
        # 1n-FMM approximation error in surface nodes compared to 1n-iter
        epsFMM[i,j] = errorCalc(paFMM,pa,M)  

        close(f);

    end
end

# export
writedlm( "epsDiscr.csv",  epsDiscr, ',');
writedlm( "epsFMM.csv",  epsFMM, ',');


#==========================================================================================
                            Export discretization and FMM approximation error wrt pa for linear physics discretization
==========================================================================================#
i = 0;
epsDiscrLin = zeros(Float64, 2,9);
epsFMMLin = zeros(Float64, 2,9);

for freq in [100,1000]
    i+=1;
    j = 0;
    for M in [125,234,422,662,1399,2255,3267,17261,38496] # linear meshes
        
        j +=1;

        data_file = "resultsFMM1x1_$(M)DOFs_$(Int(freq))Hz.JLD2"
        f = jldopen(data_file);

        paFMM = f["paFMM"]
        pa = f["pa"]
        pasAN = f["pasAN"]

        # 1n-iter discretization error in surface nodes compared to analytical solution
        epsDiscrLin[i,j] =  errorCalc(pa,pasAN,M)
        # 1n-FMM approximation error in surface nodes compared to 1n-iter
        epsFMMLin[i,j] = errorCalc(paFMM,pa,M)  

        close(f);

    end
end

# export
writedlm( "epsDiscrLin.csv",  epsDiscrLin, ',');
writedlm( "epsFMMLin.csv",  epsFMMLin, ',');

#==========================================================================================
                            Export discretization error wrt pa and vv
==========================================================================================#
# RESULTS rel error total
i = 0;
pa_error = zeros(Float64, 12,7);
vn_error = zeros(Float64, 12,7);
vt_error = zeros(Float64, 12,7);

for freq in [100,500,1000]
    i+=1;
    j = 0;
    for M in [494,930,1682,2642,5590,9014,13062]
        j +=1;
        data_file = "runtimes1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        paDense = f["paDense"]
        pa = f["pa"]
        hist_pa = f["hist_pa"]
        pasAN = f["pasAN"]
        v_rAN_V = f["v_rAN_V"]
        v_thetaAN_V = f["v_thetaAN_V"]
        v_n0 = f["v_n0"]
        vt_sum = f["vt_sum"]

        # 1n iter
        pa_error[i,j] = errorCalc(pa,pasAN,M); # discretization error in surface nodes compared to analytical solution
        vn_error[i,j] = errorCalc(v_n0,v_rAN_V,M);
        vt_error[i,j] = errorCalc(vt_sum,v_thetaAN_V,M);
    
        # 1n dense
        pa_error[i+3,j] = errorCalc(paDense,pasAN,M);

        close(f);

        data_file = "results1x1_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);
        vn_error[i+3,j] = errorCalc(v_n0Dense,v_rAN_V,M);
        vt_error[i+3,j] = errorCalc(vt_sumDense,v_thetaAN_V,M);
        close(f);

        data_file = "runtimes4x4_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        pa = f["pa"]
        pasAN = f["pasAN"]
        v_rAN_V = f["v_rAN_V"]
        v_thetaAN_V = f["v_thetaAN_V"]
        v_n0 = f["v_n0"]
        vt_sum = f["vt_sum"]

        # 4n dense 
        pa_error[i+6,j] = errorCalc(pa,pasAN,M);
        vn_error[i+6,j] = errorCalc(v_n0,v_rAN_V,M);
        vt_error[i+6,j] = errorCalc(vt_sum,v_thetaAN_V,M);

        close(f);

        data_file = "runtimes10x10_$(M)DOFs_$(Int(freq))Hz.JLD2";
        f = jldopen(data_file);

        pa = f["pa"]
        pasAN = f["pasAN"]
        v_rAN_V = f["v_rAN_V"]
        v_thetaAN_V = f["v_thetaAN_V"]
        v_n0 = f["v_n0"]
        vt_sum = f["vt_sum"]

        # 10n dense 
        pa_error[i+9,j] = errorCalc(pa,pasAN,M);
        vn_error[i+9,j] = errorCalc(v_n0,v_rAN_V,M);
        vt_error[i+9,j] = errorCalc(vt_sum,v_thetaAN_V,M);
 
        close(f);

    end
end

# export
writedlm( "pa_error.csv",  pa_error, ',');
writedlm( "vn_error.csv",  vn_error, ',');
writedlm( "vt_error.csv",  vt_error, ',');



