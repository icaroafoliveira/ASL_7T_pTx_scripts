# scripts to run ASL analysis
# Icaro Oliveira, 2023
# ===========

# Packages
using MriResearchTools
using Statistics


function savePerf(foldername, dir, filename)
    
    datc = readmag(joinpath(foldername, dir, filename));
    dimx, dimy, dimz, dimT= size(datc.raw);
    control_idx = Array(1:2:dimT);
    label_idx = Array(2:2:dimT);
        
    control_dat = datc.raw[:,:,:,control_idx];
    label_dat = datc.raw[:,:,:,label_idx];

    # Mean ASL control
    meanControl = mean(control_dat, dims=4);
    replace!(meanControl, Inf=>0);
    replace!(meanControl, NaN=>0);
    meanControl[meanControl.<=0].=0;

    # Mean ASL label
    meanLabel = mean(label_dat, dims=4);
    replace!(meanLabel, Inf=>0);
    replace!(meanLabel, NaN=>0);
    meanLabel[meanLabel.<=0].=0;


    meanPerf = meanControl.-meanLabel;
    replace!(meanPerf, Inf=>0);
    replace!(meanPerf, NaN=>0);
    meanPerf[meanPerf.<=0].=0;

    
    # saving the data

    #changing filename
    control_fname = string(splitext(filename)[1], "_meanControl");
    meanPerf_fname = string(splitext(filename)[1], "_meanPerf");
    

    savenii(meanControl, control_fname, joinpath(foldername, dir), header(datc))
    savenii(meanPerf, meanPerf_fname, joinpath(foldername, dir), header(datc))
    #savenii(meanLabel, "meanLabel", foldername, header(datc))

end


function seriesASL(foldername, dir, filename)
    
    datc = readmag(joinpath(foldername, dir, filename));
    dimx, dimy, dimz, dimT= size(datc.raw);
    control_idx = Array(1:2:dimT);
    label_idx = Array(2:2:dimT);
        
    control_dat = datc.raw[:,:,:,control_idx];
    label_dat = datc.raw[:,:,:,label_idx];

    dyn_perf = zeros(dimx, dimy, dimz, Int(dimT/2));

    for i =1:size(dyn_perf, 4)
        dyn_perf[:,:,:,i] = control_dat[:,:,:,i] - label_dat[:,:,:,i];
    end

    replace!(dyn_perf, Inf=>0);
    replace!(dyn_perf, NaN=>0);
    
    # saving the data
    #changing filename
    dyn_fname = string(splitext(filename)[1], "_dyn_perf");

    savenii(dyn_perf, dyn_fname, joinpath(foldername, dir), header(datc))
end

### for fMRI

function savefMRIASL(foldername, dir, filename)
    
    datc = readmag(joinpath(foldername, dir, filename));
    dimx, dimy, dimz, dimT= size(datc.raw);
    control_idx = Array(1:2:dimT);
    label_idx = Array(2:2:dimT);
        
    control_dat = datc.raw[:,:,:,control_idx];
    label_dat = datc.raw[:,:,:,label_idx];

    dyn_perf = zeros(dimx, dimy, dimz, Int(dimT/2));

    for i =1:size(dyn_perf, 4)
        dyn_perf[:,:,:,i] = control_dat[:,:,:,i] - label_dat[:,:,:,i];
    end

    replace!(dyn_perf, Inf=>0);
    replace!(dyn_perf, NaN=>0);

    # Mean Perf
    meanPerf = mean(dyn_perf, dims=4);
    replace!(meanPerf, Inf=>0);
    replace!(meanPerf, NaN=>0);
    meanPerf[meanPerf.<=0].=0;

    # Std ASL label
    stdPerf = std(dyn_perf, dims=4);
    replace!(stdPerf, Inf=>0);
    replace!(stdPerf, NaN=>0);
    stdPerf[stdPerf.<=0].=0;

    # Tsnr Map

    PtSNR = meanPerf./stdPerf;
    replace!(PtSNR, Inf=>0);
    replace!(PtSNR, NaN=>0);
      

    # saving the data
  

    #changing filename
    dyn_fname = string(splitext(filename)[1], "_dyn_perf");
    meanPerf_fname = string(splitext(filename)[1], "_meanPerf");
    tsnr_filename = string(splitext(filename)[1], "_tSNR");

    savenii(dyn_perf, dyn_fname, joinpath(foldername, dir), header(datc))
    savenii(meanPerf, meanPerf_fname, joinpath(foldername, dir), header(datc))
    savenii(PtSNR, tsnr_filename, joinpath(foldername, dir), header(datc))

end

### ====== Perfusion extraction function
    
    """
    This function is used specific to the this project to extract perfusion values using 2 different masks

    Perf_region(foldername, dir, filename, mask1_foldername, mask2_foldername, dir_anat)

TBW
"""
function Perf_region(foldername, dir, filename, mask1_foldername, mask2_foldername, dir_anat)
    # mask1 GM LEFT hemisphere
    # mask2 GM RIGHT hemisphere

    datc_vol = readmag(joinpath(foldername, dir, filename));
    datc = datc_vol.raw;
    replace!(datc, NaN=>0);
    #datc = datc .+ 40; # non-zeros

    datc1_vol = readmag(joinpath(foldername, dir_anat, mask1_foldername));
    datc1 = datc1_vol.raw;

    datc2_vol = readmag(joinpath(foldername, dir_anat, mask2_foldername));
    datc2 = datc2_vol.raw;

    # creating masks
    GM_mask1 = datc1 .> 0.8;
    GM_mask2 = datc2 .> 0.8;

    GM_L_perf = datc.*GM_mask1;     
    replace!(GM_L_perf, Inf=>0);
    replace!(GM_L_perf, NaN=>0);
    GM_L_perf[GM_L_perf.<=0].=0;
    mean_GM_L_perf = mean(GM_L_perf[GM_L_perf.>0]);
    sd_GM_L_perf = std(GM_L_perf[GM_L_perf.>0]);

    sSNR_GM_L = mean_GM_L_perf/sd_GM_L_perf;

    #GM_L_values = [mean_GM_L_perf, sd_GM_L_perf];


    GM_R_perf = datc.*GM_mask2;
    replace!(GM_R_perf, Inf=>0);
    replace!(GM_R_perf, NaN=>0);
    GM_R_perf[GM_R_perf.<=0].=0;
    mean_GM_R_perf = mean(GM_R_perf[GM_R_perf.>0]);
    sd_GM_R_perf = std(GM_R_perf[GM_R_perf.>0]);

    sSNR_GM_R = mean_GM_R_perf/sd_GM_R_perf;
    #GM_R_values = [mean_GM_R_perf, sd_GM_R_perf];

    # global_perf = mean(datc[datc.>0]);

    perf_values = [mean_GM_L_perf, mean_GM_R_perf];

    return(perf_values)

end

function cov(foldername, dir, filename, dir_anat, mask1_foldername, mask2_foldername)

    datc = readmag(joinpath(foldername, dir, filename));

    cov_map = std(datc,dims=4)./mean(datc,dims=4);
    replace!(cov_map, Inf=>0);
    replace!(cov_map, NaN=>0);
    cov_map[cov_map.<=0].=0;
    cov_map[cov_map.>=10].=0;

    datc1 = readmag(joinpath(foldername, dir_anat, mask1_foldername));
    GM_mask1 = datc1.raw; 

    datc2 = readmag(joinpath(foldername, dir_anat, mask2_foldername));
    GM_mask2 = datc2.raw;

    masked_cov_c1 = cov_map .* GM_mask1;
    masked_cov_c2 = cov_map .* GM_mask2;

    cov_values_c1 = mean(masked_cov_c1[masked_cov_c1.>0]);
    cov_values_c2 = mean(masked_cov_c2[masked_cov_c2.>0]);

    cov_values = [cov_values_c1, cov_values_c2];

    return cov_values
end

### ======= CoV Maps

function cov_map(foldername, dir, filename )
    datc = readmag(joinpath(foldername, dir, filename));
    dimx, dimy, dimz, dimT= size(datc.raw);
    control_idx = Array(1:2:dimT);
    label_idx = Array(2:2:dimT);
        
    control_dat = datc.raw[:,:,:,control_idx];
    label_dat = datc.raw[:,:,:,label_idx];

    control_dat = datc.raw[:,:,:,control_idx];
    label_dat = datc.raw[:,:,:,label_idx];

    perf_dyn=zeros(dimx, dimy, dimz, Int64(dimT/2));

    for i=1:Int64(dimT/2)
        perf_dyn[:,:,:,i] = control_dat[:,:,:,i] - label_dat[:,:,:,i];
    end

    cov_map = std(perf_dyn,dims=4)./mean(perf_dyn,dims=4);
    replace!(cov_map, Inf=>0);
    replace!(cov_map, NaN=>0);
    cov_map[cov_map.<=0].=0;

    cov_fname = string(splitext(filename)[1], "_CoV");

    savenii(cov_map, cov_fname, joinpath(foldername, dir), header(datc))

end


### ======= Create masks from GRAY/WHITE matter segmentation

function create_mask(foldername, fn_mask_classic, fn_mask_multix)
    c_mask_classic = readmag(joinpath(foldername, fn_mask_classic));
    c_mask_multix = readmag(joinpath(foldername, fn_mask_multix));

    datc_mask_classic = c_mask_classic.raw;
    datc_mask_multix = c_mask_multix.raw;

    mask_classic = datc_mask_classic .>30;
    mask_multix = datc_mask_multix .>30;

    final_mask = mask_classic .* mask_multix;

    outputdir = foldername;
    savenii(final_mask,"final_mask",outputdir,header(c_mask_classic))

end

### ======== sSNR

function sSNR_perf(main_folder, dir_asl, fn_perf, dir_anat, fn_mask1, fn_mask2)
    # load nifti
    c_perf = readmag(joinpath(main_folder, dir_asl, fn_perf));
    c_mask1 = readmag(joinpath(main_folder, dir_anat, fn_mask1));
    c_mask2 = readmag(joinpath(main_folder, dir_anat, fn_mask2));

    datc_perf = c_perf.raw;
    datc_mask1 = c_mask1.raw;
    datc_mask2 = c_mask2.raw;

    perf_mask1 = datc_perf .* datc_mask1;
    perf_mask2 = datc_perf .* datc_mask2;
    
    mean_perf_mask1 = mean(perf_mask1[perf_mask1.>0]);
    mean_perf_mask2 = mean(perf_mask2[perf_mask2.>0]);

    sd_perf_mask1 = std(perf_mask1[perf_mask1.>0]);
    sd_perf_mask2 = std(perf_mask2[perf_mask2.>0]);

    # N_mask1 = size(perf_mask1[perf_mask1.>0],1);
    # N_mask2 = size(perf_mask2[perf_mask2.>0],1);

    sSNR_mask1 = mean_perf_mask1/(sd_perf_mask1);#* sqrt(N_mask1/(N_mask1-1)));
    sSNR_mask2 = mean_perf_mask2/(sd_perf_mask2);#* sqrt(N_mask2/(N_mask2-1)));
   
    return sSNR_mask1, sSNR_mask2
    
end

function cnr_calc(fn_img, fn_mask1, fn_mask2)

    datc = readmag(fn_img);

    datc1 = readmag(fn_mask1);
    datc2 = readmag(fn_mask2);

    masked_datc1 = datc[:,:,:,1] .* datc1;
    masked_datc2 = datc[:,:,:,1] .* datc2;

    μGM = mean(masked_datc1[datc1.>0]);
    μWM = mean(masked_datc2[datc2.>0]);

    sigmaGM = std(masked_datc1[datc1.>0]);
    sigmaWM = std(masked_datc2[datc2.>0]);

    return CNR = abs(μGM - μWM) / √(sigmaGM^2 + sigmaWM^2);

end