#!/usr/bin/env python3
import subprocess
from pathlib import Path

# === USER CONFIGURATION ===
input_dir = Path("/volatile/clas12/users/tjhellst/ExclusiveRhoPlus_RGA_processedData/pippi0_rgaMC_in_fa18_pass2")
Asymmetry_types = ["A10","Alin_t"] # the types of asymmetries you want to inject
bkg_asym = "A01"


config_file = "config/pippi0_RGAinbending_tbinning.txt"
fit_types = ["MhMLM","MhChi2","MxMLM","MxChi2"]

dep_list = []  # Holds hadd job IDs for the final step


def sbatch_submit(args,dry_run=False):
    """Run sbatch and return the job ID."""
    print(f"Submitting: {' '.join(args)}")
    
    if dry_run:
        fake_job_id = str(abs(hash(' '.join(args))) % 100000)
        print(f"  → Fake JobID: {fake_job_id}\n")
        return fake_job_id
    #when it is not a dry_run...
    try:
        result = subprocess.run(args, capture_output=True, text=True, check=True)
        if result.stderr:
            print(f"  STDERR: {result.stderr}")
        
        # Typical sbatch output: "Submitted batch job 12345"
        job_id = result.stdout.strip().split()[-1]
        print(f"  → JobID: {job_id}")
        return job_id
    except subprocess.CalledProcessError as e:
        print(f"  ERROR: Command failed with exit code {e.returncode}")
        print(f"  STDOUT: {e.stdout}")
        print(f"  STDERR: {e.stderr}")
        raise

for i,A_type in enumerate(Asymmetry_types):
    tag = input_dir.name
    output_name = tag + "_injAsym_" + A_type
    output_dir = input_dir.parent / output_name
    


    n_files = len(list(input_dir.glob("*.root")))
    if n_files == 0:
        print(f"No files found in {input_dir}, skipping…")
        continue
    array_spec = f"0-{n_files-1}"

    # 1️⃣ Submit array job for this directory
    convert_args = [
        "sbatch", 
        f"--array={array_spec}",
        "submit_injectionStudy.sbatch",
        str(input_dir), A_type, bkg_asym
    ]
    jid_convert = sbatch_submit(convert_args)

    # 2️⃣ Submit hadd job, dependent on the array completion
    hadd_args = [
        "sbatch", f"--dependency=afterany:{jid_convert}",
        "pipeline1_run_hadd.sbatch", str(output_dir),f"{output_name}.root"
    ]
    jid_hadd = sbatch_submit(hadd_args)

    for typ in fit_types:
        fit_args = [
            "sbatch",
            f"--dependency=afterok:{jid_hadd}",
            "pipeline3_run_AsymmetryFitting.sbatch", f"{output_name}.root", f"{output_name}/"
            ,typ,config_file]
        
        sbatch_submit(fit_args)
