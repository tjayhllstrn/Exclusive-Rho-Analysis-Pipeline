#!/usr/bin/env python3
import subprocess
from pathlib import Path

# === USER CONFIGURATION ===
input_dirs = [Path("/lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis"),Path("/lustre24/expphy/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/nSidis")]#,Path("/lustre24/expphy/cache/clas12/rg-a/production/montecarlo/clasdis_pass2/fa18_inb")
input_file_base = 'nSidis_*.hipo' #corresponds with the file names that you want to run in input_dirs

output_base = Path("/w/hallb-scshelf2102/clas12/users/tjhellst/Exclusive-Rho-Analysis-Pipeline/out")
output_names = ["pippi0_fall2018_in_pass2","pippi0_spring2019_in_pass2"] # must correspond with the order of the entries in input_dirs
final_root_output_name = "pippi0_merged_in_pass2"

config_file = "config/pippi0_RGAinbending_tbinning.txt"

dep_list = []  # Holds hadd job IDs for the final step

fit_types = ["MhMLM","MxMLM","MhChi2","MxChi2"]

def sbatch_submit(args,dry_run=True):
    """Run sbatch and return the job ID."""
    print(f"Submitting: {' '.join(args)}")
    
    if dry_run:
        # Extract the script name and its arguments
        script_idx = next(i for i, arg in enumerate(args) if arg.endswith('.sbatch'))
        script_name = args[script_idx]
        script_args = args[script_idx + 1:]  # Arguments after the script name
        
        # Run the script directly with bash for testing
        bash_cmd = ['bash', script_name] + script_args
        print(f"  → Dry run: executing {' '.join(bash_cmd)}")
        
#        try:
#            result = subprocess.run(bash_cmd, capture_output=True, text=True, check=True)
#            print(f"  STDOUT:\n{result.stdout}")
#            if result.stderr:
#                print(f"  STDERR:\n{result.stderr}")
#            
            # Return fake job ID
        fake_job_id = str(abs(hash(' '.join(args))) % 100000)
        print(f"  → Fake JobID: {fake_job_id}\n")
        return fake_job_id
#        except subprocess.CalledProcessError as e:
#            print(f"  ERROR: Script failed with exit code {e.returncode}")
#            print(f"  STDOUT: {e.stdout}")
#            print(f"  STDERR: {e.stderr}")
#            raise
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

for i,input_dir in enumerate(input_dirs):
    tag = input_dir.name
    output_dir = output_base / output_names[i]
    output_dir.mkdir(parents=True, exist_ok=True)
    n_files = len(list(input_dir.glob(input_file_base)))
    if n_files == 0:
        print(f"No files found in {input_dir}, skipping…")
        continue
    array_spec = f"0-{n_files-1}"

    # 1️⃣ Submit array job for this directory
    convert_args = [
        "sbatch", 
        f"--array={array_spec}",
        "pipeline0_run_single_hipo.sbatch",
        str(input_dir), str(output_dir),input_file_base
    ]
    jid_convert = sbatch_submit(convert_args)

    # 2️⃣ Submit hadd job, dependent on the array completion
    hadd_args = [
        "sbatch", f"--dependency=afterany:{jid_convert}",
        "pipeline1_run_hadd.sbatch", str(output_dir),f"{output_names[i]}.root"
    ]
    jid_hadd = sbatch_submit(hadd_args)

    dep_list.append(jid_hadd)

# 3️⃣ Submit final AsymmetryFitting job after all hadds finish
dependency_str = ":".join(dep_list)
merged_outputs=[]
for i,d in enumerate(input_dirs):
    merged_outputs.append(f"{output_base}/{output_names[i]}/{output_names[i]}.root")


if len(input_dirs) > 1:
    merged_root_output_dir = f"{output_base}/{final_root_output_name}"
    hadd_merge_args = [
        "sbatch",
        f"--dependency=afterok:{dependency_str}",
        "pipeline2_run_hadd.sbatch"
    ] + merged_outputs + [merged_root_output_dir]
    
    jid_final_hadd = sbatch_submit(hadd_merge_args)
    
    # Update dependency to use the final merged file
    dependency_str = jid_final_hadd
    merged_outputs = [f"{merged_root_output_dir}/{final_root_output_name}.root"]


for typ in fit_types:
    fit_args = [
        "sbatch",
        f"--dependency=afterok:{dependency_str}",
        "pipeline3_run_AsymmetryFitting.sbatch", f"{final_root_output_name}.root", f"{final_root_output_name}/"
        ,typ,config_file]
    
    sbatch_submit(fit_args)
