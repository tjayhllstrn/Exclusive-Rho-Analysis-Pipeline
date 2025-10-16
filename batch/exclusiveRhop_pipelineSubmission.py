#!/usr/bin/env python3
import subprocess
from pathlib import Path

# === USER CONFIGURATION ===
input_dirs = [
    Path("/lustre24/expphy/cache/clas12/rg-a/production/recon/spring2019/torus-1/pass2/dst/train/nSidis"),
    Path("/lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass2/main/train/nSidis"),
]

output_base = Path("/w/hallb-scshelf2102/clas12/users/tjhellst/clas-ana-scaffold-tyler/out")
output_names = ["pippi0_spring2019_in_pass2","pippi0_fall2018_in_pass2"]

config_file = "pippi0_RGAinbending_zbinning"

dep_list = []  # Holds hadd job IDs for the final step



def sbatch_submit(args,dry_run=False):
    """Run sbatch and return the job ID."""
    print(f"Submitting: {' '.join(args)}")
    
    if dry_run:
        cmd_str = " ".join(args)
        fake_job_id = str(abs(hash(cmd_str)) % 100000)  # fake unique id
        print(f"  → Fake JobID: {fake_job_id}")
        print("")
        return fake_job_id
        
    result = subprocess.run(args, capture_output=True, text=True, check=True)
    # Typical sbatch output: "Submitted batch job 12345"
    job_id = result.stdout.strip().split()[-1]
    print(f"  → JobID: {job_id}")
    return job_id

for i,input_dir in enumerate(input_dirs):
    tag = input_dir.name
    output_dir = output_base / output_names[i]
    output_dir.mkdir(parents=True, exist_ok=True)
    n_files = len(list(input_dir.glob("nSidis_*")))
    if n_files == 0:
        print(f"No files found in {input_dir}, skipping…")
        continue
    array_spec = f"0-{n_files-1}"

    # 1️⃣ Submit array job for this directory
    convert_args = [
        "sbatch", 
        f"--array={array_spec}",
        "pipeline0_run_single_hipo.sbatch",
        str(input_dir), str(output_dir)
    ]
    jid_convert = sbatch_submit(convert_args)

    # 2️⃣ Submit hadd job, dependent on the array completion
    hadd_args = [
        "sbatch", f"--dependency=afterok:{jid_convert}",
        "pipeline1_run_hadd.sbatch", str(output_dir),f"{output_names[i]}.root"
    ]
    jid_hadd = sbatch_submit(hadd_args)

    dep_list.append(jid_hadd)

# 3️⃣ Submit final AsymmetryFitting job after all hadds finish
dependency_str = ":".join(dep_list)
merged_outputs=[]
for i,d in enumerate(input_dirs):
    merged_outputs.append(f"{output_base}/{output_names[i]}/{output_names[i]}.root")

fit_args = [
    "sbatch",
    f"--dependency=afterok:{dependency_str}",
    "pipeline2_run_AsymmetryFitting.sbatch", config_file
] + merged_outputs

sbatch_submit(fit_args)