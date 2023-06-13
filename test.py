import biolib
import os
import subprocess

signalp_6 = biolib.load('DTU/SignalP_6')
subprocess.run('biolib run DTU/SignalP_6 --fastafile vacsol_ai/sequences.fasta --output_dir output')
signalp6_job = signalp_6.cli(args='--fasta /content/sequences.fasta --output_dir output')
signalp6_job.save_files



